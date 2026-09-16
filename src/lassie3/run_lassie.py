"""Drive Lassie v1 by building its config in memory and calling its search."""

from __future__ import annotations

import json
import logging
import os
import shutil
from pathlib import Path

import numpy as np

from pyrocko import model as pyrocko_model
from pyrocko import util
from pyrocko.gf import TPDef
from pyrocko.gf import meta as gf_meta

from lassie import common as lassie_common
from lassie import config as lassie_config
from lassie import core as lassie_core
from lassie import grid as lassie_grid
from lassie import ifc as lassie_ifc
from lassie import receiver as lassie_receiver
from lassie import shifter as lassie_shifter
from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

APPROACH = "lassie"

#: Placed far above any attainable image-function level, so the calibration pass
#: stacks the whole window without ever triggering a detection.
CALIBRATION_THRESHOLD = 1e12


def build_config(
    settings: RunSettings,
    threshold: float,
    save_figures: bool = True,
) -> lassie_config.Config:
    """Assemble the Lassie v1 `Config` object the CLI would otherwise parse from YAML.

    Lassie stacks a characteristic function per phase: an STA/LTA onset for P and
    an envelope wave-packet detector for S, each shifted to the grid with
    Cake travel times through the shared 1D model.
    """
    grid = settings.grid
    run_dir = settings.run_dir(APPROACH)

    earthmodel = lassie_common.CakeEarthmodel(
        id="west_bohemia",
        earthmodel_1d=_load_earthmodel(settings.velocity_model_file),
    )

    fmin, fmax = settings.bandpass

    # Lassie's pile reads every channel in the staged files, so the BH? (20 Hz)
    # and LH? (1 Hz) streams reach the IFCs and their 30 Hz lowpass exceeds the
    # Nyquist. Select by NSLC here, mirroring Qseek's `channel_selector`.
    channel_filter = lassie_ifc.TraceSelector(
        white_list_regex=[f"*.*.*.{channel}" for channel in settings.channels]
    )

    onset_p = lassie_ifc.OnsetIFC(
        name="P",
        weight=30.0,
        fmin=fmin,
        fmax=fmax,
        short_window=0.2,
        window_ratio=8.0,
        fsmooth=1.0,
        fnormalize=0.1,
        shifter=lassie_shifter.CakePhaseShifter(
            timing=gf_meta.Timing("{stored:p}"),
            earthmodel_id=earthmodel.id,
        ),
        trace_selector=channel_filter,
    )

    # S energy in this band is broader and slower, so the packet detector runs
    # narrower and smoother than the P onset.
    packet_s = lassie_ifc.WavePacketIFC(
        name="S",
        weight=1.0,
        fmin=fmin,
        fmax=min(fmax, 15.0),
        fsmooth=0.5,
        shifter=lassie_shifter.CakePhaseShifter(
            timing=gf_meta.Timing("{stored:s}"),
            earthmodel_id=earthmodel.id,
        ),
        trace_selector=channel_filter,
    )

    return lassie_config.Config(
        receivers=_receivers(settings),
        data_paths=[str(settings.stage_dir)],
        run_path=str(run_dir / f"{settings.run_name}.turd"),
        tmin=settings.tmin.timestamp(),
        tmax=settings.tmax.timestamp(),
        grid=lassie_grid.Carthesian3DGrid(
            lat=grid.lat,
            lon=grid.lon,
            # Lassie's x is north and y is east, the transpose of Qseek's
            # (east, north) ordering.
            xmin=grid.north_bounds[0],
            xmax=grid.north_bounds[1],
            ymin=grid.east_bounds[0],
            ymax=grid.east_bounds[1],
            zmin=grid.depth_bounds[0],
            zmax=grid.depth_bounds[1],
            dx=grid.lassie_spacing,
            dy=grid.lassie_spacing,
            dz=grid.lassie_spacing,
        ),
        image_function_contributions=[onset_p, packet_s],
        sharpness_normalization=False,
        detector_threshold=threshold,
        save_figures=save_figures,
        # Gaps would otherwise abort the run; the selected day has none, so this
        # only guards against a partially staged archive.
        fill_incomplete_with_zeros=True,
        tabulated_phases=[
            TPDef(id="p", definition="P,p"),
            TPDef(id="s", definition="S,s"),
        ],
        earthmodels=[earthmodel],
        cache_path=str(run_dir / "cache"),
        # Lassie holds nodes x samples x 8 bytes per window unless told to
        # chunk; a 1 km grid (~55k nodes) would need ~6.5 GB per window.
        # Chunking is exact, it only bounds memory.
        stacking_blocksize=1000 if _node_count(grid) > 20_000 else None,
    )


def _node_count(grid) -> int:
    east, north, depth = grid.extents()
    step = grid.lassie_spacing
    return int((east / step + 1) * (north / step + 1) * (depth / step + 1))


def _receivers(settings: RunSettings) -> list:
    """Receivers at their true elevation, built from the shared stations file.

    Lassie's own `stations_path` loader sets each receiver's `z` from
    `Station.depth`, which is the sensor's depth below ground (0 to 34 m here),
    not its elevation, so it would place every station at sea level while
    Qseek traces from the real 510-870 m. `z` is down-positive from sea
    level, the convention Lassie's Cake shifter passes to `zstop`.
    """
    return [
        lassie_receiver.Receiver(
            codes=station.nsl(),
            lat=station.lat,
            lon=station.lon,
            z=station.depth - station.elevation,
        )
        for station in pyrocko_model.load_stations(str(settings.stations_file))
    ]


def _load_earthmodel(path: Path) -> object:
    """Read an `nd` file into the layered model Lassie's Cake shifter expects."""
    from pyrocko import cake

    return cake.load_model(str(path), format="nd")


def stacked_maximum(rundir: Path) -> np.ndarray:
    """Concatenate the per-window stack maxima Lassie writes to `ifm/`.

    Each window file holds two traces: the maximum over the grid (empty location
    code) and the index of the grid node it came from (location `i`). Only the
    former describes the detector level.
    """
    from pyrocko import io

    values = []
    for path in sorted(rundir.glob("ifm/*.mseed")):
        values.extend(
            trace.ydata for trace in io.load(str(path)) if trace.nslc_id[2] == ""
        )
    if not values:
        raise RuntimeError(f"no image function traces under {rundir / 'ifm'}")
    return np.concatenate(values)


def calibrate_threshold(rundir: Path, mad_multiple: float) -> float:
    """Place the detector threshold `mad_multiple` MADs above the typical level.

    Median and MAD are used rather than mean and standard deviation because the
    image function is dominated by a stable noise floor with a long upper tail;
    the earthquakes we are trying to detect sit in that tail and would otherwise
    inflate the very statistic meant to describe the background.
    """
    values = stacked_maximum(rundir)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    threshold = median + mad_multiple * mad
    logger.info(
        "image function: median %.2f, MAD %.3f, max %.2f -> threshold %.2f",
        median,
        mad,
        values.max(),
        threshold,
    )
    return threshold


def run(settings: RunSettings, force: bool = False) -> Path:
    """Execute the Lassie v1 search and return its run directory.

    With an "auto" threshold this stacks twice: once with detection effectively
    disabled to measure the image function, then again at the calibrated level.
    Lassie recomputes the stack rather than reusing the saved one, and stacking
    is most of each window's runtime, so pass an explicit threshold to skip the
    calibration pass once a run has established the level.
    """
    run_dir = settings.run_dir(APPROACH)
    run_dir.mkdir(parents=True, exist_ok=True)
    turd = run_dir / f"{settings.run_name}.turd"

    util.setup_logging("lassie3", "info")

    threshold = settings.lassie_detector_threshold
    if threshold == "auto":
        logger.info("calibration pass: measuring the image function")
        _search(settings, turd, CALIBRATION_THRESHOLD, save_figures=False, force=True)
        threshold = calibrate_threshold(turd, settings.lassie_mad_multiple)

    logger.info("detection pass at threshold %.2f", threshold)
    config = _search(
        settings,
        turd,
        float(threshold),
        save_figures=settings.lassie_save_figures,
        force=True,
    )

    config_path = run_dir / f"{settings.run_name}.config.yaml"
    lassie_config.write_config(config, str(config_path))
    logger.info("wrote Lassie config to %s", config_path)

    from lassie3.prepare import input_manifest

    (turd / "inputs.json").write_text(json.dumps(input_manifest(settings), indent=2))
    return turd


def _search(
    settings: RunSettings,
    turd: Path,
    threshold: float,
    save_figures: bool,
    force: bool,
) -> lassie_config.Config:
    """Run one Lassie pass, returning the config it was driven with."""
    config = build_config(settings, threshold=threshold, save_figures=save_figures)
    config.set_basepath(str(settings.run_dir(APPROACH)))
    config.set_config_name(settings.run_name)

    if turd.exists() and force:
        shutil.rmtree(turd)

    # `nparallel` goes straight to pyrocko's OpenMP stacker; it only helps when
    # pyrocko was built with OpenMP, which `pyproject.toml` arranges.
    lassie_core.search(config, force=force, nparallel=os.cpu_count() or 1)
    return config
