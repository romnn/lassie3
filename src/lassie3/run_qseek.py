"""Drive Qseek by writing its JSON config and calling its async search."""

from __future__ import annotations

import asyncio
import json
import logging
import os
import shutil
from pathlib import Path

from lassie3.prepare import station_xml_files
from lassie3.settings import RunSettings

# Qseek fixes the set of accepted `station_corrections` classes when
# `qseek.search` is first imported, so the plugin has to be registered before
# any function here touches Qseek; a lazy import inside `build_search` came
# too late because `run` imports `qseek.search` first.
from lassie3.qseek_ssst import SSSTCorrections  # noqa: E402

logger = logging.getLogger(__name__)

APPROACH = "qseek"


def build_search(settings: RunSettings):
    """Assemble the Qseek `Search` model.

    Qseek replaces Lassie's hand-tuned characteristic functions with PhaseNet
    P/S probabilities, and its regular grid with a refining octree, but stacks
    the same way: shift by Cake travel times through the shared 1D model, then
    detect on the stacked semblance.
    """
    from qseek.images.seisbench import SeisBench
    from qseek.models.station import StationInventory
    from qseek.octree import Octree
    from qseek.pre_processing.frequency_filters import Bandpass
    from qseek.search import Search
    from qseek.models.layered_model import LayeredEarthModel1D
    from qseek.tracers.cake import CakeTracer
    from qseek.waveforms.squirrel import PyrockoSquirrel

    grid = settings.grid
    fmin, fmax = settings.bandpass

    return Search(
        project_dir=settings.run_dir(APPROACH),
        stations=StationInventory(
            station_xmls=station_xml_files(settings),
            blacklist=list(settings.stations_blacklist),
        ),
        data_provider=PyrockoSquirrel(
            environment=settings.squirrel_dir,
            waveform_dirs=[settings.stage_dir],
            start_time=settings.tmin,
            end_time=settings.tmax,
            # Band codes only, e.g. "HH". Restricting here is what keeps the
            # duplicate BH?/LH? streams out of Qseek's stack, matching the
            # channel filter applied to Lassie's stations file.
            channel_selector=sorted({channel[:2] for channel in settings.channels}),
        ),
        pre_processing=[Bandpass(bandpass=(fmin, fmax))],
        octree=Octree(
            location=_reference_location(grid),
            root_node_size=grid.root_node_size,
            n_levels=grid.n_levels,
            east_bounds=grid.east_bounds,
            north_bounds=grid.north_bounds,
            depth_bounds=grid.depth_bounds,
        ),
        image_functions=[
            SeisBench(
                model="PhaseNet",
                pretrained="original",
                # No CUDA on this machine. Apple GPU support is applied
                # separately in `run`, since Qseek knows only CUDA and CPU.
                torch_use_cuda=False,
                torch_cpu_threads=settings.qseek_torch_threads or os.cpu_count() or 4,
            )
        ],
        ray_tracers=[
            CakeTracer(
                earthmodel=LayeredEarthModel1D(filename=settings.velocity_model_file),
            )
        ],
        detection_threshold=settings.qseek_detection_threshold,
        station_corrections=(
            SSSTCorrections(terms_path=settings.station_terms.resolve())
            if settings.station_terms is not None
            else None
        ),
    )


def _reference_location(grid):
    from qseek.models.location import Location

    return Location(lat=grid.lat, lon=grid.lon)


def run(settings: RunSettings, force: bool = False) -> Path:
    """Execute the Qseek search and return its run directory."""
    import nest_asyncio

    from qseek.search import Search

    run_dir = settings.run_dir(APPROACH)
    run_dir.mkdir(parents=True, exist_ok=True)

    # Clear the rundir here rather than through `force_rundir`: Qseek's
    # no-backup branch calls `Path.move`, which pathlib does not provide.
    rundir = run_dir / settings.run_name
    if rundir.exists() and force:
        shutil.rmtree(rundir)

    search = build_search(settings)
    config_path = run_dir / f"{settings.run_name}.json"
    config_path.write_text(search.model_dump_json(by_alias=False, indent=2))
    logger.info("wrote Qseek config to %s", config_path)

    # Qseek derives its rundir name from the config file stem, so reload through
    # the supported entry point rather than running the in-memory model.
    search = Search.from_config(config_path)

    # Qseek's own CLI applies this; its search nests event loops internally.
    nest_asyncio.apply()
    asyncio.run(search.start())

    from lassie3.prepare import input_manifest

    (rundir / "inputs.json").write_text(json.dumps(input_manifest(settings), indent=2))
    return rundir
