"""Build the inputs both detectors consume: waveforms, stations, velocity model."""

from __future__ import annotations

import logging
from pathlib import Path

from pyrocko import io as pyrocko_io
from pyrocko import model as pyrocko_model
from pyrocko.io import stationxml

from lassie3.settings import WEST_BOHEMIA_ND, RunSettings

logger = logging.getLogger(__name__)


def stage_waveforms(settings: RunSettings) -> list[Path]:
    """Extract the selected day and channels into a flat per-day directory.

    This rewrites the data rather than symlinking it because both detectors
    derive global properties from every stream they can see, not just the ones
    they end up stacking. Lassie sets its image-function sample rate from the
    coarsest channel in the pile, so a 1 Hz LH? stream would drop detection
    timing to 1 s even though a trace selector excludes it from the stack.

    Staging also keeps each run off the full 2.2 GB, 14-day archive, which both
    tools would otherwise index in full.
    """
    stage = settings.stage_dir
    stage.mkdir(parents=True, exist_ok=True)
    settings.squirrel_dir.mkdir(parents=True, exist_ok=True)

    # Clear stale output so a re-run after changing day or channels cannot mix
    # them with the previous selection.
    for old in stage.glob("*.mseed"):
        old.unlink()

    sources = sorted(settings.data_dir.glob(f"mseed/*/*/*.{settings.day_str}.mseed"))
    if not sources:
        raise FileNotFoundError(
            f"no miniSEED for {settings.day_str} under {settings.data_dir / 'mseed'}"
        )

    staged = []
    for src in sources:
        traces = [
            trace
            for trace in pyrocko_io.load(str(src))
            if trace.channel in settings.channels
        ]
        if not traces:
            logger.warning("%s has none of %s, skipping", src.name, settings.channels)
            continue
        out = stage / src.name
        pyrocko_io.save(traces, str(out))
        staged.append(out)

    logger.info("staged %d waveform files for %s", len(staged), settings.day_str)
    return staged


def staged_station_codes(settings: RunSettings) -> set[tuple[str, str]]:
    """Return (network, station) pairs that actually have data for the day.

    Derived from the staged filenames rather than from the StationXML: metadata
    is present for all 21 stations, but only the nine open ones (SX, CZ.NKC,
    GQ.LNDWU) have waveforms, and feeding the detectors phantom stations skews
    the stack normalisation and the azimuthal-coverage statistics.
    """
    codes = set()
    for path in sorted(settings.stage_dir.glob("*.mseed")):
        net, sta, _ = path.name.split(".", 2)
        codes.add((net, sta))
    return codes


def write_stations(settings: RunSettings) -> Path:
    """Convert the StationXML inventory to a Pyrocko stations file.

    Lassie v1 reads Pyrocko stations; Qseek reads the StationXML directly. Both
    are filtered to the same set, so the two runs see an identical network.
    """
    available = staged_station_codes(settings)
    blacklist = set(settings.stations_blacklist)

    # `get_pyrocko_stations` yields one entry per channel-epoch group, so a
    # station recording BH?, HH? and LH? appears three times. Keep only the
    # entry that still has channels after filtering, keyed by NET.STA.
    by_code: dict[tuple[str, str], pyrocko_model.Station] = {}
    for xml_file in sorted(settings.data_dir.glob("stationxml/*.xml")):
        inventory = stationxml.load_xml(filename=str(xml_file))
        for station in inventory.get_pyrocko_stations(timespan=(
            settings.tmin.timestamp(),
            settings.tmax.timestamp(),
        )):
            key = (station.network, station.station)
            if key not in available or f"{key[0]}.{key[1]}" in blacklist:
                continue
            channels = [
                channel
                for channel in station.get_channels()
                if channel.name in settings.channels
            ]
            if not channels or key in by_code:
                continue
            station.set_channels(channels)
            by_code[key] = station

    stations = [by_code[key] for key in sorted(by_code)]

    if not stations:
        raise RuntimeError("no stations left after filtering to staged waveform data")

    out = settings.stations_file
    out.parent.mkdir(parents=True, exist_ok=True)
    pyrocko_model.dump_stations(stations, str(out))
    logger.info("wrote %d stations to %s", len(stations), out)
    return out


def station_xml_files(settings: RunSettings) -> list[Path]:
    """StationXML files for the stations that have data, for Qseek's inventory."""
    available = staged_station_codes(settings)
    blacklist = set(settings.stations_blacklist)
    return [
        path
        for path in sorted(settings.data_dir.glob("stationxml/*.xml"))
        if tuple(path.stem.split(".", 1)) in available and path.stem not in blacklist
    ]


def write_velocity_model(settings: RunSettings) -> Path:
    """Materialise the shared 1D model unless the caller supplied their own."""
    out = settings.velocity_model_file
    if settings.velocity_model is not None:
        if not out.exists():
            raise FileNotFoundError(f"velocity model not found: {out}")
        return out

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(WEST_BOHEMIA_ND)
    logger.info("wrote velocity model to %s", out)
    return out


def input_manifest(settings: RunSettings) -> dict:
    """Hashes and settings of everything a run reads, captured at run time.

    Written into each run directory by the runners, so that two runs can be
    shown to have consumed identical inputs from records made when they ran,
    not from hashing whatever is on disk later.
    """
    import datetime as dt
    import hashlib

    files = sorted(settings.stage_dir.glob("*.mseed")) + [settings.stations_file, settings.velocity_model_file]
    files += station_xml_files(settings)
    return {
        "captured": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "day": settings.day_str,
        "window_hours": settings.hours,
        "channels": list(settings.channels),
        "bandpass_hz": list(settings.bandpass),
        "grid": {
            "lat": settings.grid.lat, "lon": settings.grid.lon,
            "east_bounds": list(settings.grid.east_bounds),
            "north_bounds": list(settings.grid.north_bounds),
            "depth_bounds": list(settings.grid.depth_bounds),
        },
        "blacklist": list(settings.stations_blacklist),
        "sha1": {
            str(path.relative_to(settings.results_dir.parent)): hashlib.sha1(path.read_bytes()).hexdigest()
            for path in files
        },
    }


def prepare(settings: RunSettings) -> None:
    """Run every preparation step in dependency order."""
    stage_waveforms(settings)
    write_stations(settings)
    write_velocity_model(settings)
