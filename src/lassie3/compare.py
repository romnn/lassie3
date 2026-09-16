"""Read both catalogs into a common shape and plot the comparison."""

from __future__ import annotations

import csv
import datetime as dt
import json
from dataclasses import replace
import logging
from dataclasses import dataclass
from pathlib import Path

from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

#: Minimum picks for a Qseek location to count as well constrained. Six is two
#: phases on three stations, the point at which the depth
#: distribution stops piling up at the surface.
MIN_PICKS = 6

#: What "supported" means per detector, for legends and reports.
SUPPORT_CRITERIA = {
    "lassie": "interior of the search volume",
    "qseek": f"≥{MIN_PICKS} picks, azimuthal coverage > 0",
}


@dataclass(frozen=True)
class Event:
    """One detection, normalised across the two detectors.

    `strength` is the detector's own confidence measure and is *not* comparable
    between approaches: Lassie reports a stacked image-function amplitude on an
    arbitrary scale, Qseek a semblance that sums its P and S stacks and so can
    exceed 1 (1.58 on this day).
    """

    time: dt.datetime
    lat: float
    lon: float
    depth: float
    strength: float

    #: Qseek reports how many phase picks and how much azimuthal coverage
    #: support each location. Lassie's `detections.list` carries no equivalent,
    #: so these stay `None` for its events and quality filtering applies to the
    #: Qseek catalog only.
    n_picks: int | None = None
    azimuthal_coverage: float | None = None

    #: Lassie has no pick count, so its equivalent of a poorly supported
    #: location is one the stack pushed to the edge of the search volume, which
    #: is where noise ends up when nothing coherent is inside it. Set at read
    #: time from the grid; Qseek rejects such hits itself.
    boundary: bool = False

    def is_well_constrained(self) -> bool:
        """Whether the detector's own evidence supports this location.

        For Qseek: at least `MIN_PICKS` picks and non-zero azimuthal coverage.
        A 7-station array north of the target zone produces many detections
        that clear the semblance threshold on one or two picks alone; those
        pile up at zero depth. For Lassie: not within one grid node of the
        search boundary, where the stack drifts when nothing coherent is
        inside the volume.

        Neither is an uncertainty estimate. Matched Qseek events still have
        azimuthal gaps of ~240 degrees on this array; passing the flag means
        "supported by the detector's own criteria", nothing more. The exact
        criteria are spelled out in `SUPPORT_CRITERIA` for figure labels.
        """
        if self.n_picks is None or self.azimuthal_coverage is None:
            return not self.boundary
        return self.n_picks >= MIN_PICKS and self.azimuthal_coverage > 0.0


def read_lassie(settings: RunSettings) -> list[Event]:
    """Parse Lassie's `detections.list`.

    Columns are whitespace-separated: id, date, time, image-function value,
    reference latitude, reference longitude, north offset, east offset, depth.
    The date and time are two fields, so a record has ten columns, not nine.

    The two coordinate columns are the *grid reference*, identical for every
    detection; the location is the north/east offset from it in metres. Lassie's
    own `geo.point_coords` performs that conversion, so use it rather than
    reporting the grid origin for all 973 events.
    """
    from lassie import geo

    path = _lassie_rundir(settings) / "detections.list"
    if not path.exists():
        logger.warning("no Lassie detections at %s", path)
        return []

    events = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 9:
            continue
        time = dt.datetime.strptime(
            f"{fields[1]} {fields[2]}", "%Y-%m-%d %H:%M:%S.%f"
        ).replace(tzinfo=dt.timezone.utc)
        point = geo.Point(
            lat=float(fields[4]),
            lon=float(fields[5]),
            x=float(fields[6]),
            y=float(fields[7]),
            z=float(fields[8]),
        )
        lat, lon = geo.point_coords(point, system="latlon")
        events.append(
            Event(
                time=time,
                strength=float(fields[3]),
                lat=float(lat),
                lon=float(lon),
                depth=point.z,
                boundary=_at_boundary(point.x, point.y, point.z, settings),
            )
        )
    return sorted(events, key=lambda e: e.time)


def _at_boundary(north: float, east: float, depth: float, settings: RunSettings) -> bool:
    """Within one Lassie grid node of any face of the search volume."""
    grid = settings.grid
    node = lassie_grid_spacing(settings)
    return (
        north <= grid.north_bounds[0] + node or north >= grid.north_bounds[1] - node
        or east <= grid.east_bounds[0] + node or east >= grid.east_bounds[1] - node
        or depth <= grid.depth_bounds[0] + node or depth >= grid.depth_bounds[1] - node
    )


def read_qseek(settings: RunSettings) -> list[Event]:
    """Parse Qseek's `csv/detections.csv`."""
    path = _qseek_rundir(settings) / "csv" / "detections.csv"
    if not path.exists():
        logger.warning("no Qseek detections at %s", path)
        return []

    events = []
    with path.open() as handle:
        for row in csv.DictReader(handle):
            events.append(
                Event(
                    time=dt.datetime.fromisoformat(row["time"]),
                    lat=float(row["lat"]),
                    lon=float(row["lon"]),
                    depth=float(row["depth"]),
                    strength=float(row["semblance"]),
                    n_picks=int(row["n_picks"]),
                    azimuthal_coverage=float(row["azimuthal_coverage"]),
                )
            )
    return sorted(events, key=lambda e: e.time)


def thresholds_used(settings: RunSettings) -> dict[str, str]:
    """Read back the thresholds each run actually used.

    Both are read from the config each tool persisted into its run directory,
    not from `RunSettings`, so the figure reports what produced the catalogs on
    disk even if the defaults have since changed.
    """
    used = {}

    lassie_config = _lassie_rundir(settings) / "config.yaml"
    if lassie_config.exists():
        for line in lassie_config.read_text().splitlines():
            if line.startswith("detector_threshold:"):
                used["lassie"] = f"image function > {float(line.split(':')[1]):.1f}"

    qseek_config = _qseek_rundir(settings) / "search.json"
    if qseek_config.exists():
        threshold = json.loads(qseek_config.read_text()).get("detection_threshold")
        # Qseek's "MAD" is 10x the median absolute deviation of each window's
        # semblance, required as both peak height and prominence.
        used["qseek"] = (
            "semblance peak > 10×MAD of the window, with equal prominence (Qseek default)"
            if threshold == "MAD"
            else f"semblance > {threshold}"
        )

    return used


def _lassie_rundir(settings: RunSettings) -> Path:
    return settings.run_dir("lassie") / f"{settings.run_name}.turd"


def _qseek_rundir(settings: RunSettings) -> Path:
    """Qseek run directory for these settings.

    A tagged variant normally changes only Lassie (a finer grid, another
    threshold), so when no Qseek run carries the tag the comparison falls back
    to the untagged baseline and says so, rather than reporting an empty
    catalog as zero error.
    """
    tagged = settings.run_dir("qseek") / settings.run_name
    if settings.tag and not tagged.exists():
        base = replace(settings, tag="")
        fallback = base.run_dir("qseek") / base.run_name
        if fallback.exists():
            logger.info("no Qseek run tagged %r; comparing against baseline %s", settings.tag, fallback.name)
            return fallback
    return tagged


def lassie_grid_spacing(settings: RunSettings) -> float:
    """Node spacing of the Lassie run on disk, in metres.

    Read from the config Lassie persisted into its run directory, so a tagged
    1 km run is drawn and reviewed with a 1 km lattice even when the caller
    did not repeat `--lassie-spacing`. Falls back to the settings when no run
    exists yet.
    """
    config = _lassie_rundir(settings) / "config.yaml"
    if config.exists():
        for line in config.read_text().splitlines():
            if line.strip().startswith("dx:"):
                return float(line.split(":")[1])
    return settings.grid.lassie_spacing


def match(
    reference: list[Event], other: list[Event], tolerance: float = 3.0
) -> list[tuple[Event, Event | None]]:
    """Pair events by origin time within `tolerance` seconds.

    Time is the only robust key here: the two detectors locate on different
    grids, so their epicentres differ by more than any sensible distance
    tolerance even for the same earthquake. Each `other` event is consumed at
    most once, nearest-in-time first.
    """
    unused = sorted(other, key=lambda e: e.time)
    pairs: list[tuple[Event, Event | None]] = []

    for event in reference:
        best, best_delta = None, tolerance
        for candidate in unused:
            delta = abs((candidate.time - event.time).total_seconds())
            if delta <= best_delta:
                best, best_delta = candidate, delta
        if best is not None:
            unused.remove(best)
        pairs.append((event, best))

    return pairs


def summarise(settings: RunSettings) -> dict:
    """Collect the headline numbers for both runs."""
    lassie_events = read_lassie(settings)
    qseek_events = read_qseek(settings)
    pairs = match(lassie_events, qseek_events)
    matched = sum(1 for _, other in pairs if other is not None)

    well_constrained = [e for e in qseek_events if e.is_well_constrained()]

    from lassie3.reference import accuracy, load as load_reference, recall

    reference = load_reference(settings)

    return {
        "day": settings.day_str,
        "lassie": _stats(lassie_events),
        "qseek": _stats(qseek_events),
        "qseek_well_constrained": _stats(well_constrained),
        "matched_within_3s": matched,
        "lassie_only": len(lassie_events) - matched,
        "qseek_only": len(qseek_events) - matched,
        "reference_catalog": {
            "source": "ISC (WBNET + PRU)",
            "lassie": {**recall(lassie_events, reference),
                       **accuracy(lassie_events, reference)},
            "qseek": {**recall(qseek_events, reference),
                      **accuracy(qseek_events, reference)},
        },
    }


def _stats(events: list[Event]) -> dict:
    if not events:
        return {"n": 0}
    depths = sorted(e.depth for e in events)
    return {
        "n": len(events),
        "depth_median_km": depths[len(depths) // 2] / 1e3,
        "depth_range_km": (depths[0] / 1e3, depths[-1] / 1e3),
        "strength_min": min(e.strength for e in events),
        "strength_max": max(e.strength for e in events),
    }
