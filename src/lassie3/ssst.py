"""Source-specific station terms (SSST) derived from a reference catalog.

Neither Lassie v1 nor Qseek corrects its travel times for what the 1D velocity
model gets wrong along particular source-receiver paths. This module is the
first of the two passes that add such corrections:

1. **Measure** (here): for every reference event with a known hypocentre, pick
   P and S at every station with PhaseNet and take the residual against the
   Cake travel time through the shared 1D model. Residuals are demeaned per
   event, which is what re-solving the origin time would do, so only the part
   that can move a location survives.
2. **Apply** (`lassie_ssst`, `qseek_ssst`): at every search node, each detector
   adds the distance-weighted mean residual of the reference events around
   that node to the modelled travel time before stacking.

The weighting follows SCOTER (Nooshiri et al. 2019, doi:10.5880/GFZ.2.1.2019.002):
a bicube taper over hypocentral separation inside a cutoff radius, with a
static per-station term as the fallback where no reference events are close.
SCOTER shrinks the radius over relocation iterations; the reference hypocentres
here are fixed, so a single pass is the whole computation.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import logging
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np
from pyrocko import cake
from pyrocko import model as pyrocko_model
from pyrocko import orthodrome as od

from lassie3 import reference
from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

#: Phase definitions shared with both detectors' tabulated phases
#: ("P,p" and "S,s"); the first arrival of either branch is used.
PHASE_DEFINITIONS = {"P": ("P", "p"), "S": ("S", "s")}

#: Scale factor turning a median absolute deviation into a robust standard
#: deviation for normally distributed residuals, as in SCOTER's `smad_normal`.
MAD_TO_SIGMA = 1.482602218505602


@dataclass(frozen=True)
class TermsSettings:
    """Knobs of the term derivation, persisted into `terms.json`."""

    #: Hypocentral separation beyond which a reference event stops informing a
    #: node's term, in metres. SCOTER's `cutoff_dist`.
    radius: float = 3000.0

    #: Weight of the static term relative to one reference event at zero
    #: distance. A node with fewer than about this many events nearby falls
    #: back toward the static term, and a node with none gets exactly it.
    prior_weight: float = 5.0

    #: Residuals farther than this many robust standard deviations from zero
    #: are rejected before any averaging. SCOTER's dynamic default is 6.
    outlier_level: float = 6.0

    #: Minimum PhaseNet probability for a pick, SeisBench's default threshold.
    min_probability: float = 0.3

    #: Half-widths of the search window around the predicted P and S arrival
    #: in which the PhaseNet maximum counts as the pick, in seconds.
    search_window: tuple[float, float] = (2.0, 3.0)

    #: Events with fewer accepted picks cannot be demeaned meaningfully.
    min_picks_per_event: int = 4

    #: Station-phase pairs with fewer residuals get no static term (0 s).
    min_residuals: int = 5

    #: Band applied before PhaseNet, zero-phase, matching Qseek's pre-processing.
    bandpass: tuple[float, float] = (2.0, 30.0)


@dataclass
class Pick:
    """One PhaseNet pick against its prediction; `residual` is observed minus modelled."""

    event: str
    station: str
    phase: str
    predicted: float
    picked: float
    probability: float
    residual: float
    #: Set after demeaning and outlier rejection; only `kept` picks enter the terms.
    demeaned: float | None = None
    kept: bool = False


@dataclass
class StationTerms:
    """Static and source-specific terms, evaluable at any point.

    Reference-event positions are stored in a local north/east/depth frame
    (metres) around `reference_lat`/`reference_lon`; callers pass geographic
    coordinates so that Lassie's and Qseek's differently anchored grids need
    no bookkeeping.
    """

    reference_lat: float
    reference_lon: float
    settings: TermsSettings
    #: "NET.STA" -> {"P": seconds, "S": seconds}
    static: dict[str, dict[str, float]]
    #: (n_events, 3): north, east, depth in metres.
    events: np.ndarray
    #: (station, phase) -> (event indices, demeaned residuals)
    residuals: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]
    meta: dict = field(default_factory=dict)

    def stations(self) -> list[str]:
        return sorted(self.static)

    def evaluate(
        self,
        lat: np.ndarray,
        lon: np.ndarray,
        depth: np.ndarray,
        station: str,
        phase: str,
    ) -> np.ndarray:
        """Correction in seconds to add to the modelled travel time at each point.

        Weighted mean of the residuals of reference events within `radius`,
        each weighted by SCOTER's bicube taper of hypocentral separation, with
        the static term entering as `prior_weight` pseudo-events. `depth` is
        metres below sea level, matching the detectors' node depths.
        """
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        lon = np.atleast_1d(np.asarray(lon, dtype=float))
        depth = np.atleast_1d(np.asarray(depth, dtype=float))
        static = self.static.get(station, {}).get(phase, 0.0)

        indices, values = self.residuals.get((station, phase), (np.array([], int), np.array([])))
        if indices.size == 0:
            return np.full(lat.shape, static)

        north, east = od.latlon_to_ne_numpy(self.reference_lat, self.reference_lon, lat, lon)
        events = self.events[indices]
        separation = np.sqrt(
            (north[:, None] - events[None, :, 0]) ** 2
            + (east[:, None] - events[None, :, 1]) ** 2
            + (depth[:, None] - events[None, :, 2]) ** 2
        )
        weights = bicube(separation, self.settings.radius)
        prior = self.settings.prior_weight
        return (weights @ values + prior * static) / (weights.sum(axis=1) + prior)

    def save(self, path: Path) -> None:
        payload = {
            "reference_lat": self.reference_lat,
            "reference_lon": self.reference_lon,
            "settings": asdict(self.settings),
            "meta": self.meta,
            "static": self.static,
            "events": self.events.tolist(),
            "residuals": [
                {"station": station, "phase": phase, "events": idx.tolist(), "values": val.tolist()}
                for (station, phase), (idx, val) in sorted(self.residuals.items())
            ],
        }
        path.write_text(json.dumps(payload, indent=1))

    @classmethod
    def load(cls, path: Path) -> StationTerms:
        payload = json.loads(Path(path).read_text())
        settings = payload["settings"]
        settings["search_window"] = tuple(settings["search_window"])
        settings["bandpass"] = tuple(settings["bandpass"])
        return cls(
            reference_lat=payload["reference_lat"],
            reference_lon=payload["reference_lon"],
            settings=TermsSettings(**settings),
            static=payload["static"],
            events=np.array(payload["events"], dtype=float).reshape(-1, 3),
            residuals={
                (row["station"], row["phase"]): (
                    np.array(row["events"], dtype=int),
                    np.array(row["values"], dtype=float),
                )
                for row in payload["residuals"]
            },
            meta=payload["meta"],
        )


def bicube(distance: np.ndarray, cutoff: float) -> np.ndarray:
    """SCOTER's distance weight: (1 - (d / cutoff)^3)^3 inside the cutoff, 0 beyond."""
    return np.maximum(1.0 - (distance / cutoff) ** 3, 0.0) ** 3


def terms_dir(settings: RunSettings) -> Path:
    """Where the terms for the run's day live; independent of the run tag."""
    return settings.results_dir / "ssst" / f"eger-{settings.day_str}"


def terms_path(settings: RunSettings) -> Path:
    return terms_dir(settings) / "terms.json"


def training_days(settings: RunSettings) -> list[dt.date]:
    """Every day with waveform data except the evaluation day.

    Holding the evaluation day out keeps the second pass honest: its reference
    events must not have informed the corrections it is scored against.
    """
    days = set()
    for path in settings.data_dir.glob("mseed/*/*/*.mseed"):
        try:
            days.add(dt.date.fromisoformat(path.stem.split(".")[-1]))
        except ValueError:
            continue
    return sorted(day for day in days if day != settings.day)


def training_events(settings: RunSettings, days: list[dt.date]) -> list[dict]:
    """WBNET hypocentres on the training days, fetched per day and cached.

    PRU entries are dropped: they duplicate WBNET events with a regional
    location and a fixed 10 km depth, which would put a bias straight into the
    residuals.
    """
    events = []
    for day in days:
        day_settings = replace(settings, day=day, hours=24.0)
        reference.fetch(day_settings)
        for index, event in enumerate(reference.load(day_settings)):
            if event["author"] != "WBNET":
                continue
            # Checked on the event's own time as well as on the day list, so
            # the evaluation window can never leak in through a cached
            # catalog or an odd `hours` setting.
            if settings.tmin <= event["time"] < settings.tmax:
                continue
            event = dict(event, id=f"{day.isoformat()}-{index:03d}")
            events.append(event)
    logger.info("%d WBNET reference events on %d training days", len(events), len(days))
    return events


def travel_time(earthmodel, phase: str, distance: float, source_depth: float, receiver_depth: float) -> float | None:
    """First arrival of the phase in seconds, as Lassie's Cake shifter computes it."""
    phases = [cake.PhaseDef(name) for name in PHASE_DEFINITIONS[phase]]
    rays = earthmodel.arrivals(
        phases=phases,
        distances=[distance * cake.m2d],
        zstart=source_depth,
        zstop=receiver_depth,
    )
    times = [ray.t for ray in rays]
    return min(times) if times else None


def _station_key(station) -> str:
    return f"{station.network}.{station.station}"


def measure_picks(
    settings: RunSettings,
    terms_settings: TermsSettings,
    events: list[dict],
    limit: int | None = None,
) -> list[Pick]:
    """Pick P and S of every reference event with PhaseNet and take residuals.

    Uses the same PhaseNet weights ("original") and the same zero-phase
    band-pass as the Qseek run, so the picks describe the data, not a third
    picker's habits. The pick is the probability maximum inside a window
    around the modelled arrival; a maximum below `min_probability` is no pick.
    """
    import obspy
    import seisbench.models as sbm
    from pyrocko import io as pyrocko_io
    from pyrocko.obspy_compat.base import to_obspy_trace

    model = sbm.PhaseNet.from_pretrained("original")
    earthmodel = cake.load_model(str(settings.velocity_model_file), format="nd")
    stations = pyrocko_model.load_stations(str(settings.stations_file))
    fmin, fmax = terms_settings.bandpass
    half_p, half_s = terms_settings.search_window

    picks: list[Pick] = []
    events = sorted(events, key=lambda e: e["time"])
    if limit:
        events = events[:limit]
    by_day: dict[dt.date, list[dict]] = {}
    for event in events:
        by_day.setdefault(event["time"].date(), []).append(event)

    for day, day_events in sorted(by_day.items()):
        # One load per station and day; the daily files hold nothing but the
        # selected day, and chopping windows out of memory is far cheaper than
        # re-reading a file per event.
        traces_by_station = {}
        for station in stations:
            key = _station_key(station)
            path = settings.data_dir / "mseed" / station.network / station.station / f"{key}.{day.isoformat()}.mseed"
            if not path.exists():
                logger.warning("%s has no data on %s", key, day)
                continue
            traces_by_station[key] = [
                trace for trace in pyrocko_io.load(str(path)) if trace.channel in settings.channels
            ]
        logger.info("%s: %d events, %d stations", day, len(day_events), len(traces_by_station))

        for event in day_events:
            origin = event["time"].timestamp()
            stream = obspy.Stream()
            predicted = {}
            for station in stations:
                key = _station_key(station)
                if key not in traces_by_station:
                    continue
                distance = od.distance_accurate50m(event["lat"], event["lon"], station.lat, station.lon)
                receiver_depth = station.depth - station.elevation
                times = {
                    phase: travel_time(earthmodel, phase, distance, event["depth"], receiver_depth)
                    for phase in PHASE_DEFINITIONS
                }
                if times["P"] is None or times["S"] is None:
                    continue
                predicted[key] = times

                # PhaseNet discards 2.5 s at each end of what it annotates, so
                # the window reaches well past both search ranges.
                tmin = origin + times["P"] - half_p - 15.0
                tmax = origin + times["S"] + half_s + 20.0
                for trace in traces_by_station[key]:
                    try:
                        piece = trace.chop(tmin, tmax, inplace=False)
                    except Exception:
                        continue
                    if piece.tmax - piece.tmin < 0.9 * (tmax - tmin):
                        continue
                    obspy_trace = to_obspy_trace(piece)
                    obspy_trace.detrend("demean")
                    obspy_trace.taper(0.05)
                    obspy_trace.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
                    stream += obspy_trace

            if not predicted or len(stream) == 0:
                continue
            annotated = model.annotate(stream)

            for key, times in predicted.items():
                for phase, half in (("P", half_p), ("S", half_s)):
                    probability = annotated.select(id=f"{key}..PhaseNet_{phase}")
                    if not probability:
                        continue
                    trace = probability[0]
                    centre = origin + times[phase]
                    window = trace.slice(
                        obspy.UTCDateTime(centre - half), obspy.UTCDateTime(centre + half)
                    )
                    if window.stats.npts == 0:
                        continue
                    index = int(np.argmax(window.data))
                    value = float(window.data[index])
                    if value < terms_settings.min_probability:
                        continue
                    picked = float(window.stats.starttime.timestamp) + index * window.stats.delta
                    picks.append(
                        Pick(
                            event=event["id"],
                            station=key,
                            phase=phase,
                            predicted=centre,
                            picked=picked,
                            probability=value,
                            residual=picked - centre,
                        )
                    )
    logger.info("%d picks on %d events", len(picks), len(events))
    return picks


def demean_and_reject(picks: list[Pick], terms_settings: TermsSettings) -> list[Pick]:
    """Remove each event's origin-time freedom, then drop outliers, and re-centre.

    The per-event mean over all its P and S residuals is what a location
    routine would absorb into the origin time; subtracting it leaves the
    station-differential part that station terms are meant to capture.
    Outliers are judged per phase against `outlier_level` robust standard
    deviations of the demeaned residuals, as SCOTER's dynamic rejection does.
    Rejection changes the means, so the sequence is demean, reject, demean,
    reject, demean: the kept picks of every event average to zero at the end.
    """
    by_event: dict[str, list[Pick]] = {}
    for pick in picks:
        pick.kept = True
        by_event.setdefault(pick.event, []).append(pick)

    for _ in range(2):
        _demean(by_event, terms_settings.min_picks_per_event)
        _reject_outliers(picks, terms_settings.outlier_level)
    _demean(by_event, terms_settings.min_picks_per_event)
    return picks


def _demean(by_event: dict[str, list[Pick]], min_picks: int) -> None:
    """Centre every event's residuals on the mean of its kept picks."""
    for event_picks in by_event.values():
        kept = [p for p in event_picks if p.kept]
        if len(kept) < min_picks:
            for p in event_picks:
                p.kept = False
                p.demeaned = None
            continue
        mean = float(np.mean([p.residual for p in kept]))
        for p in event_picks:
            p.demeaned = p.residual - mean


def _reject_outliers(picks: list[Pick], level: float) -> None:
    """Drop kept picks whose demeaned residual is beyond `level` robust sigmas of their phase."""
    for phase in PHASE_DEFINITIONS:
        values = np.array([p.demeaned for p in picks if p.kept and p.phase == phase and p.demeaned is not None])
        if values.size == 0:
            continue
        sigma = MAD_TO_SIGMA * float(np.median(np.abs(values - np.median(values))))
        cutoff = level * max(sigma, 1e-3)
        for p in picks:
            if p.phase == phase and p.kept and p.demeaned is not None and abs(p.demeaned) > cutoff:
                p.kept = False


def build_terms(
    settings: RunSettings,
    terms_settings: TermsSettings,
    events: list[dict],
    picks: list[Pick],
    meta: dict,
) -> StationTerms:
    """Assemble static terms and the per-event residual store from kept picks."""
    index_of = {event["id"]: i for i, event in enumerate(events)}
    north, east = od.latlon_to_ne_numpy(
        settings.grid.lat,
        settings.grid.lon,
        np.array([e["lat"] for e in events]),
        np.array([e["lon"] for e in events]),
    )
    positions = np.column_stack([north, east, [e["depth"] for e in events]])

    grouped: dict[tuple[str, str], list[Pick]] = {}
    for pick in picks:
        if pick.kept and pick.demeaned is not None:
            grouped.setdefault((pick.station, pick.phase), []).append(pick)

    static: dict[str, dict[str, float]] = {}
    residuals = {}
    counts = {}
    for (station, phase), group in sorted(grouped.items()):
        values = np.array([p.demeaned for p in group])
        indices = np.array([index_of[p.event] for p in group], dtype=int)
        residuals[(station, phase)] = (indices, values)
        term = float(values.mean()) if values.size >= terms_settings.min_residuals else 0.0
        static.setdefault(station, {})[phase] = term
        counts[f"{station}:{phase}"] = {
            "n": int(values.size),
            "mean": float(values.mean()),
            "std": float(values.std()),
            "median": float(np.median(values)),
        }

    return StationTerms(
        reference_lat=settings.grid.lat,
        reference_lon=settings.grid.lon,
        settings=terms_settings,
        static=static,
        events=positions,
        residuals=residuals,
        meta=dict(meta, counts=counts),
    )


def write_picks(picks: list[Pick], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["event", "station", "phase", "predicted", "picked", "probability", "residual", "demeaned", "kept"])
        for p in picks:
            writer.writerow([
                p.event, p.station, p.phase,
                dt.datetime.fromtimestamp(p.predicted, dt.timezone.utc).isoformat(timespec="milliseconds"),
                dt.datetime.fromtimestamp(p.picked, dt.timezone.utc).isoformat(timespec="milliseconds"),
                f"{p.probability:.3f}", f"{p.residual:.3f}",
                "" if p.demeaned is None else f"{p.demeaned:.3f}", int(p.kept),
            ])


def derive(
    settings: RunSettings,
    terms_settings: TermsSettings | None = None,
    days: list[dt.date] | None = None,
    limit: int | None = None,
    force: bool = False,
) -> Path:
    """Run the measurement pass and write `terms.json`, `picks.csv` and a figure."""
    terms_settings = terms_settings or TermsSettings()
    out_dir = terms_dir(settings)
    out = terms_path(settings)
    if out.exists() and not force:
        logger.info("station terms exist at %s (use --force to redo)", out)
        return out
    out_dir.mkdir(parents=True, exist_ok=True)

    days = days or training_days(settings)
    if settings.day in days:
        # The evaluation day never informs the terms it is scored against,
        # whichever way the days were chosen.
        logger.warning("dropping %s from the training days: it is the evaluation day", settings.day_str)
        days = [day for day in days if day != settings.day]
    events = training_events(settings, days)
    if not events:
        raise RuntimeError("no reference events on the training days")

    picks = measure_picks(settings, terms_settings, events, limit=limit)
    picks = demean_and_reject(picks, terms_settings)
    write_picks(picks, out_dir / "picks.csv")

    used_ids = {p.event for p in picks if p.kept}
    meta = {
        "evaluation_day": settings.day_str,
        "training_days": [d.isoformat() for d in days],
        "n_reference_events": len(events),
        "n_events_with_terms": len(used_ids),
        "n_picks": len(picks),
        "n_picks_kept": sum(p.kept for p in picks),
        "picker": "SeisBench PhaseNet 'original', zero-phase %g-%g Hz" % terms_settings.bandpass,
        "velocity_model": str(settings.velocity_model_file),
        "stations_file": str(settings.stations_file),
        "reference_catalog": "ISC event service, author WBNET",
    }
    terms = build_terms(settings, terms_settings, events, picks, meta)
    terms.save(out)
    logger.info("wrote station terms to %s", out)

    from lassie3.ssst_plot import plot_terms

    plot_terms(settings, terms, picks, events, out_dir / "terms.png")

    # Kept apart from terms.json, which the detector runs hash into their
    # input manifests.
    held_out, held_out_picks = held_out_check(settings, terms, terms_settings)
    (out_dir / "held-out.json").write_text(json.dumps(held_out, indent=2))
    write_picks(held_out_picks, out_dir / "held-out-picks.csv")
    logger.info("held-out check on %s: %s", settings.day_str, held_out)
    return out


def held_out_check(
    settings: RunSettings,
    terms: StationTerms,
    terms_settings: TermsSettings,
) -> tuple[dict, list[Pick]]:
    """Do the terms predict the evaluation day's own residuals?

    The evaluation day's WBNET events never entered the terms, so measuring
    their residuals the same way and subtracting the term at each hypocentre
    tests the terms directly, independent of either detector: a residual
    scatter that shrinks says the terms carry structure that generalises, one
    that grows says they are fitting noise or the wrong day.
    """
    events = [
        dict(event, id=f"{settings.day_str}-{index:03d}")
        for index, event in enumerate(reference.load(settings))
        if event["author"] == "WBNET"
    ]
    if not events:
        return {"n_events": 0}, []
    picks = demean_and_reject(measure_picks(settings, terms_settings, events), terms_settings)
    by_id = {event["id"]: event for event in events}

    result = {"n_events": len(events), "n_picks_kept": sum(p.kept for p in picks)}
    for phase in PHASE_DEFINITIONS:
        kept = [p for p in picks if p.kept and p.phase == phase and p.demeaned is not None]
        if not kept:
            continue
        residual = np.array([p.demeaned for p in kept])
        term = np.array([
            terms.evaluate(
                [by_id[p.event]["lat"]], [by_id[p.event]["lon"]], [by_id[p.event]["depth"]], p.station, phase
            )[0]
            for p in kept
        ])
        static = np.array([terms.static.get(p.station, {}).get(phase, 0.0) for p in kept])
        rms = lambda values: float(np.sqrt(np.mean(np.square(values))))  # noqa: E731
        result[phase] = {
            "n": len(kept),
            "rms_residual_s": round(rms(residual), 4),
            "rms_after_static_s": round(rms(residual - static), 4),
            "rms_after_terms_s": round(rms(residual - term), 4),
            "correlation_residual_term": round(float(np.corrcoef(residual, term)[0, 1]), 3) if len(kept) > 2 else None,
        }
    return result, picks


def render(terms: StationTerms) -> str:
    """Markdown table of the static terms, one row per station."""
    counts = terms.meta.get("counts", {})
    lines = ["| station | P term (s) | n | σ (s) | S term (s) | n | σ (s) |", "|---|---:|---:|---:|---:|---:|---:|"]
    for station in terms.stations():
        cells = [station]
        for phase in ("P", "S"):
            stats = counts.get(f"{station}:{phase}")
            if stats:
                cells += [f"{terms.static[station].get(phase, 0.0):+.3f}", str(stats["n"]), f"{stats['std']:.3f}"]
            else:
                cells += ["—", "0", "—"]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)
