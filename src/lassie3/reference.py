"""Fetch an independent catalog to validate the detectors against.

Neither detector's output means much on its own: both produce hundreds of
detections per day on this array, and whether those are earthquakes can only be
settled against a catalog nobody in this pipeline produced.

ISC is used because it redistributes the **WBNET** (West Bohemia Network)
catalog, the dense local network operated over this exact swarm zone. WBNET
reports events down to about ML 0.0, which matters enormously here: on
2024-03-20 the EMSC/BGR catalog contains zero events in the region while WBNET
contains thirteen, all between ML 0.13 and 0.69. A regional catalog with an
ML ~1.5 completeness limit would have wrongly implied that every detection on
this day was a false alarm.
"""

from __future__ import annotations

import datetime as dt
import logging
import urllib.parse
import urllib.request
from pathlib import Path

from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

ISC_EVENT_SERVICE = "http://www.isc.ac.uk/fdsnws/event/1/query"

#: Region covering the SX array and the Vogtland/NW-Bohemia source zones.
SEARCH_REGION = {
    "minlatitude": 50.0,
    "maxlatitude": 50.6,
    "minlongitude": 12.0,
    "maxlongitude": 12.8,
}


def catalog_path(settings: RunSettings) -> Path:
    return settings.results_dir / "_stage" / f"reference-{settings.day_str}.txt"


def fetch(settings: RunSettings, force: bool = False) -> Path:
    """Download the reference catalog for the run's day, caching it locally."""
    out = catalog_path(settings)
    if out.exists() and not force:
        return out

    # ISC rejects percent-encoded colons in timestamps, so pass whole days.
    query = urllib.parse.urlencode({
        "starttime": settings.tmin.strftime("%Y-%m-%d"),
        "endtime": (settings.tmax + dt.timedelta(days=1)).strftime("%Y-%m-%d"),
        "format": "text",
        **SEARCH_REGION,
    })
    with urllib.request.urlopen(f"{ISC_EVENT_SERVICE}?{query}", timeout=90) as response:
        # An FDSN event service answers 204 with no body when nothing matches,
        # which is a valid empty catalog rather than a failure.
        body = response.read().decode() if response.status == 200 else ""

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(body)
    logger.info("reference catalog: %d events -> %s", len(load(settings)), out)
    return out


def load(settings: RunSettings) -> list[dict]:
    """Parse the cached catalog into time/lat/lon/depth/magnitude records."""
    path = catalog_path(settings)
    if not path.exists():
        return []

    events = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("|")
        time = dt.datetime.fromisoformat(fields[1].replace("Z", "+00:00"))
        if time.tzinfo is None:
            time = time.replace(tzinfo=dt.timezone.utc)
        # Whole days are requested, so trim to the run window.
        if not settings.tmin <= time < settings.tmax:
            continue
        events.append({
            "time": time,
            "author": fields[5],
            "lat": float(fields[2]),
            "lon": float(fields[3]),
            # FDSN text reports depth in km; the detectors work in metres.
            "depth": float(fields[4]) * 1e3,
            "magnitude": float(fields[10]) if fields[10] else None,
        })
    return sorted(events, key=lambda e: e["time"])


def nearest(detections, event: dict, tolerance: float = 5.0, max_distance_km: float | None = None):
    """The detection closest in origin time to `event`, or None beyond `tolerance`.

    `max_distance_km` adds an epicentral gate: a detection 15 km away that
    happens to fall inside the time tolerance is an association, not a
    location, and the gate lets recall be reported both ways.
    """
    from pyrocko import orthodrome as od

    if max_distance_km is not None:
        detections = [
            d for d in detections
            if od.distance_accurate50m(event["lat"], event["lon"], d.lat, d.lon) / 1e3
            <= max_distance_km
        ]
    candidate = min(
        detections,
        key=lambda d: abs((d.time - event["time"]).total_seconds()),
        default=None,
    )
    if candidate is None:
        return None
    if abs((candidate.time - event["time"]).total_seconds()) > tolerance:
        return None
    return candidate


def recall(detections, reference: list[dict], tolerance: float = 5.0,
           max_distance_km: float | None = None) -> dict:
    """How many reference events a detector recovered, matching on origin time.

    Time alone is the key: a detection displaced by kilometres is still a
    recovery of the same earthquake, whereas one displaced by seconds is not.
    The tolerance is wide enough to absorb both the reference catalog's own
    origin-time error and the difference between the two velocity models.
    """
    if not reference:
        return {"reference_events": 0, "recovered": 0, "recall": None}

    recovered = sum(
        nearest(detections, event, tolerance, max_distance_km) is not None
        for event in reference
    )
    return {
        "reference_events": len(reference),
        "recovered": recovered,
        "recall": round(recovered / len(reference), 3),
    }


def accuracy(detections, reference: list[dict], tolerance: float = 5.0) -> dict:
    """Median origin-time, epicentre and depth error against the reference.

    This is the comparison that means something. Detection counts depend on
    where each threshold was placed, but the error against thirteen
    independently located earthquakes is a property of the detector.
    """
    import numpy as np
    from pyrocko import orthodrome as od

    errors, offsets = [], []
    for event in reference:
        match = nearest(detections, event, tolerance)
        if match is None:
            continue
        north, east = od.latlon_to_ne(event["lat"], event["lon"], match.lat, match.lon)
        errors.append((
            abs((match.time - event["time"]).total_seconds()),
            od.distance_accurate50m(
                event["lat"], event["lon"], match.lat, match.lon
            ) / 1e3,
            abs(match.depth - event["depth"]) / 1e3,
        ))
        offsets.append((
            (match.time - event["time"]).total_seconds(),
            north / 1e3,
            east / 1e3,
            (match.depth - event["depth"]) / 1e3,
        ))

    if not errors:
        return {"matched": 0}

    matrix = np.array(errors)
    bias = np.array(offsets).mean(axis=0)
    return {
        "matched": len(errors),
        "origin_time_median_s": round(float(np.median(matrix[:, 0])), 3),
        "epicentre_median_km": round(float(np.median(matrix[:, 1])), 2),
        "epicentre_max_km": round(float(matrix[:, 1].max()), 2),
        "depth_median_km": round(float(np.median(matrix[:, 2])), 2),
        # Means expose what medians hide: a bimodal depth error averages to
        # several kilometres even when its median is small.
        "origin_time_mean_s": round(float(matrix[:, 0].mean()), 3),
        "epicentre_mean_km": round(float(matrix[:, 1].mean()), 2),
        "depth_mean_km": round(float(matrix[:, 2].mean()), 2),
        # Signed mean offsets (detection minus reference): a systematic
        # displacement shows up here and not in the absolute errors.
        "bias_time_s": round(float(bias[0]), 3),
        "bias_north_km": round(float(bias[1]), 2),
        "bias_east_km": round(float(bias[2]), 2),
        "bias_depth_km": round(float(bias[3]), 2),
    }
