"""Adversarial checks on the baseline results.

A recall of 13/13 for both detectors reads as conclusive, yet it is also what a
detector that fires constantly would score. Each check here quantifies one way
the headline numbers could mislead: recovery by chance, boundary artefacts,
scores that do not separate signal from noise, duplicated reference entries,
and the two runs quietly disagreeing on their shared inputs.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import logging

import numpy as np
from pyrocko import orthodrome as od

from lassie3.compare import lassie_grid_spacing, read_lassie, read_qseek, thresholds_used
from lassie3.reference import accuracy, load as load_reference, nearest, recall
from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

#: Seconds within which a detection counts as recovering a reference event.
#: Wide enough to absorb the reference catalog's own origin-time error and the
#: difference between its velocity model and ours.
TOLERANCE = 5.0

#: Reference entries closer than this may be one earthquake listed twice.
DUPLICATE_WINDOW = 60.0


def chance_matches(detections, reference, settings: RunSettings) -> dict:
    """Recoveries expected if the reference events were placed at random.

    Every detection claims a window of `2 * TOLERANCE` seconds; the union of
    those windows, as a fraction of the day, is the probability that a
    randomly placed reference event lands inside one. Overlapping windows are
    merged so a burst of detections is not counted several times over.

    The circular-shift null is the stronger test: the whole reference pattern
    is slid around the day by a random offset (keeping its own clustering)
    and recovery is recounted. If the true recall is not clearly above the
    largest value any shift produces, temporal association is not shown.
    """
    if not detections:
        return {"n_detections": 0}
    window = (settings.tmax - settings.tmin).total_seconds()
    times = np.sort([(d.time - settings.tmin).total_seconds() for d in detections])

    merged = []
    # Windows are clipped to the day: a detection seconds before midnight
    # cannot claim time the reference catalog does not cover.
    for start, end in zip(np.clip(times - TOLERANCE, 0, window), np.clip(times + TOLERANCE, 0, window)):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    coverage = min(1.0, sum(end - start for start, end in merged) / window)

    ref = np.array([(e["time"] - settings.tmin).total_seconds() for e in reference])
    rng = np.random.default_rng(0)
    shifted_hits = []
    # Shifts closer than a minute to zero would recover the real matches.
    for shift in rng.uniform(60.0, window - 60.0, size=10_000):
        moved = (ref + shift) % window
        index = np.searchsorted(times, moved)
        left = np.abs(moved - times[np.clip(index - 1, 0, len(times) - 1)])
        right = np.abs(moved - times[np.clip(index, 0, len(times) - 1)])
        shifted_hits.append(int((np.minimum(left, right) <= TOLERANCE).sum()))
    shifted_hits = np.array(shifted_hits)

    return {
        "n_detections": len(detections),
        "fraction_of_window_covered": round(coverage, 4),
        "expected_chance_recoveries": round(coverage * len(reference), 2),
        "actual_recoveries": recall(detections, reference, TOLERANCE)["recovered"],
        "circular_shift_null": {
            "shifts": len(shifted_hits),
            "max_recoveries": int(shifted_hits.max()),
            "p99_recoveries": int(np.percentile(shifted_hits, 99)),
            "mean_recoveries": round(float(shifted_hits.mean()), 2),
        },
    }


def per_event_table(catalogs: dict, reference: list[dict]) -> list[dict]:
    """One row per reference event with each detector's nearest detection.

    Medians hide a single wild miss and a lucky noise hit alike; the rows make
    both visible. `n_picks` is reported so a Qseek match can be checked for
    being a real pick-supported location rather than a low-semblance blip that
    happened to fall inside the tolerance.
    """
    rows = []
    for number, event in enumerate(reference, start=1):
        row = {
            "n": number,
            "time": event["time"].strftime("%H:%M:%S.%f")[:-3],
            "author": event["author"],
            "ML": event["magnitude"],
            "depth_km": round(event["depth"] / 1e3, 1),
        }
        for name, detections in catalogs.items():
            match = nearest(detections, event, TOLERANCE)
            if match is None:
                row[name] = None
                continue
            row[name] = {
                "dt_s": round((match.time - event["time"]).total_seconds(), 2),
                "epi_km": round(
                    od.distance_accurate50m(event["lat"], event["lon"], match.lat, match.lon)
                    / 1e3,
                    2,
                ),
                "ddepth_km": round((match.depth - event["depth"]) / 1e3, 2),
                "strength": round(match.strength, 3),
                "n_picks": match.n_picks,
            }
        rows.append(row)
    return rows


def boundary_artefacts(detections, settings: RunSettings, node: float) -> dict:
    """Detections within one node of the search-volume edge.

    When nothing coherent is inside the volume, a stacker's argmax drifts to
    wherever the travel-time moveout happens to align noise, which is
    disproportionately the boundary. Edge detections are therefore the
    signature of noise triggers. Qseek rejects lateral and bottom boundary
    hits itself (`ignore_boundary = "without_surface"`), so for it only the
    surface count should be non-zero; Lassie has no such rejection at all.
    """
    if not detections:
        return {"n": 0}

    grid = settings.grid
    lats = np.array([e.lat for e in detections])
    lons = np.array([e.lon for e in detections])
    depths = np.array([e.depth for e in detections])
    north, east = od.latlon_to_ne_numpy(grid.lat, grid.lon, lats, lons)

    lateral = (
        (east <= grid.east_bounds[0] + node)
        | (east >= grid.east_bounds[1] - node)
        | (north <= grid.north_bounds[0] + node)
        | (north >= grid.north_bounds[1] - node)
    )
    surface = depths <= grid.depth_bounds[0] + node
    bottom = depths >= grid.depth_bounds[1] - node
    return {
        "n": len(detections),
        "node_m": node,
        "lateral_edge": int(lateral.sum()),
        "surface": int(surface.sum()),
        "bottom": int(bottom.sum()),
        "clean_interior": int((~lateral & ~surface & ~bottom).sum()),
    }


def strength_separation(detections, reference: list[dict]) -> dict:
    """Where the recovered earthquakes rank in the detector's own score.

    If the true events are the strongest detections, the surplus is a tuning
    problem: raising the threshold would drop noise before it drops any
    earthquake. If they sit mid-pack, the score does not separate signal from
    noise and no threshold fixes it.
    """
    if not detections or not reference:
        return {}

    matched = {
        id(match)
        for event in reference
        if (match := nearest(detections, event, TOLERANCE)) is not None
    }
    ranked = sorted(detections, key=lambda e: e.strength, reverse=True)
    ranks = [index + 1 for index, event in enumerate(ranked) if id(event) in matched]
    top = len(reference)
    return {
        "n_detections": len(detections),
        "matched_unique_detections": len(ranks),
        "median_rank": int(np.median(ranks)) if ranks else None,
        "worst_rank": max(ranks) if ranks else None,
        f"true_events_in_top_{top}": sum(1 for rank in ranks if rank <= top),
    }


def temporal_association(detections, reference: list[dict], settings: RunSettings,
                         half_window: float = 1800.0) -> dict:
    """Fraction of detections falling within `half_window` s of any reference event.

    Swarm seismicity is clustered in time and noise is not, so if the detections
    beyond the thirteen catalogued events are real sub-catalog micro-earthquakes
    they should crowd the same hours. The comparison is against the fraction of
    the day those windows cover: a ratio near 1 means no association.
    """
    if not detections or not reference:
        return {}
    starts = sorted(max(settings.tmin, e["time"] - dt.timedelta(seconds=half_window))
                    for e in reference)
    ends = sorted(min(settings.tmax, e["time"] + dt.timedelta(seconds=half_window))
                  for e in reference)
    # Merge overlapping windows so the covered fraction is not double counted.
    merged = []
    for start, end in sorted(zip(starts, ends)):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    covered = sum((end - start).total_seconds() for start, end in merged)
    fraction_of_day = covered / (settings.tmax - settings.tmin).total_seconds()

    inside = sum(
        any(start <= d.time <= end for start, end in merged) for d in detections
    )
    fraction_inside = inside / len(detections)
    return {
        "half_window_s": half_window,
        "fraction_of_day_near_reference": round(fraction_of_day, 3),
        "fraction_of_detections_near_reference": round(fraction_inside, 3),
        "enrichment": round(fraction_inside / fraction_of_day, 2) if fraction_of_day else None,
    }


def constant_depth_null(detections, reference: list[dict]) -> dict:
    """Depth error a detector-free guess would score.

    The reference depths span less than a kilometre, so predicting their own
    median for every event already scores a ~0.1 km median error. A detector's
    depth error must be read against that, as agreement with a narrow catalog
    rather than demonstrated depth resolution.
    """
    if not reference:
        return {}
    depths = np.array([e["depth"] for e in reference]) / 1e3
    guess = float(np.median(depths))
    return {
        "constant_guess_km": round(guess, 2),
        "median_error_km": round(float(np.median(np.abs(depths - guess))), 3),
        "mean_error_km": round(float(np.mean(np.abs(depths - guess))), 3),
    }


def reference_composition(reference: list[dict]) -> dict:
    """Who located the reference events, and which carry magnitudes."""
    authors = {}
    for event in reference:
        authors.setdefault(event["author"], []).append(event["magnitude"])
    return {
        author: {
            "n": len(mags),
            "with_magnitude": sum(m is not None for m in mags),
            "ML_range": (
                [min(m for m in mags if m is not None), max(m for m in mags if m is not None)]
                if any(m is not None for m in mags) else None
            ),
        }
        for author, mags in authors.items()
    }


def recall_vs_tolerance(detections, reference: list[dict]) -> dict:
    """Recall at several matching tolerances.

    A recall that only holds at a wide tolerance is a different claim from one
    that holds at a tight one: Lassie's origin times run ~2.5 s late whenever
    it trades depth for time, so its recall at ±1 s is the honest measure of
    how often it times an event correctly, and ±5 s of whether it noticed it.
    """
    if not detections or not reference:
        return {}
    return {
        f"{tolerance:g}s": recall(detections, reference, tolerance)["recovered"]
        for tolerance in (0.5, 1.0, 2.0, 3.0, 5.0)
    }


def reference_duplicates(reference: list[dict]) -> list[dict]:
    """Reference pairs close enough in time to be one earthquake listed twice.

    ISC keeps the WBNET and PRU solutions of the same earthquake as separate
    events when their origin times disagree. Swarms also genuinely produce
    events seconds apart, so a pair here is a caveat, not a verdict; recall is
    reported both with and without the second member of each pair.
    """
    pairs = []
    for first, second in zip(reference, reference[1:]):
        gap = (second["time"] - first["time"]).total_seconds()
        if gap <= DUPLICATE_WINDOW:
            pairs.append({
                "first": first["time"].strftime("%H:%M:%S"),
                "second": second["time"].strftime("%H:%M:%S"),
                "gap_s": round(gap, 1),
                "authors": [first["author"], second["author"]],
            })
    return pairs


def setup_consistency(settings: RunSettings) -> dict:
    """Confirm from the persisted configs that both runs shared their inputs.

    Everything is read back from what each tool wrote into its run directory,
    not from `RunSettings`, because the point is to catch the two runs having
    drifted apart without anyone noticing.
    """
    from lassie3.compare import _lassie_rundir, _qseek_rundir

    report: dict = {}

    model_file = settings.velocity_model_file
    report["velocity_model_sha1"] = hashlib.sha1(model_file.read_bytes()).hexdigest()[:12]

    lassie_config = _lassie_rundir(settings) / "config.yaml"
    if lassie_config.exists():
        text = lassie_config.read_text()
        bands = [
            float(line.split(":")[1])
            for line in text.splitlines()
            if line.strip().startswith(("fmin:", "fmax:"))
        ]
        report["lassie"] = {
            "bandpass_hz": sorted(set(bands)),
            "stations_file": next(
                (line.split(":", 1)[1].strip().strip("'") for line in text.splitlines()
                 if line.startswith("stations_path:")),
                None,
            ),
            "earthmodel_ids": [
                line.split(":")[1].strip().strip("'")
                for line in text.splitlines()
                if line.strip().startswith("id:")
            ],
        }

    qseek_config = _qseek_rundir(settings) / "search.json"
    if qseek_config.exists():
        config = json.loads(qseek_config.read_text())
        bandpass = next(
            (step["bandpass"] for step in config["pre_processing"]
             if step["process"] == "bandpass"),
            None,
        )
        tracer = config["ray_tracers"][0]
        report["qseek"] = {
            "bandpass_hz": bandpass,
            "velocity_model_is_shared_file": tracer["earthmodel"]["filename"] == str(model_file),
            "n_station_xmls": len(config["stations"]["station_xmls"]),
            "channel_selector": config["data_provider"]["channel_selector"],
        }

    report["thresholds"] = thresholds_used(settings)

    # Manifests each runner wrote when it ran are the provenance; hashing the
    # staged files now only shows what is on disk today.
    manifests = {}
    for name, rundir in (("lassie", _lassie_rundir(settings)), ("qseek", _qseek_rundir(settings))):
        path = rundir / "inputs.json"
        manifests[name] = json.loads(path.read_text()) if path.exists() else None
    report["run_manifests"] = {
        name: (None if m is None else {"captured": m["captured"], "n_files": len(m["sha1"])})
        for name, m in manifests.items()
    }
    if all(manifests.values()):
        a, b = manifests["lassie"]["sha1"], manifests["qseek"]["sha1"]
        report["runs_consumed_identical_inputs"] = a == b
        report["input_differences"] = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
    else:
        report["runs_consumed_identical_inputs"] = None
    return report


def review(settings: RunSettings) -> dict:
    """Run every check and persist the report next to the figure."""
    lassie_events = read_lassie(settings)
    qseek_events = read_qseek(settings)
    qseek_strong = [e for e in qseek_events if e.is_well_constrained()]
    reference = load_reference(settings)

    catalogs = {
        "lassie": lassie_events,
        "qseek": qseek_events,
        "qseek_well_constrained": qseek_strong,
    }
    grid = settings.grid
    nodes = {
        "lassie": lassie_grid_spacing(settings),
        "qseek": grid.root_node_size / 2 ** (grid.n_levels - 1),
        "qseek_well_constrained": grid.root_node_size / 2 ** (grid.n_levels - 1),
    }

    duplicates = reference_duplicates(reference)
    deduplicated = [
        event for index, event in enumerate(reference)
        if index == 0
        or (event["time"] - reference[index - 1]["time"]).total_seconds() > DUPLICATE_WINDOW
    ]

    report = {
        "run": settings.run_name,
        "reference_events": len(reference),
        "reference_composition": reference_composition(reference),
        "constant_depth_null": constant_depth_null(None, reference),
        "recall_with_5km_gate": {
            name: recall(events, reference, TOLERANCE, max_distance_km=5.0)
            for name, events in catalogs.items()
        },
        # The two PRU solutions have no magnitude and a coarser origin; the
        # bias against the WBNET locations alone is the cleaner statement.
        "bias_vs_wbnet_only": {
            name: {k: v for k, v in accuracy(
                events, [e for e in reference if e["author"] == "WBNET"], TOLERANCE
            ).items() if k.startswith("bias_") or k == "matched"}
            for name, events in catalogs.items()
        },
        "setup": setup_consistency(settings),
        "chance": {
            name: chance_matches(events, reference, settings)
            for name, events in catalogs.items()
        },
        "recall_without_possible_duplicates": {
            name: recall(events, deduplicated, TOLERANCE)
            for name, events in catalogs.items()
        },
        "possible_reference_duplicates": duplicates,
        "boundary": {
            name: boundary_artefacts(events, settings, nodes[name])
            for name, events in catalogs.items()
        },
        "separation": {
            name: strength_separation(events, reference)
            for name, events in catalogs.items()
        },
        "temporal": {
            name: temporal_association(events, reference, settings)
            for name, events in catalogs.items()
        },
        "recall_vs_tolerance": {
            name: recall_vs_tolerance(events, reference)
            for name, events in catalogs.items()
        },
        "per_event": per_event_table(catalogs, reference),
    }

    out = settings.results_dir / f"review-{settings.run_name}.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    logger.info("wrote review to %s", out)
    return report


def render(report: dict) -> str:
    """Compact plain-text rendering of the review for the terminal."""
    lines = [f"Adversarial review — {report['run']} — {report['reference_events']} reference events", ""]

    comp = ", ".join(
        f"{n} {author}" + (f" (ML {c['ML_range'][0]}–{c['ML_range'][1]})" if c["ML_range"] else " (no ML)")
        for author, c in report["reference_composition"].items() for n in [c["n"]]
    )
    lines.append(f"Reference: {comp}")
    null = report["constant_depth_null"]
    lines.append(
        f"Depth null: guessing {null['constant_guess_km']} km for every event scores "
        f"{null['median_error_km']} km median / {null['mean_error_km']} km mean error"
    )
    lines.append("")
    lines.append("Recovery by chance (merged-window expectation, and the best of 10,000 circular shifts):")
    for name, chance in report["chance"].items():
        if not chance.get("n_detections"):
            continue
        null = chance["circular_shift_null"]
        lines.append(
            f"  {name:24s} {chance['n_detections']:5d} detections cover "
            f"{100 * chance['fraction_of_window_covered']:5.1f}% -> {chance['expected_chance_recoveries']:5.2f} expected; "
            f"shifted pattern: mean {null['mean_recoveries']}, p99 {null['p99_recoveries']}, max {null['max_recoveries']}; "
            f"actual {chance['actual_recoveries']:2d}"
        )
    lines.append("  with a 5 km epicentral gate: " + ", ".join(
        f"{name} {r['recovered']}/{r['reference_events']}"
        for name, r in report["recall_with_5km_gate"].items()
    ))

    if report["possible_reference_duplicates"]:
        lines.append("")
        lines.append("Possible duplicate reference entries:")
        for pair in report["possible_reference_duplicates"]:
            lines.append(f"  {pair['first']} / {pair['second']}  gap {pair['gap_s']} s  {pair['authors']}")
        lines.append("  recall without them: " + ", ".join(
            f"{name} {r['recovered']}/{r['reference_events']}"
            for name, r in report["recall_without_possible_duplicates"].items()
        ))

    lines.append("")
    lines.append("Boundary artefacts (within one node of the search-volume edge):")
    for name, boundary in report["boundary"].items():
        if boundary.get("n"):
            lines.append(
                f"  {name:24s} lateral {boundary['lateral_edge']:4d}  surface {boundary['surface']:4d}  "
                f"bottom {boundary['bottom']:4d}  clean interior {boundary['clean_interior']:4d} / {boundary['n']}"
            )

    lines.append("")
    lines.append("Do the true events rank as the strongest detections?")
    for name, separation in report["separation"].items():
        if separation:
            top = next(k for k in separation if k.startswith("true_events_in_top_"))
            lines.append(
                f"  {name:24s} median rank {separation['median_rank']:4d}, worst {separation['worst_rank']:4d} "
                f"of {separation['n_detections']};  {separation[top]} of the top {top.rsplit('_', 1)[1]} are true events"
            )

    lines.append("")
    lines.append(f"Recall vs matching tolerance (of {report['reference_events']}):")
    for name, table in report["recall_vs_tolerance"].items():
        if table:
            lines.append(f"  {name:24s} " + "  ".join(f"±{k}: {v:2d}" for k, v in table.items()))

    lines.append("")
    lines.append("Are the surplus detections clustered in time with the reference sequence?")
    for name, temporal in report["temporal"].items():
        if temporal:
            lines.append(
                f"  {name:24s} {100 * temporal['fraction_of_detections_near_reference']:5.1f}% of detections "
                f"lie within ±{temporal['half_window_s'] / 60:.0f} min of a reference event, "
                f"vs {100 * temporal['fraction_of_day_near_reference']:5.1f}% of the day -> "
                f"enrichment x{temporal['enrichment']}"
            )
    lines.append("")
    lines.append("Per reference event (dt s / epicentre km / depth error km / strength):")
    header = f"  {'#':>2s} {'time':12s} {'auth':6s} {'ML':>4s} {'z km':>5s} | {'lassie':30s} | {'qseek':30s} | qseek picks"
    lines.append(header)
    for row in report["per_event"]:
        def cell(match):
            if match is None:
                return f"{'-- missed --':30s}"
            return f"{match['dt_s']:+6.2f} {match['epi_km']:6.2f} {match['ddepth_km']:+6.2f} {match['strength']:8.3f}"
        picks = row["qseek"]["n_picks"] if row["qseek"] else "-"
        ml = f"{row['ML']:.2f}" if row["ML"] is not None else "  --"
        lines.append(
            f"  {row['n']:2d} {row['time']:12s} {row['author']:6s} {ml:>4s} {row['depth_km']:5.1f} | "
            f"{cell(row['lassie'])} | {cell(row['qseek'])} | {picks}"
        )

    lines.append("")
    lines.append("Mean signed offset vs WBNET events only (detection − reference):")
    for name, bias in report["bias_vs_wbnet_only"].items():
        if bias.get("matched"):
            lines.append(
                f"  {name:24s} N {bias['bias_north_km']:+.2f} km  E {bias['bias_east_km']:+.2f} km  "
                f"Z {bias['bias_depth_km']:+.2f} km  t {bias['bias_time_s']:+.2f} s  ({bias['matched']} matched)"
            )
    lines.append("")
    setup = dict(report["setup"])
    identical = setup.pop("runs_consumed_identical_inputs", None)
    lines.append("Run-time input manifests: " + json.dumps(setup.pop("run_manifests", None))
                 + f" -> identical inputs: {identical}"
                 + (f", differences {setup.pop('input_differences')}" if setup.get("input_differences") else ""))
    setup.pop("input_differences", None)
    lines.append("Shared-setup check: " + json.dumps(setup, default=str))
    return "\n".join(lines)
