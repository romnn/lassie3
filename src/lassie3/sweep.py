"""Threshold sweep from real Lassie runs, tabulated against the reference.

A count of peaks in the saved image function over-states what a run at that
level would detect, because Lassie's windows overlap and its detection pass
blinds around each peak. The only trustworthy sweep is one where each level
is an actual run, so this runs them, tagged so they sit next to the baseline
instead of overwriting it.
"""

from __future__ import annotations

import logging
from dataclasses import replace

from lassie3.compare import read_lassie, read_qseek
from lassie3.reference import load as load_reference, recall
from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

#: Qseek levels shown as post-filters of the saved run. These are not fresh
#: runs: a different threshold would also change Qseek's peak selection and
#: octree refinement, so they only bound what a re-run would give.
QSEEK_LEVELS = (0.0, 0.2, 0.3, 0.5, 0.7, 1.0, 1.1)


def sweep(settings: RunSettings, mad_multiples: list[float], run: bool = True) -> list[dict]:
    """Run (or tabulate) Lassie at each `median + k * MAD` level of the baseline."""
    from lassie3.run_lassie import calibrate_threshold, run as run_lassie

    baseline_turd = settings.run_dir("lassie") / f"{settings.run_name}.turd"
    reference = load_reference(settings)
    rows = []

    baseline = read_lassie(settings)
    rows.append(_row("baseline", settings.lassie_mad_multiple, _run_threshold(settings), baseline, reference))

    for k in mad_multiples:
        threshold = calibrate_threshold(baseline_turd, k)
        variant = replace(settings, tag=f"mad{k:g}", lassie_detector_threshold=threshold)
        if run:
            logger.info("sweep: Lassie at %g MAD = %.2f", k, threshold)
            run_lassie(variant, force=True)
        events = read_lassie(variant)
        if not events and not run:
            logger.warning("no run found for %s; run without --no-run", variant.run_name)
            continue
        rows.append(_row(variant.tag, k, threshold, events, reference))

    qseek = read_qseek(settings)
    for level in QSEEK_LEVELS:
        kept = [e for e in qseek if e.strength >= level]
        rows.append({
            "detector": "qseek (post-filter)",
            "level": f"semblance ≥ {level:g}",
            "detections": len(kept),
            "supported": sum(e.is_well_constrained() for e in kept),
            "recall_5s": recall(kept, reference, 5.0)["recovered"],
            "recall_1s": recall(kept, reference, 1.0)["recovered"],
        })
    return rows


def _run_threshold(settings: RunSettings) -> float | None:
    from lassie3.compare import thresholds_used

    text = thresholds_used(settings).get("lassie", "")
    try:
        return float(text.rsplit(" ", 1)[1])
    except (IndexError, ValueError):
        return None


def _row(tag: str, k: float, threshold: float | None, events, reference) -> dict:
    return {
        "detector": f"lassie ({tag})",
        "level": f"{k:g} MAD" + (f" = {threshold:.1f}" if threshold is not None else ""),
        "detections": len(events),
        "supported": sum(e.is_well_constrained() for e in events),
        "recall_5s": recall(events, reference, 5.0)["recovered"],
        "recall_1s": recall(events, reference, 1.0)["recovered"],
    }


def render(rows: list[dict], n_reference: int) -> str:
    """Markdown table, ready for the README."""
    lines = [
        f"| detector | level | detections | of which supported | recovered ±5 s | recovered ±1 s |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['detector']} | {row['level']} | {row['detections']} | {row['supported']} | "
            f"{row['recall_5s']}/{n_reference} | {row['recall_1s']}/{n_reference} |"
        )
    return "\n".join(lines)
