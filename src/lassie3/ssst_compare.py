"""Four-way comparison: each detector with and without station terms.

The baseline runs and the `ssst`-tagged runs are read side by side and scored
against the same reference catalog, so the effect of the terms shows up as a
difference between otherwise identical runs of the same detector.
"""

from __future__ import annotations

import json
import logging
from dataclasses import replace
from pathlib import Path

import numpy as np
from pyrocko import orthodrome as od

from lassie3 import ssst
from lassie3.compare import Event, lassie_grid_spacing, read_lassie, read_qseek, thresholds_used
from lassie3.reference import accuracy, load as load_reference, nearest, recall
from lassie3.review import recall_vs_tolerance
from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

SSST_TAG = "ssst"
ORDER = ("lassie", "lassie+ssst", "qseek", "qseek+ssst")
MATCH_TOLERANCE = 5.0
TOLERANCES = (0.5, 1.0, 2.0, 3.0, 5.0)


def variant_settings(settings: RunSettings) -> dict[str, RunSettings]:
    """Settings of the baseline pair and of the station-term pair."""
    base = replace(settings, tag="") if settings.tag == SSST_TAG else settings
    tag = f"{base.tag}-{SSST_TAG}" if base.tag else SSST_TAG
    return {"": base, SSST_TAG: replace(base, tag=tag)}


def catalogs(settings: RunSettings) -> dict[str, list[Event]]:
    """The four catalogs; a missing run is an empty list, never a fallback."""
    runs = variant_settings(settings)
    result = {}
    for variant, run in runs.items():
        suffix = f"+{variant}" if variant else ""
        result[f"lassie{suffix}"] = read_lassie(run)
        qseek_dir = run.run_dir("qseek") / run.run_name
        result[f"qseek{suffix}"] = read_qseek(run) if qseek_dir.exists() else []
    return result


def _metrics(events: list[Event], reference: list[dict]) -> dict:
    supported = [e for e in events if e.is_well_constrained()]
    metrics = {"n": len(events), "n_supported": len(supported)}
    if events and reference:
        metrics["recall"] = recall_vs_tolerance(events, reference)
        metrics["recall_5km_gate"] = recall(events, reference, MATCH_TOLERANCE, max_distance_km=5.0)["recovered"]
        metrics.update(accuracy(events, reference, MATCH_TOLERANCE))
    return metrics


def compare(settings: RunSettings) -> dict:
    """Score the four catalogs, draw the figure and write the JSON report."""
    runs = variant_settings(settings)
    base = runs[""]
    reference = load_reference(base)
    found = catalogs(settings)

    report = {
        "day": base.day_str,
        "reference_events": len(reference),
        "match_tolerance_s": MATCH_TOLERANCE,
        "thresholds": {
            "baseline": thresholds_used(base),
            "ssst": thresholds_used(runs[SSST_TAG]),
        },
        "catalogs": {name: _metrics(found[name], reference) for name in ORDER},
    }

    terms_file = ssst.terms_path(base)
    if terms_file.exists():
        terms = ssst.StationTerms.load(terms_file)
        report["station_terms"] = {
            "path": str(terms_file),
            "n_reference_events": terms.meta.get("n_reference_events"),
            "n_events_with_terms": terms.meta.get("n_events_with_terms"),
            "training_days": terms.meta.get("training_days"),
            "radius_m": terms.settings.radius,
            "prior_weight": terms.settings.prior_weight,
            "static": terms.static,
        }

    report["figure"] = str(_figure(base, found, reference, report))
    out = base.results_dir / f"compare-ssst-{base.run_name}.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    report["report"] = str(out)
    return report


def _figure(settings: RunSettings, found: dict, reference: list[dict], report: dict) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from lassie3.plot import STYLE, _plot_zoom, _size_scale

    spacing = lassie_grid_spacing(settings)
    sizes = {name: _size_scale(events) for name, events in found.items()}

    fig = plt.figure(figsize=(26, 17))
    grid = fig.add_gridspec(3, 4, height_ratios=(1.0, 1.0, 0.9))
    for column, name in enumerate(ORDER):
        for row, depth_section in ((0, False), (1, True)):
            axis = fig.add_subplot(grid[row, column])
            _plot_zoom(axis, {name: found[name]}, sizes, settings, reference, spacing, depth_section)
            # The shared close-up helper writes a long generic title; the panel
            # is identified by its catalog instead, and the suptitle explains
            # the lines once for all eight panels.
            what = "full 0–20 km depth" if depth_section else "map"
            axis.set_title(f"{STYLE[name]['label']} — source zone ±4 km, {what}, exact positions", fontsize=9.5)

    _plot_agreement(fig.add_subplot(grid[2, 0:2]), found, reference)
    _plot_recall(fig.add_subplot(grid[2, 2]), found, reference)
    _plot_table(fig.add_subplot(grid[2, 3]), report)

    fig.suptitle(
        f"ICDP Eger {settings.day_str} — each detector with and without source-specific station terms "
        f"(terms learned from WBNET events on the other days; this day held out)\n"
        "close-ups: stars are the numbered ISC reference events, lines join each to the detection credited "
        "with recovering it (±5 s); ×n marks n detections on one Lassie node",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = settings.results_dir / f"comparison-ssst-{settings.run_name}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    logger.info("wrote figure to %s", out)
    return out


def _plot_agreement(axis, found, reference) -> None:
    from lassie3.plot import STYLE

    if not reference:
        axis.text(0.5, 0.5, "no reference events this day", ha="center", va="center", transform=axis.transAxes)
        axis.set_axis_off()
        return

    metrics = [
        ("origin time [s]", "origin_time_median_s", lambda m, r: abs((m.time - r["time"]).total_seconds())),
        ("epicentre [km]", "epicentre_median_km",
         lambda m, r: od.distance_accurate50m(r["lat"], r["lon"], m.lat, m.lon) / 1e3),
        ("depth [km]", "depth_median_km", lambda m, r: abs(m.depth - r["depth"]) / 1e3),
    ]
    positions = np.arange(len(metrics))
    width = 0.2
    rng = np.random.default_rng(0)
    for k, name in enumerate(ORDER):
        offset = (k - 1.5) * width
        events = found[name]
        if not events:
            for index in range(len(metrics)):
                axis.text(index + offset, 0.02, "N/A", ha="center", va="bottom",
                          transform=axis.get_xaxis_transform(), fontsize=8, color=STYLE[name]["color"])
            axis.bar([], [], color=STYLE[name]["color"], label=f"{STYLE[name]['label']} — no run")
            continue
        scores = accuracy(events, reference, MATCH_TOLERANCE)
        got = recall(events, reference, MATCH_TOLERANCE)
        medians = [scores.get(key, 0.0) for _, key, _ in metrics]
        bars = axis.bar(positions + offset, medians, width, color=STYLE[name]["color"], alpha=0.55,
                        label=f"{STYLE[name]['label']} — {got['recovered']}/{got['reference_events']} "
                              f"within ±{MATCH_TOLERANCE:g} s (bar = median, dots = events)")
        for bar, value in zip(bars, medians):
            axis.text(bar.get_x() + bar.get_width() / 2, value, f"{value:g}", ha="center", va="bottom", fontsize=8)
        for index, (_, _, error_of) in enumerate(metrics):
            errors = [error_of(match, event) for event in reference
                      if (match := nearest(events, event, MATCH_TOLERANCE)) is not None]
            jitter = rng.uniform(-width / 4, width / 4, size=len(errors))
            axis.scatter(index + offset + jitter, errors, s=12, color=STYLE[name]["color"],
                         edgecolors="black", linewidths=0.3, zorder=5)

    depths = np.array([e["depth"] for e in reference]) / 1e3
    null = float(np.median(np.abs(depths - np.median(depths))))
    axis.hlines(null, 2 - 0.45, 2 + 0.45, colors="black", linestyles="--", linewidth=1,
                label=f"constant-depth guess ({np.median(depths):.1f} km) scores {null:.2f} km")
    axis.set_xticks(list(positions))
    axis.set_xticklabels([label for label, _, _ in metrics])
    axis.set_ylabel("|detection − reference|")
    axis.set_title("Agreement with the ISC/WBNET reference (matched within ±5 s)", fontsize=10)
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3, axis="y")


def _plot_recall(axis, found, reference) -> None:
    from lassie3.plot import STYLE

    for name in ORDER:
        events = found[name]
        if not events or not reference:
            continue
        counts = [recall(events, reference, tolerance)["recovered"] for tolerance in TOLERANCES]
        gated = recall(events, reference, MATCH_TOLERANCE, max_distance_km=5.0)["recovered"]
        axis.plot(TOLERANCES, counts, marker=STYLE[name]["marker"], color=STYLE[name]["color"],
                  label=f"{STYLE[name]['label']} ({gated} also within 5 km at ±5 s)")
    if reference:
        axis.axhline(len(reference), color="black", lw=0.6, ls=":", label=f"all {len(reference)} reference events")
    axis.set_xscale("log")
    axis.set_xticks(TOLERANCES)
    axis.set_xticklabels([f"±{t:g}" for t in TOLERANCES])
    axis.set_xlabel("origin-time matching tolerance [s]")
    axis.set_ylabel("reference events recovered")
    axis.set_title("Recall vs tolerance", fontsize=10)
    axis.legend(fontsize=7, loc="lower right")
    axis.grid(alpha=0.3)


def _plot_table(axis, report: dict) -> None:
    axis.axis("off")
    rows = [("", *ORDER)]
    cats = report["catalogs"]

    def cell(name, key, fmt="{}"):
        value = cats[name].get(key)
        return "—" if value is None else fmt.format(value)

    rows.append(("detections", *[cell(n, "n") for n in ORDER]))
    rows.append(("supported", *[cell(n, "n_supported") for n in ORDER]))
    for tolerance in TOLERANCES:
        rows.append((f"recovered ±{tolerance:g} s",
                     *["—" if not cats[n].get("recall") else str(cats[n]["recall"][f"{tolerance:g}s"]) for n in ORDER]))
    rows.append(("recovered ±5 s & 5 km", *[cell(n, "recall_5km_gate") for n in ORDER]))
    for label, key, fmt in (
        ("median |Δt| s", "origin_time_median_s", "{:.2f}"),
        ("median epicentre km", "epicentre_median_km", "{:.2f}"),
        ("median depth km", "depth_median_km", "{:.2f}"),
        ("mean |Δt| s", "origin_time_mean_s", "{:.2f}"),
        ("mean epicentre km", "epicentre_mean_km", "{:.2f}"),
        ("mean depth km", "depth_mean_km", "{:.2f}"),
        ("bias north km", "bias_north_km", "{:+.2f}"),
        ("bias east km", "bias_east_km", "{:+.2f}"),
        ("bias depth km", "bias_depth_km", "{:+.2f}"),
        ("bias time s", "bias_time_s", "{:+.2f}"),
    ):
        rows.append((label, *[cell(n, key, fmt) for n in ORDER]))

    table = axis.table(cellText=[list(r) for r in rows[1:]], colLabels=list(rows[0]), loc="upper center",
                       cellLoc="center", colWidths=[0.3, 0.175, 0.175, 0.175, 0.175])
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1.0, 1.25)
    terms = report.get("station_terms")
    note = (f"terms: {terms['n_events_with_terms']} of {terms['n_reference_events']} WBNET events on "
            f"{len(terms['training_days'])} other days, radius {terms['radius_m'] / 1e3:g} km, "
            f"prior weight {terms['prior_weight']:g}" if terms else "no station terms found")
    axis.set_title("Summary (errors on matched events; bias = detection − reference)\n" + note, fontsize=8.5)


def render(report: dict) -> str:
    cats = report["catalogs"]
    lines = [
        f"Four-way comparison for {report['day']} against {report['reference_events']} reference events "
        f"(matched within ±{report['match_tolerance_s']:g} s)",
        "",
        "| | " + " | ".join(ORDER) + " |",
        "|---|" + "---:|" * len(ORDER),
    ]

    def row(label, key, fmt="{}"):
        cells = []
        for name in ORDER:
            value = cats[name].get(key)
            cells.append("—" if value is None else fmt.format(value))
        lines.append(f"| {label} | " + " | ".join(cells) + " |")

    row("detections", "n")
    row("of which supported", "n_supported")
    for tolerance in TOLERANCES:
        cells = ["—" if not cats[n].get("recall") else str(cats[n]["recall"][f"{tolerance:g}s"]) for n in ORDER]
        lines.append(f"| recovered ±{tolerance:g} s | " + " | ".join(cells) + " |")
    row("recovered ±5 s and within 5 km", "recall_5km_gate")
    row("median origin-time error (s)", "origin_time_median_s", "{:.3f}")
    row("median epicentre error (km)", "epicentre_median_km", "{:.2f}")
    row("median depth error (km)", "depth_median_km", "{:.2f}")
    row("mean origin-time error (s)", "origin_time_mean_s", "{:.3f}")
    row("mean epicentre error (km)", "epicentre_mean_km", "{:.2f}")
    row("mean depth error (km)", "depth_mean_km", "{:.2f}")
    row("bias north (km)", "bias_north_km", "{:+.2f}")
    row("bias east (km)", "bias_east_km", "{:+.2f}")
    row("bias depth (km)", "bias_depth_km", "{:+.2f}")
    row("bias origin time (s)", "bias_time_s", "{:+.3f}")
    lines += ["", f"Figure: {report['figure']}", f"Report: {report.get('report', '')}"]
    return "\n".join(lines)
