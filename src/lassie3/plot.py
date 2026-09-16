"""Six-panel comparison figure for the two catalogs."""

from __future__ import annotations

import dataclasses
import logging
import math
from pathlib import Path

import numpy as np
from pyrocko import orthodrome as od

from lassie3.compare import (
    SUPPORT_CRITERIA,
    Event,
    lassie_grid_spacing,
    read_lassie,
    read_qseek,
    thresholds_used,
)
from lassie3.reference import accuracy, load as load_reference, nearest, recall
from lassie3.settings import RunSettings

logger = logging.getLogger(__name__)

STYLE = {
    "lassie": {"color": "#c1440e", "marker": "o", "label": "Lassie v1 (STA/LTA + wave packet)"},
    "qseek": {"color": "#1f6feb", "marker": "^", "label": "Qseek (PhaseNet + octree)"},
}

#: Half-width of the source-zone close-ups. Four kilometres comfortably holds
#: the WBNET cluster and Qseek's matches while still showing where Lassie's
#: grid puts the same events; matches outside it are listed on the panel.
ZOOM_HALF_WIDTH_M = 4_000.0

#: Marker area range, in points squared, spanned by a catalog's own strength
#: range, so that within a catalog a bigger marker is a stronger detection.
MARKER_MIN, MARKER_MAX = 8.0, 90.0

#: Tolerance used to join reference events to detections on the figure; the
#: same value the summary and review use.
MATCH_TOLERANCE = 5.0


def plot(
    settings: RunSettings,
    lassie_min: float | None = None,
    qseek_min: float | None = None,
    suffix: str = "",
    jitter: bool = True,
) -> Path:
    """Render the comparison figure and return its path.

    `lassie_min` and `qseek_min` drop detections below an absolute level on
    each detector's own scale before plotting, for a stricter view of the same
    runs without re-running them. The run thresholds themselves stay in the
    title so the figure states both what produced the catalog and what it shows.

    With `jitter`, Lassie's positions in the overview panels are displaced
    uniformly within their grid cell, with an offset fixed per detection so a
    detection sits in the same place in every view. Lassie cannot resolve
    position inside a cell, so this loses nothing and stops a cluster on one
    node from collapsing into a single marker. Qseek's positions are its own
    interpolated estimates and are never moved. The close-ups and every
    statistic use exact positions.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pyrocko import model as pyrocko_model

    catalogs = {"lassie": read_lassie(settings), "qseek": read_qseek(settings)}
    shown = {}
    for name, floor in (("lassie", lassie_min), ("qseek", qseek_min)):
        shown[name] = (
            [e for e in catalogs[name] if e.strength >= floor]
            if floor is not None
            else catalogs[name]
        )

    stations = pyrocko_model.load_stations(str(settings.stations_file))
    reference = load_reference(settings)
    spacing = lassie_grid_spacing(settings)
    # Sizes are calibrated on the full catalogs, not the shown subset, so a
    # marker keeps its size between the default and strict views.
    sizes = {name: _size_scale(events) for name, events in catalogs.items()}

    display = dict(shown)
    if jitter:
        offsets = _cell_offsets(catalogs["lassie"], spacing)
        display["lassie"] = _apply_offsets(shown["lassie"], offsets)

    figure, axes = plt.subplots(2, 3, figsize=(21, 11))
    figure.suptitle(
        _title(settings, shown, catalogs, lassie_min, qseek_min, len(stations), spacing)
        + (f"\nLassie positions in the overview panels are jittered within their "
           f"{spacing / 1e3:g} km cell (surface nodes downward only), one fixed offset per "
           "detection; Qseek positions are its interpolated estimates, drawn as is; "
           "close-ups and all statistics use exact positions" if jitter else ""),
        fontsize=10.5,
    )

    _plot_map(axes[0][0], display, sizes, stations, settings, reference, spacing)
    _plot_depth_section(axes[0][1], display, sizes, settings, reference, spacing)
    _plot_zoom(axes[0][2], shown, sizes, settings, reference, spacing, depth_section=False)
    _plot_time_series(axes[1][0], display, sizes, settings)
    _plot_accuracy(axes[1][1], shown, reference)
    _plot_zoom(axes[1][2], shown, sizes, settings, reference, spacing, depth_section=True)

    figure.tight_layout()
    out = settings.results_dir / f"comparison-{settings.run_name}{suffix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out, dpi=150)
    plt.close(figure)
    logger.info("wrote figure to %s", out)
    return out


def _title(settings, shown, catalogs, lassie_min, qseek_min, n_stations, spacing) -> str:
    used = thresholds_used(settings)
    east, north, depth = settings.grid.extents()
    lines = [
        f"ICDP Eger baseline — {settings.run_name} — {n_stations} stations, "
        f"{settings.bandpass[0]}–{settings.bandpass[1]} Hz, shared 1D model, "
        f"{east / 1e3:.0f}x{north / 1e3:.0f}x{depth / 1e3:.0f} km search volume; "
        f"Lassie grid {spacing / 1e3:g} km, Qseek octree "
        f"{settings.grid.root_node_size / 1e3:.1f} km → "
        f"{settings.grid.root_node_size / 2 ** (settings.grid.n_levels - 1):.0f} m with interpolation",
        f"run thresholds — Lassie: {used.get('lassie', 'n/a')}   ·   "
        f"Qseek: {used.get('qseek', 'n/a')}   ·   different scales, each set "
        "against its own noise floor   ·   marker size ∝ detection strength within each catalog",
    ]
    if lassie_min is not None or qseek_min is not None:
        lines.append(
            "STRICT VIEW — shown only: "
            + "   ·   ".join(
                part for part in (
                    f"Lassie image function ≥ {lassie_min:.1f} ({len(shown['lassie'])} of "
                    f"{len(catalogs['lassie'])})" if lassie_min is not None else "",
                    f"Qseek semblance ≥ {qseek_min:.2f} ({len(shown['qseek'])} of "
                    f"{len(catalogs['qseek'])})" if qseek_min is not None else "",
                ) if part
            )
        )
    return "\n".join(lines)


def _event_key(event: Event) -> tuple:
    return (event.time, event.lat, event.lon, event.depth)


def _cell_offsets(events: list[Event], cell: float, seed: int = 0) -> dict:
    """One fixed (north, east, down) offset per detection of the full catalog.

    Keyed by the detection itself, so filtering for a strict view never
    changes where a surviving detection is drawn. Surface nodes are only
    displaced downward: the part of their cell above ground is not a place an
    earthquake can be.
    """
    rng = np.random.default_rng(seed)
    offsets = {}
    for event in events:
        north, east, down = rng.uniform(-cell / 2, cell / 2, size=3)
        if event.depth <= 0.0:
            down = abs(down)
        offsets[_event_key(event)] = (north, east, down)
    return offsets


def _apply_offsets(events: list[Event], offsets: dict) -> list[Event]:
    if not events:
        return events
    moves = np.array([offsets[_event_key(e)] for e in events])
    lats = np.array([e.lat for e in events])
    lons = np.array([e.lon for e in events])
    new_lats, new_lons = od.ne_to_latlon(lats, lons, moves[:, 0], moves[:, 1])
    return [
        dataclasses.replace(event, lat=float(lat), lon=float(lon), depth=event.depth + float(down))
        for event, lat, lon, down in zip(events, new_lats, new_lons, moves[:, 2])
    ]


def _size_scale(events: list[Event]):
    """Map a catalog's strength range onto marker areas.

    The top is the 99th percentile rather than the maximum so that a single
    outlier cannot flatten every other marker to the minimum size. The
    returned callable carries `low`/`high` so legends can state the scale.
    """
    if not events:
        size = lambda e: MARKER_MIN  # noqa: E731
        size.low = size.high = None
        return size
    values = np.array([e.strength for e in events])
    low, high = float(values.min()), float(np.percentile(values, 99))
    span = max(high - low, 1e-9)

    def size(event: Event) -> float:
        return MARKER_MIN + (MARKER_MAX - MARKER_MIN) * min(max((event.strength - low) / span, 0.0), 1.0)

    size.low, size.high = low, high
    return size


def _lassie_grid_nodes(settings: RunSettings, spacing: float):
    """Lassie's fixed grid as (lat, lon) node grids and the depth nodes.

    Lassie places nodes at `xmin + i * dx`, so with symmetric bounds and an
    even spacing the nodes sit at odd-kilometre offsets and none lies on the
    reference point; drawing them explains why detections cluster on a lattice
    that is visibly off-centre.
    """
    grid = settings.grid
    north = np.arange(grid.north_bounds[0], grid.north_bounds[1] + spacing / 2, spacing)
    east = np.arange(grid.east_bounds[0], grid.east_bounds[1] + spacing / 2, spacing)
    depth = np.arange(grid.depth_bounds[0], grid.depth_bounds[1] + spacing / 2, spacing)
    nn, ee = np.meshgrid(north, east, indexing="ij")
    lats, lons = od.ne_to_latlon(grid.lat, grid.lon, nn.ravel(), ee.ravel())
    return lats.reshape(nn.shape), lons.reshape(nn.shape), depth


def _draw_lassie_grid(axis, settings, spacing, depth_section: bool, box=None) -> None:
    """Overlay Lassie's search grid as faint crosses; `box` limits to a view."""
    lats, lons, depths = _lassie_grid_nodes(settings, spacing)
    if depth_section:
        # The lattice is Cartesian, so each east column maps to a slightly
        # different longitude per latitude row; draw the row nearest the
        # reference latitude rather than smearing every row into the section.
        middle = np.argmin(np.abs(lats[:, 0] - settings.grid.lat))
        xs, ys = np.meshgrid(lons[middle], depths / 1e3, indexing="ij")
        xs, ys = xs.ravel(), ys.ravel()
    else:
        xs, ys = lons.ravel(), lats.ravel()
    if box is not None:
        x0, x1, y0, y1 = box
        keep = (xs >= x0) & (xs <= x1) & (ys >= min(y0, y1)) & (ys <= max(y0, y1))
        xs, ys = xs[keep], ys[keep]
    axis.scatter(xs, ys, marker="+", s=14, c="#8b8b8b", alpha=0.45, linewidths=0.6,
                 zorder=1, label=f"Lassie grid nodes ({spacing / 1e3:g} km)")


def _scatter_catalog(axis, name, events, size_of, x_of, y_of, total=None) -> None:
    """Draw one catalog, separating supported locations from the rest.

    Detections that fail the detector's own support criterion — Qseek
    locations with few picks, Lassie locations at the search-volume boundary —
    outnumber the rest and would otherwise hide the very events the comparison
    is about. They stay visible but recede; marker size still follows strength.
    """
    style = STYLE[name]
    strong = [e for e in events if e.is_well_constrained()]
    weak = [e for e in events if not e.is_well_constrained()]

    if weak:
        axis.scatter(
            [x_of(e) for e in weak],
            [y_of(e) for e in weak],
            s=[0.5 * size_of(e) for e in weak],
            alpha=0.18,
            color=style["color"],
            marker=style["marker"],
            linewidths=0,
        )
    if strong:
        axis.scatter(
            [x_of(e) for e in strong],
            [y_of(e) for e in strong],
            s=[size_of(e) for e in strong],
            alpha=0.75,
            color=style["color"],
            marker=style["marker"],
            edgecolors="black",
            linewidths=0.3,
            label=_catalog_label(name, len(strong), len(events), total),
        )
    elif not weak:
        axis.scatter([], [], color=style["color"], marker=style["marker"],
                     label=f"{style['label']} — no run available (N/A)")


def _catalog_label(name, n_strong, n_shown, total) -> str:
    label, criterion = STYLE[name]["label"], SUPPORT_CRITERIA[name]
    if total is not None and total != n_shown:
        return f"{label} — {n_shown} in view of {total}; {n_strong} meet: {criterion}"
    if n_strong != n_shown:
        return f"{label} — {n_shown}; {n_strong} meet: {criterion} (rest faded)"
    return f"{label} — {n_shown}"


def _reference_label(reference) -> str:
    authors = {}
    for event in reference:
        authors.setdefault(event["author"], []).append(event["magnitude"])
    parts = []
    for author, mags in authors.items():
        known = [m for m in mags if m is not None]
        parts.append(f"{len(mags)} {author}" + (f" ML {min(known):.2f}–{max(known):.2f}" if known else " no ML"))
    return "ISC reference: " + ", ".join(parts)


def _draw_reference(axis, reference, x_of, y_of, label=True, number=False) -> None:
    if reference:
        axis.scatter([x_of(e) for e in reference], [y_of(e) for e in reference],
                     marker="*", s=220, c="#d4a017", zorder=6, edgecolors="black",
                     linewidths=0.5, label=_reference_label(reference) if label else None)
        if number:
            # Numbers tie each star to its row in the review's per-event table.
            # Labels fan out around the star so a tight cluster stays legible.
            for index, event in enumerate(reference, start=1):
                angle = math.radians(137.5 * index)
                axis.annotate(str(index), (x_of(event), y_of(event)), fontsize=6,
                              xytext=(7 * math.cos(angle), 7 * math.sin(angle)),
                              textcoords="offset points", zorder=7,
                              ha="center", va="center",
                              bbox={"boxstyle": "circle,pad=0.1", "fc": "white", "ec": "none", "alpha": 0.7})
    elif label:
        # Saying so is the point: with no catalogued event that day there is
        # nothing to validate either catalog against.
        axis.scatter([], [], marker="*", s=220, c="#d4a017", edgecolors="black",
                     linewidths=0.5, label="ISC reference — none this day")


def _plot_map(axis, catalogs, sizes, stations, settings, reference, spacing) -> None:
    _draw_lassie_grid(axis, settings, spacing, depth_section=False)
    axis.scatter([s.lon for s in stations], [s.lat for s in stations], marker="v", s=90,
                 c="#333333", zorder=5, label=f"stations ({len(stations)})")
    for station in stations:
        axis.annotate(station.station, (station.lon, station.lat), fontsize=7,
                      xytext=(4, 4), textcoords="offset points")

    for name, events in catalogs.items():
        _scatter_catalog(axis, name, events, sizes[name], lambda e: e.lon, lambda e: e.lat)
    _draw_reference(axis, reference, lambda e: e["lon"], lambda e: e["lat"])

    scale = "; ".join(
        f"{STYLE[name]['label'].split(' (')[0]} {sizes[name].low:.3g}→{sizes[name].high:.3g}"
        for name in catalogs if sizes[name].low is not None
    )
    axis.text(0.99, 0.01, f"marker area spans strength {scale}", transform=axis.transAxes,
              ha="right", va="bottom", fontsize=6.5, color="#444444")
    axis.set_xlabel("longitude [°E]")
    axis.set_ylabel("latitude [°N]")
    axis.set_title("Epicentres")
    axis.legend(fontsize=7, loc="upper left", framealpha=0.85)
    axis.grid(alpha=0.3)


def _plot_depth_section(axis, catalogs, sizes, settings, reference, spacing) -> None:
    _draw_lassie_grid(axis, settings, spacing, depth_section=True)
    for name, events in catalogs.items():
        _scatter_catalog(axis, name, events, sizes[name], lambda e: e.lon, lambda e: e.depth / 1e3)
    _draw_reference(axis, reference, lambda e: e["lon"], lambda e: e["depth"] / 1e3)
    axis.invert_yaxis()
    axis.set_xlabel("longitude [°E]")
    axis.set_ylabel("depth [km]")
    axis.set_title("Depth section (E–W)")
    axis.legend(fontsize=7, loc="lower left", framealpha=0.85)
    axis.grid(alpha=0.3)


def _coincident_counts(events: list[Event]) -> dict:
    """How many detections share each exact position (a Lassie node)."""
    counts: dict = {}
    for event in events:
        key = (round(event.lat, 6), round(event.lon, 6), round(event.depth, 1))
        counts[key] = counts.get(key, 0) + 1
    return counts


def _plot_zoom(axis, catalogs, sizes, settings, reference, spacing, depth_section: bool) -> None:
    """Close-up of the reference cluster with each match joined to its event.

    This is the visual confirmation the summary numbers rest on: a line from a
    star to a marker is one catalogued earthquake and the detection credited
    with recovering it, drawn at exact positions. Matches that lie outside the
    view are listed on the panel with their distance and depth rather than
    silently clipped, and detections sharing one node carry a "×n" count.
    """
    if not reference:
        axis.text(0.5, 0.5, "no reference events for this day",
                  ha="center", va="center", transform=axis.transAxes, fontsize=11)
        axis.set_axis_off()
        return

    lat0 = float(np.mean([e["lat"] for e in reference]))
    lon0 = float(np.mean([e["lon"] for e in reference]))
    dlat = ZOOM_HALF_WIDTH_M / 111_195.0
    dlon = dlat / math.cos(math.radians(lat0))
    half_km = ZOOM_HALF_WIDTH_M / 1e3

    if depth_section:
        y_of = lambda e: e.depth / 1e3  # noqa: E731
        y_ref = lambda r: r["depth"] / 1e3  # noqa: E731
        box = (lon0 - dlon, lon0 + dlon,
               settings.grid.depth_bounds[0] / 1e3, settings.grid.depth_bounds[1] / 1e3)
    else:
        y_of = lambda e: e.lat  # noqa: E731
        y_ref = lambda r: r["lat"]  # noqa: E731
        box = (lon0 - dlon, lon0 + dlon, lat0 - dlat, lat0 + dlat)

    def in_view(event) -> bool:
        return abs(event.lon - lon0) <= dlon and abs(event.lat - lat0) <= dlat

    _draw_lassie_grid(axis, settings, spacing, depth_section=depth_section, box=box)

    for name, events in catalogs.items():
        inside = [e for e in events if in_view(e)]
        _scatter_catalog(axis, name, inside, sizes[name], lambda e: e.lon, y_of, total=len(events))
        if name == "lassie":
            for (lat, lon, depth), count in _coincident_counts(inside).items():
                if count > 1:
                    axis.annotate(f"×{count}", (lon, depth / 1e3 if depth_section else lat),
                                  fontsize=6.5, color=STYLE[name]["color"], xytext=(5, -3),
                                  textcoords="offset points", zorder=7)

    off_panel = []
    for index, event in enumerate(reference, start=1):
        for name, events in catalogs.items():
            match = nearest(events, event, tolerance=MATCH_TOLERANCE)
            if match is None:
                continue
            axis.plot([event["lon"], match.lon], [y_ref(event), y_of(match)],
                      color=STYLE[name]["color"], linewidth=0.9, alpha=0.8, zorder=4)
            if not in_view(match):
                distance = od.distance_accurate50m(event["lat"], event["lon"], match.lat, match.lon) / 1e3
                off_panel.append(f"{name} #{index}: {distance:.1f} km away, {match.depth / 1e3:.0f} km deep")

    _draw_reference(axis, reference, lambda e: e["lon"], y_ref, number=True)

    if off_panel:
        # Opposite corner from the legend, so the two never overlap.
        corner = (0.99, 0.99, "right", "top") if depth_section else (0.99, 0.01, "right", "bottom")
        axis.text(corner[0], corner[1], "matches outside this view:\n" + "\n".join(off_panel),
                  transform=axis.transAxes, ha=corner[2], va=corner[3], fontsize=6.5,
                  bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#888888", "alpha": 0.9},
                  zorder=8)

    axis.set_xlim(box[0], box[1])
    if depth_section:
        axis.set_ylim(box[3], box[2])
        axis.set_ylabel("depth [km]")
        axis.set_title(f"Source zone ±{half_km:.0f} km, full 0–{box[3]:.0f} km depth, exact positions",
                       fontsize=10)
    else:
        axis.set_ylim(box[2], box[3])
        axis.set_aspect(1 / math.cos(math.radians(lat0)))
        axis.set_ylabel("latitude [°N]")
        axis.set_title(f"Source zone ±{half_km:.0f} km, exact positions — lines join each "
                       "catalogued event to its recovering detection", fontsize=10)
    axis.set_xlabel("longitude [°E]")
    axis.legend(fontsize=6.5, loc="lower left" if depth_section else "upper left", framealpha=0.85)
    axis.grid(alpha=0.3)


def _plot_time_series(axis, catalogs, sizes, settings) -> None:
    for name, events in catalogs.items():
        if not events:
            continue
        hours = [(e.time - settings.tmin).total_seconds() / 3600 for e in events]
        # Each detector's strength is on its own scale, so normalise to its own
        # maximum; only the relative pattern within a catalog is meaningful.
        peak = max(e.strength for e in events) or 1.0
        levels = [e.strength / peak for e in events]
        axis.vlines(hours, 0, levels, color=STYLE[name]["color"], linewidth=0.6, alpha=0.5)
        axis.scatter(hours, levels, s=[0.6 * sizes[name](e) for e in events],
                     color=STYLE[name]["color"], marker=STYLE[name]["marker"],
                     alpha=0.8, label=f"{STYLE[name]['label']} — {len(events)}")
    axis.set_xlim(0, settings.hours)
    axis.set_xlabel(f"hours after {settings.day_str} 00:00 UTC")
    axis.set_ylabel("detector strength / catalog maximum")
    axis.set_title("Detections over the day (all shown detections)")
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)


def _plot_accuracy(axis, catalogs, reference) -> None:
    """Compare both detectors against the independent reference catalog.

    Bars are medians; the dots are the individual events, because a median
    hides a bimodal error — Lassie's depth error, for instance, is either
    about +2 km or about -10 km and its median describes neither. A catalog
    with no run is drawn as N/A rather than as zero error.
    """
    if not reference:
        axis.text(0.5, 0.5, "no reference events for this day\nnothing to validate against",
                  ha="center", va="center", transform=axis.transAxes, fontsize=11)
        axis.set_axis_off()
        return

    metrics = [
        ("origin time\n[s]", "origin_time_median_s", lambda m, r: abs((m.time - r["time"]).total_seconds())),
        ("epicentre\n[km]", "epicentre_median_km",
         lambda m, r: od.distance_accurate50m(r["lat"], r["lon"], m.lat, m.lon) / 1e3),
        ("depth\n[km]", "depth_median_km", lambda m, r: abs(m.depth - r["depth"]) / 1e3),
    ]
    positions = np.arange(len(metrics))
    width = 0.38
    rng = np.random.default_rng(0)

    for offset, name in ((-width / 2, "lassie"), (width / 2, "qseek")):
        events = catalogs[name]
        if not events:
            for index in range(len(metrics)):
                axis.text(index + offset, 0.02, "N/A", ha="center", va="bottom",
                          transform=axis.get_xaxis_transform(), fontsize=8, color=STYLE[name]["color"])
            axis.bar([], [], color=STYLE[name]["color"], label=f"{STYLE[name]['label']} — no run available")
            continue
        scores = accuracy(events, reference, MATCH_TOLERANCE)
        found = recall(events, reference, MATCH_TOLERANCE)
        medians = [scores.get(key, 0.0) for _, key, _ in metrics]
        bars = axis.bar(positions + offset, medians, width, color=STYLE[name]["color"], alpha=0.55,
                        label=f"{STYLE[name]['label']} — recovered "
                              f"{found['recovered']}/{found['reference_events']} within ±{MATCH_TOLERANCE:g} s "
                              "(bar = median, dots = events)")
        for bar, value in zip(bars, medians):
            axis.text(bar.get_x() + bar.get_width() / 2, value, f"{value:g}",
                      ha="center", va="bottom", fontsize=9)
        for index, (_, _, error_of) in enumerate(metrics):
            errors = [error_of(match, event) for event in reference
                      if (match := nearest(events, event, MATCH_TOLERANCE)) is not None]
            jitter = rng.uniform(-width / 4, width / 4, size=len(errors))
            axis.scatter(index + offset + jitter, errors, s=14, color=STYLE[name]["color"],
                         edgecolors="black", linewidths=0.3, zorder=5)

    # A detector-free guess of the reference's own median depth scores this;
    # depth agreement below it says nothing about depth resolution.
    depths = np.array([e["depth"] for e in reference]) / 1e3
    null = float(np.median(np.abs(depths - np.median(depths))))
    axis.hlines(null, 2 - 0.45, 2 + 0.45, colors="black", linestyles="--", linewidth=1,
                label=f"constant-depth guess ({np.median(depths):.1f} km) scores {null:.2f} km")

    axis.set_xticks(list(positions))
    axis.set_xticklabels([label for label, _, _ in metrics])
    axis.set_ylabel("|detection − reference| (agreement with the catalog, not absolute accuracy)")
    axis.set_title(f"Agreement with the ISC reference — {_reference_label(reference)[15:]}", fontsize=10)
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3, axis="y")
