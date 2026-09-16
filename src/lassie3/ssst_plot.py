"""Diagnostic figure for derived station terms."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from pyrocko import model as pyrocko_model
from pyrocko import orthodrome as od

from lassie3 import reference
from lassie3.settings import RunSettings
from lassie3.ssst import Pick, StationTerms

logger = logging.getLogger(__name__)

PHASE_COLOURS = {"P": "#1f77b4", "S": "#d62728"}


def plot_terms(
    settings: RunSettings,
    terms: StationTerms,
    picks: list[Pick],
    events: list[dict],
    path: Path,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stations = pyrocko_model.load_stations(str(settings.stations_file))
    fig, axes = plt.subplots(2, 3, figsize=(21, 12.5))
    (ax_map, ax_bars, ax_hist), (ax_p, ax_s, ax_text) = axes

    _plot_map(ax_map, settings, terms, picks, events, stations)
    _plot_static(ax_bars, terms)
    _plot_histograms(ax_hist, picks, terms)

    station = _busiest_station(terms)
    depth = float(np.median(terms.events[:, 2])) if terms.events.size else 10_000.0
    for axis, phase in ((ax_p, "P"), (ax_s, "S")):
        _plot_field(axis, settings, terms, stations, station, phase, depth)
    _plot_text(ax_text, terms)

    fig.suptitle(
        f"Station terms for {settings.day_str}, learned from WBNET events on the other days — "
        f"{terms.meta.get('n_events_with_terms', 0)} events, "
        f"{terms.meta.get('n_picks_kept', 0)} of {terms.meta.get('n_picks', 0)} PhaseNet picks kept",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    logger.info("wrote %s", path)
    return path


def _search_box(settings: RunSettings):
    grid = settings.grid
    north = np.array([grid.north_bounds[0], grid.north_bounds[1], grid.north_bounds[1], grid.north_bounds[0], grid.north_bounds[0]])
    east = np.array([grid.east_bounds[0], grid.east_bounds[0], grid.east_bounds[1], grid.east_bounds[1], grid.east_bounds[0]])
    return od.ne_to_latlon(grid.lat, grid.lon, north, east)


def _plot_map(axis, settings, terms, picks, events, stations) -> None:
    used = {p.event for p in picks if p.kept}
    lat = np.array([e["lat"] for e in events])
    lon = np.array([e["lon"] for e in events])
    day = np.array([e["time"].day for e in events])
    mask = np.array([e["id"] in used for e in events])

    if (~mask).any():
        axis.scatter(lon[~mask], lat[~mask], s=12, c="#bbbbbb", label=f"reference events without terms ({(~mask).sum()})")
    if mask.any():
        handle = axis.scatter(lon[mask], lat[mask], s=14, c=day[mask], cmap="viridis", label=f"reference events used ({mask.sum()})")
        bar = axis.figure.colorbar(handle, ax=axis, fraction=0.04, pad=0.02)
        bar.set_label("day of March 2024")

    held_out = reference.load(settings)
    if held_out:
        axis.scatter([e["lon"] for e in held_out], [e["lat"] for e in held_out], marker="*", s=110,
                     c="gold", edgecolors="black", linewidths=0.5, zorder=6,
                     label=f"{settings.day_str} reference events (held out, {len(held_out)})")

    axis.scatter([s.lon for s in stations], [s.lat for s in stations], marker="v", s=80, c="#333333", zorder=5,
                 label=f"stations ({len(stations)})")
    for station in stations:
        axis.annotate(station.station, (station.lon, station.lat), fontsize=7, xytext=(4, 4), textcoords="offset points")

    box_lat, box_lon = _search_box(settings)
    axis.plot(box_lon, box_lat, color="#888888", lw=0.8, ls="--", label="search volume")
    axis.set_aspect(1.0 / np.cos(np.radians(settings.grid.lat)))
    axis.set_xlabel("longitude")
    axis.set_ylabel("latitude")
    axis.set_title("Reference events the terms are learned from")
    axis.legend(fontsize=7, loc="lower left")


def _plot_static(axis, terms) -> None:
    stations = terms.stations()
    counts = terms.meta.get("counts", {})
    x = np.arange(len(stations))
    width = 0.38
    for offset, phase in ((-width / 2, "P"), (width / 2, "S")):
        values = [terms.static[s].get(phase, 0.0) for s in stations]
        errors = []
        for s in stations:
            stats = counts.get(f"{s}:{phase}")
            errors.append(stats["std"] / np.sqrt(stats["n"]) if stats and stats["n"] else 0.0)
        axis.bar(x + offset, values, width, yerr=errors, capsize=2, color=PHASE_COLOURS[phase],
                 alpha=0.85, label=f"{phase} static term ± standard error")
        for xi, s, v in zip(x + offset, stations, values):
            stats = counts.get(f"{s}:{phase}")
            if stats:
                axis.annotate(f"n={stats['n']}", (xi, v), fontsize=6, ha="center",
                              xytext=(0, 4 if v >= 0 else -10), textcoords="offset points")
    axis.axhline(0, color="black", lw=0.6)
    axis.set_xticks(x)
    axis.set_xticklabels(stations, rotation=35, ha="right", fontsize=8)
    axis.set_ylabel("term (s), observed − modelled, after per-event demeaning")
    axis.set_title("Static station terms (the fallback far from any reference event)")
    axis.legend(fontsize=8)


def _plot_histograms(axis, picks, terms) -> None:
    bins = np.arange(-2.0, 2.01, 0.05)
    for phase in ("P", "S"):
        kept = [p.demeaned for p in picks if p.phase == phase and p.kept and p.demeaned is not None]
        dropped = [p.demeaned for p in picks if p.phase == phase and not p.kept and p.demeaned is not None]
        axis.hist(np.clip(kept, bins[0], bins[-1]), bins=bins, color=PHASE_COLOURS[phase], alpha=0.5,
                  label=f"{phase} kept ({len(kept)})")
        if dropped:
            axis.hist(np.clip(dropped, bins[0], bins[-1]), bins=bins, histtype="step",
                      color=PHASE_COLOURS[phase], lw=1.2, label=f"{phase} rejected ({len(dropped)})")
    axis.axvline(0, color="black", lw=0.6)
    axis.set_xlabel("demeaned residual (s), clipped to ±2 s")
    axis.set_ylabel("picks")
    axis.set_title(f"Residuals entering the terms (outliers beyond {terms.settings.outlier_level:g} robust σ rejected)")
    axis.legend(fontsize=8)


def _busiest_station(terms) -> str:
    best, best_n = "", -1
    for (station, phase), (idx, _) in terms.residuals.items():
        if phase == "P" and idx.size > best_n:
            best, best_n = station, idx.size
    return best or (terms.stations()[0] if terms.stations() else "")


def _plot_field(axis, settings, terms, stations, station, phase, depth) -> None:
    grid = settings.grid
    north = np.linspace(grid.north_bounds[0], grid.north_bounds[1], 121)
    east = np.linspace(grid.east_bounds[0], grid.east_bounds[1], 121)
    nn, ee = np.meshgrid(north, east, indexing="ij")
    lat, lon = od.ne_to_latlon(grid.lat, grid.lon, nn.ravel(), ee.ravel())
    # Shown relative to the static term: the absolute term is one number per
    # station, the map is about where the reference events pull it away.
    static = terms.static.get(station, {}).get(phase, 0.0)
    field = terms.evaluate(lat, lon, np.full(lat.shape, depth), station, phase).reshape(nn.shape) - static
    limit = max(0.02, float(np.abs(field).max()))
    mesh = axis.pcolormesh(lon.reshape(nn.shape), lat.reshape(nn.shape), field, cmap="RdBu_r",
                           vmin=-limit, vmax=limit, shading="auto")
    axis.figure.colorbar(mesh, ax=axis, fraction=0.04, pad=0.02).set_label(f"{phase} term − static term (s)")

    ev_lat, ev_lon = od.ne_to_latlon(terms.reference_lat, terms.reference_lon, terms.events[:, 0], terms.events[:, 1])
    axis.scatter(ev_lon, ev_lat, s=4, c="black", alpha=0.6, label="reference events")
    for s in stations:
        marker = "v"
        colour = "red" if f"{s.network}.{s.station}" == station else "#333333"
        axis.scatter(s.lon, s.lat, marker=marker, s=70, c=colour, zorder=5)
        axis.annotate(s.station, (s.lon, s.lat), fontsize=7, xytext=(4, 4), textcoords="offset points")
    axis.set_aspect(1.0 / np.cos(np.radians(grid.lat)))
    axis.set_title(f"{phase} term of {station} (red) at {depth / 1e3:.1f} km depth, static {static:+.3f} s — "
                   f"radius {terms.settings.radius / 1e3:g} km, prior weight {terms.settings.prior_weight:g}", fontsize=9)
    axis.set_xlabel("longitude")
    axis.set_ylabel("latitude")


def _plot_text(axis, terms) -> None:
    meta = terms.meta
    lines = [
        "How the terms are made",
        "",
        f"reference catalog: {meta.get('reference_catalog', '')}",
        f"training days: {', '.join(d[5:] for d in meta.get('training_days', []))}",
        f"evaluation day held out: {meta.get('evaluation_day', '')}",
        f"picker: {meta.get('picker', '')}",
        f"search window: ±{terms.settings.search_window[0]:g} s (P), ±{terms.settings.search_window[1]:g} s (S) around the Cake arrival",
        f"minimum pick probability: {terms.settings.min_probability:g}",
        "",
        "residual = pick − (origin time + Cake travel time from the WBNET hypocentre)",
        "per event, the mean over all its picks is removed (origin-time freedom)",
        f"outliers beyond {terms.settings.outlier_level:g} robust σ per phase are rejected, then demeaned again",
        "",
        "term at a node = Σ w·residual + prior·static over events within the radius, / (Σ w + prior)",
        f"w = (1 − (d/{terms.settings.radius / 1e3:g} km)³)³ (SCOTER's bicube), prior = {terms.settings.prior_weight:g} pseudo-events",
        f"static term = mean demeaned residual per station and phase (≥{terms.settings.min_residuals} residuals, else 0)",
        "",
        f"events: {meta.get('n_reference_events', 0)} referenced, {meta.get('n_events_with_terms', 0)} with ≥{terms.settings.min_picks_per_event} kept picks",
        f"picks: {meta.get('n_picks_kept', 0)} kept of {meta.get('n_picks', 0)}",
    ]
    axis.axis("off")
    axis.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=9, family="monospace", transform=axis.transAxes)
