"""Qseek travel-time corrections carrying source-specific station terms.

Qseek accepts any `TravelTimeCorrections` subclass as `station_corrections`
and adds its delays to the ray tracer's travel times before stacking. It
collects the subclasses into its config union when `qseek.search` is first
imported, so this module must be imported before that — `run_qseek` does so.
Like the Lassie shifter this is a plugin: Qseek's search, octree refinement
and detection code run unmodified.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Literal, Sequence

import numpy as np
from pydantic import PrivateAttr
from pyrocko import orthodrome as od
from qseek.corrections.base import TravelTimeCorrections

from lassie3.ssst import StationTerms

if TYPE_CHECKING:
    from qseek.models.station import StationInventory
    from qseek.octree import Node, Octree
    from qseek.utils import NSL, PhaseDescription

logger = logging.getLogger(__name__)


class SSSTCorrections(TravelTimeCorrections):
    """Per-node, per-station delays evaluated from `terms.json`."""

    corrections: Literal["SSSTCorrections"] = "SSSTCorrections"

    terms_path: Path

    _terms: StationTerms | None = PrivateAttr(default=None)
    _octree: Octree | None = PrivateAttr(default=None)

    @property
    def terms(self) -> StationTerms:
        if self._terms is None:
            self._terms = StationTerms.load(self.terms_path)
        return self._terms

    @property
    def n_stations(self) -> int:
        return len(self.terms.stations())

    async def prepare(
        self,
        stations: StationInventory,
        octree: Octree,
        phases: Iterable[PhaseDescription],
        rundir: Path,
    ) -> None:
        self._octree = octree
        known = set(self.terms.stations())
        missing = [sta.nsl.pretty for sta in stations if f"{sta.nsl.network}.{sta.nsl.station}" not in known]
        if missing:
            logger.warning("no station terms for %s; they get 0 s", ", ".join(missing))
        # Keep a copy next to the run so the catalog stays explainable after
        # the terms are re-derived.
        (rundir / "station-terms.json").write_text(Path(self.terms_path).read_text())
        logger.info(
            "station terms from %s: %d stations, %d reference events, radius %.0f m",
            self.terms_path, self.n_stations, self.terms.events.shape[0], self.terms.settings.radius,
        )

    @staticmethod
    def _phase_name(phase: str) -> str:
        """Qseek names phases like "cake:P"; the terms are keyed "P" and "S"."""
        return phase.split(":")[-1].upper()

    @staticmethod
    def _station_key(nsl: NSL) -> str:
        return f"{nsl.network}.{nsl.station}"

    def _node_coordinates(self, nodes: Sequence[Node]):
        origin = (self._octree or nodes[0].tree).location
        lat0, lon0 = origin.effective_lat_lon
        north = np.array([node.north for node in nodes], dtype=float)
        east = np.array([node.east for node in nodes], dtype=float)
        # Node depths are measured from the octree's reference location, whose
        # effective depth is its depth below sea level.
        depth = np.array([node.depth for node in nodes], dtype=float) + origin.effective_depth
        lat, lon = od.ne_to_latlon(lat0, lon0, north, east)
        return lat, lon, depth

    def get_delay(self, station_nsl: NSL, phase: PhaseDescription, node: Node | None = None) -> float:
        station = self._station_key(station_nsl)
        name = self._phase_name(phase)
        if node is None:
            return self.terms.static.get(station, {}).get(name, 0.0)
        lat, lon, depth = self._node_coordinates([node])
        return float(self.terms.evaluate(lat, lon, depth, station, name)[0])

    async def get_delays(
        self,
        station_nsls: Sequence[NSL],
        phase: PhaseDescription,
        nodes: Sequence[Node],
    ) -> np.ndarray:
        lat, lon, depth = self._node_coordinates(nodes)
        name = self._phase_name(phase)
        return np.column_stack([
            self.terms.evaluate(lat, lon, depth, self._station_key(nsl), name)
            for nsl in station_nsls
        ])
