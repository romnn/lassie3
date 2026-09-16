"""Lassie v1 shifter that adds station terms to its Cake travel times.

Lassie's `Shifter` is a guts extension point: an image-function contribution
takes any `Shifter` subclass, and `core.search` only ever asks it for the
travel-time table. This subclass is therefore a plugin rather than a patch —
Lassie's stacking, peak search and detection code run unmodified and merely
see a different table. The class is registered under the `lassie3` guts
prefix, so a persisted config names it as `!lassie3.CorrectedCakePhaseShifter`
and re-loading that config needs this module imported first.
"""

from __future__ import annotations

import numpy as np
from lassie import shifter as lassie_shifter
from pyrocko.guts import String

from lassie3.ssst import StationTerms

guts_prefix = "lassie3"


class CorrectedCakePhaseShifter(lassie_shifter.CakePhaseShifter):
    """Cake travel times plus the source-specific term of every node/receiver pair."""

    terms_path = String.T(help="terms.json written by `lassie3 ssst-terms`")
    phase = String.T(help='which term family to add, "P" or "S"')

    def get_table(self, grid, receivers):
        table = super().get_table(grid, receivers)
        terms = StationTerms.load(self.terms_path)
        lat, lon, depth = grid_node_coordinates(grid)
        for index, receiver in enumerate(receivers):
            station = ".".join(receiver.codes[:2])
            table[:, index] += terms.evaluate(lat, lon, depth, station, self.phase)
        return table


def grid_node_coordinates(grid):
    """Latitude, longitude and depth of every node, in the table's row order.

    Lassie lays the table out depth-major: `depths()` repeats each level over
    the whole surface, and `surface_points()` walks the surface in the same
    order `lateral_distances()` does, so tiling the surface once per level
    reproduces the row order without touching the grid's private state.
    """
    surface_lat, surface_lon = grid.surface_points(system="latlon")
    depth = grid.depths()
    levels = depth.size // surface_lat.size
    return np.tile(surface_lat, levels), np.tile(surface_lon, levels), depth
