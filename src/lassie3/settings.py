"""Shared run parameters for the Lassie-v1 / Qseek baseline comparison.

Both detectors are driven from the same values here so that differences in the
resulting catalogs come from the detectors, not from the setup.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_DIR / "data_icdp_eger"
RESULTS_DIR = PROJECT_DIR / "results"

#: Day selected as the baseline: all nine open stations have uninterrupted
#: 100 Hz HH? coverage, verified against the miniSEED headers.
DEFAULT_DAY = dt.date(2024, 3, 20)

#: Approximate layered model for the West Bohemia / Vogtland swarm region,
#: as depth[km] vp[km/s] vs[km/s] rho[g/cm^3] rows in Pyrocko's `nd` format.
#: Vs follows Vp/Vs = 1.70, the value commonly reported for the Novy Kostel
#: focal zone. This is a representative model, not a verbatim published one --
#: swap in your preferred model with `--velocity-model` and both detectors pick
#: it up, since they read the same file.
#:
#: The profile starts above sea level because Qseek traces from true station
#: elevations, and Cake refuses to extract a layer outside the model. SX.TANN
#: sits at 870 m, so -1.5 km leaves margin for any station added later.
WEST_BOHEMIA_ND = """\
-1.5 5.20 3.06 2.55
 0.0 5.20 3.06 2.55
 1.0 5.40 3.18 2.60
 2.0 5.70 3.35 2.65
 4.0 5.90 3.47 2.68
 6.0 6.05 3.56 2.70
 8.0 6.15 3.62 2.72
10.0 6.25 3.68 2.74
12.0 6.35 3.74 2.76
15.0 6.50 3.82 2.79
20.0 6.70 3.94 2.83
28.0 7.00 4.12 2.90
mantle
31.0 8.00 4.71 3.30
50.0 8.10 4.76 3.35
"""


@dataclass(frozen=True)
class Grid:
    """Search volume shared by both detectors.

    Offsets are metres from (`lat`, `lon`); north and east are positive, and
    `depth` grows downwards from sea level.
    """

    lat: float
    lon: float
    east_bounds: tuple[float, float]
    north_bounds: tuple[float, float]
    depth_bounds: tuple[float, float]

    #: Edge length of an unrefined Qseek octree node. Every extent must be an
    #: integer multiple of it or Qseek rejects the octree.
    root_node_size: float

    #: Qseek halves nodes `n_levels - 1` times, so the finest node is
    #: `root_node_size / 2 ** (n_levels - 1)`.
    n_levels: int

    #: Lassie v1 stacks a fixed regular grid instead of refining an octree, so
    #: it needs its own, necessarily coarser, node spacing.
    lassie_spacing: float

    def extents(self) -> tuple[float, float, float]:
        return (
            self.east_bounds[1] - self.east_bounds[0],
            self.north_bounds[1] - self.north_bounds[0],
            self.depth_bounds[1] - self.depth_bounds[0],
        )


#: Centred so that the SX network footprint and the Novy Kostel focal zone
#: (~10 km south, at the very edge of the array) both fall well inside.
FULL_APERTURE_GRID = Grid(
    lat=50.32,
    lon=12.35,
    east_bounds=(-25_000.0, 25_000.0),
    north_bounds=(-25_000.0, 25_000.0),
    depth_bounds=(0.0, 20_000.0),
    root_node_size=2_500.0,
    n_levels=4,
    lassie_spacing=2_000.0,
)


@dataclass(frozen=True)
class RunSettings:
    """Everything a single baseline run needs."""

    day: dt.date = DEFAULT_DAY
    grid: Grid = FULL_APERTURE_GRID
    data_dir: Path = DATA_DIR
    results_dir: Path = RESULTS_DIR

    #: Restricts both runs to the broadband 100 Hz channels. SX.GUNZ, SX.MULD,
    #: SX.TANN and SX.WERN also carry BH? and LH?, which would otherwise enter
    #: the stack as duplicate, lower-rate copies of the same ground motion.
    #: GQ.LNDWU adds a 200 Hz rotational HJ? triplet, which is rotation rate
    #: rather than translational ground motion and is excluded the same way.
    channels: tuple[str, ...] = ("HHE", "HHN", "HHZ")

    #: Band common to both detectors, above the regional microseism peak and
    #: below the 50 Hz Nyquist.
    bandpass: tuple[float, float] = (2.0, 30.0)

    #: Lassie v1 only supports an absolute level on the stacked image
    #: function, whose scale depends on the station count, weights and noise of
    #: the particular run. "auto" derives it from the run's own image function
    #: instead, which is what makes it comparable to Qseek's adaptive default.
    lassie_detector_threshold: float | str = "auto"

    #: Multiples of the median absolute deviation above the median image-function
    #: level at which "auto" places the Lassie threshold.
    #:
    #: The image function has a long upper tail, so detection count falls very
    #: steeply here: on 2024-03-20 this array yields ~1750 peaks/day at 12 MAD,
    #: 224 at 20, 74 at 30 and 42 at 40. Below roughly 20 the catalog is
    #: dominated by noise, visible as detections scattered evenly across the
    #: whole search grid rather than clustered near the array.
    lassie_mad_multiple: float = 30.0

    #: Qseek defaults to "MAD", a noise-adaptive threshold on the semblance.
    qseek_detection_threshold: str | float = "MAD"

    #: Threads PhaseNet inference may use. Qseek's own default is 4, which
    #: leaves most of an 11-core machine idle; None means "all cores".
    qseek_torch_threads: int | None = None

    #: Torch device for PhaseNet. Qseek only knows CUDA and CPU, so "mps" is
    #: applied by this project after Qseek has built the model.
    qseek_device: str = "cpu"

    #: Length of the processed window from 00:00 UTC. Shorten it to smoke-test
    #: a configuration change without paying for a full day of inference.
    hours: float = 24.0

    #: Lassie renders one matplotlib figure per detection, which dominates the
    #: runtime once a day yields hundreds of them. Off by default; enable it
    #: when inspecting individual detections rather than comparing catalogs.
    lassie_save_figures: bool = False

    velocity_model: Path | None = None
    stations_blacklist: tuple[str, ...] = field(default_factory=tuple)

    #: Suffix for variant runs (a finer Lassie grid, a threshold sweep) so they
    #: land next to the baseline instead of overwriting it.
    tag: str = ""

    #: `terms.json` from `lassie3 ssst-terms`. When set, both detectors add the
    #: source-specific station terms to their Cake travel times at every search
    #: node; the baseline runs leave it unset.
    station_terms: Path | None = None

    @property
    def tmin(self) -> dt.datetime:
        return dt.datetime.combine(self.day, dt.time.min, tzinfo=dt.timezone.utc)

    @property
    def tmax(self) -> dt.datetime:
        return self.tmin + dt.timedelta(hours=self.hours)

    @property
    def day_str(self) -> str:
        return self.day.isoformat()

    @property
    def stage_dir(self) -> Path:
        """Waveforms for `day` only, symlinked out of the 14-day archive.

        Both Lassie and Qseek index every file they are pointed at, so handing
        them the whole archive would cost a 2.2 GB scan per run.
        """
        return self.results_dir / "_stage" / self.day_str

    @property
    def squirrel_dir(self) -> Path:
        """Squirrel environment kept inside the project so runs stay self-contained.

        Without it Squirrel falls back to a cache in the user's home directory.
        """
        return self.results_dir / "_stage" / "squirrel"

    @property
    def stations_file(self) -> Path:
        return self.results_dir / "_stage" / "stations.yaml"

    @property
    def velocity_model_file(self) -> Path:
        return self.velocity_model or (self.results_dir / "_stage" / "velocity.nd")

    @property
    def run_name(self) -> str:
        name = f"eger-{self.day_str}"
        if self.hours != 24.0:
            name += f"-{self.hours:g}h"
        if self.tag:
            name += f"-{self.tag}"
        return name

    def run_dir(self, approach: str) -> Path:
        return self.results_dir / approach
