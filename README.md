# lassie3 — Lassie v1 vs. Qseek baseline on ICDP Eger

A reproducible side-by-side run of the two Pyrocko earthquake detectors over one
day of the ICDP-EGER data set, so their catalogs can be compared on identical
inputs.

Both detectors stack shifted characteristic functions over a search volume, but
differ in both halves of that:

| | Lassie v1 | Qseek |
|---|---|---|
| characteristic function | STA/LTA onset (P) + envelope wave packet (S) | PhaseNet P/S probabilities (SeisBench) |
| search volume | fixed regular grid | octree, refined toward the peak |
| threshold | absolute level on the image function | adaptive: 10x the MAD of each window's semblance, as both peak height and prominence |

Shared on purpose: the same day, the same stations at their true
elevations, the same `HH?` channels, the same 1D velocity model, the same
50x50x20 km search volume and reference point, and the same 2-30 Hz band for
the waveform/P side. **Not** shared, because each is a stock capability of one
tool only, and worth knowing before attributing a difference to "the
algorithm":

| | Lassie v1 | Qseek |
|---|---|---|
| location resolution | 2 km fixed nodes (a 1 km variant is run for comparison) | 312 m finest octree node, plus interpolation between nodes |
| boundary handling | none | rejects detections within 2.5 km of the lateral/bottom faces |
| station weighting | none | distance taper (`mean_interstation`), at least 4 closest stations |
| filtering | causal 4th-order high/low-pass per characteristic function | zero-phase 4th-order bandpass before PhaseNet |
| S band | 2-15 Hz wave-packet detector | 2-30 Hz waveform band into PhaseNet |

## Quick start

```sh
task sync        # install from uv.lock, then repair the OpenMP linkage
task all         # prepare -> lassie -> qseek -> figure -> summary
```

`task smoke` runs the same pipeline over one hour, for checking a config change
without paying for a full day.

Run `task --list` for everything else. Any task takes `DAY` and `HOURS`:

```sh
task all DAY=2024-03-21 HOURS=6
```

## Outputs

```
results/
  _stage/              inputs both detectors read (waveforms, stations, velocity model)
  lassie/<run>.turd/   detections.list, events.list, ifm/ stacks, figures/
  qseek/<run>/         csv/detections.csv, search.json, cake/ plots
  comparison-<run>.png        six-panel figure, including source-zone close-ups
  comparison-<run>-strict.png the same runs shown above stricter levels (`task plot:strict`)
  review-<run>.json           adversarial checks (see below)
  archive-sx7/                seven-station variant runs and their figures (see Data)
```

## Reading the figure

- **Grey crosses** are Lassie's fixed 2 km grid nodes. They sit at odd-kilometre
  offsets from the reference point (`xmin + i * dx`), so the lattice is
  visibly off-centre; Lassie detections can only ever fall on these crosses,
  which is why its epicentres look quantised.
- **Marker size** grows with the detector's own strength (image function for
  Lassie, semblance for Qseek), scaled within each catalog, so the threshold's
  effect can be judged by eye without a re-run.
- **Faded markers** are detections whose own evidence does not support the
  location: Qseek ones with fewer than 6 picks or no azimuthal spread, Lassie
  ones within one node of the search-volume boundary (Qseek rejects those
  itself; Lassie has no such rejection, so it is applied when reading its
  catalog, never inside Lassie). Legends spell out the criterion; it is a
  support flag, not an uncertainty.
- **Lassie positions in the two overview panels are jittered** uniformly
  within their grid cell (±1 km on a 2 km grid; surface nodes downward only),
  with one fixed offset per detection so a detection sits in the same place
  in every view. Lassie cannot resolve position inside a cell, so this loses
  nothing, and it stops a cluster of detections on one node from collapsing
  into a single marker. Qseek's positions are its own interpolated estimates
  and are never moved. The close-ups draw exact positions (with a "×n" where
  Lassie detections share a node) and every statistic uses exact positions;
  `lassie3 plot --no-jitter` draws the raw nodes everywhere.
- **Gold stars** are the ISC reference events, numbered as in the review's
  per-event table; in the two close-ups a line joins each star to the
  detection credited with recovering it. The depth close-up spans the full
  0-20 km; matches outside a close-up's ±4 km window are listed on the panel
  with their distance and depth rather than silently clipped.
- **Agreement bars** are medians with every event drawn as a dot, because a
  median hides Lassie's bimodal depth error; the dashed line is what a
  constant 9.9 km guess would score.

`task plot:strict` renders the same runs showing only detections above
`median + 40 MAD` for Lassie (84.2 here) and semblance 0.5 for Qseek; override
with `LASSIE_MIN_MAD=… QSEEK_MIN_SEMBLANCE=…`. On 2024-03-20 that keeps 25 of
Lassie's 47 (11 of the 13 reference events survive) and 94 of Qseek's 1105
(85 pick-supported, all 13 survive). The run threshold is left at 30 MAD
because the two reference events Lassie loses at 40 MAD sit in the middle of
its catalog by strength: Lassie's score does not separate them from the noise,
so any stricter level trades recall for purity.

`lassie3 summary` prints catalog counts, depth distributions, the overlap
between the two catalogs, and validation against an independent reference.

## Validation

Detection counts alone say nothing — they follow from where each threshold was
placed. `task catalog` fetches an independent catalog from **ISC**, and
`task summary` scores both detectors against it.

ISC is used specifically because it redistributes **WBNET**, the dense local
network over this swarm zone, which reports down to about ML 0.0. This matters:
on 2024-03-20, EMSC/BGR lists **zero** events in the region while ISC lists
**thirteen** (11 WBNET with ML 0.13–0.69, 2 PRU without magnitude). Validating
against a regional catalog with an
ML ~1.5 completeness limit would have wrongly suggested every detection that day
was a false alarm.

The reference is **11 WBNET events with ML 0.13-0.69 and 2 PRU events without
magnitude**; its depths span 9.4-10.2 km and carry no uncertainty, so a
detector's depth error against it measures *agreement with a narrow catalog*,
not depth resolution: guessing 9.9 km for every event already scores a 0.10 km
median error.

Results for 2024-03-20, against those 13 events:

| | recovered ±5 s (±1 s) | origin time med / mean | epicentre med / mean | depth med / mean | detections |
|---|---|---|---|---|---|
| Lassie v1 | **13/13** (7/13) | 0.63 / 1.59 s | 2.05 / 3.35 km | 2.60 / 5.44 km | 47 |
| Qseek | **13/13** (13/13) | **0.04 / 0.18 s** | **0.96 / 1.21 km** | **0.09 / 0.22 km** | 1105 (109 pick-supported) |

Both notice every catalogued earthquake. Qseek's locations agree with WBNET
about an order of magnitude better; Lassie's medians flatter it, because its
errors are bimodal (see the review below), which is why the means are given
too. Both detectors sit systematically **north-east of the reference**
(mean offset against the 11 WBNET events: Lassie N +1.2 / E +2.3 km, Qseek
N +0.8 / E +0.5 km). The common displacement is consistent with shared
modelling or geometry effects — the same waveforms, one-sided array, 1D model
and travel-time engine feed both — but these runs do not isolate its cause.

The counts are not comparable as-is: Lassie was run at a strict operating point
(47/day) and Qseek at its permissive adaptive default (1105/day, of which 109
have at least 6 picks and non-zero azimuthal coverage). All 13 reference
events fall inside that pick-supported subset. Neither flag is an uncertainty
estimate — matched Qseek events still have ~240 degree azimuthal gaps on this
array — so the surplus detections in both catalogs are **unverified
candidates**, not a validated catalog.

Receiver elevations matter and were initially wrong on the Lassie side: its
`stations_path` loader sets receiver depth from the sensor's burial depth
(0-34 m), not its elevation, so every station sat at sea level while Qseek
traced from the true 510-870 m. Explicit receivers (a stock config feature)
fix that; the fix changed Lassie's catalog from 49 to 47 detections, left every
recall figure unchanged, and moved its shallow location mode from the 2 km node
to the surface.

### Adversarial review

`task review` re-examines those numbers for the ways they could mislead, and
writes `results/review-<run>.json`. Findings for 2024-03-20:

- **Recovery is not chance.** With a ±5 s tolerance, Lassie's 47 detections
  cover 0.5 % of the day (0.07 recoveries expected at random); Qseek's full
  1105 cover 12.0 % (1.56 expected), its pick-supported 109 cover 1.0 %
  (0.14). Sliding the whole reference pattern around the day 10,000 times
  never recovers more than 3 (Lassie), 8 (Qseek, all) or 3 (Qseek,
  pick-supported) events; the real pattern recovers 13. Every Qseek match
  carries 12–14 picks and a semblance of 1.10–1.58, near the top of its
  catalog (worst rank 27 of 1105).
- **"13/13" is temporal association.** The matcher takes the nearest
  detection in time within ±5 s, with no spatial condition. With a 5 km
  epicentral gate Lassie recovers 12/13 (the 06:24:55 match is 14.7 km away
  at 18 km depth), Qseek 13/13.
- **Lassie's depth median hides a bimodal failure.** Per event, its depth error
  is either about +2 km (the 12 km node for events at ~9.9 km: six of the
  thirteen, one at 10 km) or about **−10 km** (five placed at the surface,
  each with an origin time ~2.5 s late), plus one at 18 km. Before the
  receiver-elevation fix the shallow mode sat at the 2 km node instead; the
  trade-off itself is the depth/origin-time ambiguity of a one-sided array,
  which a 2 km grid cannot resolve. Qseek's depth errors are all within
  ±0.8 km.
- **Lassie's 13/13 only holds at a wide tolerance.** Recall against the 13
  events at ±0.5 / 1 / 2 / 3 / 5 s is 5 / 7 / 7 / 11 / 13 for Lassie and
  11 / 13 / 13 / 13 / 13 for Qseek. Lassie *notices* every event; it times
  barely half of them to within a second, because the late origin time is
  the other half of its depth error — and the sweep shows that this does not
  change with the threshold.
- **One Lassie recovery is dubious.** The 06:24:55 event is matched by a
  detection 14.7 km away at 18 km depth, 3.8 s early, ranked 27th of 47 by
  strength. Honest reading: 12/13 recovered plus one ambiguous association
  inside the tolerance.
- **Lassie's catalog is mostly boundary.** 34 of 47 detections lie within one
  node of the surface, 3 at the lateral edge, 2 at the bottom; 10 are clean
  interior. Qseek rejects lateral and bottom boundary hits itself and leaves
  only the surface (485 of 1105; 3 of the 109 pick-supported).
- **The true events are not simply the strongest.** 8 of Lassie's and 8 of
  Qseek's top-13 detections are catalogued events (Lassie's worst true event
  ranks 31st of 47, Qseek's 27th of 1105); the rest of Qseek's top ranks sit
  in the same 10 km cluster and may be uncatalogued micro-events, but nothing
  here can confirm that.
- **Possible duplicate reference entry.** ISC lists a PRU solution at 04:04:40
  and a WBNET one at 04:05:11, 31.5 s apart. Both detectors report separate,
  pick-supported detections for each, so they are probably two events; recall
  without the pair is 12/12 either way.
- **The surplus detections cluster in time with the swarm.** 66 % of Lassie's
  detections and 58 % of Qseek's pick-supported ones fall within ±30 min of
  a catalogued event, against 28.5 % of the day covered by such windows
  (enrichment x2.3 and x2.0). Swarm activity clusters in time and random noise
  does not, so this is *consistent with* much of the clutter near the source
  zone being sub-catalog micro-seismicity — but it is not proof: coda,
  aftershock-triggered noise and shared time-of-day noise conditions cluster
  too, and nothing here inspects waveforms. Qseek's full 1105 show no such
  enrichment (x1.15): outside the pick-supported subset the catalog is
  dominated by one-or-two-pick detections with no temporal association,
  including a persistent blob at 12.24°E / 50.27°N, 0–3 km deep, south-west of
  ROHR and WERN.
- **Lassie mostly avoids its 10 km node.** The grid has a node ~0.6 km from
  the reference cluster at exactly 10 km depth; only one recovered event
  landed on it, the rest went to 12 km or to the surface. A 2 km grid cannot
  by itself explain choosing 0 km over an available 10 km node. The nominal
  30:1 P:S weighting is a suspect, but nominal weights are not effective
  contributions (P's smoothing and normalisation rescale it), so the cause
  stays open until the P and S characteristic functions of the failed events
  are inspected directly.
- **Qseek's epicentres are biased, not scattered.** Its matches form a tight
  cluster displaced ~1 km north-east of the WBNET locations — the 0.96 km
  median error is systematic, most likely velocity-model and geometry bias.
- **The passband is not perfectly shared.** Qseek band-passes waveforms at
  2–30 Hz before PhaseNet; Lassie's P onset uses 2–30 Hz but its S wave-packet
  detector 2–15 Hz. The velocity model, stations, channels, day and search
  volume are verified identical from the persisted run configs.

## Data

`download.sh` fetches StationXML and miniSEED for 2024-03-16..30.

**Nine stations are open data:** the seven SX stations plus two found later,
**CZ.NKC** (the Nový Kostel site, whose broadband STS-2 is published openly in
the Czech Regional Seismic Network on GEOFON while the same site's WB.NKC
stream stays restricted) and **GQ.LNDWU** (BGR's Landwüst borehole site).
LNDWU also records a 200 Hz rotational sensor (`HJ1-3`); only its `HH?`
channels are used, which the `HH?` channel selection in `RunSettings.channels`
enforces for every station. The WB (7), 6A (4) and 1D (1) stations remain
restricted at GEOFON and need an EIDA token, without which the download
silently yields metadata but no waveforms:

```sh
EIDA_TOKEN=~/.eidatoken task data:download
```

Geometry matters for the result, not just the station count. The reference
events of 2024-03-20 lie at 50.37°N 12.49°E, ~10 km deep, near Klingenthal on
the German–Czech border, about 15 km north of Nový Kostel. The SX stations see
that source from the north-west quadrant only, at 5–25 km (TANN 5 km N, MULD
7 km NW, GUNZ 11 km W, WERN 12 km SW; WERD, ROHR and TRIB further out), so the
seven-station array has no southern or eastern coverage and depth is poorly
constrained. NKC (16 km S) and LNDWU (17 km SW) add the southern azimuths; the
east stays empty, and no open station sits above the source.

`task data:inventory` reports which stations actually have data for a given day.

2024-03-20 is the default because all nine open stations have uninterrupted
100 Hz `HH?` coverage that day.

**The results in this README were computed with the seven SX stations**, before
CZ.NKC and GQ.LNDWU were downloaded. Those runs are kept under the tag `sx7`
(`results/lassie/eger-2024-03-20-sx7.turd`, `results/qseek/eger-2024-03-20-sx7`,
`results/comparison-eger-2024-03-20-sx7.png`; the tagged variants and sweeps
under `results/archive-sx7/`), and `lassie3 --tag sx7 summary|review|plot`
reads them. Every untagged run from now on includes all nine stations.

## Velocity model

`RunSettings.WEST_BOHEMIA_ND` is an approximate layered model for the West
Bohemia / Vogtland region, not a verbatim published one. Both detectors read the
same file, so substituting your own keeps the comparison fair:

```sh
uv run lassie3 --velocity-model my_model.nd run-qseek
```

## Thresholds

The two detectors do not share a threshold scale, so they are calibrated
separately against their own noise floor:

- Qseek uses its built-in `MAD` semblance threshold.
- Lassie has only an absolute level, whose scale depends on station count and
  weights. `--lassie-threshold auto` (the default) therefore stacks twice: once
  with detection disabled to measure the image function, then again at
  `median + 30 x MAD`. Tune with `--lassie-mad-multiple`, or pass an absolute
  number to skip calibration.

Detection counts are consequently **not** a like-for-like sensitivity
comparison; they depend on where each threshold was placed, and they are very
sensitive to it. The Lassie sweep below comes from **actual runs** at each
level (tagged `mad20`, `mad40`, `mad60` next to the baseline), not from
counting peaks in the saved image function, which over-counts across
overlapping windows. The Qseek rows are post-filters of the one saved run, and
a post-filter is not a fresh run: a different threshold would also change
Qseek's peak selection and octree refinement.

Sweep on 2024-03-20 (`task sweep`; "supported" = interior of the search
volume for Lassie, ≥6 picks with azimuthal coverage for Qseek):

| detector | level | detections | of which supported | recovered ±5 s | recovered ±1 s |
|---|---|---:|---:|---:|---:|
| lassie (baseline) | 30 MAD = 71.4 | 47 | 10 | 13/13 | 7/13 |
| lassie (mad20) | 20 MAD = 58.6 | 140 | 24 | 13/13 | 7/13 |
| lassie (mad40) | 40 MAD = 84.2 | 25 | 8 | 11/13 | 7/13 |
| lassie (mad60) | 60 MAD = 109.8 | 5 | 5 | 5/13 | 5/13 |
| qseek (post-filter) | semblance ≥ 0 | 1105 | 109 | 13/13 | 13/13 |
| qseek (post-filter) | semblance ≥ 0.2 | 317 | 109 | 13/13 | 13/13 |
| qseek (post-filter) | semblance ≥ 0.3 | 165 | 106 | 13/13 | 13/13 |
| qseek (post-filter) | semblance ≥ 0.5 | 94 | 85 | 13/13 | 13/13 |
| qseek (post-filter) | semblance ≥ 0.7 | 56 | 56 | 13/13 | 13/13 |
| qseek (post-filter) | semblance ≥ 1 | 34 | 34 | 13/13 | 13/13 |
| qseek (post-filter) | semblance ≥ 1.1 | 27 | 27 | 13/13 | 13/13 |

Below roughly 20 MAD, Lassie's catalog is visibly noise: detections spread
evenly across the whole search grid instead of clustering near the array. Qseek
degrades more gracefully because its per-event pick count and azimuthal coverage
let weak locations be filtered afterwards; `lassie3 summary` reports both the
raw and the pick-supported counts for that reason.

## Variants

Tagged runs live next to the baseline (`results/lassie/eger-<day>-<tag>.turd`)
and are compared against the baseline Qseek run. The runs below were made
with the seven SX stations and are archived under `results/archive-sx7/`.

**1 km Lassie grid** (`task lassie:grid1km`, ~35 min): the obvious objection to
Lassie's location errors is its 2 km node spacing — nearest-node assignment
alone would cost ~0.9 km horizontally and ~0.2 km in depth. Halving the
spacing (54,621 nodes, chunked stacking) changes essentially nothing:

| | 2 km baseline | 1 km variant |
|---|---|---|
| detections at 30 MAD of its own image function | 47 (threshold 71.4) | 42 (threshold 72.4) |
| recovered at ±0.5 / 1 / 3 / 5 s | 5 / 7 / 11 / 13 | 4 / 7 / 11 / 13 |
| median origin time / epicentre / depth error | 0.63 s / 2.05 km / 2.60 km | 0.63 s / 2.05 km / 2.60 km |
| mean origin time / epicentre / depth error | 1.59 s / 3.35 km / 5.44 km | 1.60 s / 3.38 km / 5.21 km |
| matched depth nodes | 0 km ×5, 10 ×1, 12 ×6, 18 ×1 | 1 km ×5, 11 ×1, 12 ×6, 19 ×1 |
| within one node of the surface | 34 of 47 | 29 of 42 |

The same five events still go shallow with a ~2.5 s late origin time and the
same seven go to 11–12 km, on a grid that now has nodes at 9 and 10 km. The
epicentre error does not move either. Lassie's errors on this array are a
property of its characteristic functions and the one-sided geometry, not of
its grid. (The finer grid also halves Lassie's travel-time interpolation
tolerance, which is tied to the node spacing, so this run is not *only* a
resolution change; it is nonetheless the closest stock Lassie gets to Qseek's
312 m octree.)

## Second opinion

An independent adversarial review by a different model family (OpenAI
`gpt-6-astra`, via agentmux, read-only, transcript under
`~/Library/Application Support/com.romnn.agentmux/runs/20260916T201259-fb2fbd17/`)
was run on the v1 figure, the raw catalogs and this code. It confirmed the
stock-code claim (476/476 installed files match their wheel records), the
Lassie coordinate conversion, the day/station/channel equivalence, and that
chance cannot produce the recoveries. It found, and this revision fixes or
documents: the receiver-elevation mismatch; two mis-stated numbers in the
earlier review (five, not six, shallow Lassie matches; the boundary match ranks
28th, not 45th); the "everything else shared" over-claim; that 13/13 is
temporal association rather than 13 credible hypocentres; an over-counted
chance estimate; the reference composition; that the depth agreement must be
read against a constant-depth null; that a peak-count sweep contradicted the
actual run; and that "well constrained" and "mostly real" were unsupported
labels. Its verdict on v1 — useful exploratory evidence favouring Qseek, not
yet a controlled baseline — stands until the finer-grid and sweep runs below
are in.

## Stock code

Lassie, Qseek and pyrocko run unmodified. Every installed `.py` file was
compared byte-for-byte against a clean checkout of the locked commit (Lassie
14/14, Qseek 115/115, pyrocko 347/347 identical), nothing in `src/lassie3`
monkeypatches them, and the C extensions are compiled from the upstream
sources. The only changes are build flags (below) and, from `doctor`, which
copy of `libomp.dylib` the compiled extensions load — linkage, not code. The
boundary de-emphasis for Lassie and the pick filter for Qseek are applied when
reading the catalogs, after both detectors have finished.

## Upstream workarounds

Three upstream problems are worked around in this project rather than patched
in the dependencies:

1. **OpenMP conflict (crashes Qseek).** Qseek's C extensions link Homebrew's
   `libomp`; PyTorch loads its own. Two OpenMP runtimes in one process segfault
   as soon as PhaseNet touches a tensor. `lassie3 doctor` repoints Qseek's
   extensions at PyTorch's copy, and `run-qseek` applies it automatically.
   **`uv sync` reverts it**, which is why `task sync` runs `doctor` after.
2. **`Search.init_rundir` calls `Path.move`**, which pathlib does not provide,
   so `--force` cannot use Qseek's own no-backup path. `run_qseek.run` removes
   the run directory itself instead.
3. **pyrocko builds its stacker single-threaded on macOS.** Its `setup.py`
   enables OpenMP only if `cc -fopenmp` compiles, which Apple clang refuses, so
   `parstack` — most of Lassie's runtime — ignored `nparallel` entirely. The
   `pyrocko` entry in `[tool.uv.extra-build-variables]` passes the Apple-clang
   spelling (`-Xpreprocessor -fopenmp`, `-lomp`) instead. Two traps inside that
   fix: a `CFLAGS` variable *replaces* the interpreter's flags on macOS, silently
   dropping `-O3` (8x slower per thread) unless they are repeated; and `uv`
   reuses a previously built wheel after the variables change, so a real
   rebuild needs `uv sync --reinstall-package pyrocko --no-cache`. `doctor`
   then repoints all twelve pyrocko extensions at PyTorch's `libomp`, the same
   way it handles Qseek's.

Three further pins live in `pyproject.toml` with comments: `setuptools<81` for
Lassie's build, `matplotlib<3.11` for `pyrocko.gui`, and the Homebrew libomp
include/lib paths Qseek needs at build time (Apple-Silicon-specific).
