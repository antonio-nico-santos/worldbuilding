# Tappa 4 — Hydrology (flow routing, stream network, basins)

Status: **closed**. Depression-filled DEM, D8 flow direction, contributing
area, precipitation-weighted discharge proxy, a stream network by
contributing-area threshold, a lake/depression mask, and major drainage
basins — all computed at the DEM's native 30 m resolution (5,334x4,334) by
`run_tappa4_hydrology.py` from the v3 eroded DEM and Tappa 2's annual
precipitation field.

This document is the decision record for the stage, per the project's
per-stage workflow. As with the three prior documents it is written to work
either as a log of what was decided, or as a reusable recipe.

## 0. Two decisions locked in before any code, and why they had to come first

Tappa 0 scoped this stage as "`r.watershed` + `r.stream.extract` on the
DEM" — but GRASS is not available in this sandbox (same constraint as Tappa
1-3: no PyPI/apt egress for `grass-session`, `rasterio`, `osgeo`, or
`numba`; confirmed directly this session, not assumed). Everything here is a
from-scratch reimplementation on top of numpy/scipy, same posture as
`noise.py`/`erosion.py`. Two forks had to be settled before writing anything,
both decided with the user up front (same pattern as Tappa 2's wind
direction — a "what world is this" choice, not an implementation detail):

1. **Resolution: native 30 m DEM, not the 120 m climate grid.** Flow
   routing determines the drainage network's own topology (branching
   density, tributary count) directly from grid resolution — unlike
   temperature/precipitation, which Tappa 2 S2 established are smooth
   enough at 120 m to lose nothing. A drainage network computed at 120 m
   would be visibly blockier and thinner in tributaries. The cost: 23.1M
   cells instead of 1.4M, and no GRASS/numba to lean on for the
   inherently-sequential flow-accumulation traversal. Timed and validated
   at increasing crop sizes before committing (300x300 through 4M cells,
   ~7.4 us/cell and scaling near-linearly) — full run: **442 s**.
2. **Precipitation-weighted flow accumulation, not plain cell-count.**
   Couples this stage to Tappa 2's climate stack the same way Tappa 2 itself
   derived `Cw`/`Hw` from temperature rather than hand-tuning them (S4c) —
   consistent cross-stage physical coupling rather than a topology-only
   textbook run. Concretely: windward cells contribute their actual modelled
   annual precipitation (mm) to the flow-accumulation sum instead of `1`, so
   a river's size reflects the catchment's real modelled rainfall, not just
   its cell count. Both a plain-area accumulation (for stream-network
   topology) and a precipitation-weighted one (for a discharge proxy) are
   computed — see S2 and S4.

## 1. Algorithm: Priority-Flood+Epsilon, filling and D8 direction in one pass

r.watershed's job is really two things GRASS keeps as separate internal
steps (depression-fill the DEM, then compute D8 direction and accumulate)
which the **Priority-Flood+Epsilon** algorithm (Barnes, Lehman & Mulla 2014,
*Computers & Geosciences* 62) computes together, in one pass:

A min-heap is seeded with every drainage sink — every ocean cell
(`dem <= 0`) plus every domain-boundary cell (this DEM's land reaches all
four edges, per `02_tappa2_climate.md` S4f, so the boundary itself is a
legitimate outlet, not just the coastline). Repeatedly pop the lowest open
cell `c` and flood into its unclosed neighbors `n`: if `dem[n] <= dem[c]`,
raise it to `dem[c] + epsilon` before closing it. Two things fall out for
free:

- **Depression filling**, with the classic "then a separate flat-resolution
  pass" problem solved by construction: the `+ epsilon` guarantees every
  closed cell is *strictly* higher than the cell that flooded into it, so
  there are no dead-flat plateaus left over for D8 to choke on, and no
  second pass is needed (`epsilon = 1e-4` m — negligible against this DEM's
  ~4,540 m relief, verified: filled elevation never drops below the
  original anywhere).
- **D8 flow direction**, as a literal byproduct: whichever cell was popped
  when `n` got closed IS `n`'s receiver (the only way `n` could have been
  reached was by flooding in from there). No separate steepest-descent pass
  over the filled DEM.

The pop order is also a ready-made reverse-topological order of the flow
DAG (every cell pops strictly after its receiver), which is what makes flow
accumulation a single **O(n) reverse pass** (`accumulate_flow` in
`src/hydrology/flow.py`) rather than a sort-and-iterate — sum each cell's
current weight into its receiver's, processing from most-recently-popped
(most upstream) to least. No priority queue, no convergence loop.

Correctness checked on synthetic cases before touching the real DEM: a
single pit surrounded by a rim (filled correctly, strictly-decreasing path
to the sea confirmed by walking receivers), a dead-flat plateau (every cell
still reaches a sink in a bounded number of hops — no cycles, the epsilon
tie-break working as intended), and a mass-conservation check (total
weight in == total weight summed at every outlet, exact to float
precision). See `scripts/test_flow_small.py`.

**Performance**, since the one part of this algorithm that cannot be
vectorized is the heap traversal itself: benchmarked at 300x300 through
2000x2000 crops (90K to 4M cells) before the full run, ~6-7.5 us/cell for
fill+direction and ~0.6 us/cell for the accumulation reverse pass, both
scaling close to linearly (the `log n` factor in the heap barely moves
across this range). Extrapolated the full 23.1M-cell domain at ~3 minutes;
actual run was 442 s total (fill+direction 335 s, both accumulation passes
combined ~55 s, the rest basin labeling and I/O) — same order as Tappa 1's
erosion pass (1,324 s), not out of line for a "final run" in this project.

## 2. Contributing area and the (attempted, rejected) slope-dependent threshold

Plain D8 contributing area (`accumulate_flow` with weight = 1 on land, 0 on
ocean) gives area in cells, converted to km² by the fixed 30 m cell size.
Real channel-initiation literature (Montgomery & Dietrich, *Where do
channels begin?*, *Nature* 336, 1988; *Source areas, drainage density, and
channel initiation*, *Water Resources Research* 25(8), 1989) establishes
that the contributing area needed to initiate a channel is **slope-
dependent** — steep terrain channelizes with much less drainage area than
gentle terrain, because less water is needed to overcome a smaller
threshold for erosion. A combined criterion (area x slope², in the spirit
of Montgomery & Foufoula-Georgiou 1993) was tried first for exactly this
reason.

**It was tried and rejected.** On this DEM's steepest terrain (95th/99th
percentile slope 0.60/0.94 m/m — genuinely alpine), `slope²` alone is large
enough that even single-digit-cell contributing areas clear any threshold
worth using, so the criterion doesn't pick out a thin channel network at
all on steep ground — it floods whole hillslopes solid blue
(`docs/decisions/assets/04_tappa4_hydrology/
04_threshold_sweep_rejected_areaslope.png`, all four candidate thresholds).
This is a real, documented failure mode of the naive area x slope^n
formulation without the local calibration Montgomery & Dietrich's own
compiled dataset uses (their published thresholds are fit per-study-area
against a real channel-head survey, not a fixed constant applied
worldwide) — reproducing that calibration would need a per-terrain-type
fitting exercise this sandbox has no channel-head ground truth to check
against, unlike, say, Tappa 3's NIWA station data.

**Adopted instead: plain contributing-area threshold**, exactly
`r.stream.extract`'s simplest and most common real-world mode. Calibrated
by a visual/cartographic sweep (`04_method_crop_test.png`,
0.1/0.5/2.0 km² panel) plus a drainage-density cross-check against Horton's
long-standing classification (Horton 1945; low <1, moderate 1-3, high 3-5,
very high 5-8 km/km²) — the general literature range for channel-initiation
area itself spans roughly 0.001-1+ km² across climates and slopes (same
Montgomery & Dietrich body of work), consistent with this project's honesty
standard elsewhere: **no Southern-Alps-specific channel-initiation or
drainage-density figure was found** in the sources checked this session
(same kind of gap as Tappa 3 S1's `sigma_day_c`) — the number below is
picked from within a defensible literature range and a plausible drainage-
density outcome, not independently verified against a named station.

```
stream_threshold_km2 = 0.3
```

Result: 377,281 stream cells, ~11,318 km of channel (cell-count x 30 m, a
rough proxy — undercounts true length on diagonal-heavy reaches by up to
~30%), drainage density **1.14 km/km²** — solidly Horton's "moderate" band,
defensible for a wet, steep, young landscape without overclaiming a
specific real-world match.

## 3. Precipitation-weighted accumulation and the discharge proxy

`src/hydrology/weighting.upsample_precip_to_dem` bilinear-resamples Tappa
2's 120 m annual precipitation field onto the 30 m DEM grid (not a
block-repeat: Tappa 2 S2 already established precipitation has no signal
below ~1 km, so interpolating is free of any real information loss, and
avoids stamping visible 120 m block edges onto the accumulated field). Land
cells contribute their local annual precipitation (mm) as accumulation
weight; ocean cells contribute zero (rain falling directly on the sea does
not feed a river's discharge).

The resulting `accum_precip_mm` is converted to a **discharge proxy**:

```
discharge_proxy_m3s = accum_precip_mm * (30 m)^2 / 1000 (mm->m) / seconds_per_year
```

Stated plainly as an **upper bound, not a discharge model**: this treats
100% of a catchment's annual precipitation as becoming streamflow
instantaneously, with no infiltration, evapotranspiration, or channel
routing/lag — the same category of honest scope-cut as Tappa 3 S6's
mass-balance model (no ice flow, no refreezing). A real discharge estimate
would need a runoff coefficient this pipeline has no independent way to
derive (no soil/vegetation model). Max value: **103.1 m³/s**, for the
751 km² catchment draining the Spine's wet flank — for comparison, a rough
specific-discharge rule of thumb for wet temperate catchments
(0.03-0.15 m³/s per km²) puts 103.1/751 = 0.137 m³/s/km² at the upper end
of that range but not outside it, consistent with treating this number as
an upper bound rather than a calibrated estimate.

**Internal consistency check** (not an independent validation, but a sanity
check that the weighting propagated through the drainage network the way
it should): windward cells' mean discharge proxy is 2.21 m³/s against the
leeward side's 0.38 m³/s — a **5.9:1** ratio, closely tracking Tappa 2 S5's
established 5.6:1 windward:leeward *precipitation* ratio. The two maxima are
closer (103.1 vs 44.5 m³/s, 2.3:1) because the largest rivers' catchments
straddle both precipitation zones on their way from the Spine to the coast,
diluting the ratio at the outlet — an expected consequence of catchment
geometry, not a modelling inconsistency.

## 4. Lakes / depressions

Any cell the priority-flood raised more than `lake_fill_threshold_m = 2.0`
above its original elevation is flagged as a lake/depression cell — this is
not a bug-fix, it is exactly what filling a genuine closed basin is
supposed to do: the interior becomes a flat "virtual lake" surface that
drains out through a single spillway, which is real behavior for an actual
endorheic basin (confirmed by checking the underlying elevation: the
clearest example, visible as a fan/spoke pattern converging on one outlet
in early crop tests, sits in genuinely rough terrain — 131 m local std
dev, ~790 m local relief — not a flat plateau artifact). **176.5 km²**
(1.8% of land) came out flagged this way, the largest a plausible alpine
lake basin near the domain center (`04_zoom_streams_lakes.png`). Worth a
visual QA pass in QGIS before treating every flagged cell as a "real" named
lake — 2.0 m is a floor against tiny numerical fill noise, not a curated
lake catalogue, and some flagged area is likely valley-floor filling too
shallow/small to read as a distinct feature at map scale.

## 5. Basins

`label_basins` (one forward pass over the same pop order, mirroring
`accumulate_flow`'s reverse pass) tags every cell with the flat index of
the outlet it ultimately drains to — exact by construction, but with
~96,751 distinct outlets, because every one of the thousands of individual
ocean-boundary cells is technically its own trivial sink, most of them
capturing a handful of land cells with no real catchment. Cartographically
meaningless at that count, so basins are aggregated by contributing land
area: the **54 basins >= 20 km²** are kept and relabeled 1-54 by descending
size (6,818 of 9,913 km² total land, 69%); everything smaller collapses to
`0` ("minor/direct coastal drainage") rather than being individually
labeled. This is a display simplification, not a hydrological one — no
data is discarded, `contributing_area_km2` and `stream_mask` still carry
the full-resolution truth for anyone who wants to re-derive a different
basin count.

## 6. Validation against South Island, NZ

| metric | this world | reference |
|---|---|---|
| land area | 9,912.5 km² | matches Tappa 2's independently-computed 9,913 km² (`02_tappa2_climate.md` S3b) — internal cross-check, not a new measurement |
| drainage density | 1.14 km/km² | Horton (1945) "moderate" band (1-3 km/km²); no Southern-Alps-specific figure found this session |
| max discharge proxy / catchment area | 0.137 m³/s per km² | within the 0.03-0.15 m³/s/km² rule-of-thumb range for wet temperate catchments (upper end, consistent with an unrouted upper-bound proxy) |
| windward:leeward discharge ratio (mean) | 5.9:1 | tracks Tappa 2's validated 5.6:1 precipitation ratio (S5) — internal consistency, not independent |

Honest reading: this stage leans more on **internal consistency** (does the
hydrology agree with what Tappa 2/3 already established and validated)
than on fresh external validation — unlike Tappa 2/3, no NZ station-level
streamflow or drainage-density dataset was located and checked this
session. The land-area match and the windward:leeward ratio match are
genuine, meaningful checks (both are independent computations from the same
DEM/climate stack arriving at the same answer two different ways); the
drainage-density and specific-discharge comparisons are against general
literature ranges, not a named reference site the way Tappa 2/3's NIWA
station comparisons were.

## 7. Locked-in parameters

```
resolution_m               = 30          # native DEM grid, not 120 m climate grid
priority_flood_epsilon_m   = 1e-4
stream_threshold_km2       = 0.3          # visual/drainage-density calibrated, S2
lake_fill_threshold_m      = 2.0          # floor against fill noise, not a lake catalogue
major_basin_min_km2        = 20.0         # display aggregation only, S5
```

## 8. Outputs

`run_tappa4_hydrology.py` reads the v3 eroded DEM and Tappa 2's monthly
precipitation, writes to `data/processed/hydrology/` (gitignored,
regenerates in ~440 s):

| file | contents |
|---|---|
| `filled_dem_30m.npy` | float32, depression-filled DEM |
| `flow_direction_code.npy` / `.bin+.hdr+.prj` | int16, D8 direction, ESRI convention (E1 SE2 S4 SW8 W16 NW32 N64 NE128, 0=outlet) |
| `contributing_area_km2.npy` / `.bin+.hdr+.prj` | float32, plain D8 contributing area |
| `discharge_proxy_m3s.npy` / `.bin+.hdr+.prj` | float32, precipitation-weighted discharge upper-bound proxy (S3) |
| `stream_mask.npy` / `.bin+.hdr+.prj` | int16 (0/1), area >= 0.3 km² |
| `lake_mask.npy` / `.bin+.hdr+.prj` | int16 (0/1), fill raise > 2 m |
| `basin_labeled.npy` / `.bin+.hdr+.prj` | int16, 1-54 major basins, 0 = minor/coastal |
| `tappa4_hydrology_meta.json` | parameters, summary stats, windward/leeward asymmetry check |

Same CRS caveat as Tappa 1-3: `.hdr`'s `map info` carries the affine
georeferencing, a `.prj` sidecar carries the PROJ4 string, and QGIS needs
"Fictional World LCC" assigned manually if it loads with an unknown CRS.

Summary of the locked run: land area 9,912.5 km², max contributing area
751.1 km², max discharge proxy 103.1 m³/s, 377,281 stream cells (~11,318 km,
drainage density 1.14 km/km²), 196,094 lake cells (176.5 km²), 96,751 total
basins / 54 major (>= 20 km²).

## 9. Consequences for the site's visual plan

`00_pre_project_planning.md` lists rivers as an "always-on context layer"
in the interactive map. `stream_mask` at 30 m is almost certainly too dense
to render directly as a web vector layer at domain-overview zoom (377K
raster cells) — this needs a raster-to-vector + simplification pass
(e.g. skeletonize then Douglas-Peucker) before it is a `InteractiveMap.astro`
layer, not attempted here since it is a cartography/export step rather than
a hydrology one. `discharge_proxy_m3s` is a natural candidate for stream
line-width scaling once vectorized (bigger rivers draw thicker), and
`lake_mask` for a simple polygon fill layer after the same vectorization.
`basin_labeled` (54 major basins) could work as its own qualitative-palette
layer if Tappa 5/6 want a "which catchment" context, though nothing in this
project's plan currently calls for it.

## 10. Open follow-ups (not done in this stage, deliberately left open)

- ~~No Strahler stream order~~ — **done, see S11.** Built as a side effect
  of reach segmentation rather than the CSR contributor-inversion originally
  anticipated here.
- ~~Stream-to-vector export not attempted~~ — **done, see S11.**
- **`lake_mask` is unreviewed** (S4) — worth a QGIS pass to confirm the
  176.5 km² flagged area reads as sensible discrete lakes rather than
  diffuse valley-floor noise before treating any of it as narratively
  "real" (naming a lake, etc.).
- **No area x slope channel-initiation criterion** (S2) — tried, rejected
  for this terrain's slope range, documented rather than silently dropped.
  A properly locally-calibrated version (fitting the constant against this
  DEM's own slope distribution rather than importing one) is possible
  future work if the plain-area network turns out to need finer headwater
  detail for Tappa 7's zoom stage.
- **No streamflow/drainage-density validation against a named NZ site**
  (S6) — this stage leaned on internal consistency with Tappa 2/3 instead;
  a real station comparison (e.g. a specific West Coast gauged catchment of
  comparable area) would strengthen this if one turns up later.
- `config/parameters.yml`'s `hydrology:` section (added this stage) should
  be checked against this document if either changes independently.

## 11. Addendum — reach segmentation, Strahler order, and a smoothed vector export

Added after this stage was first closed, prompted by a direct visual
question: the raster `stream_mask`, viewed at zoom (`04_zoom_streams_lakes.png`),
makes visibly angular turns rather than natural curves. Root cause,
confirmed by tracing an actual reach through `flow_direction_code` rather
than assumed: D8 restricts every flow step to one of 8 directions (45°
increments), so a true bearing between two of those (e.g. ~20°) can only be
represented as an alternating zigzag between the nearest two — one traced
reach was 49% direction-changes, bouncing NE/N/NW rather than holding a
diagonal. This is cosmetic (drainage topology, accumulation, and discharge
are all unaffected by how jagged the geometric bearing looks) but real, and
purely a rendering/geometry problem — not something fixing the flow-routing
algorithm (e.g. switching to D-infinity) would actually solve, since a
cartographic line still needs one committed path per reach regardless.
**No re-run of the priority-flood was needed** — `flow_direction_code`
already encodes the full receiver graph losslessly, so this is a pure
export step over already-computed rasters (`src/hydrology/vectorize.py`,
driven by `scripts/export_streams_vector.py`, deliberately NOT folded into
`run_tappa4_hydrology.py` — that script stays focused on producing the
hydrologically-correct rasters, matching how plotting/export lives in
separate scripts for Tappa 1-3).

**Reach segmentation.** `stream_mask` alone has no notion of individual
reaches — smoothing needs actual polylines, and naively walking every
stream cell to its receiver would tangle every tributary into one
self-overlapping path. `segment_reaches` classifies every stream cell by
in-degree (contributors within the stream network only): 0 = channel head,
1 = ordinary mid-reach cell, >=2 = confluence. Reaches are the maximal
chains between two nodes (head, confluence, or outlet), each reach
including both of its endpoint nodes — so adjacent reaches share exactly
one coincident vertex at a confluence, the same topology a real
hydrographic vector layer (e.g. NHDPlus) uses. Full domain: **17,678
reaches** (9,622 heads, 8,056 confluences, 1,408 outlets), longest reach
260 cells — segmentation naturally keeps reaches short by breaking at every
confluence, which matters for S11's RDP note below. Runtime: 8 s.

**Strahler order fell out almost for free.** The follow-up above anticipated
needing a CSR-style inversion of the full per-cell receiver array (millions
of cells) to get contributor lists; reach segmentation already produces
that adjacency, but over ~17,700 reaches instead of 11M land cells, cheap
enough for a plain recursive resolve with memoization (`strahler_order`,
<0.1 s). Standard rule (a head is order 1; two order-N reaches meeting at a
confluence produce order N+1, otherwise the max order carries through
unchanged). Result: order 1 (9,622 reaches, expected — equal to the head
count) down to order 6 (179 reaches, the handful of major trunk rivers) —
a smooth, plausible distribution with no gaps or spikes.

**Smoothing and the file-size trap.** `chaikin_smooth` (corner-cutting,
4 iterations, endpoints preserved exactly so confluences still connect)
reproduces natural-looking curves from the jagged D8 path — verified on
one reach before the full run (`04_smooth_demo` in the delivered assets:
same path, same endpoints, 46 vertices -> 736). Applied to the whole
network, though, this is where the first full export went wrong: Chaikin
roughly doubles vertex count per iteration, so 4 iterations inflated
~395K reach cells to a **133 MB** GeoJSON — not "lightweight web-ready" by
any definition. Fixed with `simplify_rdp` (Ramer-Douglas-Peucker,
15 m tolerance — half a cell) applied AFTER smoothing: Chaikin needs the
original jagged vertices to know where a curve should bend, but its own
output packs in far more near-collinear points than the resulting smooth
curve needs to render, and RDP strips exactly those back out without
visibly changing the shape (`04_vector_zoom.png` — same crop as the raster
zoom that prompted this, now genuinely smooth). Result: **6.9 MB**, 17,678
features, 6.8 vertices/reach on average.

One implementation note worth keeping: `simplify_rdp` is iterative (an
explicit stack), not the textbook recursive formulation. A near-straight
reach is RDP's worst case (O(n) recursion depth instead of the typical
O(log n)) — and this DEM has genuinely long straight stretches (visible
again in this stage's own crops; see `01_tappa1_terrain.md`'s ridge-shelf
geometry for why), which would have hit Python's default recursion limit
on this project's longer post-Chaikin reaches (~4,200 points before
simplification) had it been left recursive.

### Output

`scripts/export_streams_vector.py` writes `data/exports/streams.geojson`
(6.9 MB — `data/exports/` is meant to be committed, unlike
`data/processed/`) as a `FeatureCollection` of `LineString`s in the
project's local "Fictional World LCC" CRS (metres, same convention as
`data/input/*.geojson` — NOT WGS84 lon/lat; assign the CRS manually in
QGIS after loading). Per-feature properties: `reach_id`, `strahler_order`,
`n_cells`, `max_contributing_area_km2`, `mean_discharge_proxy_m3s`,
`max_discharge_proxy_m3s` — `strahler_order` or the discharge proxy are
both ready to drive line-width scaling directly in
`InteractiveMap.astro` without further processing.

### Locked-in parameters (this addendum)

```
chaikin_iterations   = 4
rdp_tolerance_m       = 15.0   # half a cell; removes smoothing's own
                                # redundant points, not real shape detail
```

### Still open after this addendum

- `lake_mask` still has no equivalent vector export (S10) — the same
  segment-then-simplify approach doesn't directly apply (lakes are
  polygons from a connected-component fill, not a graph of reaches), but
  the RDP piece would carry over directly to whatever boundary-tracing
  step produces the polygon vertices.
- The dead-straight diagonal segment visible in both `04_zoom_streams_lakes.png`
  and this addendum's `04_vector_zoom.png` (upper-center) is a real feature
  of the underlying terrain (confirmed: flow simply follows whatever
  gradient the DEM already has), not introduced by vectorization — worth a
  QGIS look at the DEM itself if it reads as too regular, but out of scope
  for a hydrology-stage fix.

## 12. Addendum — seasonal / intermittent flow

Prompted by a direct fieldwork observation: years of updating river maps on
the ground showed a lot of mapped "streams" as dry or near-dry channels
outside the wet season, with a sharp wet/dry contrast — something a single
annual `discharge_proxy_m3s` grid (S3) cannot represent at all, since it
collapses twelve months of hydrology into one number. This addendum asks,
per stream cell, how many months of the year it plausibly carries flow.

### Why this needs real snow timing, not monthly precipitation directly

Routing each month's *precipitation* through the flow network was
considered and rejected as too crude for this terrain. Tappa 3 already
established that a large fraction of this world's precipitation at
mid-to-high elevation falls as snow and is held rather than immediately
contributing to runoff (`03_tappa3_snow.md`); treating precipitation as
instantaneous runoff would show winter snow accumulation as "flow" the
month it falls and miss the actual melt-season pulse entirely — exactly
backwards for a snow-fed catchment, and exactly the kind of error that
would corrupt the wet/dry classification this addendum exists to produce.

`src/climate/snow.py` gained a new function, `monthly_water_input`, built
by extending the state Tappa 3's `annual_mass_balance` already discards:
Tappa 3 only needed a single annual accumulate-vs-melt balance (permanent
snow: yes/no), but seasonal routing needs the actual month-by-month
sequence of *when* meltwater is released. `monthly_water_input` carries an
explicit snowpack state (mm w.e.) forward across months — accumulating
each month's modelled snowfall (`snow_fraction`, unchanged from Tappa 3),
releasing melt via the same degree-day model (`expected_positive_degree_days`,
`degree_day_factor_mm_per_c_day` — both unchanged, same constants, same
citations) — and returns `rain + melt` as that month's actual water input
to the ground, which is what gets routed, not raw precipitation.

One addition beyond Tappa 3: the simulation runs 3 full annual cycles
("spin-up") and only keeps the last one. A cold-start January with an
empty snowpack is artificially dry at any cell that should be carrying
over a multi-month accumulation from the previous winter; 3 cycles is
enough for the snowpack to settle into a repeating annual pattern rather
than the arbitrary zero it started from (checked directly, not assumed —
see validation below).

### Two bugs caught during this build (documented, not silently fixed)

**Float32/epsilon precision collision.** The first attempt tried to avoid
re-running the ~5.5-minute priority-flood pass by reconstructing
`pop_order` from the already-saved `filled_dem_30m.npy` via
`np.argsort`. A sanity check against the original `discharge_proxy_m3s.npy`
came back wrong (correlation 0.87, not ~1.0) — traced to
`filled_dem_30m.npy` having been saved as **float32** in
`run_tappa4_hydrology.py`. At this DEM's elevation range (1000-4000 m),
float32's representable precision (~0.0001-0.001 m) is the same order of
magnitude as the priority-flood's own `epsilon=1e-4` m tie-breaker, so
21,944 cells that were correctly given distinct parent/child elevations
during the original float64 pass silently collapsed to identical float32
values on save — breaking the strict-descent invariant `receiver_from_codes`
and any reconstruction from disk depends on. Fix: `run_tappa4_seasonal.py`
re-runs `priority_flood_d8` fresh, in-memory, float64 throughout, rather
than reconstructing anything from the float32 export. Confirmed fixed: the
same sanity check (reproduce the annual discharge with the fresh
receiver/pop_order and the original weighting) now matches the saved
annual run at correlation 1.000000, max abs diff 3.8e-6 m³/s.

**Wrong seconds-per-month divisor.** The first working version of the
monthly routing loop divided every month's accumulated water volume by a
full year of seconds (`SECONDS_PER_YEAR`) instead of that month's own
seconds (`days_in_month * 86400`). This is a straightforward mean-rate
unit error, not a hydrology bug, but caught the same way any unit bug
should be: a magnitude check against a value already trusted. Monthly
maxima came out 6.7-11.9 m³/s, roughly an order of magnitude below the
annual max of 103.1 m³/s — not physically possible for a mean-monthly rate
compared against a mean-annual one; if anything, wet-season monthly maxima
should exceed the annual mean rate. Root cause and fix both in
`run_tappa4_seasonal.py` line 108-109 (now divides by `DAYS[m] * 86400.0`
per month). The already-computed (expensive) flow-routed arrays were
corrected in place with a scalar per month rather than re-running the
~9-minute routing pass, since the bug was purely in the final unit
conversion, not the routing itself; `months_flowing` was then
recomputed from the corrected values (perennial fraction moved from
84.8% to 84.6% — the bug was *nearly* but not exactly uniform across
months, since `days_in_month` varies, so the correction is not a pure
rescale).

### Classification and results

A stream cell counts as "flowing" in a given month if that month's
discharge is >= 10% of that same cell's own mean annual discharge — a
relative, self-scaling threshold (`FLOWING_FRACTION = 0.10`) rather than
one absolute discharge constant applied everywhere, because this domain
spans headwater trickles to trunk rivers differing by orders of magnitude
in absolute discharge; a single absolute cutoff would misclassify one end
or the other. `months_flowing` (0-12) is that count, per stream cell,
over the full year.

Full domain (377,281 stream cells):

```
perennial (12/12 months):     319,305 cells (84.6%)
intermittent (1-11 months):    57,976 cells (15.4%)
zero months flowing:                0 cells (0.0%)

months-flowing distribution (intermittent cells only):
  7:   5,600      9:  14,986     11:  5,539
  8:  23,780     10:   8,071     12: 319,305 (perennial, shown for scale)
```

No cell falls below 7 months flowing — i.e. nothing in this model reads as
a truly ephemeral, flash-flood-only channel; the driest classified cells
still flow more than half the year. That ceiling is a direct consequence
of what this model does and doesn't include (see limitations below), not
a claim that this world has no flashier channels than these.

Monthly discharge maxima (post-fix, m³/s), showing the expected
southern-hemisphere seasonal shape — low in austral winter (Jun-Aug),
peaking with spring melt:

```
Jan 117.0  Apr  95.3  Jul  81.4  Oct 105.5
Feb 105.0  May  93.3  Aug  79.1  Nov 139.5
Mar  99.4  Jun  87.1  Sep  82.5  Dec 139.7
```

### Validation

**Elevation trend** (`04_seasonal_overview.png`, left panel): mean months
flowing is a clean, monotonic function of elevation — 12.0/12 (fully
perennial) from sea level through ~1,600 m, then a steady decline to ~7/12
at the highest sampled band (~3,300 m). This is the expected shape for a
snow-timed system: low elevations get rain input close to year-round,
while high elevations increasingly gate flow behind a short melt window,
and there is no discontinuity or reversal anywhere in the curve.

**Windward/leeward, and a confound caught before it was written up.** The
first comparison — wet-tercile vs. dry-tercile cells by annual
precipitation — gave a result that would have been reported backwards if
taken at face value: wetter cells were *less* often perennial (75.3%) than
drier cells (92.0%). Before accepting that, this was checked against the
same elevation/continentality confound already documented in
`02_tappa2_climate.md` S3b: the "wet" tercile disproportionately captures
high-elevation orographic cells (more precipitation *and* colder, snowier,
more melt-gated), so the raw comparison was mixing two effects together
rather than isolating precipitation's own effect. Re-run controlled for
elevation band:

```
elevation band     wet (months)   dry (months)   diff
0-600 m                 12.00          12.00      0.00
2400-3600 m               8.21           8.60     -0.39
```

At matched elevation the gap collapses almost entirely. The small residual
at high elevation runs the *opposite* direction from the raw comparison
(wet slightly less perennial than dry, even controlled) and has a
physically sensible explanation: at the same elevation and temperature,
more precipitation means more snow, which means a bigger, more
concentrated single melt pulse rather than a steadier trickle — consistent
with, not contradicting, the snow-timed model this addendum is built on.
The headline "wetter is less perennial" claim from the raw comparison does
not survive controlling for elevation and would have been a
confound-driven false conclusion if reported as-is.

### Honest limitations (stated before this was built, restated here)

- **No groundwater or baseflow model.** Real intermittent streams often
  keep trickling between rain/melt events on stored subsurface water this
  model has no representation of at all; every cell here goes fully to
  zero the instant modelled surface water input drops below threshold,
  with nothing smoothing the transition.
- **No channel-bed transmission losses.** Real ephemeral channels
  (especially on alluvial or gravel beds) lose flow into the bed over
  distance, which can dry out a channel independent of upstream discharge
  — not modelled; discharge here is a pure upstream-accumulation proxy,
  identical in kind to S3's caveat for the annual figure.
- **These two gaps are almost certainly why nothing classifies below 7
  months flowing.** A groundwater/baseflow term would likely push some
  currently-perennial low-order headwater cells down into a true
  intermittent or ephemeral category; a transmission-loss term would do
  the same for long channels crossing permeable valley floors. This
  model's 84.6%/15.4% split should be read as an upper bound on
  perenniality, not a precise forecast — directionally consistent with
  the fieldwork observation that prompted this addendum (dry season
  really does reduce flow substantially here), but almost certainly
  understating how much of the network goes fully dry in the driest
  months.
- Same as S3/S6: no independent streamflow-gauge validation for this
  world exists to check the *absolute* monthly discharge numbers against;
  what's validated here is internal consistency (mass conservation in the
  snowpack model, the elevation trend, the confound check) — not an
  external ground truth.

### Locked-in parameters (this addendum)

```
flowing_fraction_of_annual_mean = 0.10   # relative, per-cell threshold
spinup_cycles                   = 3      # snowpack spin-up before the
                                          # kept simulation year
```

### Output

`run_tappa4_seasonal.py` writes to `data/processed/hydrology/seasonal/`
(gitignored, regenerate locally — ~9 minute runtime, dominated by the
in-memory priority-flood re-run): `monthly_discharge_proxy_m3s.npy/.bin`
(float32, 12 bands), `months_flowing.npy/.bin` (int16, 0-12),
`tappa4_seasonal_meta.json`. `04_seasonal_overview.png` (delivered asset)
shows the elevation trend and a spatial view of `months_flowing` over the
stream network.

### Still open after this addendum

- No vector export of `months_flowing` (S11's smoothing/simplification
  pipeline would carry over directly — attach `months_flowing` as a
  per-reach property the same way `strahler_order` and the discharge
  fields already are, using each reach's own mean).
- A groundwater/baseflow term and a transmission-loss term (see
  limitations above) are the two most likely next improvements if the
  perennial/intermittent split needs to be more realistic than an upper
  bound.
