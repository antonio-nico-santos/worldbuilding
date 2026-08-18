# Tappa 8 — Geomorphology (lithology, caves, jade/pounamu resource)

First of the three domains scoped into Tappa 8 by `07_tappa7_regional_scenario.md` S9
(geology/lithology, caves, resources). Closed and committed this session. Caves were
already built and approved (v2) before this session started working on lithology;
lithology went through five full methodological iterations before landing; caves' lava
tube layer was then re-run once more to stay consistent with the final lithology.

Dangerous creatures/conflict zones, dangerous seas, and transportation (the other three
Tappa 7 S9 domains folded into "Tappa 8") are **not** covered here — still open.

## 0. Why lithology took five iterations, in one paragraph

S9's original plan (07 doc S1) was distance-based: schist = within each ridge's own
`falloff_km` of its axis, greywacke = the flank beyond that out to shelf reach, basin
fill = the authored plains/plateau zone polygons, volcanic = the SW Island landmass.
Versions 1-4 all iterated on *how that distance geometry was derived* (the authored line
directly, then a DEM-snapped crest cross-section, then an actual DEM ridge network, then
noise-warped) without ever fixing the structural problem: distance-to-a-line cannot
distinguish a flat plateau sitting near a mountain range from the mountain range itself.
Nico caught this directly against a real NZ geology reference map — the authored SE
plains zone, geometrically close to the Spine ridge, was reading as schist/greywacke
under every distance-based version, when the real analog (Canterbury Plains next to the
Southern Alps) is sedimentary regardless of proximity. v5 replaced the entire geometric
approach with a two-field DEM-native classification (elevation + local relief) that
makes flatness itself the deciding signal, not distance to anything authored.

## 1. Lithology v1-v4 — distance-based, all superseded

All four kept for their code/history (`src/geomorphology/real_crest.py`,
`boundary_noise.py`, `run_tappa8_lithology.py` through `_v4.py`), none of their data
outputs re-committed here (large, and methodologically superseded — see S0). Summary:

- **v1** — literal S1 plan: ridge `falloff_km` radius = schist, shelf reach = greywacke,
  authored zone polygons = basin fill, SW Island landmass = volcanic. Visually "too
  smooth, too authored-looking."
- **v2 (Option B)** — schist radius sourced from the nearest cell of the actual DEM
  ridge-crest, found via a radius-limited search on inverted-flow-accumulation cells
  (Priority-Flood+D8 run on `dem.max()-dem`, i.e. treating divides as "streams" of the
  inverted surface). Still resembled v1 too closely per Nico's read.
- **v3 (Option A)** — perpendicular cross-section DEM-elevation argmax per arclength
  station along the authored line, plus height-normalized falloff (1.0x base at the
  ridge's highest point, 0.6x at its lowest — Nico's explicit spec).
  `height_normalized_falloff_m()` is the only v3 piece later reused as-is by nothing
  downstream (v4 replaced the source geometry, keeping this same falloff function).
- **v4** — replaced the authored line as *source* geometry with an actual DEM ridge
  NETWORK: same inverted-flow-accumulation field, thresholded per-ridge. First attempt
  used a single percentile over the whole corridor and put 84% of the Spine's 64.4 km
  length with zero qualifying crest cells — self-caught by checking arclength coverage
  before this was shown to Nico. Fixed with a windowed/local percentile instead
  (`extract_real_crest_network_local`: 1000m station spacing, 750m window half-length,
  97th percentile *within each window* — raised from an initial 90th after checking that
  90th let overlapping windows swallow 47-88% of the entire corridor as "crest"). Added a
  synthetic multi-octave value-noise domain-warp on the ridge-distance query points only
  (`boundary_noise.py` — an independent implementation for this stage, **not** a reuse of
  Tappa 1's actual `domain_warp`, which lives in `terrain.generate`/`terrain.erosion` and
  was never needed here since v5 dropped noise entirely). Visually the best of the
  distance-based versions, but still shared v1-v3's structural flaw (S0).

## 2. Lithology v5 — DEM-native, locked

`src/geomorphology/terrain_relief.py`, driven by `run_tappa8_lithology_v5.py`
(the file kept as-is, versioned name and all — no `_v5` was stripped at commit time,
see S6 for why).

**Method.** Two DEM-derived raster fields only, no authored line/zone/corridor geometry
enters the classification decision at all:

- **Elevation** — the DEM itself.
- **Local relief** — windowed max−min elevation, `window_m=2000.0`, via
  `scipy.ndimage.maximum_filter`/`minimum_filter` (separable, O(n) regardless of window
  size). A flat high plateau reads *low* on this even though it's high — the exact
  distinction "distance to a ridge" could never make.

**Mainland split** (thresholds calibrated as percentiles of mainland land's own
elevation/relief distributions, same "calibrate against this world's own data"
convention already used for jade grade percentiles and the biome moisture terciles):

1. `schist = (elev >= p60) & (relief >= p75)`, **OR** `elev >= p90` regardless of relief.
2. `greywacke = relief >= p50`, for everything not already schist.
3. `basin_fill` = everything else classifiable (the default/leftover class, not a
   competing polygon — mirrors S1's original precedence rule, just decided by threshold
   now instead of paint order over authored zones).
4. Applied to all non-volcanic land, including ~200+ minor islets — **not** restricted to
   the mainland connected-component (an early bug: restricting both calibration *and*
   application to the mainland component left 268,101 islet cells / ~241 km² silently
   defaulted to the ocean class code; caught by an ocean-cell-count mismatch before
   anything was shown).

**The `elev >= p90` OR-branch** exists because of a second, real terrain property: this
world's ridges (Tappa 1's Gaussian falloff decay) have smoothly rounded summits, not
knife edges. A fixed 2 km relief window centered exactly on a rounded high point often
measures *less* relief than a window centered on the steep flank just below it — verified
directly by sampling: 33-47% of the actual authored crest lines came out `basin_fill` at
elevations up to 3350m, with as little as 159-265m of measured local relief at those
points. No amount of relief-window tuning fixes this, because the terrain really is
locally smooth exactly there. The OR-branch (elevation alone above mainland's own p90,
~2612m) recovers these cells regardless of relief. Checked before adding it that this
can't reopen the SE-plains fix: the highest authored plateau (Central plateau) tops out
at 2463m, below the 2612m bar, so no authored flat zone can trigger it.

**Island split** (Nico's explicit follow-up request — extend the same logic to the SW
Island, previously unconditionally 100% volcanic since Tappa 7 S9/Tappa 8 v1 through v4,
a scope lock the relief redesign hadn't touched on its own). Real Banks Peninsula (this
project's own volcanic citation) is mostly-but-not-entirely volcanic bedrock — it has a
real flat isthmus/harbour-margin apron. Calibrated against the **island's own**
elevation/relief population (its whole elevation range, 0-737m, sits below every
mainland threshold, so reusing mainland numbers would have painted the entire island
`basin_fill`, wrong in the other direction):

- `island_volcanic = relief >= p25(island land)`; `island_basin_fill = relief < p25`.
- No schist/greywacke tier on the island — a young monogenetic shield doesn't have that
  metamorphic-flank story.

**Precedence, final**: volcanic (island only, by relief split) > schist > greywacke >
basin fill (default/leftover), unchanged from S1's rule in spirit.

## 3. Lithology v5 — resolved thresholds and validation

| threshold | value |
|---|---|
| elev (schist AND-branch floor, mainland p60) | 1061 m |
| relief, schist AND-branch floor (mainland p75) | 469 m |
| relief, greywacke floor (mainland p50) | 272 m |
| elev, schist OR-branch floor (mainland p90) | 2612 m |
| relief, island volcanic floor (island p25) | 202 m |

Cells recovered by the OR-branch alone (i.e. failed the AND test but qualify on
elevation alone): 728,746 — meaningfully large, consistent with how much of the crest
network the rounded-summit artifact affects.

Class areas (final, mainland + island combined):

| class | km² |
|---|---|
| sedimentary basin fill | 4437.3 |
| greywacke/argillite | 3098.5 |
| schist | 1780.9 |
| volcanic | 595.7 |

Island split: 595.7 km² volcanic (75.0%) / 198.6 km² basin fill (25.0%) — exact 75/25 by
construction of the p25 cut, not itself a meaningful number; what was checked instead was
*where* the margin landed (the low-lying coastal/isthmus fringe, not a random speckle —
confirmed visually).

**Authored-zone validation** (QA only — these zone polygons are no longer a
classification input, just a check against them):

| zone | basin_fill % | greywacke % | schist % | volcanic % |
|---|---|---|---|---|
| NW Plateau | 79.6 | 18.9 | 1.6 | 0.0 |
| North plains | 87.0 | 11.8 | 1.2 | 0.0 |
| Central plateau | 56.6 | 40.8 | 2.6 | 0.0 |
| South plains | 90.8 | 9.1 | 0.1 | 0.0 |
| SE plains | 78.5 | 21.5 | 0.0 | 0.0 |
| Island plateau | 79.8 | 0.0 | 0.0 | 20.2 |

SE plains at 78.5% basin fill / 0% schist is the direct fix for the complaint that
started v5. "Island plateau" flipping from 99.99% volcanic (trivial under v1-v4, when
the whole island was volcanic by construction) to 79.8% basin fill is a large, real
change, not a bug: sampled directly, the zone's own elevation (median 123m, max 397m)
and relief (median 119m, 75th percentile 185m — below the island's own 202m threshold)
are genuinely low relative to the rest of the island. **Confirmed by Nico**: "Island
plateau" was authored to represent the flat isthmus/apron, not the volcanic highland, so
this result is correct, not a regression.

Basin fill grounded against the real stream network (S1's hydrology check, re-verified
against v5's larger basin-fill footprint): 77.0% of basin-fill area sits inside one of
the six authored plains/plateau zones (the rest is basin-fill elsewhere on the DEM,
which is expected now that the class isn't zone-authored).

## 4. Caves — re-grounded against lithology v5

Caves themselves were built and approved (v2) *before* lithology's v5 redesign; only the
lava tube layer consumes a lithology raster (`lithology == CLASS_VOLCANIC`), so it's the
only one of the four cave types this touches. Checked directly against every candidate
function in `caves.py`, not assumed: `talus_pseudokarst_candidates` (slope + stream
proximity + landmass ID), `sea_cave_candidates` (slope + coastal distance), and
`glacier_moulin_candidates_v2` (snow mask + slope + margin depth) take no lithology
input at all, and `mainland_mask` (used by talus) comes from plain connected-component
landmass labeling, never lithology. Re-running only the lava tube layer is therefore
complete, not partial.

**v1 → v2** (already approved, unchanged this session): narrowed v1's "entire volcanic
class" candidate (794.3 km², literally the whole island) to `volcanic & vent_weight>=0.1
& slope<=p25(within volcanic zone)` — real vent proximity (Gaussian falloff from the 3
`geothermal.geojson` points, computed in v1 but never actually used as a filter until
v2) and a gentle-slope preference (inverse of talus/sea's steep preference — lava tubes
form in channelized flows down gentle-to-moderate flanks, not on cliff faces or crater
rims). Glacier/moulin got the same v1→v2 narrowing treatment (snow-mask-local slope p75
AND ice-margin proximity); talus/sea caves are identical logic since v1, re-approved
as-is.

**v2 → v3** (this session, triggered by lithology v5's island split): re-masked the
*same* v2 test against `lithology_v5`'s `CLASS_VOLCANIC` instead of v2's. Two effects,
both verified numerically before reporting:

1. **Footprint loss** — the volcanic-eligible zone shrank from the full island (794.3
   km²) to just the core (595.7 km²); 9,463 cells (8.5 km²) that were candidates under v2
   are now outside the volcanic zone and dropped.
2. **Threshold recalibration, counterintuitive direction** — the gentle-slope cutoff is
   a local p25 computed *within* the volcanic zone. Removing the flat margin from that
   population makes the remaining population steeper on average, so its own p25 rises:
   6.17% → 8.94%. A higher "gentle enough" bar means *more* of the smaller, steeper core
   qualifies: 16,904 cells (15.2 km²) newly eligible.

Net: total lava tube candidate area went **up**, 26.3 km² → 33.0 km², despite the
eligible zone getting smaller, because the threshold is relative to its own population,
not absolute. Vent proximity weighting and its threshold (0.1) are unchanged placeholders
— checked with Nico and left as-is.

## 5. Jade/pounamu + Vértice materials

Mechanics unchanged throughout all lithology versions — only the schist mask and grade
field feeding it changed. `place_jade_pods()` (stochastic, weighted by grade, minimum
5 km separation, 10 pods, seed=13) was never touched this session.

`schist_grade` (v5): rank-based, `0.5*elevation_rank + 0.5*relief_rank`, both ranked only
among schist cells themselves (not all mainland land), replacing v2-v4's ridge-falloff
Gaussian grade — the ranking method changed because the underlying geometry it would have
graded against (distance to a ridge) no longer exists in v5's classification.

Final placement: 10 pods, 6.00 km² pod area, 356.2 km² suitability zone (grade ≥ p80).
Full pod centers/radii are in `tappa8_lithology_v5_meta.json`.

`src/geomorphology/resources.py` ships the Vértice material lookup table (schist →
gold-bearing quartz veins/Onda + muscovite/Energia; greywacke → laumontite/Matéria;
basin fill → vivianite/Bios + weak reworked magnetite/Campo; volcanic → magnetite/Campo +
native silver/Mente) as structured data, reproducing what was already locked in
`07_tappa7_regional_scenario.md` S3 rather than leaving it only in prose. Non-spatial —
a lookup table keyed by lithology class, not a raster.

## 6. Committed data — what's in and what's left out

Every version of the code is committed (all of `src/geomorphology/*.py`, all
`run_tappa8_lithology*.py` / `run_tappa8_caves*.py` / `run_tappa8_ridge_extraction.py`)
— cheap, and it's the actual record of what was tried and why each earlier attempt was
dropped (S0-S1 above). Filenames keep their in-session version suffixes (`_v2`, `_v5`,
etc.) rather than being renamed to a single canonical name at commit time — this project
has no other Tappa 8 commit to be consistent with yet, and renaming without re-running
would have meant re-verifying results under new names for no functional benefit.

Data is **not** kept for every version, unlike the code. Only the final round's rasters
are committed:

- Lithology: `lithology_v5.*`, `basin_fill_grounded_v5.npy`, `jade_suitable_v5.npy`,
  `jade_pods_v5.*`, `tappa8_lithology_v5_meta.json`.
- Caves: `cave_lava_tube_v3.*` (final), `cave_glacier_moulin_v2.*` (unchanged since v2,
  still current), `cave_talus_pseudokarst.*` and `cave_sea_cave.*` (unchanged since v1,
  still current/approved), `tappa8_caves_v2_meta.json` (carries v1 numbers forward for
  comparison) and `tappa8_caves_v3_meta.json` (final lava tube numbers).

**Mechanical constraint discovered mid-commit, and how it was handled**: the file
delivery path used to reach the device turned out to have two stacked limits, found the
hard way (both rejected a real upload attempt before being worked around) --
30 MiB/file to get a file into the conversation at all, then a *stricter* 20 MiB/file
(100 MiB/call total) on the actual device-write step. Three fixes were applied
depending on what the file actually was:

1. **Classification/mask rasters, dtype** (`lithology_v5`, `jade_pods_v5`, and all four
   `cave_*` layers) are small-integer data (0-4 class codes or 0/1 booleans) that were
   being exported to ENVI `.bin` at int16 (2 bytes/cell, ~44-46 MiB) purely because
   `write_envi_raw` only offered `f4`/`i2` before now. Added a third dtype, `u1`
   (unsigned byte, ENVI data type code 1 — standard, not invented) to
   `src/terrain/raster_io.py`'s `_ENVI_DTYPE` table — additive only, existing `f4`/`i2`
   callers elsewhere in the project are untouched. Re-exported these six layers at
   `u1` (22.0 MiB each) and verified byte-for-byte round-trip against the source `.npy`
   before shipping (`np.fromfile(...).reshape(...)` compared with `np.array_equal`
   against `np.load(...)`, all six matched exactly).
2. **Same rasters, still over the stricter 20 MiB device-write limit even as uint8**
   (22.0 MiB > 20 MiB, and the `.npy` files were the same ~22 MiB) — gzipped instead of
   shipped raw. This data is overwhelmingly boolean/5-class categorical over huge
   contiguous regions, so gzip -9 crushes it: all fourteen files (the eight `.npy` +
   six `.bin`, `.hdr`/`.prj`/meta JSON were already small) came down to 20 KB-644 KB
   each, 2.7 MB combined. Shipped as `<name>.gz` at the same repo-relative path, plus a
   new committed script, `scripts/decompress_tappa8_data.py` (stdlib `gzip` only, no
   new dependency), that decompresses every `*.gz` under
   `data/processed/geomorphology/` back to its exact original filename and deletes the
   `.gz`. Verified end-to-end in a throwaway directory before shipping: compressed all
   fourteen, decompressed with the exact script being committed, byte-compared every
   output against the pre-compression original — all fourteen matched exactly.
   **Action needed once, locally**: run `python scripts/decompress_tappa8_data.py`
   from the repo root after pulling this commit to materialize the real `.npy`/`.bin`
   files; the pipeline scripts and QGIS both expect the decompressed form, not `.gz`.
3. **Continuous float32 fields** (`schist_grade_v5.npy`, `relief_2km.npy`,
   `slope_pct_30m.npy`, `dist_to_stream_km_30m.npy`, `dist_to_ocean_km_30m.npy`,
   `lava_tube_vent_weight.npy`, all ~88-89 MiB) — no lossless dtype narrowing gets these
   under 30 MiB without either compressing them (which the device write path can't
   decompress on its own, since this session has no remote-execution access on the
   device, only file read/write) or quantizing them (a real precision loss this doc
   isn't deciding unilaterally). **Dropped from the data commit entirely.** All six are
   cheap to regenerate (each is seconds of compute: a filter, a distance transform, or a
   Gaussian-falloff sum over 3 points) from already-committed inputs (the DEM,
   `stream_mask.npy`, `land_mask`) via already-committed code — `relief_2km` and
   `schist_grade_v5` come out of `run_tappa8_lithology_v5.py`,
   `slope_pct_30m`/`dist_to_stream_km_30m`/`dist_to_ocean_km_30m`/`lava_tube_vent_weight`
   out of `run_tappa8_caves.py`/`_v2.py`. Same "regenerate with, in order" convention
   already used for Tappa 6 in `parameters.yml` — nothing here is unrecoverable, just not
   worth shipping as a static binary when the code to remake it in seconds is right
   there.

**Explicitly dropped**: `lithology.npy`/`_v2`/`_v3`/`_v4` and their `schist_grade`/
`jade_pods`/`jade_suitable`/`basin_fill_grounded` companions, `cave_lava_tube.npy`(v1)/
`_v2.npy` (v2's own raster, superseded by v3 — the v2 *numbers* survive in
`tappa8_caves_v2_meta.json`), `cave_glacier_moulin.npy` (v1), `ridge_accum_cells.npy`,
and the v1-v4 lithology/v1 caves meta JSONs. These represent genuinely superseded
methodology (S0-S1), not intermediate pipeline steps of the winning approach — unlike,
say, Tappa 4's `flow_direction_code.npy` or `contributing_area_km2.npy`, which are
committed intermediates of hydrology's one (correct) pipeline. Nothing here is lost:
the code that produced them is committed and the numbers that matter are carried
forward in prose (this doc) or in the surviving meta JSONs. Flagging this as a real
policy choice made without asking first — if full binary history across all five
lithology rounds is wanted after all, it still exists in this session's workspace and
can be committed separately.

Rendering scripts committed to `scripts/`: `make_tappa8_lithology_plots.py`
(v4-vs-v5 comparison, authored-zone overlay, jade overview, island zoom panel),
`make_tappa8_caves_plots.py` (v1-vs-v2 and v2-vs-v3 comparisons, island lava-tube delta
panel). Also `decompress_tappa8_data.py` — **run this once** after pulling this commit,
see the gzip note above; it's not optional, the shipped rasters are `.gz` until it runs.

## 7. Locked-in parameters (recap)

```
lithology (v5, terrain_relief.classify_from_terrain):
  relief_window_m: 2000.0
  elev_percentile: 60.0
  relief_schist_percentile: 75.0
  relief_greywacke_percentile: 50.0
  high_elev_percentile: 90.0
  island_relief_percentile: 25.0

jade:
  n_pods: 10
  min_separation_km: 5.0
  radius_range_m: [300.0, 800.0]
  grade_percentile_threshold: 80.0
  seed: 13

caves (v2/v3):
  steep_threshold_pct (talus/sea, mainland p75): 20.14
  glacier_moulin: snow-mask-local slope p75 = 22.55%, margin depth threshold = median = 2.48 km
  lava_tube: vent_weight_threshold = 0.1, gentle-slope threshold = local p25 within
             the volcanic zone (6.17% under lithology v2, 8.94% under lithology v5)
  talus stream_buffer_km: 0.5
  sea coastal_buffer_km: 0.5
```

**Every threshold above is a first-pass placeholder**, calibrated against this world's
own data distribution (percentiles), not against any external number — same status as
every other percentile-based threshold already locked elsewhere in this project (jade
grade, biome moisture terciles, Tappa 6's slope/water limits). None have been
independently recalibrated beyond the checks narrated in S2-S4 above.

## 8. Open follow-ups (not done this stage, deliberately left open)

- Schist fracture caves — still explicitly parked (07 doc S2), not built.
- Dangerous creatures/conflict zones, dangerous seas, transportation — the other three
  Tappa 7 S9 domains folded into "Tappa 8," not started.
- No independent calibration exists yet for any of the percentile thresholds in S7 above
  beyond the checks already narrated (SE-plains fix, crest-line recovery, island split,
  lava tube footprint/threshold effect).
