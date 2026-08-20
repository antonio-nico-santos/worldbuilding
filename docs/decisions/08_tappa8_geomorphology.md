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

*(Superseded by a later addition, not by anything wrong here: once `lithology_v6`
existed — S8 — a fifth cave type, `karst_cave_candidates`, was added, and it DOES
consume lithology, marble/sedimentary_limestone specifically. See S8c for that one;
this section's "only lava tube" claim was true for v2/v3 caves and stays true for the
original four, it just no longer describes the full cave-type roster.)*

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
  `jade_pods_v5.*`, `tappa8_lithology_v5_meta.json` — v5 is kept even though v6
  supersedes it as "current," since v6 is a composite ON TOP of v5 (see S8), not an
  independent reclassification; `basin_fill_grounded_v5`/jade both still key off v5's
  own schist mask, unaffected by the authored overlay.
- Lithology authored overlay (v6, new this session — S8): `lithology_v6.*`,
  `tappa8_lithology_v6_meta.json` (per-zone composition table, S8b).
  `data/input/lithology_authoral.geojson` (the hand-drawn source) lives under
  `data/input/`, not `data/processed/`, alongside `terrain_ridges.geojson`/
  `terrain_zones.geojson` — same convention, authored input vs. derived output.
- Caves: `cave_lava_tube_v3.*` (final), `cave_glacier_moulin_v2.*` (unchanged since v2,
  still current), `cave_talus_pseudokarst.*` and `cave_sea_cave.*` (unchanged since v1,
  still current/approved), `tappa8_caves_v2_meta.json` (carries v1 numbers forward for
  comparison) and `tappa8_caves_v3_meta.json` (final lava tube numbers).
- Caves, karst (v4, new this session — S8c): `cave_karst.*`, `tappa8_caves_v4_meta.json`.
- Caves, blended bitmask (new this session — S8d): `cave_blend.*`,
  `tappa8_cave_blend_meta.json`. Additional convenience layer; the five individual cave
  rasters above are still committed and still the authoritative per-type source.
- Resource pods (new this session — S8e): `resource_laumontite.*`,
  `resource_vivianite.*`, `resource_placer_magnetite.*`, `resource_silver_copper.*`,
  `tappa8_resource_pods_meta.json` (per-material weight basis, pod centers/radii, the
  placer-magnetite citation-correction flag). Schist's gold/mica and volcanic's
  magnetite (primary) deliberately have no new file — see S8e for why.
- Transport lithology-cost multiplier (new this session — S8f): lives under
  `data/processed/suitability/` (a Tappa 6 output directory, not geomorphology's own —
  this layer belongs to the cost-distance graph it modifies, not to lithology itself):
  `transport_friction_multiplier_120m.*` (the multiplier raster, 5.8 MiB float32,
  under both size limits, shipped uncompressed) and
  `tappa8_transport_friction_meta.json` (friction table + citations, class areas,
  the baseline-vs-friction pairwise comparison for all 17 already-placed Circulos).
  `circulo_candidate_sites.geojson`/`tappa6_site_selection_meta.json` are unchanged —
  this is a new, separate layer, not a Tappa 6 re-run (see S8f).
- Excavation effort (new this session — S8g): `excavation_effort_multiplier_120m.*`
  (5.8 MiB float32, uncompressed) and `tappa8_excavation_effort_meta.json`, under
  `data/processed/geomorphology/` (unlike S8f's friction layer, this has no cost-graph
  dependency, so it lives alongside lithology's other outputs, not Tappa 6's).
- Basin fill zonation (new this session — S9): `basin_zonation_30m.*` (native 30m, u1 —
  gzipped before shipping, same size reason as every other 30m categorical raster this
  session, see the mechanical-constraint note below) and
  `tappa8_basin_zonation_meta.json` (calibration record, area table, and the separate
  basin_fill/Grassland overlap check).
- Iron + aluminium (new this session — S10): `resource_vivianite.*` **overwritten**
  (re-placed into `wetland_backswamp`), `resource_bog_iron.*` (new, identical pod
  geometry to vivianite by design), `resource_bauxite.*` (new), all gzipped before
  shipping (same 20 MiB device-write constraint as every other native-30m raster this
  session), plus `tappa8_iron_aluminium_meta.json` (both decisions' full reasoning,
  pod counts/areas, the placer-magnetite consistency flag).
- Resource blend (new this session — S11): `resource_blend.*` (native 30m, u1, all 7
  materials OR-packed, gzipped before shipping — 27 KB, 0.12%, even smaller than the
  individual pod rasters' own gzip since the same bit pattern repeats across most of the
  grid) and `tappa8_resource_blend_meta.json` (bit table, area/overlap diagnostics, the
  bog_iron/vivianite redundancy flag, round-trip verification result). The seven
  individual pod rasters (plus `jade_pods_v5.npy`, re-staged from the device's S5
  commit) stay committed too — this is a convenience layer, not a replacement, same
  status as `cave_blend.npy` (S8d).
- Placer magnetite re-placement + resource blend update (new this session — S12):
  `resource_placer_magnetite.*` **overwritten** (re-placed into `estuarine_coastal`) and
  `resource_blend.*` **overwritten** (bit 16 updated to match), both gzipped before
  shipping, plus `tappa8_placer_magnetite_restrict_meta.json` (the structural
  non-overlap check against vivianite/bog_iron, area before/after).
- Quartz/mica/gold + resource blend re-run at int16 (new this session — S13):
  `resource_quartz.*`, `resource_mica.*`, `resource_gold.*` (all new), `resource_blend.*`
  **overwritten again** (now ENVI `"i2"`/int16, 10 bits, ~46 MiB raw, gzipped to ~52 KB),
  plus `tappa8_schist_vein_materials_meta.json`. `schist_grade_v5.npy`, regenerated to
  make this possible, is NOT re-committed — same "cheap to regenerate, drop the static
  binary" policy S6 already applied to it and five other continuous float32 fields.
- Sedimentary-limestone/granite pods + resource blend re-run again, still int16 (new
  this session — S14): `resource_calcite.*`, `resource_mica_granite.*`,
  `resource_quartz_granite.*` (all new), `resource_blend.*` **overwritten a third time**
  (13 bits now, same `"i2"` dtype — 13 bits' max value 8191 still fits, no further dtype
  change needed), plus `tappa8_limestone_granite_materials_meta.json`. Closes the
  marble/sedimentary_limestone/granite item S13's own follow-ups had flagged as blocked
  on the Scenario chat.

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
   **Reused as-is for this session's three new categorical rasters** (`lithology_v6`,
   `cave_karst`, `cave_blend` — all `u1`, all ~22.0 MiB raw, same shape/reasons as
   above): no script changes needed, since `decompress_tappa8_data.py` already globs
   every `*.gz` under `data/processed/geomorphology/` rather than a fixed filename
   list. Compressed sizes: `lithology_v6` 244 KB (1.06%, five real classes now instead
   of four), `cave_karst` 32 KB (0.14%, small isolated footprint), `cave_blend` 900 KB
   (3.90%, higher entropy since it's a 5-bit combination, not a single boolean/5-class
   field) — all three round-trip-verified byte-for-byte the same way as before.
   **Same again for the four S8e resource-pod rasters**: tiny sparse boolean masks (a
   handful of small disks in an otherwise-empty 5334×4334 grid), so gzip crushes them
   even harder — 23-24 KB each (0.10%), all four round-trip-verified.
   **Same again for S9's `basin_zonation_30m`**: a 4-class categorical raster (mostly
   `not_basin_fill`, the other 3 classes partition basin_fill's own ~4335 km² footprint)
   — compressed to 150 KB / 150 KB (`.npy`/`.bin`, 0.66%), round-trip-verified.
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

caves (v4, karst -- new this session, see S8):
  stream_buffer_km: 0.5 (same number as talus/sea, reused not re-derived)
  no slope condition

lithology authoral zones (v6 -- new this session, see S8):
  priority_rank tiers: 1 beats {schist, greywacke, basin_fill}
                       2 beats {greywacke, basin_fill}
                       3 beats {basin_fill}
  (volcanic is never overwritable by any tier; there is no tier below basin_fill)
```

**Every threshold above is a first-pass placeholder**, calibrated against this world's
own data distribution (percentiles), not against any external number — same status as
every other percentile-based threshold already locked elsewhere in this project (jade
grade, biome moisture terciles, Tappa 6's slope/water limits). None have been
independently recalibrated beyond the checks narrated in S2-S4 above.

## 8. Authored lithology additions — marble, sedimentary limestone, granite (v6, implemented)

Raised in a follow-up conversation about Urban Scale construction materials, continuing
into a parallel scenario-discussion chat (see `scenario_reference.md` S22). Geometry was
hand-drawn by Nico (`data/input/lithology_authoral.geojson`, six polygons, iterated
several rounds against real-time QA feedback — geometry validity, feature_type spelling,
`priority_rank` semantics, and composition-against-`lithology_v5` checks, all before any
of this was authored into code); this session wrote the code that consumes it, ran it,
and is recording the result here as locked, not a candidate anymore.

**Limestone — real NZ precedent, structurally different from every existing class.**
Real NZ limestone (Oamaru stone, Waitomo/Te Kuiti Group, Punakaiki) is Oligocene (~34-24
Ma), deposited during a period when Zealandia was tectonically quiet, low-relief, and
partially submerged — a carbonate "cover" laid unconformably over the older basement,
specifically in areas that were *not* undergoing the uplift/erosion that produces
schist/greywacke exhumation. That's why it can't be a sixth threshold in
`terrain_relief.classify_from_terrain`: elevation and local relief describe current
terrain shape, and limestone's presence is a fact about geological age/depositional
history, orthogonal to shape. Two distinct flavors, both drawn, not mutually exclusive:
plain sedimentary limestone (Oamaru-style — dimension stone, and the fix for both the
mortar gap and the plaster gap below via lime plaster/lime mortar) and a metamorphosed
marble variant (Takaka/Arthur Marble-style — a separate, much older Ordovician terrane
story, with its own karst caves — Harwoods Hole, Oparara).

**Granite — confirmed real, but a minor/monumental-stone addition, not structural.**
Coromandel granite is a real, citable NZ dimension-stone tradition (Parliament House,
Auckland Museum, war memorials), but genuinely niche relative to greywacke/schist/basalt
volume — decorative/monumental use, not everyday construction. Geologically unrelated to
the schist/greywacke basement (a Miocene volcanic-arc intrusion, not a relief signature),
same as limestone/marble — confirmed by its drawn zone landing 100% inside basin_fill
(see S8b), not touching either mountain-core class.

### 8a. Implementation — class codes, priority tiers, and the compositing step

Added directly to `lithology.py` (not the legacy pre-v5 `classify()`/`build_zone_fields`
path, which is still dead code from v2-v4 — see S1): three new class codes,
`CLASS_MARBLE=5`, `CLASS_SEDIMENTARY_LIMESTONE=6`, `CLASS_GRANITE=7`, and two new
functions. `load_authoral_zones()` parses `lithology_authoral.geojson` — deliberately
NOT `terrain.skeleton.build_zone_fields`, since that loader's `ZoneField` carries
Tappa-1-specific fields (`target_elevation_m`, `amplitude_scale`, `edge_transition_km`)
meaningless for a rock-type layer, and its softened `blend_weight` edges are wrong here
too (a rock-type boundary is a hard structural fact, same reasoning the legacy
`classify()` already used). Handles both `Polygon` and `MultiPolygon` geometry types —
QGIS re-saved the file as single-part `MultiPolygon` partway through Nico's editing
session, exterior ring only per part, holes not yet exercised by any of the six zones.

`apply_authoral_zones()` composites the authored zones onto the v5 DEM-native array
using each zone's `priority_rank` as a tier in the chain
`volcanic > 1 > schist > 2 > greywacke > 3 > basin_fill` (Nico's own scheme, chat this
session — not a citation). A zone at rank *N* can only overwrite base classes strictly
weaker than its own slot; volcanic is never overwritable by any rank, and there's no
tier below basin_fill (a rank-4 zone would be a permanent no-op — exactly the bug caught
and fixed in chat, mid-session, before this code existed: `Granite South` briefly sat at
rank 4, which would have made it invisible in the output despite being a well-cited,
cleanly-drawn zone). Zones are painted weakest-rank-first so a stronger rank wins any
future overlap between two authored zones — not currently exercised (zero pairwise
overlap among all six, verified before this code was written).

### 8b. Result — `lithology_v6`, per-zone composition, real numbers

```
class_areas_km2 (v6):
  sedimentary_basin_fill: 4334.5   (v5: 4437.3 -- lost 102.8 km2 to the three rank-3 zones)
  greywacke_argillite:    3030.1   (v5: 3098.5 -- lost 68.5 km2 to the three rank-1 marble zones)
  schist:                 1772.8   (v5: 1780.9 -- lost   8.1 km2 to the same marble zones)
  volcanic:                595.7   (unchanged -- no authored zone touches it)
  marble:                   77.8
  sedimentary_limestone:    87.7
  granite:                  13.8
```

Per-zone, `priority_rank`, and the pre-authoral composition each polygon actually landed
on (the real grounded-check evidence, computed by `apply_authoral_zones`, not asserted
by the geojson's own `grounded` property):

| zone | type | rank | nominal km² | painted km² | pre-authoral composition |
|---|---|---|---|---|---|
| Marble North | marble | 1 | 42.40 | 42.40 | 89.7% greywacke, 10.3% schist |
| Marble Wall | marble | 1 | 24.08 | 24.08 | 93.9% greywacke, 6.1% schist |
| Forest Marble | marble | 1 | 11.37 | 11.37 | 69.0% greywacke, 19.9% schist, 11.1% basin_fill |
| North Coast Limestone | sedimentary_limestone | 3 | 63.60 | 48.01 | 75.5% basin_fill, 21.5% ocean, 3.0% greywacke |
| Sedimentary Bay | sedimentary_limestone | 3 | 55.80 | 39.71 | 71.2% basin_fill, 28.8% ocean |
| Granite South | granite | 3 | 13.84 | 13.84 | 100% basin_fill |

Two things worth being explicit about, not just implicit in the table: (1) "painted" is
smaller than "nominal" for `North Coast Limestone` and `Sedimentary Bay` because roughly
a fifth to a third of each polygon was deliberately drawn over ocean (Nico's own call —
easier to draw that way — confirmed in chat, not a bug); ocean cells never paint,
regardless of rank. (2) All three marble zones sit majority-greywacke with a real but
minority bite into schist (6-20%) — read as "at the margin of the schist belt," matching
Takaka Marble's real siting relative to the Otago Schist basement, not "replacing the
core of it." All six zones' `grounded: true` claim is now backed by this table, not just
asserted — this is the first time that flag has meant anything checkable.

### 8c. Karst caves — a fifth, mechanistically distinct cave type, now unlocked

`caves.py` gains `karst_cave_candidates()`: soluble lithology (marble OR
sedimentary_limestone, from `lithology_v6`) AND within 0.5 km of a stream (reusing
talus/pseudokarst's own `stream_buffer_km` convention — Tappa 6's `water_gentle_km` —
not a new number). Deliberately no slope condition: unlike talus/sea/lava-tube, real NZ
karst doesn't consistently correlate with steep relief (Waitomo is rolling hill country,
Harwoods Hole is alpine, Punakaiki is a low coastal platform), so adding one would be an
invented constraint, not a citation. Result: **144.1 km² total** (75.2 km² inside marble
— 96.6% of that class — and 68.9 km² inside sedimentary_limestone — 78.5% of that
class), committed as `cave_karst.*` / `tappa8_caves_v4_meta.json`.

### 8d. All five cave types, blended into one bitmask raster

Nico asked whether the five separate cave rasters could be one file with the pixel value
carrying type info. Checked before picking an encoding: the five masks are **not**
mutually exclusive — `talus_pseudokarst` alone shares 111.0 km² with `glacier_moulin`,
142.6 km² with `sea_cave`, and 55.2 km² with `karst` (all geologically sensible: steep
terrain near a stream can simultaneously sit in the alpine snow zone, the coastal
buffer, or a soluble-rock zone). 308.4 km² carries 2+ types on the same cell, one small
patch (0.4 km²) carries 3. A single mutually-exclusive category code would have silently
discarded that overlap, so this is a **bitmask**, not a category code: `1=lava_tube,
2=talus_pseudokarst, 4=glacier_moulin, 8=sea_cave, 16=karst`, values 0-26 observed
(0-31 possible), one byte/cell (`u1`), decoded per-type with `(value & bit) > 0`.
Round-trip verified: decoding all five bits back out of `cave_blend.npy` reproduces
every one of the five source rasters exactly. Committed as `cave_blend.*` /
`tappa8_cave_blend_meta.json`; the five individual rasters are also still committed
unchanged — this is an additional convenience layer, not a replacement.

### 8e. Resource pod placement — generalizing jade's method to the other four materials

Nico asked to map the other Vertice materials' possible locations "the same way we did
with jade." `resources.py`'s `VERTICE_MATERIALS` table had six material entries, all
non-spatial (a per-class lookup) except jade itself, which got real stochastic pod
placement back in S5 (`place_jade_pods`, weighted by `schist_grade`, min 5 km
separation). `lithology.py` gains `place_material_pods()` — the same technique,
factored out to take an arbitrary eligibility mask + optional weight field instead of
being hardcoded to schist. `place_jade_pods` itself is untouched (already
locked/committed).

**Only four of the six materials get a new raster, not six — two are excluded on
purpose:**
- Schist's gold-bearing quartz veins / muscovite mica: `resources.py` already said
  these co-locate with the jade/pounamu high-grade band (`spatial_note`, pre-existing).
  Giving them independent pods would contradict that note — they reuse `jade_pods_v5.npy`
  directly.
- Volcanic's magnetite (primary): titanomagnetite in basalt is a disseminated bulk
  mineral (a few percent through the rock), not vein/joint/bog-localized like the other
  five. Treating it as discrete pods would be LESS geologically honest than what
  `resources.py` already did — a class-level fact. Stays that way.

**The four that do get pods, each weighted by a different, deliberately-chosen field**
(`run_tappa8_resource_pods.py`, 8 target pods/material, same 5 km min-separation and
300-800 m radius range as jade — first-pass placeholders, not independently calibrated):

| material | class | weight field | pods placed | pod area km² |
|---|---|---|---|---|
| laumontite | greywacke | **uniform** (no citable within-class gradient exists) | 8 | 8.3 |
| vivianite | basin_fill | wetness proxy, 1/(1+dist_to_stream_km) | 8 | 8.7 |
| placer magnetite | basin_fill | coastal proximity × volcanic-landmass proximity | 8 | 4.7 |
| native silver/copper | volcanic | `vent_weight` (same field lava tubes use) | 6 | 5.1 |

Laumontite's uniform weighting is a deliberate choice, not a gap: schist has a real
metamorphic-grade-rank citation to weight jade/gold by; greywacke has no equivalent
within-class gradient in this project's sources, so inventing one just to have
*something* to weight by would be worse than admitting none exists. Vivianite's
wetness-proxy weighting is a **lightweight stand-in**, not the real thing — it's a
placeholder for the fuller basin-fill Arable/Wetland-backswamp/Estuarine zonation
`scenario_reference.md` S22.4 already proposes but hasn't built as a raster; revisit
vivianite's placement once that lands. Silver/copper only placed 6 of the 8 targeted —
the volcanic zone is small (595.7 km²) and the 5 km separation constraint genuinely
limits how many non-overlapping pods fit; each placed pod was independently rolled
silver (rarer) vs copper (common) at an unreviewed 25%/75% placeholder split, no
citation for that specific ratio (1 silver, 5 copper this run).

**A real citation problem caught while doing this, not yet fixed in `resources.py`
itself:** placer magnetite's existing citation says titanomagnetite is "reconcentrated
by rivers" from the volcanic zone into basin fill. That doesn't hold up geographically —
mainland basin_fill and the volcanic zone are two separate landmasses, 16.44 km of open
water apart (Tappa 7 S1), with no shared river catchment for fluvial reworking to
happen through. Real NZ ironsand beaches (Taranaki, Westport) actually form by
**coastal/longshore redistribution** of eroded volcanic material, a mechanism that DOES
work across a water gap. This pass's spatial weighting already uses that corrected
mechanism (coastal × volcanic-landmass proximity, not any river field) — but
`resources.py`'s own citation *text* still says "rivers," flagged in both files
pending Nico's sign-off on the wording change, not silently rewritten.

**Plaster, resolved without new geology — worth recording since it constrains the
limestone decision above.** This world has no other landmass to trade with, so gypsum
(which real NZ doesn't produce domestically either — it's an evaporite mineral needing an
arid enclosed marine basin this world's wet, stream-threaded climate never had) can't be
imported the way real NZ does. Resolution reached: clay plaster (basin-fill clay + plant
fiber, zero new geology, real and ancient — predates mud brick in some regions) as the
baseline, with lime plaster as a second option *if* limestone is added (independently
ancient tradition, not a gypsum substitute). A third, more speculative option surfaced:
laumontite (already greywacke's Matéria-domain crystal) is a real zeolite, and zeolites
are documented pozzolanic lime additives — but laumontite specifically is a genuinely
unstable mineral (dehydrates, crumbles to powder at surface conditions), so raw laumontite
is a poor binder as-is; this would need to be a processed/calcined refinement, not a raw
resource, if used at all.

**Iron — currently confined to a small footprint, real fix identified but not NZ-specific.**
Titanomagnetite is currently the only iron source, and it's tied entirely to the volcanic
zone (595.7 km², ~6% of total land) plus a "weak secondary" reworked placer form in
basin fill. Real fix under consideration: bog iron (goethite/limonite, precipitated by
iron-oxidizing processes in anoxic, organic-rich wetland/floodplain sediment) is a
well-documented global pre-industrial ore source (Iron Age Scandinavia through colonial
North America) and is environmentally the same setting vivianite already occupies in this
project (basin fill's own primary Vértice material, also an authigenic anoxic-floodplain
mineral) — so it would extend iron into basin fill's much larger footprint (4437 km²)
without contradicting anything already locked. Flagging honestly: no NZ-specific citation
found for bog iron — Te Ara's NZ iron history is exclusively the ironsand story. This
would be cited as a general real-world process, not an NZ-attested one, unlike nearly
everything else in this project's resource layer.

**Aluminium — genuinely unresolved, two different problems that look like one.** Checked
directly: New Zealand does have real bauxite, small relict Pliocene-Pleistocene lateritic
deposits (Northland, largest ~20 Mt at Otoroa/Matauri Bay), but they're explicitly
subeconomic — too small/low-grade, never mined. So "no accessible aluminium ore" is
defensible either way: a tiny "known but unworkable" pocket is a real option if a
narrative hook is wanted (something characters know exists but can't practically use),
or omitting it entirely is equally defensible. Separately — and this may be the more
useful framing — aluminium metal was harder to extract than gold or silver for all of
human history until the 1886 Hall-Héroult electrolytic process (Napoleon III reserved
aluminium tableware for his most honored guests over gold/silver; the 1884 Washington
Monument capstone was aluminium specifically because it was the era's most precious
displayable metal). If this society doesn't have cheap industrial-scale electrolysis,
aluminium's absence may be a **technology gate, not a resource gate** — worth deciding
which problem is actually being solved (no ore exists vs. no one can extract it yet)
before picking a geological fix, since they have different answers.

### 8f. Transport lithology-cost multiplier — a friction layer on Tappa 6's cost graph

Nico's own framing for this item was "a weight cost to calculate how this lithology
could help/difficult the **construction** of the transport system." Read literally that's
excavation cost, foundation stability, tunnel/bridge risk — none of which is modelled
here, or modellable from anything currently in this project (no unit-cost data for any
of it exists anywhere in the pipeline). Scoped down with Nico beforehand (AskUserQuestion,
"just a lithology cost multiplier") to the nearest thing the EXISTING cost graph can
actually represent: a **travel-TIME friction weight** on Tappa 6's Tobler-hiking-function
+ Dijkstra cost-distance graph (`cost_distance.py`, used for Circulo isochrones/tier
separation, not path/network design — no predecessor extraction, no biome-differentiated
costs, no rail grade ceiling; that larger scope was explicitly declined as "starting
Tappa 9"). Worth being explicit about the gap between what was asked and what got built —
this answers "how much slower is travel across this rock," not "how much would a road
here cost to build."

**Implementation, additive only, nothing already-locked touched:**
- `geomorphology/lithology.py` gains `LAND_TRAVEL_FRICTION` (a per-class multiplier
  dict) and `travel_friction_multiplier(lithology, ...)` (maps a class-code raster to a
  float multiplier field, unmapped codes — i.e. ocean/coastline mismatch — default to
  1.0 neutral, not an error).
- `suitability/cost_distance.py`'s `build_cost_graph()` gains an optional
  `friction_multiplier` parameter (a per-cell array; edge friction = mean of the edge's
  two endpoint cells). Passing `None` (the default) reproduces the exact prior behaviour
  — every already-committed Tappa 6 result stays reproducible untouched. Applied to
  LAND-LAND edges only; sea/boat edges are unaffected (lithology isn't defined offshore).
- `suitability/terrain_metrics.py` gains `block_mode()` — this pipeline's first
  CATEGORICAL block-reduce (existing `block_mean`/`block_max`/`block_any` are for
  continuous/boolean data). Majority-vote per 4×4 block via `n_classes` one-hot
  `block_mean` calls + `argmax`, needed to bring `lithology_v6` (native 30 m, 8 classes)
  down to Tappa 6's 120 m grid without averaging class CODES into meaningless
  fractional values.

**The friction table itself (`LAND_TRAVEL_FRICTION`) — UNREVIEWED first-pass estimates,
pending Nico's sign-off, not written to `config/parameters.yml`:**

| class | multiplier | grounding |
|---|---|---|
| sedimentary_basin_fill | 1.00 | baseline/anchor — Tobler's function is itself calibrated on ordinary, unobstructed ground |
| greywacke_argillite | 0.85 | NZ high-country off-track travel (scree, tussock, jointed rock) — real accounts describe a material time penalty vs. formed track, no single canonical number |
| schist | 0.85 | same value as greywacke on purpose — nothing differentiates the two specifically |
| volcanic | 0.70 | rugged basalt; deliberately NOT worst-case aa-lava (real accounts put that under 1 km/h) since this pipeline doesn't distinguish aa/pahoehoe |
| marble | 0.60 | karst dissolution terrain — same mechanism `karst_cave_candidates()` already uses |
| sedimentary_limestone | 0.60 | same as marble, same reasoning |
| granite | 0.90 | mild penalty — real granite slab/tor terrain is often GOOD footing; tiny footprint (13.5 km²) anyway, low-leverage |

None of these are a citation of a specific per-lithology hiking-speed figure — no such
dataset was found or is claimed to exist. They're directional judgement calls sized so
the two classes covering most of the map (greywacke+schist, ~4770 km² of ~9880 km² land)
get the mildest, not the harshest, penalty — an aggressive value there would have an
outsized, poorly-justified effect on the whole network.

**Result, run via `run_tappa8_transport_friction.py`** (block-mode 120 m class areas:
basin_fill 4339.0, greywacke 3007.9, schist 1763.7, volcanic 593.1, sedimentary_limestone
86.9, marble 77.2, granite 13.5 km² — minor differences from the native-30 m S8b table are
the expected majority-vote effect at block boundaries, not a bug; 0.19% of effective-land
120 m cells had no matching class due to `land_mask`/`lithology_v6` coastline disagreement
and fell back to neutral):

Rather than silently re-running Tappa 6's actual site placement, this built a SEPARATE
friction-adjusted cost graph alongside the original baseline graph and directly compared
pairwise cost-distance (hours) between the 17 **already-placed, locked** Circulo sites
under both — `circulo_candidate_sites.geojson` and `tappa6_site_selection_meta.json` were
read-only inputs, never modified. Sanity check passed first: this script's own baseline
graph reproduced Tappa 6's committed smallest tier margin exactly (0.1665 h both ways).

Of the 64 tier-constrained pairs re-checked: mean travel-time increase +3.06%, median
+2.28%, max +13.74% (`Circulo_B_35k` ↔ `Circulo_E2_2k`, a large-medium pair); 8 of 64
pairs showed exactly zero change — their shortest path is sea/boat-dominated, which this
friction layer never touches. **0 pairs would now violate their tier's hour threshold.**

**Important honest caveat about that last number — it is not real evidence the friction
values are "safely small."** Every multiplier in `LAND_TRAVEL_FRICTION` is ≤ 1.0 (slower
only, nothing sped up), so cost-distance HOURS between any two points can only increase
or stay the same once friction is applied — it can mathematically never decrease. That
means an already-satisfied `>= threshold` constraint can only get MORE satisfied, never
less: "0 new violations" was guaranteed by construction the moment every value was chosen
≤ 1.0, regardless of which specific numbers were picked. It genuinely does confirm the
graph-building code is wired correctly (friction is doing real, measurable work — up to
~14% on some pairs) and that this experiment doesn't retroactively break anything, but it
is not a validation that the multiplier magnitudes themselves are reasonable — that
judgement is still entirely Nico's to make from the table above, independent of this
check.

**What this does NOT do, on purpose:** re-place the 17 Circulos under this friction (a
materially bigger, separate decision — `place_circulos()` re-run from scratch with
`friction_multiplier` wired through `cost_distance_fn` could move sites, not just change
measured margins, since it would also affect which candidate cells satisfy `valid` at
placement time); model actual construction/excavation cost (see the framing note above);
or touch the sea/boat portion of the cost graph at all.

### 8g. Excavation effort — a second, deliberately separate lithology-cost axis

Direct follow-up to S8f. Nico clarified his original "construction" framing wasn't
asking for a full engineering-cost model (excavation crew-hours, blasting budgets,
haulage) — just a relative, per-class number, "how easier it is to work on basins than
schist or granite." That's a genuinely different physical property from S8f's travel
friction (how fast can you cross this ground vs. how hard is it to dig into/build on it),
so this is a **separate** table and raster, not an extension of `LAND_TRAVEL_FRICTION`.

**Implementation, additive only, no cost-graph involvement at all** (unlike S8f, this
needed no Dijkstra/DEM/lake/boat-speed machinery — excavation effort is a pure per-cell
function of lithology class, not a property of an edge between two cells):
`geomorphology/lithology.py` gains `EXCAVATION_EFFORT_MULTIPLIER` (a per-class ratio,
`>=1.0`, relative to `basin_fill=1.0`) and `excavation_effort_multiplier(lithology, ...)`.
`run_tappa8_excavation_effort.py` computes it at the same 120 m grid as S8f purely for
easy QGIS overlay (block-mode downsample of `lithology_v6`, reusing the S8f `block_mode`
function) — nothing about excavation effort actually required that resolution.

**The table itself — UNREVIEWED first-pass estimates, pending Nico's sign-off, not
written to `config/parameters.yml`:**

| class | multiplier | grounding |
|---|---|---|
| sedimentary_basin_fill | 1.0 | baseline/anchor — unconsolidated alluvium, diggable by hand/basic tools |
| sedimentary_limestone | 1.3 | soft carbonate rock — this project's own Oamaru-stone citation (S8) is literally famous for being soft/easy to carve |
| marble | 1.6 | same calcite mineral as limestone, but DIFFERENT value here (unlike S8f, where they share one) — metamorphic recrystallization makes marble's fabric denser and measurably harder to carve than ordinary limestone |
| greywacke_argillite | 2.3 | quartz-rich, well-indurated bedrock |
| schist | 2.3 | same value as greywacke on purpose — nothing differentiates them specifically; real asymmetry NOT modeled: schist splits easily along foliation planes (Otago schist flagstone is a real tradition) but is harder to cut across them — a single scalar can't represent that directionality |
| volcanic | 2.6 | basalt — dense, tough, the standard "needs blasting" rippability-chart case |
| granite | 3.0 | hardest value — coarse crystalline, no foliation weakness at all, the canonical hardest-to-excavate common rock type |

Grounding is relative mineral hardness (calcite ~Mohs 3 for limestone/marble vs.
quartz/feldspar ~Mohs 6-7 for greywacke/schist/granite/basalt) plus a genuinely
well-documented, if non-NZ-specific, historical fact: ancient Egyptian/Mediterranean
quarrying cut limestone with simple copper/bronze tools, while granite needed much
harder abrasives/dolerite pounders — not a specific numeric engineering source, an
ordinal judgement call same as S8f's.

**Important caveat, same spirit as the placer-magnetite citation flag in S8e — surfaced,
not glossed over:** this index does NOT capture karst's foundation-**stability** risk
(sinkholes/voids under marble and sedimentary_limestone, see S8c). Taken alone, it makes
marble/limestone look like easy building ground — true for the quarrying/cutting half of
the question, not for "is this safe ground to found a structure on." A fuller model
would need those as two separate fields; this stage only built the one Nico asked for.

**The interesting part — this table's ranking deliberately INVERTS S8f's in two places,
because the two axes measure different things:**
- **Granite**: mild travel penalty (0.90, good footing on slab) but the *hardest*
  excavation value here (3.0, no foliation weakness) — easy to walk across, brutal to
  quarry into.
- **Marble/limestone**: the *worst* travel penalty in S8f (0.60, karst hazard) but among
  the *easiest* excavation values here (1.3-1.6, soft calcite) — hard to walk across,
  comparatively easy to cut once you're there.

Area-weighted mean excavation effort across all land: 1.734 — pulled toward the low end
by basin_fill+greywacke+schist covering ~92% of land; granite/volcanic's harsher values
(2.6-3.0) only cover a small footprint (607 km² combined).

**Not wired into anything yet, on purpose** — no Vértice-materials linkage, no Urban
Scale build-time estimate, no interaction with S8f's cost graph. A pure reference layer
for now, matching the "just a relative value" scope Nico asked for.

## 9. Basin fill internal zonation — Arable / Wetland-backswamp / Estuarine-coastal

`scenario_reference.md` S22.4 authored a three-way sub-zonation of basin_fill (flagged
there as doing too many jobs at once: presumed farmland, vivianite's source, and home to
two new candidate resources) and explicitly left "exact area/threshold percentiles...
for Tappa 8 to calibrate against the actual rasters" — that section's own words. This is
that calibration, run via `run_tappa8_basin_zonation.py`. **The logic (which fields,
coastal-claims-priority ordering) is S22.4's, not authored here** — only the specific
threshold values are new.

**Implementation — derived from a DEM/distance field, not hand-authored, on purpose**,
the same "DEM-native beats line/polygon-distance" reasoning S0-S2 already established for
lithology v5 itself, applied here to resolve `scenario_reference.md` S22.19's explicit
open question ("whether Tappa 8 actually authors it as real geometry"): answer is no,
it's classified from two fields Tappa 8's pipeline already computes elsewhere, no new
spatial primitive.
- `geomorphology/terrain_relief.py`'s existing `compute_local_relief` (the SAME
  `relief_2km` field lithology v5's `classify_from_terrain` already uses for
  schist/greywacke/basin_fill) — drainage proxy: higher relief = better-drained = Arable,
  lower = poorly-drained = Wetland-backswamp.
- `geomorphology/caves.py`'s existing `distance_to_ocean_km` (the SAME field sea caves
  and placer magnetite's coastal weighting already use) — coastal proximity: below the
  threshold, regardless of relief, = Estuarine/coastal margin.
- New: `geomorphology/basin_zonation.py`'s `classify_basin_fill_zones()`. Classification
  order matters and mirrors S22.4's own prose exactly: the coastal check runs FIRST and
  claims cells unconditionally (S22.4 defines wetland as "poorly-drained, low local
  relief, INLAND (not coastal-adjacent)" — that phrasing presupposes the coastal band is
  already carved out); the relief split (Arable vs. Wetland) only applies to whatever
  basin_fill remains after that.

**Thresholds actually used — UNREVIEWED first-pass, calibrated against basin_fill's own
data distribution, pending Nico's sign-off, not written to `config/parameters.yml`:**
basin_fill's own `dist_to_ocean_km`/`relief_2km` percentiles were checked directly before
picking anything (p10/p25/p50/p75/p90 = 0.98/2.72/6.41/11.73/16.47 km and
42.0/57.5/96.8/182.5/235.4 m respectively). `coastal_threshold_km=1.5` (an ABSOLUTE
distance — unlike relief, "how far does estuarine/tidal influence plausibly reach" is a
physically meaningful km quantity, not a percentile; scaled up from the 0.5 km
point-feature buffers karst/talus/sea-caves use, since this defines a broad zone, not a
cave-formation trigger — a first attempt at 3.0 km left Arable at only 54.6%, not the
"clear majority" S22.4 specifies, so it was tightened). `wetland_relief_percentile=25.0`
(bottom quartile of the NON-coastal population's own relief distribution — not
mainland's, since basin_fill was already excluded from schist/greywacke's higher-relief
thresholds, so re-percentiling within basin_fill itself is the locally-consistent choice).

**Result**: Arable 2770.3 km² (63.9%), Wetland-backswamp 923.4 km² (21.3%), Estuarine-
coastal 640.8 km² (14.8%) — Arable is a clear majority, both other classes stay
substantial rather than vanishing, matching S22.4's framing without either number being
independently derived beyond that one sanity check.

**Bonus, unrelated open item answered while the inputs were already in hand**:
`scenario_reference.md` S22.3/S22.19 separately flagged "basin fill / Grassland-biome
spatial overlap... a hypothesis, not yet checked against the rasters." Checked directly
here (120 m biome grid, basin_fill block-mode downsampled to match): 38.3% of basin_fill's
area falls inside the Grassland biome, and 56.2% of Grassland falls inside basin_fill —
a real, substantial, but partial overlap. Confirms the hypothesis is directionally right
without being a near-total match — worth the Scenario chat updating S22.19 to reflect
this is now checked, not still open.

**Not done in this script, done in S10 instead**: vivianite's spatial placement (S8e)
still used its original wetness-proxy stand-in (`1/(1+dist_to_stream_km)`) over the whole
basin_fill mask at the time this script ran, explicitly called a "lightweight stand-in"
pending this zonation landing. Re-pointing it at `wetland_backswamp` specifically would
change `resource_vivianite.npy`'s already-committed output, so it wasn't done unasked
inside this script — it's the first thing S10 does, on Nico's explicit go-ahead.

## 10. Iron and aluminium — resolving S8e's two open items (Nico's explicit direction)

S8e flagged both as genuinely open ("Iron — currently confined to a small footprint, real
fix identified but not NZ-specific" / "Aluminium — genuinely unresolved, two different
problems that look like one"). Nico picked concrete resolutions for both; this section
implements them (`run_tappa8_iron_aluminium.py`) — **the decisions are Nico's, only the
implementation is this stage's.**

**Iron — bog iron, co-located with vivianite, vivianite re-placed into `wetland_backswamp`
first.** Nico's call: implement fully, including moving vivianite. Two steps:
1. **Vivianite re-placement**: S8e's original placement used the whole `basin_fill` mask
   (4334.5 km²) because S9's sub-zonation didn't exist yet. Re-run with the identical
   seed and weight field (`1/(1+dist_to_stream_km)`), only the eligible mask changed, now
   restricted to S9's `wetland_backswamp` (923.4 km², 21.3% of basin_fill). Still placed
   8/8 target pods (the shrunk mask remained large enough) — pod area 7.57 km²,
   `resource_vivianite.npy` **overwritten**.
2. **Bog iron**: goethite/limonite, precipitated by iron-oxidizing processes in the same
   anoxic, organic-rich floodplain/bog setting vivianite occupies. Per
   `scenario_reference.md` S22.4's own framing — "a single wetland patch plausibly yields
   both, worked by different specialists" — bog iron's pod footprint is the **identical**
   geometry to vivianite's, not an independent stochastic draw: `resource_bog_iron.npy` ==
   `resource_vivianite.npy` by construction. `resources.py`: `mundane_only` (no citable
   electrical/magnetic/optical gating property, same category as schist's gold/jade).
   **Citation-honesty flag, carried over from S8e verbatim**: this is a well-documented
   GENERAL pre-industrial process (Iron Age Scandinavia through colonial North America),
   but NOT NZ-specific — Te Ara's NZ iron history is exclusively the ironsand/
   titanomagnetite story. First exception in this project's resource layer to the
   NZ-specific-citation norm every other material has met.

**Aluminium — bauxite added as an active resource, technology gate closed narratively.**
Nico's exact framing: *"Adicionar recurso. Eletrólise está em desenvolvimento com auxílio
de Vértices"* — treat the real (subeconomic-in-reality) Northland bauxite as workable
here, with the historical Hall-Héroult (1886) electrolysis gate closed by Vértice-assisted
extraction rather than left as an unresolved ore-vs-technology ambiguity.
- **Class**: `volcanic`, not a new class — the real citation (Otoroa/Matauri Bay,
  Northland) is lateritic weathering of the Tangihua Complex's basaltic rock specifically,
  matching this world's volcanic zone, not granite or any other class.
- **Spatial weighting**: flatness (`1/(1+slope_pct)`) — real lateritic bauxite caps form
  and survive on low-relief ground and get stripped off steep slopes by erosion. Reuses
  the existing slope field, same direction of logic `caves.py` already uses for lava
  tubes (gentle-slope preference), not a new spatial primitive.
- **Deliberately small footprint**: 2 pods, 200–500 m radius range (vs the other five
  materials' 8 pods / 300–800 m) — matches the real citation's own framing (small, ~20 Mt
  at the largest real deposit, historically subeconomic, never mined). Pod area 0.69 km²
  within 595.7 km² eligible volcanic zone. Vértice-assisted electrolysis closes the
  **technology** gate; it does not make the ore itself abundant — kept as two separate
  constraints, not conflated.
- **Domain**: `mundane_only` — the ore itself has no citable Vértice-domain-gating
  property. The Vértice involvement Nico described is at the **extraction/process**
  level (assisted electrolysis), not a property of the mineral. That process-level
  mechanic — how, exactly, Vértices assist an electrolytic process — is a Scenario-chat
  question if further mechanical detail is wanted; recorded here only as citation context
  for why the resource is narratively viable, not elaborated.
- **Not resolved here, flagged**: how RARE/hard-won this bauxite should feel in actual
  play (the real deposit was subeconomic — does that framing carry over narratively, or
  does Vértice-assisted extraction make it straightforwardly workable?) is a tone
  question, not a geology one.

**Flagged for consistency, not changed in this pass** — `resources.py`'s placer magnetite
(basin_fill `secondary_weak`) still uses the whole `basin_fill` mask weighted by coastal ×
volcanic-landmass proximity (S8e), not restricted to S9's `estuarine_coastal` sub-class,
even though that sub-class's own definition is coastal-proximity-driven and would be a
natural fit. Out of this task's explicit scope (iron, aluminium only) — a real follow-up,
not forgotten. **Resolved in S12** (same session, later still) — see that section.

## 11. All seven resource-pod materials, blended into one bitmask raster

Nico asked to blend the resource layers the same way S8d blended the five cave types —
one file, pixel value carries which material(s) are present, instead of seven files.
Same convention exactly (`scripts/make_tappa8_resource_blend.py`, uint8, powers-of-two
bits, OR-packed): **1=jade, 2=laumontite, 4=vivianite, 8=bog_iron, 16=placer_magnetite,
32=silver_copper, 64=bauxite** — seven of the eight named Vértice/mundane materials that
exist in `resources.py`. Excluded, same reasoning as S8e/S10: volcanic's magnetite
(primary) is a disseminated bulk mineral through the whole basalt class, not
vein/joint/pod-localized like the other seven — no raster exists for it, so no bit.
Schist's gold-bearing quartz veins and mica also get no separate bit — `resources.py`'s
own `spatial_note` says they co-locate with jade, so a separate bit would double-count
jade's own footprint. `jade_pods_v5.npy` (S5, closed stage) was staged fresh from the
device repo for this — it hadn't been touched since, unlike the six S8e/S10 rasters.

**Checked directly rather than assumed: unlike the five cave types, overlap here is
structurally CONSTRAINED, not open.** Each material's eligible mask is drawn from
exactly one lithology class (jade/schist, laumontite/greywacke,
vivianite+bog_iron+placer_magnetite/basin_fill, silver_copper+bauxite/volcanic), and
lithology classes are mutually exclusive by construction (`lithology_v6` assigns one
code per cell) — so materials from *different* classes cannot spatially overlap, full
stop, not just empirically. The script asserts this directly (`cross_class_overlap_found
== False`) rather than taking it on faith, and it holds. Within a shared class, overlap
is only *possible*, not guaranteed — checked, not assumed: **the only nonzero pairwise
overlap found is `vivianite+bog_iron`, at 7.5663 km², exactly equal to both materials'
own individual area** — i.e. bog_iron's bit is a byte-for-byte structural duplicate of
vivianite's (S10's co-location design working exactly as intended), asserted directly in
the script (`np.array_equal(masks["vivianite"], masks["bog_iron"])`). Vivianite and
placer_magnetite (different weight fields, different seeds, same basin_fill population)
turned out **not** to overlap at all, and neither did silver_copper and bauxite within
volcanic — genuinely incidental non-overlap, not designed. Max simultaneous materials on
any one cell: 2 (only from the vivianite/bog_iron pair).

**Per-material area, as of this run**: jade 6.00, laumontite 8.32, vivianite 7.57, bog_iron
7.57, placer_magnetite 4.74, silver_copper 5.10, bauxite 0.69 km² — any-resource-material
union 32.41 km² (less than the simple sum, 39.98 km², solely because of the vivianite/
bog_iron duplication, not because of any other pod overlap). **Superseded by S12**:
placer_magnetite was re-placed into `estuarine_coastal` immediately after this ran (same
session), changing its own area to 5.16 km² and the union to 32.83 km² — see S12 for the
re-run and why the overlap picture itself doesn't change.

**Round-trip verified in-script** (an assertion, not just a claim): decoding all seven
bits back out of `resource_blend.npy` reproduces every one of the seven source masks
exactly (`np.array_equal` per material, all seven pass). Committed as `resource_blend.*`
/ `tappa8_resource_blend_meta.json` — the seven individual rasters (plus
`jade_pods_v5.npy`, unchanged from S5) stay committed too; this is an additional
convenience layer, same status as `cave_blend.npy`, not a replacement.

**Honest limitation worth flagging, not smoothed over**: bog_iron's bit adds zero new
spatial information over vivianite's — anyone decoding this raster bit-by-bit should
know bit 8 is redundant with bit 4 before building anything on top of it (e.g. counting
"how many distinct material types" per cell from `n_types_per_cell` would silently
double-count the same physical deposit as two "materials" unless this is accounted for).
Kept as its own bit anyway for catalogue completeness — `resources.py` lists bog_iron as
its own named material, and dropping its bit here would make this raster disagree with
that source of truth — but it's flagged here and in the meta JSON, not hidden.

## 12. Placer magnetite re-placement + resource blend update (closing S10/S11's last flag)

Nico's request: fix placer magnetite (the one consistency item flagged, unresolved, in
both S10 and S11 — "still uses the whole `basin_fill` mask... even though [S9's
`estuarine_coastal` sub-class] would be a natural fit"), then update `resource_blend.npy`
(S11) to reflect it. Both done together (`run_tappa8_placer_magnetite_restrict.py`, then
re-running `scripts/make_tappa8_resource_blend.py`).

**Re-placement**: same treatment S10 already gave vivianite — SAME seed (133) and SAME
weight field (coastal × volcanic-landmass proximity, unchanged from S8e), only the
eligible mask changed: from S8e's whole `basin_fill` (4334.5 km²) to S9's
`estuarine_coastal` sub-class (640.8 km²). This isn't introducing a new criterion —
placer magnetite's weight field was already a coastal-proximity signal; restricting it to
the sub-class that's *authoritatively* coastal (S9's own `coastal_threshold_km` test)
tightens an existing bias rather than changing what the placement is trying to represent.
Still placed 8/8 target pods despite the much smaller eligible ground — pod area grew
slightly, 5.16 km² (was 4.74 km² under S8e's placement), because the same weight field
now concentrates entirely within a zone where it was already scoring highest.
`resources.py`'s citation text is unchanged by this — only the `spatial` note was
updated to record the restriction; the `citation_flagged_for_review` note about
"reconcentrated by rivers" still stands, unrelated to this fix.

**Structural non-overlap, checked directly, not assumed**: S9's three `basin_fill`
sub-zones (arable / wetland_backswamp / estuarine_coastal) are mutually exclusive by
construction (`classify_basin_fill_zones` assigns exactly one code per cell). Vivianite/
bog_iron are restricted to `wetland_backswamp` (S10); placer magnetite is now restricted
to `estuarine_coastal`. These eligible masks cannot overlap, so the pods drawn from them
can't either — asserted directly in the script (0 overlapping cells found), not inferred
from the sub-zone logic alone.

**Resource blend (S11) re-run**: `resource_blend.npy`'s bit 16 (placer_magnetite) updated
to the new geometry; the other six bits are unchanged (their own source rasters weren't
touched). Round-trip re-verified for all seven bits. Diagnostics after the update: still
only one nonzero pairwise overlap (`vivianite+bog_iron`, unchanged, 7.5663 km²) —
placer_magnetite doesn't gain or lose any overlap with anything, consistent with the
structural non-overlap argument above (it had zero overlap under the old mask too, so
this wasn't expected to change anything about the blend's overlap picture, only its own
area/shape). Any-resource-material union grew slightly to 32.83 km² (was 32.41 km²),
tracking placer_magnetite's own small area increase.

## 13. Quartz/mica/gold separated from jade; resource blend re-run at int16

Feedback relayed from the Scenario chat (Nico: "mica (Energia), quartz (Onda), and gold
(mundane) are all currently bundled as a single 'co-locates with jade' note with no
independent raster — separating them into their own `place_material_pods()`").

**Why the bundling existed and why separating it isn't a contradiction**: S8e explicitly
excluded schist's gold-quartz-mica assemblage from the pod-placement generalization,
reasoning that giving them independent pods "would contradict" the existing note that
they co-locate with jade's high-grade band. That reasoning conflates two different
claims — "these minerals are found within the same high-grade schist zone as jade" (a
real geological fact, unchanged by this work) and "their discoverable footprints must be
pixel-identical to jade's" (not implied by the first claim, and arguably LESS honest —
jade's own `place_jade_pods` docstring already argues real deposits aren't a smooth
function of a suitability zone: "if [deposits] were [smooth], prospecting would be
trivial"). A real vein system is heterogeneous at the pod scale — quartz, mica, and gold
concentrate unevenly within the same overall high-grade zone, same as jade's gold-vein
neighbors do in the real Otago Schist. Three independent stochastic draws within the SAME
eligible zone is the more geologically honest representation, not a break from the
existing note.

**Implementation** (`run_tappa8_schist_vein_materials.py`): eligible mask and weight
field are UNCHANGED from jade's own established machinery — `jade_eligible_mask`
(80th-percentile `schist_grade` threshold) + `schist_grade` itself, not a new criterion.
Generalizes `place_material_pods()` to quartz/mica/gold the same way S8e generalized it
to laumontite/vivianite/etc: same 8-pod/5km-separation/300-800m-radius placeholders as
that family, three new distinct seeds (130+10/+11/+12, chosen to avoid colliding with
S8e/S10/S12's existing +1..+6 offsets). Eligible ground: 354.6 km². Results: quartz 8/8
pods (4.63 km²), mica 8/8 pods (5.67 km²), gold 8/8 pods (5.23 km²). Checked directly:
zero pairwise overlap among the three this run — incidental (all three draw from the same
zone/weight field with only the seed differing), not designed that way.

**`schist_grade_v5.npy` regenerated, verified before trusting it**: this field wasn't in
the data commit (S6's "continuous float32 fields... dropped, cheap to regenerate" list).
Regenerated by re-running `run_tappa8_lithology_v5.py` in full — entirely deterministic
(connected-component landmass ID, windowed max-min relief, percentile-threshold
classification; no RNG until jade's own seeded pod placement at the very end) — and
verified byte-for-byte identical to the already-committed `lithology_v5.npy` AND
`jade_pods_v5.npy` before trusting the regenerated `schist_grade_v5.npy` for anything.
Both matched exactly. Not re-committed (same "cheap to regenerate" policy as before).

**One deliberate deviation from jade's own placement, flagged explicitly, not
reconciled**: `jade_eligible_mask` was recomputed here against `lithology_v6` (current,
authoritative) rather than `lithology_v5` (what `jade_pods_v5.npy` itself still uses,
unchanged since S5). Checked directly: 9,013 cells (8.11 km², ~0.46% of v5's schist area)
that were schist in v5 became marble in v6's priority-tier compositing (marble's
`priority_rank=1` outranks schist). Using v6 avoids seeding a new mineral pod on ground
that, per the CURRENT lithology, is no longer schist at all — but it means quartz/mica/
gold's eligible zone and jade's own aren't computed identically anymore. `jade_pods_v5.npy`
carries this staleness unchanged; not fixed here (would mean re-running jade's own locked
placement, out of this task's scope, a decision for Nico if it matters enough to revisit).
Real, small consequence surfaced by the resource-blend re-run below: jade's pods overlap
each of quartz/mica/gold's by a small amount (0.12/0.08/0.35 km²) — not zero, unlike
quartz/mica/gold's mutual non-overlap — directly attributable to this v5-vs-v6 divergence
plus the independent seeding.

**`resources.py` updated**: schist's `primary`/`secondary` entries gained `spatial`
fields (`resource_quartz.npy`, `resource_mica.npy`); `mundane_only`'s gold entry gained
one (`resource_gold.npy`); the old `spatial_note` (which said all three "reuse
jade_pods_v5.npy directly") was rewritten to reflect the new independent-raster status
and to carry the v5-vs-v6 flag forward into the source-of-truth file, not just this doc.

**Resource blend (S11) re-run, dtype changed int16 from uint8, flagged not silent**:
adding three more bits pushed the bitmask past uint8's 8-bit/255-value capacity (bit 512
alone exceeds it) — switched `resource_blend.npy`/`.bin` to ENVI `"i2"` (signed int16,
already-supported dtype, no new code needed in `raster_io.py`). Ten bits now: `1=jade,
2=laumontite, 4=vivianite, 8=bog_iron, 16=placer_magnetite, 32=silver_copper, 64=bauxite,
128=quartz, 256=mica, 512=gold` — max possible value 1023, comfortably inside int16's
positive range. Round-trip re-verified for all ten. Two structural checks re-asserted and
still holding: cross-lithology-class overlap still impossible (confirmed False), and
placer_magnetite/vivianite+bog_iron overlap still exactly 0.0 km² (S12's fix holds).
Any-resource-material union: 47.81 km² (up from S12's 32.83 km², reflecting the three new
materials' own footprints plus jade's small overlaps with them).

## 14. Sedimentary limestone (calcite) and granite (mica/quartz) resource pods; marble's explicit null

Closes the item S13's own follow-ups flagged as blocked: Nico relayed the Scenario chat's
Vértice material assignments for the three v6 rock classes, pulled verbatim from the now-
locked `scenario_reference.md` §21. This is a **sync, not a new decision** — the content
was already settled in the Scenario chat; the gap was only that it hadn't been carried
into this chat/`resources.py` yet.

**Assignments (verbatim from Nico's relay, `scenario_reference.md` §21):**

- **`sedimentary_limestone` → Onda, secondary.** Material: optical-grade calcite
  ("Iceland spar" — clear calcite in limestone cavities, not the bulk rock). Does not
  replace quartz/schist as Onda's primary; a rarer upgrade alongside it. Citation:
  calcite's birefringence (Δn≈0.17) is ~20x stronger than quartz's (Δn≈0.009) — the real
  physical basis of the Iceland-spar/Viking-sunstone effect — but clear, intact crystals
  large enough to use are much rarer in limestone cavities than vein quartz is in schist.
- **`granite` → Energia, secondary (mica) + Onda, tertiary (quartz).** Granite is felsic
  igneous rock naturally composed of quartz+feldspar+mica — the SAME two minerals already
  gating Energia/Onda via schist, just a second, unrelated host rock. No new physical-
  property claim; diversified sourcing only, so a Círculo doesn't need to sit on schist
  specifically to reach either domain.
- **`marble` → NO Vértice domain. Explicit null, not a gap.** Marble stays purely mundane
  (construction/monumental stone). Reasoning: recrystallization destroys the parallel-
  lattice clarity limestone's calcite needs for its birefringence — marble is optically
  cloudy/light-scattering, not transparent. Separately, calcite's crystal structure is
  centrosymmetric, so it was never piezoelectric either — it wouldn't have qualified for
  Energia by that route regardless. `VERTICE_MATERIALS['marble']` is zeroed out
  deliberately, per Nico's explicit instruction not to treat it as a placeholder.

**Implementation (`run_tappa8_limestone_granite_materials.py`, new)**: same
`place_material_pods()` machinery as S13, applied to `lithology_v6 ==
CLASS_SEDIMENTARY_LIMESTONE` (87.72 km², spanning the authored "North Coast Limestone" +
"Sedimentary Bay" zones — the raster doesn't distinguish which authored zone a cell came
from, same as every other class) and `lithology_v6 == CLASS_GRANITE` (13.84 km²,
"Granite South" alone). UNIFORM weight for all three (calcite, mica_granite,
quartz_granite) — no citable within-class spatial gradient exists for either host rock,
same reasoning laumontite's placement already used for greywacke.

**Deliberately NOT rescaling pod count/radius for calcite's "much rarer" framing**: the
handoff text explains WHY calcite is secondary rather than primary (a tier decision,
already encoded), not an instruction to shrink its footprint. Schist's own mica (S13,
"rarer than the quartz itself") got the identical `n_pods`/`radius_range_m` as quartz and
gold in that class — tier already carries relative scarcity in-fiction; pod count/
footprint has not been conflated with tier elsewhere in this project (the same separation
S8g draws between excavation effort and S8f's friction). Same standard placeholders used
here: `n_pods=8, min_separation_km=5.0, radius_range_m=(300, 800)`.

**Real results, checked directly, not assumed**:

| material | host | eligible km² | pods placed | pod area km² |
|---|---|---|---|---|
| calcite | sedimentary_limestone | 87.72 | 7/8 | 6.38 |
| mica_granite | granite | 13.84 | 3/8 | 1.33 |
| quartz_granite | granite | 13.84 | 2/8 | 1.45 |

**Flag worth surfacing, not silently absorbed**: granite's eligible ground (13.84 km²) is
~43x smaller than volcanic's (595.7 km², which itself only fit 6/8 silver_copper pods
under this same rule, S8e). At granite's scale, the project's standard 5 km
`min_separation_km` leaves room for only 2–3 pods rather than 8 — this is a geometric
consequence of applying one constant across wildly different zone sizes, not a bug. Worth
deciding whether granite specifically warrants a smaller separation constant if having a
comparable pod count to the other materials matters for play; not changed unilaterally
here. mica_granite vs quartz_granite overlap checked directly: 0.0 km² (same tiny zone,
independent seeds — could plausibly have been nonzero given how little room 5 km leaves,
so this was worth checking rather than assuming either outcome).

**`resources.py` updated**: new `sedimentary_limestone` entry (`secondary`, calcite);
new `granite` entry (`secondary` mica_granite + `tertiary` quartz_granite — `tertiary` is
a new tier key, not used elsewhere in this catalogue before now); new `marble` entry
(`no_vertice_domain: True`, no `spatial` raster, citation for the null recorded in full).
Module docstring updated to reflect thirteen materials with rasters, two deliberate
exceptions (volcanic magnetite, marble).

**Resource blend (§11) re-run a third time, three more bits, dtype UNCHANGED (int16
already had headroom)**: `1024=calcite, 2048=mica_granite, 4096=quartz_granite` — 13 bits
total, max possible value 8191, comfortably inside int16's 32767 ceiling, so no dtype
change needed this round (unlike S13's uint8→int16 change). All structural checks
re-asserted and still holding: cross-lithology-class overlap still False (confirmed —
`sedimentary_limestone`/`granite` are two more mutually-exclusive classes added to the
existing set), `vivianite`+`bog_iron` still exact match, `placer_magnetite` vs
`vivianite`/`bog_iron` still 0.0 km². All 13 materials round-trip exactly. Any-resource-
material union: 56.98 km² (up from S13's 47.81 km², reflecting calcite/mica_granite/
quartz_granite's own footprints — no new overlap introduced, since sedimentary_limestone
and granite are lithology classes distinct from every other material's host class).

**Shipped**: `run_tappa8_limestone_granite_materials.py` (new), `scripts/
make_tappa8_resource_blend.py` (updated), `resources.py` (updated), this decision doc,
`resource_calcite.*`, `resource_mica_granite.*`, `resource_quartz_granite.*`,
`resource_blend.*` (overwritten again), `tappa8_limestone_granite_materials_meta.json`,
`tappa8_resource_blend_meta.json` (updated).

**Closes the last item that was open specifically because of a missing decision from
another chat.** What's left for Tappa 8 (§15 below) is now purely Tappa-8-internal.

## 15. Open follow-ups (not done this stage, deliberately left open)

- Schist fracture caves — still explicitly parked (07 doc S2), not built.
- Dangerous creatures/conflict zones, dangerous seas — two of the three Tappa 7 S9 domains
  folded into "Tappa 8," still not started. **Transportation, the third, is DECOUPLED from
  this list (2026-08-20, Nico's explicit call)**: it will NOT be a Tappa 8 add-on — it's
  Tappa 9's full scope (path/network design, biome-differentiated cost, rail grade
  ceiling), a separate future chat, not folded in here.
- No independent calibration exists yet for any of the percentile thresholds in S7 above
  beyond the checks already narrated (SE-plains fix, crest-line recovery, island split,
  lava tube footprint/threshold effect).
- **RESOLVED, not open anymore (2026-08-20, Nico's explicit call): the 17 Círculos will
  NOT be re-placed under S8f's transport friction multiplier.** The friction layer stays
  exactly what S8f already said it was — a measured-but-unapplied reference layer on
  `cost_distance.py`'s graph (opt-in via `friction_multiplier`, `None` by default) — not a
  half-finished decision anymore. Site coordinates from Tappa 6 remain final and current
  for all downstream Tappas, same status `scenario_reference.md` §22.1 already gave them
  for narrative relocation.
- Basin-fill area zonation is now IMPLEMENTED (S9); vivianite's re-placement into
  `wetland_backswamp`, bog iron's co-located addition (S10), and placer magnetite's
  re-placement into `estuarine_coastal` (S12) are now ALL IMPLEMENTED — every material
  S8e originally placed on the whole `basin_fill` mask is now restricted to its correct
  S9 sub-zone. Nothing left open in this family.
- Aluminium is now IMPLEMENTED (S10, bauxite/volcanic/`mundane_only`) — the remaining
  follow-up is a tone question, not a geology one: how rare/hard-won this bauxite should
  feel in play, given the real-world deposit it's based on was subeconomic and never
  mined. Also open: the Vértice-assisted-electrolysis mechanic itself (how, mechanically)
  is only recorded as citation context here — a Scenario-chat question if more detail is
  wanted.
- **RESOLVED (2026-08-20, S14): Vértice material assignments for marble /
  sedimentary_limestone / granite.** `scenario_reference.md` §21 locked `sedimentary_
  limestone → Onda secondary (calcite)`, `granite → Energia secondary (mica) + Onda
  tertiary (quartz)`, `marble → no Vértice domain (explicit null)` — relayed verbatim by
  Nico from the Scenario chat and implemented the same session: `resources.py` updated,
  `resource_calcite.npy`/`resource_mica_granite.npy`/`resource_quartz_granite.npy` placed
  (7/8, 3/8, 2/8 pods respectively — granite's small 13.84 km² zone limits pod count under
  the standard 5 km separation rule, flagged not silently accepted), resource blend
  re-run at 13 bits. Nothing left open in this family. The only genuinely new follow-up
  from this closure: whether granite specifically warrants a smaller
  `min_separation_km` than the project's standard 5 km, given how few pods fit there —
  a tuning question, not a blocked decision.
- Standing-rule suggestion from the Scenario chat (2026-08-20, relayed by Nico, not yet
  acted on): whenever a `VERTICE_MATERIALS`-relevant assignment gets locked in
  `scenario_reference.md`, add a one-line note to the Tappa 8 handoff doc
  (`docs/reference/tappa8_new_rocks_vertice_handoff.md`) at the same time, rather than
  relying on Nico to manually relay it into this chat later — would have caught S14's
  sync gap before it needed a question round-trip. Worth Nico's call whether to adopt this
  as a standing convention for future chat-boundary handoffs, not decided unilaterally
  here.
