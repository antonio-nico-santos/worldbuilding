# Tappa 6 — Solarpunk settlement suitability

Status: **closed**. A weighted multi-criteria suitability model at the 120 m
climate/biome working resolution, producing two separate composite surfaces
(Círculo / open-field-agriculture economy vs. Povo Livre / forestry-foraging
economy), then a greedy, settlement-size-aware site-selection pass on the
Círculo surface that picks and validates 17 candidate locations. Criteria
were deliberately left undefined at the end of pre-project planning
(`00_pre_project_planning.md`) so they could be set from the user's own
scenario/narrative context rather than generic defaults — this document
records that context as it was actually supplied, session by session.

As with the five prior documents, this is written to work either as a log
of what was decided or as a reusable recipe.

## 0. Architecture, locked before any layer was built

Two kinds of layer, kept conceptually separate:

* **núcleo** layers — continuous 0-1 desirability, combined by a per-
  population WEIGHTED SUM into a composite index. Six were built: slope,
  water (distance-to-stream), agriculture (TWI proxy), solar exposure,
  biome, and (informationally, not weighted into the sum) settlement size
  itself, read off downstream rather than as its own criterion.
* **exclusão** layers — 0-1 MULTIPLIERS applied to the composite, not
  additive terms. Only one was built this stage: Povo Silencioso
  archipelago exclusion. A second (dangerous-creature territory) was
  planned but reclassified before being built — see S9.

Architectural style (Terrapedra/Trançada/Salina/Caldária/Modular) and
settlement-size footprint fit are explicitly POST-HOC lookups against
whatever site gets chosen, not their own weighted criteria — this was
locked in during pre-project planning and held throughout the stage.

Two separate composite surfaces were built from the same five núcleo
layers, not one shared index — see S5 for why.

## 1. Núcleo layer 1-4: slope, water, agriculture, Povo Silencioso distance

Built together in `src/suitability/terrain_metrics.py` +
`run_tappa6_core_metrics.py`, at 120 m (matching the climate/biome grid,
not the DEM/hydrology native 30 m — this is a siting decision, not a
per-building layout).

**Slope**: `np.gradient` at native 30 m, then block-**MAX** (not block-mean)
downsampled to 120 m — block-mean underestimated land mean slope by ~44%
(16.76% vs. 29.69%), the wrong direction of error for a hazard criterion; a
locally steep sub-cell must not disappear inside an otherwise-gentle 120 m
cell. `slope_suitability`: 1.0 at/below 5% grade, smootherstep decay to 0.0
at 30% — placeholder generic land-use thresholds, not calibrated to this
world, flagged as too permissive for the largest Círculos and too strict for
the smallest (the downstream settlement-size filter, S6-8, is meant to
tighten/loosen the *effective* requirement per size class).

**Water**: Euclidean distance-to-stream (km), then `water_suitability` — 1.0
at/below 0.5 km, decay to 0.0 by 5 km, thresholds set from this world's own
percentiles (0.5 km ≈ where "trivially close" stops covering nearly
everyone; 5 km ≈ p99.5). **Comes out nearly flat in practice** (land mean
0.9907, p5 0.9981) because this world's stream network is dense (82.5% of
land within 0.5 km of a stream, 97.7% within 1 km) — by design, not a bug,
but this criterion barely discriminates except at the margin.

**Agriculture (TWI proxy)**: `TWI = ln(specific catchment area / tan(slope))`
from Tappa 4's contributing-area raster (block-mean downsampled) and the
same slope field, normalized to the 2nd/98th land percentiles. Explicitly a
geomorphic proxy — no real soil/pedology data (texture, depth, drainage
class, pH) exists anywhere in this project; it ranks wet, low-relief,
well-drained-looking land higher, nothing more.

**Povo Silencioso distance**: Euclidean distance (km) to the nearest cell of
the NE archipelago (5 connected-component labels identified once by hand
against `land_mask`, not re-derived geometrically each run — a bounding-box
heuristic over-collects unrelated islets elsewhere). `povo_silencioso_
exclusion_factor`: 0.0 within a 5 km hard buffer, smootherstep ramp to 1.0
(no penalty) by 15 km. A localized corner-of-the-map effect by construction
— only 1.02% of non-archipelago land (~100 km²) falls within 10 km of the
archipelago at all. Buffer sizes are placeholders; no treaty/lore distance
was ever specified for this stage.

## 2. Núcleo layer 5: solar exposure

`src/suitability/solar.py` + `run_tappa6_solar.py` — built from scratch (no
`r.sun`/GRASS available in this sandbox), FAO-56 elevation-corrected
clear-sky Rso scaled by a numerically-integrated tilted/horizontal ratio
(16-direction, 10 km horizon shading, 0.5 hr time steps, 12 mid-month
representative days), diffuse = 15% clear-sky fraction × slope+terrain
sky-view-factor.

Annual land-mean insolation 7475.4 MJ/m²/yr (p5-p95: 6219.8-8595.2).
**Validated finding**: at this world's -44° latitude, north-facing slopes
receive ~37% more annual insolation than south-facing slopes at matched
moderate grade (8578 vs. 6257 MJ/m²/yr on 8-25° slopes) — the
southern-hemisphere-correct inverse of the northern-hemisphere intuition,
confirmed rather than assumed. Honest limitations, documented in the module
and not fixed this stage: clear-sky only (no cloud climatology), the 15%
diffuse fraction is an uncalibrated placeholder, and 120 m/16-direction
horizon shading is coarse relative to real terrain silhouettes.

## 3. Núcleo layer 6: biome suitability

`src/suitability/biome_lookup.py` + `run_tappa6_biome_suitability.py` — a
**hand-scored** 0-1 lookup by `biome_id`, not a derived physical field like
layers 1-5, preceded by a conservative 3×3 majority filter on `biome_id`
itself (ties keep the original class; changed only 0.10% of land cells in
practice, since most fragmented patches are compact ~4-cell blobs, not
single-pixel speckle).

| biome | score |
|---|---|
| Permanent Snow & Ice | 0.00 |
| Alpine Fellfield | 0.05 |
| Alpine Tundra | 0.20 |
| Subalpine Dry Scrub | 0.30 |
| Subalpine Woodland | 0.45 |
| Subalpine Wet Forest | 0.55 |
| Temperate Forest | 0.85 |
| Woodland / Shrubland | 0.90 |
| Lowland Steppe / Grassland | 1.00 |

Two explicit value decisions, both the user's, both checked against real
data before being finalized rather than taken on faith:

* **Economic base = open-field agriculture** → Grassland > Woodland/
  Shrubland > Temperate Forest. These three classes cover 58.5% of all
  land, so this single ordering dominates the layer's effect more than the
  other six classes combined.
* **Permanent Snow & Ice / Alpine Fellfield stay núcleo-only**, not also an
  exclusão multiplier. I initially assumed slope would zero these out on
  its own; checked directly and found that wrong — 41% of Permanent Snow &
  Ice cells and 21% of Alpine Fellfield cells have `slope_suitability_120m`
  above 0.3 (flat ice-cap/plateau terrain exists). Told the user the
  assumption was wrong before they made the final call; they chose to
  accept the trade-off (a flat, sunny, well-watered patch of permanent ice
  can still pull a moderate composite score) rather than add a second
  exclusion.

## 4. A second biome variant, for the Povo Livre

Forestry/foraging economy, not open-field agriculture — mirrors only the
Forest/Grassland endpoints (Temperate Forest 1.00, Woodland/Shrubland
unchanged at 0.90, Lowland Steppe/Grassland 0.85). Before building this, I
checked whether the fork was geographically real or cosmetic: Temperate
Forest averages 4936 mm/yr precipitation vs. Grassland's 531 mm/yr (9.3×),
the two classes sit on opposite sides of the domain (mean X -17.5 km vs.
+19.4 km) with **0% spatial overlap**, and 36.9% of land changes score by up
to ±0.15 between the two versions. Confirmed real; built as
`BIOME_SUITABILITY_POVO_LIVRE` in the same lookup module +
`run_tappa6_biome_suitability_povo_livre.py` (reuses the same
majority-filtered `biome_id`, does not recompute smoothing). The other six
biome classes (ice/rock/tundra/subalpine) are left unexamined and identical
between the two variants — not verified to fit a foraging economy, flagged
as such rather than silently assumed correct.

## 5. Composite: two separate weighted indices

`src/suitability/composite.py` + `run_tappa6_composite.py`. Two composites,
not one shared index, because the two populations' economies are different
enough that a single weighting would misrepresent one of them:

| criterion | Círculo | Povo Livre |
|---|---|---|
| slope | 0.20 | 0.10 |
| water | 0.05 | 0.05 |
| agriculture (TWI) | 0.25 | 0.10 |
| solar | 0.20 | 0.15 |
| biome | 0.30 | 0.60 |

Povo Livre reasoning: slope reduced not zeroed (temporary camps tolerate
more relief, but rugged terrain still matters for mobility/safety);
agriculture reduced sharply (foraging/hunting economy, wetter land still
means more biomass/game, so not zero); solar reduced not zeroed (portable
panels remove the fixed-orientation dependency, but horizon shading still
applies regardless of portability, and the annual-average formulation is a
known, undocumented-fix temporal mismatch for a seasonally-mobile
population — would need an actual migration-route model); biome dominant
at 0.60 (vegetation type is close to the whole story for a foraging
economy). Water stays low (0.05) for both and barely matters either way,
per S1.

Povo Silencioso exclusion is multiplied into **both** composites — a
default choice (the archipelago is framed as the Povo Silencioso's own
territory to respect generally, not a Círculo-specific concern), not
explicitly confirmed by the user. Flagged as revisit-if-wrong, not
revisited this stage.

Validated the two composites actually favor different regions, not just
different scores in the same place: top-decile mean X +3.1 km (Círculo,
central/east) vs. -22.4 km (Povo Livre, west/wetter side), 50.5% overlap
between the two top-10% sets. Land means: Círculo 0.560 (p5/p50/p95 =
0.227/0.565/0.829), Povo Livre 0.611 (0.178/0.696/0.851).

## 6. Site selection: method

`src/suitability/site_selection.py` + `run_tappa6_site_selection.py` — the
settlement-size downstream filter promised in S0. Population → required
footprint area via an assumed density (**4000 ppl/km²**, an adjustable
placeholder — "dense but green, walkable, with internal agriculture/parks",
denser than car-dependent suburbia but far short of a packed historic
core) → candidate cells screened by MEAN suitability over a **square**
window sized to that footprint (not a true circular disk — Tappa 7/8's
"urban zoom" does the actual layout; window area ≈1.27× the true disk
area, an accepted coarseness) → placement **greedy, largest population
first** (most area-constrained), irrevocable once chosen — not
jointly-optimal, standard practice for this class of problem, transparent
about which choice happened first.

17 Círculos: 4 large (40k/35k/25k/20k), 5 medium (2k each), 8 small. The 8
smallest only had a combined 5,000-person total specified in the scenario,
not individual sizes — split evenly (625 each), a simplifying assumption.
For those smallest Círculos the required window is only ~4×4 120 m cells
(~0.23 km² against a ~0.156 km² target) — the grid is genuinely coarse
relative to their size; their reported sites are indicative, not precise.

A candidate window must be ≥98% land to be considered at all (keeps sites
off the coast where the window would otherwise average in ocean cells).

## 7. Site selection: minimum-distance mechanisms, three rounds

This is where most of the stage's iteration happened — three real design
corrections, each caught either by this project's before-delivery
clean-mirror testing discipline or by the user inspecting the actual
output, in that order.

**Round 1 — straight-line, uniform, buggy.** First attempt: each Círculo's
own `min_distance_km` (initially planned to mirror the user's "greek
city-state" idea, imagined at 50-150 km between large/medium Círculos)
checked against EVERY already-placed site regardless of tier. This
collapsed placement from a working 17/17 (no minimum at all) to 5/17: a
60 km-radius exclusion disk is ~11,310 km², **larger than this world's
entire land area** (~9,904 km²). Caught by this project's standard
before-delivery clean-mirror test (one site landed at the domain's far NE
corner with ≈0 suitability, the rest failed to place).

**Round 2 — straight-line, tier-pair-scoped, working.** Fixed by scoping
minimum distance to explicit same-tier PAIRS only (`frozenset({tier}) ->
km`), with cross-tier or untiered pairs falling back to ordinary
footprint-buffer spacing. Feasibility checked directly against this map's
actual geometry (domain diagonal only 206 km, not the scale the "greek
city-state" idea originally imagined): 60 km alone comfortably feasible for
the 4 large Círculos against each other (up to 70 km still works); 25 km
comfortably feasible for the 5 medium Círculos (tested to 40 km); a
UNIFORM minimum across all 9 large+medium Círculos together is NOT feasible
beyond ~40 km — the 50-150 km original vision only becomes achievable here
by scoping to same-tier pairs. Result: 17/17 placed, 0 violations, smallest
margin 0.02 km. **But**: all 17 sites landed in the same biome (Lowland
Steppe/Grassland) — pure suitability-maximization funnels every settlement
into the single highest-scoring class, flattening the scenario's intended
per-region architecture-style diversity. Flagged, not silently fixed;
offered two paths (accept it, or add an explicit biome-diversity
constraint). The user chose neither — see Round 3.

**Round 3 — cost-distance (walking/boat TIME, hours), the user's own
alternative fix.** Rather than force biome diversity directly, the user
proposed that mountains and open sea should count for MORE effective
separation than the same straight-line km over flat terrain, and open flat
ground for LESS — testing whether that alone would break up the
single-biome clustering. Built `src/suitability/cost_distance.py`: an
8-connected DIRECTED graph over the whole domain (land AND sea, ~1.45M
nodes, ~11.6M edges), land-land edges costed by **Tobler's hiking function**
(`speed_kmh = 6·exp(-3.5·|slope_ratio+0.05|)`, ~5.04 km/h flat, signed
slope from a 120 m block-MEAN DEM — mean, not the block-max used for
`slope_pct_120m`, because this wants a realistic average edge slope, not a
hazard-conservative worst case), any edge touching non-land at a flat
**BOAT_SPEED_KMH = 6.0** (the user's explicit choice — boat travel, not an
impassable barrier; a placeholder for "simple rowing/small-sail boats, no
motors"). Solved via `scipy.sparse.csgraph.dijkstra` (fast in practice:
graph build ~5 s, ~0.2-0.5 s per already-placed site). Two independent
minimum-distance mechanisms now coexist: straight-line footprint
non-overlap (always enforced, a physical constraint) and cost-distance
tier rules in HOURS (same-tier pairs, feasibility-swept the same way as
Round 2: large-large 12h, ceiling ~15h; medium-medium 8h, ceiling ~12-14h;
small-small 2.5h, the user's own midpoint of "2-3 horas de caminhada",
never found a ceiling in the tested range).

**Important side effect, flagged to the user, not fixed by their choice**:
`BOAT_SPEED_KMH=6.0` is slightly FASTER than Tobler's flat-ground pace
(~5.04 km/h), and open water carries no slope penalty — so, as
parameterised, sea is not really a barrier in this model, closer to a
slight shortcut vs. walking real terrain. Only mountains still act as a
genuine separator under these settings. Result: biome distribution improved
to 12 Grassland / 4 Woodland-Shrubland / 1 Temperate Forest (from the
Round 2 result of effectively 14/17-Grassland) — the user's hypothesis was
right, worth testing. One large Círculo (Circulo_D) landed on the SW
"Island" (label 135, ~794 km², used for the Caldária volcanic zone),
directly confirming islands are live, competitive candidates, not token
inclusions, per the user's explicit ask to check.

## 8. Site selection: round 4, two problems the user found by inspection

The user caught two more real problems by looking at the Round 3 output
directly, not by re-running an abstract feasibility sweep:

**Lakes.** The user spotted Círculos sitting inside lakes and estimated ~4.
Checked precisely: **7 of 17**, not 4 — 2 of them large Círculos, with
53-95% of their own site window covered by lake water. Root cause:
`land_mask` only excludes ocean, never inland lakes; `lake_mask.npy`
already existed in the repo from Tappa 4 hydrology but had never been wired
into Tappa 6. Fixed by building an EFFECTIVE land mask
(`land_mask & ~lake_mask_120m`, `lake_mask` block-ANY downsampled to 120 m,
the same hazard-conservative convention `stream_mask` already uses) and
substituting it everywhere in site selection — window screening, the cost
graph, everything, not just a post-hoc filter. Result: 0/17 on a lake.

**Large-medium proximity.** The user noticed a medium Círculo only ~4 km
from a large one — cross-tier pairs had no minimum at all beyond ordinary
footprint spacing. Added `frozenset({"large","medium"}): 6.0h` to
`TIER_MIN_HOURS`, feasibility-swept together with everything else (17/17 up
to 10h; 11h+ drops to 13/17) — picked 6h, deliberately between small-small
(2.5h) and medium-medium (8h): a medium town should feel more independent
from a nearby metropolis than from another village, but not as separated as
two mediums need from each other. Nearest large-medium pair moved from
~4 km / <1h to **29 km / 6.19h**.

Checked, but explicitly NOT acted on: the user's separate observation that
increasing `BOAT_SPEED_KMH` (to represent launch/dock preparation time, not
just transit) probably would not have changed much, since few of the 17
sites sit directly on the coast. Verified directly: closest site to the
ocean is <1 km, most sit 2-17 km inland. Confirmed correct; parameter left
unchanged.

## 9. Dangerous-creature territory: reclassified, not built

Originally scoped (pre-project planning, and again mid-Tappa-6) as a
possible exclusão multiplier. The user corrected this after Tappa 6 closed:
creature territories are a **conflict layer**, informational, used to help
define architectural/cultural character per site (the same role
`biome_at_site` already plays in the site-selection GeoJSON output), not a
suitability penalty requiring a calibrated weight. This is a simpler
engineering commitment than an exclusion multiplier — a lookup, not a
number that needs justifying — but still blocked on the user defining
creature ecology/distribution logic. Not built this stage; deferred to a
future stage (see S11).

## 10. Final locked results

17/17 Círculos placed, 0 minimum-distance violations, 0 on a lake.

| name | tier | population | x (km) | y (km) | mean suitability | biome | location |
|---|---|---|---|---|---|---|---|
| Circulo_A_40k | large | 40,000 | 38.0 | -45.6 | 0.837 | Lowland Steppe/Grassland | mainland |
| Circulo_B_35k | large | 35,000 | 19.3 | 35.3 | 0.829 | Lowland Steppe/Grassland | mainland |
| Circulo_C_25k | large | 25,000 | -40.8 | 49.5 | 0.797 | Temperate Forest | mainland |
| Circulo_D_20k | large | 20,000 | -31.1 | -56.2 | 0.779 | Lowland Steppe/Grassland | **SW Island** |
| Circulo_E1_2k | medium | 2,000 | -12.2 | 54.8 | 0.851 | Lowland Steppe/Grassland | mainland |
| Circulo_E2_2k | medium | 2,000 | -36.2 | 15.9 | 0.845 | Woodland/Shrubland | mainland |
| Circulo_E3_2k | medium | 2,000 | 7.0 | -15.3 | 0.806 | Temperate Forest | mainland |
| Circulo_E4_2k | medium | 2,000 | 44.9 | -16.1 | 0.783 | Lowland Steppe/Grassland | mainland |
| Circulo_E5_2k | medium | 2,000 | 50.6 | 59.8 | 0.762 | Lowland Steppe/Grassland | mainland |
| Circulo_F1-F8_small | small | 625 each | (8 sites) | | 0.861-0.887 | mostly Grassland, 1 Woodland/Shrubland | mainland |

Biome distribution across all 17: 13 Lowland Steppe/Grassland, 2 Temperate
Forest, 2 Woodland/Shrubland. Not a fully even spread — grassland still
dominates, consistent with it also dominating this world's land area
generally (S3) — but a real improvement over Round 2's near-total
single-biome clustering, and a direct, verified result of the user's own
cost-distance idea plus the lake/cross-tier fixes.

## 11. Locked-in parameters

```
working_resolution_m           = 120
slope_gentle_pct                = 5.0
slope_hard_limit_pct            = 30.0
water_gentle_km                 = 0.5
water_hard_limit_km             = 5.0
povo_silencioso_hard_buffer_km  = 5.0
povo_silencioso_soft_buffer_km  = 15.0
weights_circulo    = {slope: 0.20, water: 0.05, agriculture: 0.25, solar: 0.20, biome: 0.30}
weights_povo_livre = {slope: 0.10, water: 0.05, agriculture: 0.10, solar: 0.15, biome: 0.60}
biome_suitability_circulo    = {snow_ice: 0.00, fellfield: 0.05, tundra: 0.20,
                                 dry_scrub: 0.30, woodland_sub: 0.45, wet_forest: 0.55,
                                 temperate_forest: 0.85, woodland_shrub: 0.90, grassland: 1.00}
biome_suitability_povo_livre = same as above except temperate_forest: 1.00, grassland: 0.85
density_ppl_km2       = 4000.0
buffer_factor         = 2.0     # footprint non-overlap, straight-line km
min_land_fraction     = 0.98
tobler_flat_speed_kmh = 5.04    # derived from the Tobler formula, not itself a free parameter
boat_speed_kmh        = 6.0
tier_min_hours = {large-large: 12.0, medium-medium: 8.0, small-small: 2.5, large-medium: 6.0}
lake_exclusion = "land_mask & ~lake_mask_120m (block-ANY downsample), used for ALL site-selection computations"
```

Actual run on disk (`data/processed/suitability/tappa6_site_selection_meta.json`):
17/17 placed, 0 violations, 0 on lake, smallest footprint margin 2.30 km,
smallest tier margin 0.17h.

## 12. Outputs

`run_tappa6_core_metrics.py`, `run_tappa6_solar.py`,
`run_tappa6_biome_suitability.py` (+`_povo_livre`), `run_tappa6_composite.py`,
`run_tappa6_site_selection.py` — five scripts, run in that order, writing to
`data/processed/suitability/` (gitignored, regenerate locally):

| file | contents |
|---|---|
| `slope_pct_120m` / `slope_suitability_120m` | float32 |
| `dist_to_stream_km_120m` / `water_suitability_120m` | float32 |
| `twi_120m` / `agriculture_suitability_120m` | float32 |
| `dist_to_povo_silencioso_km_120m` / `povo_silencioso_exclusion_120m` | float32 |
| `annual_insolation_MJm2_120m`, `june_insolation_MJm2_120m`, `aspect_deg_120m`, `horizon_mean_deg_120m`, `solar_suitability_annual_120m`, `solar_suitability_june_120m` | float32 |
| `biome_id_smoothed_120m` | int16 |
| `biome_suitability_120m` / `biome_suitability_povo_livre_120m` | float32 |
| `suitability_circulo_120m` / `suitability_povo_livre_120m` | float32, the final composites |
| `circulo_candidate_sites.geojson` | 17 Point features: name, population, tier, radius_km, mean_suitability, land_fraction, biome_at_site |
| `circulo_claimed_footprint_120m` | int16, sanity-check visual only (each site's straight-line footprint-buffer circle) — NOT the tier-rule exclusion zone, see S7's cost-distance isochrone note below |
| `tappa6_core_metrics_meta.json`, `tappa6_solar_meta.json`, `tappa6_biome_meta.json`, `tappa6_biome_povo_livre_meta.json`, `tappa6_composite_meta.json`, `tappa6_site_selection_meta.json` | parameters, stats, decisions, verification results, one per script |

Each `.npy` also ships as ENVI `.bin`+`.hdr`+`.prj`, same convention as
every prior stage. Same CRS caveat as Tappa 1-5: QGIS needs "Fictional
World LCC" assigned manually if a layer loads with an unknown CRS.

**Visualization note**: shading each tiered site's own cost-distance
isochrone (all cells reachable within its tier's hour threshold) into
`circulo_claimed_footprint_120m` was tried and reverted — a single large
Círculo's 12h isochrone alone covers ~24% of all land, and the union across
all 17 sites covers 99.998% of the entire grid, useless as a "what's still
open" picture. The delivered raster only shades each site's small
straight-line footprint-buffer circle. This does not weaken the actual
tier-rule enforcement, which is a per-pair lookup, verified correct
directly (S10) independent of this visualization choice.

## 13. Open follow-ups (not done this stage, deliberately left open)

* **Dangerous-creature conflict layer** (S9) — reclassified from exclusion
  to informational/cultural lookup; still blocked on the user defining
  creature ecology/distribution logic. Deferred to a future macro-scenario
  stage, not scheduled.
* **Sea-danger refinement to the cost-distance model** — a natural fix for
  the `BOAT_SPEED_KMH ≈ Tobler flat pace` side effect (S7): a spatially
  varying boat cost (safe passages cheap, dangerous ones very expensive)
  instead of one flat constant. Not built; a candidate piece of a future
  macro-scenario stage alongside road/rail routing.
* **Road/rail regional infrastructure** — not scoped this stage at all,
  raised only after Tappa 6 closed. `cost_distance.py`'s graph (Tobler
  friction + Dijkstra over the whole domain) is directly reusable for
  least-cost-PATH road routing between the 17 Círculos (needs predecessor
  extraction, not just cost values, which the current code doesn't do).
  Rail is a genuinely different problem, not a relabeling of the same
  friction function — real rail needs a hard grade ceiling (~2-4%
  sustained, without heavy tunneling/switchback investment), not a smooth
  Tobler-style exponential decay; would need its own cost function.
* **Povo Silencioso exclusion applied to both composites** (S5) — a default
  choice, never explicitly confirmed by the user. Revisit if it turns out
  the archipelago should only matter for one population.
* **`suitability_povo_livre_120m` has no site-selection pass** — only the
  Círculo composite went through S6-8's greedy placement; the Povo Livre
  surface exists but was never used to place actual nomadic-camp
  candidates. Not requested this stage.
* **Cost-distance graph does not model launch/dock overhead** — checked and
  confirmed low-impact for THIS run (S8), but the underlying simplification
  (any land↔sea edge just uses boat speed, no fixed preparation cost)
  remains in the code as-is.
* **Pipeline renumbering, SUPERSEDED — see `07_tappa7_regional_scenario.md` §9 for the current
  numbering.** (Originally: the planned "Tappa 7 (urban zoom)" was renumbered **Tappa 8**, with a
  new **Tappa 7 (regional scenario deepening)** inserted before it. Tappa 7 has since run and closed
  — see its doc — and its own roadmap planning found this Tappa-8-as-urban-zoom assignment collided
  with a later plan to use "Tappa 8" for geomorphology. Resolved by retiring urban zoom from the
  Tappa sequence entirely: it's now its own untethered **Urban Scale** track, not a numbered Tappa,
  which frees Tappa 8/9/10 for geomorphology/transports/interactions respectively.) Urban Scale work,
  when resumed, is scoped as scenario-depth for the TTRPG side of the project primarily, not
  additional portfolio material — intentionally lighter empirical-rigor
  standard than Tappa 6 used (no expectation of a feasibility sweep for
  every threshold), per the user's own steer.
