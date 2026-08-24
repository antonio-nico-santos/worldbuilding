# Tappa 10 — auxiliary settlement siting (mines, mountain huts, coastal villages, forest posts)

Fourth pass. Feature names are English-consistent (`Mine_<Class>_NN`, `Coastal_Village_NN`,
`Forest_Post_NN`, plus Tappa 7's own `Outpost_MainSpine_N`/`Outpost_SouthBranch_N` for huts,
kept as-is). Output: `auxiliary_settlements.geojson` (108 features, same CRS as
`circulo_candidate_sites.geojson`) + `auxiliary_settlements.csv`. All coordinates are real,
computed from the project's own rasters/vectors — nothing hand-placed, though several coastal
villages are now sited by explicit direction rather than pure algorithm (see below).

**This pass fixes two things you reported: an almost-empty attribute table, and
Coastal_Village_08 sitting on the wrong side of a lake.** Both are covered in their own
sections below (right after this one, and inside Coastal villages).

## Attribute schema fix (why the table looked almost empty)

Two real problems, both now fixed, same root cause: the four settlement types were each
carrying only their OWN property keys (mining posts: ~17 fields; huts: ~14 different fields;
villages: ~8; forest posts: ~8), and three of those fields — `attached_circulos`, `resources`,
`shares_structure_with` — were JSON list values, not plain strings.

Combine those two and most GIS viewers (QGIS's attribute table included) build their column
set as the UNION of every feature's keys across the whole file. Open a coastal-village row and
you'd see real values in ~8 columns and blank `NULL`s in the ~28 columns that belonged to
mining posts, huts, or forest posts — which reads exactly like "the table is almost empty,"
even though nothing was actually missing for that feature. List-valued fields make it worse:
GeoJSON arrays become OGR "StringList" fields, which several QGIS/GDAL version combinations
render blank in the attribute table widget rather than as a readable joined string.

**Fix**: every one of the 108 features now carries the exact same 37 property keys (verified —
one distinct key-set across the whole file), with `null`/blank where a key genuinely doesn't
apply to that settlement type (e.g. a coastal village has no `lithology_class`) — that's
expected and correct, not a gap. `attached_circulos`, `resources`, and `shares_structure_with`
are now semicolon-joined strings (`"Circulo_F2_small;Circulo_B_35k;Circulo_F1_small"`) instead
of JSON arrays, in both the GeoJSON and the CSV. I also rebuilt the CSV writer around a single
shared dict per feature (Python's `DictWriter`) instead of hand-positioned columns — the old
script had a real bug where mountain-hut rows silently shifted `massif`/`status_prior` text
into the `resources` column and left `note` empty; that's gone now, every column lines up with
its header for every row.

**Follow-up, same day: you reported the SAME kind of misalignment again after this fix**
(Círculo name showing under `lithology_class`, `resources` empty with its data under
`resource_note`). I checked both files this fix produced byte-for-byte, field by field, for
every feature type — printed each property key next to its value directly from the actual
GeoJSON and CSV I generated, e.g. `Mine_BasinFill_01`: `attached_circulos: "Circulo_A_40k"`,
`lithology_class: "basin_fill"`, `resources: "vivianite;bog_iron"`,
`resource_note: "bog_iron: (co-located w/ vivianite)"` — every field lines up with its own
name correctly. **The files themselves are not misaligned.** The most likely explanation is on
the viewing side: this exact file path has now been overwritten 4+ times with materially
different property SCHEMAS (different key sets and orders across passes — the very first
version used `materials` where this version uses `resources`+`resource_note`, for instance).
GDAL's GeoJSON reader (which QGIS uses under the hood) can cache a file's field schema in a
sidecar `.gfs` file the first time it's opened, and/or QGIS itself can retain a layer's
original field mapping across a simple refresh rather than re-detecting it — either one would
explain values landing under the wrong header after several schema-changing overwrites at the
same path, without the underlying file actually being wrong. **What I did**: delivered this
round's files under a new name (`auxiliary_settlements_v2.geojson`/`.csv`) specifically to
rule out any stale-cache collision with the old path. **What's worth doing on your end**: if
you still see misalignment with the `_v2` files, check the same folder for a
`auxiliary_settlements_tappa10.geojson.gfs` file (delete it if present — it's a GDAL cache
artifact, not part of the dataset) and, in QGIS, fully remove the old layer and add the new
file fresh rather than pointing the existing layer at the new path or using Refresh — a
straight reload can still reuse a layer's original field mapping in some QGIS versions.

## Summary

| type | count | population | source data |
|---|---|---|---|
| Mining posts | 71 (49 structure-bearing: 41 standalone + 8 shared hubs; 12 day-commute, no structure; 10 satellites sharing a hub's structure) | none (workers counted at home Círculo) | resource-pod rasters (siting) × real cost-distance hours over the project's own Tobler/friction/Dijkstra graph (structure tier) |
| Mountain huts | 18 | none | Tappa 7's `outpost_candidates.geojson` (fauna-suitability siting), Círculo-attached here |
| Coastal villages | 12 | own population | lithology/DEM ocean-connectivity check × slope raster, 4 sites adjusted per Nico's direct review |
| Forest posts | 7 (serving 9 Círculos) | own population | basin raster × biome raster × `streams.geojson` |

**Flagging before anything else: removing the three villages drops Circulo_E5_2k to ZERO
attached auxiliary settlements.** Circulo_D_20k keeps 13 mining posts + 1 forest post, and
Circulo_F3_small keeps 1 mining post, so losing their dedicated coastal village doesn't leave
them bare. Circulo_E5_2k has no mining pod nearby and no forest post (its basin is one of the 3
`basin_labeled==0` unresolved ones, see Forest posts below) — its coastal village was its ONLY
auxiliary settlement. If it really has no dock-worthy activity beyond what Urban Scale gives it
directly, that may be fine narratively (a small Círculo can just be a small Círculo), but it's
worth a deliberate call rather than a side effect of this edit — I did not add anything back for
it unless you say so.

## Method and three things worth flagging before this goes further

**1. The resource-pod rasters don't store pod centers as data — I decoded them.** Only
`laumontite` and `silver_copper` had `pod_centers_xy` in `tappa8_resource_pods_meta.json`; every
other material's coordinates (jade, quartz, mica, gold, vivianite, bog_iron, placer_magnetite,
bauxite, calcite, mica_granite, quartz_granite) had to be extracted from the `.npy` masks via
connected-component labeling. Raw labeling over-counts: a single placed pod's circular footprint
gets fragmented into 2-4 pixel-blobs wherever it's clipped by a class boundary (a stream, a zone
edge). I merged fragments within 2.5 km before computing centroids — after merging, every
material's pod count and minimum pairwise spacing matches its own meta JSON exactly (e.g.
laumontite: 12 raw fragments → 8 pods, min spacing 6.07 km, matching `n_pods_placed: 8` and
`min_separation_km: 5.0`). I'm confident in these coordinates, but they're derived, not read off a
field — if the pipeline's own `place_material_pods()` output ever changes, this decode needs
re-running, it won't self-update.

**2. Mining-post consolidation is my judgment call, not sourced from the handoff.** Decoding pods
material-by-material gives 78 single-material sites (many overlapping the same lithology class:
schist alone yields 34 jade/quartz/mica/gold pods separately). Treating each as its own camp reads
as absurd — adjacent veins in the same rock get worked from one camp in any real mining district.
I consolidated same-lithology-class pods within 1.5 km into one multi-material post (down to 71
total; the schist cluster shrinks from 34 sites to ~25, granite's 5 pods to 3). **This threshold is
mine, not derived from anything in the docs — if you want either finer (one post per pod, 78
sites) or coarser (one post per Círculo-class pairing, ~20-30 sites), say so and I'll redo it.**
Bog iron is reported folded into vivianite's site, not separately, per the handoff's explicit
"one site, not two."

**3. A factual correction to the handoff itself.** It states "only Circulo_D area is under 1km
from ocean per `06_tappa6_suitability.md` §8." I checked this directly against two independent
rasters (lithology_v6 at 30m and biome_id at 120m, both agree): the actual closest-to-ocean site is
**Circulo_F3_small at 0.92-0.97 km**. Circulo_D_20k is 1.26-1.32 km — just outside 1 km, essentially
tied with Circulo_E5_2k (1.27-1.29 km). The 06 doc itself never names which site is under 1 km
(I read it directly — it just says "closest site... <1 km, most sit 2-17 km inland"); the
attribution to Circulo_D appears to be an inference made somewhere downstream, and it doesn't hold
up. This matters here because it means Circulo_D's own coastal/ferry village is not coincident
with its parent Círculo — it's ~1.3 km away, same as everyone else, not a zero-distance case.

## Mining posts (71)

**Renamed to English** (`Mina_<Class>_NN` → `Mine_<Class>_NN`, e.g. `Mina_Schist_01` →
`Mine_Schist_01`, `Mina_Calcario_07` → `Mine_Limestone_07`, `Mina_Granito` → `Mine_Granite`) for
consistency with the rest of the project's descriptive labels (lithology classes, resource
materials, road-network edge types are all English already) — per Nico's direct request. `Circulo`
itself is intentionally left untranslated everywhere, including here, since it's an in-world term
(`scenario_reference.md`'s own vocabulary), not a generic English word standing in for it.

Decoded from every resource-pod raster in `data/processed/geomorphology/`, one post per
consolidated pod cluster, attached to the nearest Círculo (a second Círculo is listed when it's
within 15% of the nearest distance or within 1 km absolute — "can be shared by 2+" per the
handoff). Granite's own zone is tiny (13.84 km², vs 595-4,437 km² for the other classes) so its
two materials only place 3+2=5 raw pods total (3 posts after consolidation) — already flagged
upstream in `tappa8_limestone_granite_materials_meta.json` as a poor fit for the project's
standard 5 km pod-separation constant, not something new I'm introducing.

Granite (2 zones, no pod raster yet per the original handoff) and the two limestone zones (North
Coast Limestone, Sedimentary Bay — no, wait: sedimentary_limestone) — **this gap is actually
closed**. The handoff was written before Tappa 8's §14 landed; `resource_calcite.npy`,
`resource_mica_granite.npy`, and `resource_quartz_granite.npy` all exist now (7+3+2 = 12 pods).
I sited all of them. The only remaining true gap is **mountain huts having no resource tie** (by
design) and marble, which the Scenario chat explicitly assigned no Vértice domain and no pod
raster — correctly absent here, not an oversight.

Full breakdown by lithology class in the CSV/GeoJSON; distances to attached Círculo range from
0.9 km to 34.9 km (the far end is placer_magnetite's coastal/estuarine pods on the far side of the
island from every Círculo except D_20k and F7/F8 — worth a second look if 35 km feels too far for
a "rotating workforce" camp, though FIFO camps regularly move workers that kind of distance
in the real world).

### Revision: structure tiers by real cost-distance (addresses your granularity question)

You flagged that treating all 71 mining posts as full camps didn't make sense — some sit close
enough to a Círculo for a daily walk, and some remote posts cluster tightly enough to share one
camp rather than each getting its own. You asked for two specific methodology choices, both
implemented exactly as specified: **real cost-distance hours** (not straight-line km) for the
day-commute cut, reusing this project's own routing model byte-for-byte, not a new one; and,
for any qualifying cluster, the hub is **whichever pod is closest (by real hours) to its own
parent Círculo**, not whichever has the largest resource area.

**Method.** I rebuilt the exact graph `run_tappa9_road_network.py` builds the road network from
— same 120 m DEM (block-mean of the native 30 m DEM), same effective-land mask (`land_mask.npy`
minus `lake_mask.npy` block-OR'd to 120 m), same `combined_friction_multiplier_120m.npy`
(lithology × biome × river-crossing, already computed and locked by Tappa 9 — I loaded it
directly rather than recomputing it), same `sea_mode="impassable"` (a mining commute is a walk,
not a boat trip, so a route that would require crossing open water isn't a "commute" — it's a
sign the post needs full FIFO logistics or isn't viable from land at all). I ran Dijkstra once
from each of the 17 Círculos over this graph (9.0 s total, all 17 runs) and sampled every mining
post's real one-way travel hours to **every** Círculo, not just its straight-line nearest — this
matters: the real-hours-nearest Círculo differs from the straight-line-nearest one for several
posts (e.g. a granite pod 2.86 km straight-line from Circulo_A_40k is actually 0.60 h closer, by
real terrain cost, to Circulo_F8_small — friction and slope along the direct line to A_40k cost
more than the longer-looking route to F8_small).

**Three tiers, thresholds are my call, not sourced from anywhere — flagging exactly like the
1.5 km consolidation threshold above:**

| tier | one-way real-hours cutoff | rationale | count |
|---|---|---|---|
| `day_commute_no_structure` | ≤ 1.5 h | round trip ≤ 3 h, leaves a full working day; no quarters needed, workers stay at their home Círculo | 12 |
| `small_camp` | 1.5 h – 5 h | too far for daily but not extreme — a bunkhouse, not full facilities | 20 |
| `full_camp` | > 5 h, or land-isolated (see below) | genuine multi-day/FIFO logistics | 29 |

Distribution that motivated these cuts (68 of 71 posts reachable by land from at least one
Círculo; the other 3 are a separate finding below): min 0.22 h, p10 0.87 h, p25 2.15 h, median
5.57 h, p75 8.71 h, p90 12.34 h, max 17.46 h. 1.5 h and 5 h sit at natural breaks in that spread,
not arbitrary round numbers, but they're still a judgment call — say so if you want them moved.

**Hub consolidation.** Among the 59 posts that don't qualify for day-commute, I grouped posts
within 5 km of each other by straight-line distance (checked at 5/8/10/12 km first — 8 km and
above start *chaining* the whole schist district into one 17-25-member "cluster" via transitive
links across a district that actually spans 20-30 km end to end, which isn't a real "share one
camp" case, just single-linkage chaining; 5 km gives 8 clusters with a max internal span of 7.3 km,
which is a defensible walking-distance-between-pods-of-one-camp radius). That produced **8 hub
camps, each absorbing 1-2 satellite pods (10 satellites total, no structure of their own)** —
down from what would have been 18 separate full/small camps:

| hub (lithology, xy) | tier | real h to its Círculo | satellites absorbed |
|---|---|---|---|
| volcanic (-34235,-63985) | small_camp | 2.15 h → D_20k | 1 (basin_fill, 3.5 km away) |
| sedimentary_limestone (-39955,-998) | small_camp | 4.20 h → E2_2k | 2 (basin_fill + limestone, 4.7-5.0 km) |
| schist (12355,18545) | full_camp | 6.71 h → F1_small | 2 (schist ×2, 3.0-4.4 km) |
| schist (-9985,28524) | full_camp | 8.65 h → F2_small | 1 (schist, 3.7 km) |
| schist (37451,1563) | full_camp | 8.31 h → E4_2k | 1 (schist, 4.5 km) |
| schist (20665,5675) | full_camp | 13.02 h → F1_small | 1 (schist — the land-isolated one, see below, 3.9 km) |
| schist (8185,9965) | full_camp | 12.34 h → F1_small | 1 (schist, 4.7 km) |
| schist (-2023,14202) | full_camp | 14.54 h → F2_small | 1 (schist, 1.7 km) |

Net effect on the mining-post count that actually needs built infrastructure: **49 structure-
bearing sites** (41 standalone + 8 hubs) instead of 71 — 12 day-commute posts need nothing built,
and 10 satellite posts share their hub's structure. All 71 stay on the map as resource-tie
points (workers still walk out to the actual vein from the hub), just not all 71 get their own
building footprint anymore.

**A real finding, not a threshold artifact: 3 of the 71 posts are unreachable BY LAND from
literally every one of the 17 Círculos** (not just their nearest one — all 17 dist grids read
`inf`), the same isolation class Tappa 9 already gave Circulo_D_20k itself. One
(`basin_fill`, 9235/-74466, the same "far side of the island" pod already flagged above at
34.9 km straight-line) sits on genuinely disconnected land — a real islet/promontory finding,
consistent with what was already flagged. The other two (a `greywacke` pod and a `schist` pod)
land in a 120 m cell that `block_any`-aggregates as touching `lake_mask` — likely a lake-edge
raster artifact from resolution mismatch (the underlying 30 m pod pixel may well be dry land;
the 120 m routing cell it falls into isn't), not necessarily a real islet. All three are folded
into `full_camp` by necessity (one is a hub's satellite, sharing a camp with a reachable
neighbour 3.9 km away — plausible even without a *routed* path between them, since the routing
graph's land-only assumption is stricter than "can two adjacent points physically be walked
between"). Worth a real look at the lake-mask resolution issue if it recurs elsewhere in the
project; I did not attempt to fix it here, only to flag it and route around it defensibly.

## Mountain huts (18, rebuilt on Tappa 7's own outpost layer)

**Superseded, not merged.** My first pass sited 2 huts from scratch (road-crossing alpine
stretches from Tappa 9). Nico pointed out a real, already-existing candidate-outpost layer from
Tappa 7's fauna work — `data/processed/fauna/outpost_candidates.geojson`, 18 sites (17 on the
"main_spine" massif, 1 on "south_branch"), each carrying a real composite-suitability score
(resource access, slope, climate mildness) and even a `nacre_suitability_at_site` field (Nacre
being the alpine apex predator Tappa 7 modeled 13 iterations of — see `claude/index.md`'s "07"
row). That work was explicitly unfinished on the settlement-administration side: every one of the
18 carries `status_is_authorial_final: false` and a `status_prior` of "temporary refuge" (17) or
"abandoned" (1) — a real, still-open decision, not something I'm resolving here. Nico's framing:
this Tappa is the right place to handle the Círculo-attachment side of these sites, since Tappa 7
never touched that.

**What I did**: dropped my 2 road-crossing huts entirely (they're a reasonable but now-redundant
heuristic next to a purpose-built dataset) and took all 18 outpost candidates as the mountain-hut
set, unchanged in name/geometry/fauna properties. Added `attached_circulos` via straight-line
nearest-Círculo (same tolerance-band convention as forest posts: a second/third Círculo listed
when within 15% of the nearest distance or 1 km absolute) — real distances run 12.1-25.3 km,
notably farther than a Círculo's own basin the way mining/forest posts sit, consistent with
these being alpine-spine outposts serving a wide catchment rather than one settlement's own
backyard. `status_prior`/`status_is_authorial_final` are carried through verbatim as properties,
not acted on — **whether each of these 18 is actually active, abandoned, or narratively real at
all is still open, exactly as Tappa 7 left it.** This Tappa only added where each one
administratively belongs.

## Coastal villages (12, revised — a real bug found and fixed, plus four direct Nico edits)

**A genuine methodology bug, found while investigating Nico's report on Coastal_Village 11
(old numbering).** The original "nearest ocean cell" computation used `lithology_v6 == 0` as the
ocean mask, on the documented convention that 0 = ocean. That's true for the map's real
coastline, but I never checked whether *every* disconnected blob of `lithology==0` is actually
part of the sea. It isn't: connected-component labeling the mask found 1,231 separate blobs;
exactly two are large and touch the domain border (the real ocean, ~8,935 km² + 1,807 km²); a
handful of the rest are small but genuinely part of that coastline (border-touching too); and a
few — one 3,394-cell (3.05 km²) enclosed blob in particular, centered around (-41292, 7802),
confirmed against `lake_mask.npy` and a DEM that dips to -10.3 m there — are a below-sea-level
**enclosed lake with no border connection at all**, mislabeled as ocean by the lithology
classification. `Circulo_E2_2k`'s old nearest-coast computation (8.59 km) had locked onto
exactly this lake, not the real sea. **Fix**: rebuilt the ocean mask from only the
border-touching components ("true ocean"), and added a hard requirement that a village's final
site must itself be within 300 m of that true-ocean mask — the old refinement step only checked
"is this land and is it flat," never "is it actually on the coast," which is how a flat inland
spot near a false-ocean lake could get picked in the first place. Re-running the corrected
pipeline changed the recorded distance-to-ocean for 3 Círculos: `Circulo_F4_small` (5.27→11.0 km),
`Circulo_E2_2k` (8.59→11.5 km), `Circulo_F1_small` (16.96→17.41 km) — and reshuffled the 10 km
clustering slightly (F4_small no longer merges with F1_small; F1_small merges into the
F2_small+B_35k cluster instead). Full corrected baseline: 14 villages, all now confirmed
genuinely coastal.

**Then four direct edits from Nico's own review, on top of that corrected baseline:**

1. **Dropped the villages for Circulo_F3_small, Circulo_D_20k, and Circulo_E5_2k** — all three
   sit under 1.3 km from real ocean already, close enough that Nico will give them a dock
   directly as part of their own Urban Scale footprint rather than a separate satellite
   settlement. **Flagging plainly: this leaves Circulo_E5_2k with zero attached auxiliary
   settlements of any kind** (no mining pod nearby, no forest post — its basin is one of the 3
   unresolved `basin_labeled==0` cases). D_20k and F3_small still have mining posts. See the
   Summary section above.
2. **Relocated Circulo_E2_2k's village** off the false-ocean lake site entirely, to a real
   coastal point 1.6 km from `Mine_Limestone_07` (formerly `Mina_Calcario_07`) — Nico's explicit
   direction, prioritizing the mine having its own port over minimizing distance to the Círculo
   itself (now 20.2 km, up from the old, invalid 8.59 km reading). Slope at the new site: 0.83°,
   a genuinely good harbor.
3. **Relocated Circulo_F7_small's existing village** ~9.3 km SSW, off a 22.6° slope (the true-
   ocean-adjacency fix actually made this WORSE before Nico's edit — the old, non-coastal-checked
   refinement had found a fake-flat 5.79° spot nearby that wasn't really on the coast either) to a
   genuinely flat 1.06° real coastal point.
4. **Added a new second village for Circulo_F7_small**, near `Mine_Greywacke_03` (2.48 km away) —
   Nico's explicit addition, since that mine is itself only 0.43 km from real ocean and could ship
   ore out directly rather than hauling it inland. Honest caveat: this specific stretch of coast is
   genuinely rugged — every cell within 800 m of the mine is 15-28° slope; 6.36° is the best
   available within a 3 km search, not a flat harbor. A real rocky-cove site, same "narratively
   acceptable, not silently defective" framing as the steep sites already flagged below — not
   something I'd have picked on my own without Nico's specific direction to prioritize the mine.

**Fifth edit, this pass: Coastal_Village_08 (Circulo_F4_small) moved back near the Círculo, onto
a lake instead of the ocean.** You flagged that its relocated (true-ocean) site had actually
displaced it FROM a body of water that, while also not the real sea, is a much larger and
genuinely different feature than the small lake that caused the E2_2k relocation — and since
this village was never a ferry point, sitting on a large lake is fine for a fishing-only
settlement. I checked this directly rather than taking it on faith: the village's ORIGINAL
(pre-fix) site sat on connected-component label 383 of the false-ocean mask, a **131.5 km²**
below-sea-level basin (DEM dips to -457.6 m across it) — more than 40× the size of the 3.05 km²
lake near E2_2k, and a genuinely different feature (zero cell overlap between the two, ~370 km
apart). It sits ~201 m from the nearest real ocean at closest approach, so it's not a
near-miss/thin-landbridge case either — a real, separate, below-sea-level lake basin. I re-ran
the same flattest-land-within-1.5km refinement, targeting adjacency to THIS lake instead of true
ocean: new site **5.27 km from the lake shore / 5.72 km from Circulo_F4_small, slope 2.05°** —
better on both counts than the true-ocean-adjacent site this village had moved to (11.5 km from
the Círculo, 4.33°). Recorded on the feature: `site_type: "lake_fishing_village_no_ferry"`,
`dist_to_lake_km: 5.27`, `lake_component_area_km2: 131.52`, alongside the unchanged
`dist_to_ocean_km: 11.0` (Circulo_F4_small's own distance to the real sea, kept for reference —
this village just doesn't sit there anymore). `ferry_connection_candidate` is `false`, as it was
before. Net effect on the total: still 12 villages, no count change, just this one relocated.

Two sites still carry meaningfully steep slopes after all of this, worth a narrative read (rocky
cove villages are a real thing, not necessarily a siting problem):

| village (attached) | slope |
|---|---|
| Circulo_E4_2k | 9.19° |
| Circulo_E1_2k | 7.13° |
| Coastal_Village near Mine_Greywacke_03 (Circulo_F7_small) | 6.36° |
| Circulo_E3_2k | 6.07° |
(most others land under 4.5°, several under 1.5°)

The `Circulo_D_20k↔Circulo_E3_2k` ferry crossing (57.7 km, 10.88 h — Tappa 9's own number, not
re-derived) no longer has a dedicated village at the D_20k end (see point 1 above) — that endpoint
is now assumed to be D_20k's own future Urban Scale dock, not a separate feature in this layer.
`Circulo_E3_2k`'s own village still carries `ferry_connection_candidate: true`.

## Forest posts (7 sites, serving 9 Círculos)

The handoff's own suggestion — "the 2 Temperate Forest sites and 2 Woodland/Shrubland sites... but
check each Círculo's actual basin, don't assume only those four qualify" — undersold it. Checked
all 17 against `basin_labeled.npy` (54 major basins) × `biome_id.npy` (forest classes: Subalpine
Wet Forest, Subalpine Woodland, Temperate Forest, Woodland/Shrubland):

- **9 Círculos have forest in their own basin**: Circulo_C_25k (84% of basin forested, shares the
  basin with F6_small), Circulo_E3_2k (67%), Circulo_A_40k (53%, shares its basin with F7_small
  AND F8_small — three Círculos, one basin), Circulo_D_20k (51% — notable, since D_20k's own site
  biome is Grassland, this is a purely-upstream forest resource), Circulo_E2_2k (67%),
  Circulo_E1_2k (36%), Circulo_F5_small (29%).
- **4 Círculos confirmed zero forest in their basin**: Circulo_F4_small, Circulo_F1_small,
  Circulo_B_35k, Circulo_F2_small — all 0.0%. These four cluster tightly together geographically
  (x≈10-32 km, y≈27-35 km) — this is almost certainly the "North Plains cluster" the handoff
  named without defining. I found it from the data, not from a predefined label anywhere in the
  docs (I looked — "North Plains" isn't a formally delineated group anywhere I can find).
- **3 Círculos are undetermined, not zero**: Circulo_E4_2k, Circulo_E5_2k, Circulo_F3_small all
  sit in `basin_labeled == 0`, which isn't one basin — it's the catch-all label for every cell
  outside the 54 tracked major basins (>20 km² each). I can't tell whether these three share a
  real minor catchment or sit in three unrelated ones without running a proper pour-point
  delineation from `flow_direction_code.npy`, which I didn't do here. Flagging as genuinely open,
  not guessing.

For basins with shared Círculos (basin 1: C_25k+F6_small; basin 3: A_40k+F7_small+F8_small), I
sited **one** forest post per basin, shared — per the handoff's "could share if two Círculos draw
from the same basin, not assumed," here it's directly confirmed by the raster, not assumed.

Site selection within each qualifying basin: among `streams.geojson` reaches whose midpoint falls
in that basin's forest cover AND at higher elevation than the Círculo itself (a real "upstream"
test, using the DEM — my first pass without the elevation filter kept picking each basin's single
biggest river, which is almost always near the basin *outlet*, the opposite of upstream; C_25k's
candidate was 806 m *lower* than the Círculo before I added that filter), I picked the highest
Strahler-order reach. Distances from the post to its nearest parent Círculo: 0.6-11.0 km.

## Files

- `auxiliary_settlements.geojson` — all 108 sites, ready to merge into a new layer or into
  `circulo_candidate_sites.geojson`, same CRS. **Every feature now carries the same 37 property
  keys** (see "Attribute schema fix" above), null/blank where a key doesn't apply to that
  settlement type — check `settlement_type` to know which columns are meaningful for a given
  row. `attached_circulos`, `resources`, and `shares_structure_with` are semicolon-joined
  strings, not JSON arrays. Mining-post features carry `structure_tier`
  (`day_commute_no_structure` / `small_camp` / `full_camp` / `satellite_no_own_structure`),
  `real_cost_distance_hours_to_nearest_circulo`, `real_nearest_circulo`,
  `reachable_by_land_from_any_circulo`, `is_camp_hub`, `shares_structure_with`, `attached_hub`,
  `resources` (clean material-name list), `resource_note`. Mountain-hut features carry Tappa 7's
  own `massif`/`composite_suitability`/`resource_suit`/`slope_suit`/`climate_mildness`/
  `nacre_suitability_at_site`/`status_prior`/`status_is_authorial_final`, plus this Tappa's
  `attached_circulos`/`nearest_circulo_km_straight_line`. Coastal-village features carry
  `harbor_slope_deg`, `dist_to_ocean_km`, `site_type` (`true_ocean_village` or
  `lake_fishing_village_no_ferry` — only Coastal_Village_08 is the latter — plus
  `dist_to_lake_km`/`lake_component_area_km2` when it applies), `ferry_connection_candidate`,
  and a free-text `note` recording any manual siting rationale.
- `auxiliary_settlements.csv` — same 37-column schema as the GeoJSON (`x_m`/`y_m` pulled out of
  the geometry as their own columns), built with a dict-per-row writer so there's no risk of a
  column-position mismatch. `real_cost_distance_hours_to_nearest_circulo` (mining posts) and
  `nearest_circulo_km_straight_line` (mining posts and huts, straight-line only) are now separate
  columns — no longer sharing one overloaded column the way the previous pass did.

## Open items for you / the Scenario chat

1. ~~Mining-post consolidation threshold (1.5 km) is my call — confirm or adjust.~~ Superseded by
   #1a below — the consolidation itself (71 posts) is unchanged, only how many of those 71 get
   a built structure has changed.
1a. ~~Structure-tier thresholds (day-commute ≤1.5h, full-camp >5h one-way real hours) and the
    hub-clustering radius (5 km straight-line) are my calls, not sourced from anywhere.~~
    **Confirmed by Nico (2026-08-21).** One caveat worth keeping in mind if the narrative ever
    specifies mine-commute transport beyond walking (cart, pack animal): these hours come from
    `cost_distance.py`'s Tobler-hiking-function graph, the same foot-travel model the whole
    project's cost-distance work uses (Tappa 6 siting, Tappa 9 roads) — not a claim that miners
    literally walk. If a faster transport mode ever gets its own speed function, these same
    tiers would need re-running against it, not just relabeling.
2. ~~Whether 71 mining posts... is the right granularity for play.~~ Resolved by the tiering: all
   71 stay as resource-tie points, but only 49 now carry an actual structure (41 standalone + 8
   shared hubs), which should address the "too many full camps" concern directly.
2a. The 2 pods that read as land-isolated because their 120m routing cell touches a lake
    (`block_any` aggregation from the 30m `lake_mask`) are very likely a resolution artifact, not
    a real islet — worth checking directly against the 30m lake mask if this matters for anything
    beyond this classification (I routed around it defensively rather than fixing it).
3. E4_2k/E5_2k/F3_small's forest status needs real pour-point catchment delineation to resolve,
   not assumed either way.
4. ~~Three coastal harbor sites (E1, F7, E4) are geometrically the best available but still
   steep.~~ Re-checked with the true-ocean-adjacency fix: E1_2k (7.13°) and E4_2k (9.19°, worse
   than before — the old reading was on a false-ocean cell) are still worth the same narrative
   read as originally flagged. F7_small's original steep site (22.6°, also worse once corrected)
   was fixed per Nico's direct edit (see Coastal villages above); its NEW second village (near
   Mine_Greywacke_03, 6.36°) is a fresh instance of the same question.
5. The 06-doc "<1km from ocean" attribution to Circulo_D was wrong (it's F3_small) — worth fixing
   at the source doc if anything downstream has repeated the same claim.
6. **New: Circulo_E5_2k has zero attached auxiliary settlements** after dropping its coastal
   village (see Summary/Coastal villages above) — worth a deliberate decision, not a silent gap.
7. **New: whether each of the 18 `Outpost_*` mountain huts is actually active, abandoned, or real
   at all is still undecided** (`status_is_authorial_final: false` on all 18, inherited from Tappa
   7, not resolved here) — this Tappa only added Círculo attachment, not a narrative status call.
8. **New: the false-ocean/lake mislabeling bug** (see Coastal villages above) is a
   `lithology_v6` classification gap, not something I fixed at the source — worth checking whether
   any OTHER downstream layer (not just this one) ever treated `lithology==0` as "the ocean"
   without checking border-connectivity first.

## Network connections, v2 re-run (2026-08-21, later the same day, from the Tappa 9 chat)

This file's 108-feature v2 (the rename to English names, the attribute-schema fix, and the 18
Tappa-7-sourced huts replacing the original 2) required re-running the connection script built
for the original 94-feature file. Full implementation record lives in
`09_tappa9_transports.md`'s §14 addendum; this section covers what's specific to this file's own
data.

**Schema migration**: `attached_circulos`/`resources`/etc. moved from JSON arrays to
semicolon-joined strings in v2 (see this doc's own "Attribute schema fix" section above) —
normalized back to lists on load, no data lost.

**A real bug found in v2, not present in the original file**: all 10 `satellite_no_own_structure`
mining posts' `attached_hub` property still names their hub under the OLD `Mina_*` naming — the
rename to `Mine_*`/`Mine_Limestone*` updated every feature's own `name` but not the cross-
references other features make to it. Confirmed 10/10, not a coincidence. Worked around here
(prefix-substitution lookup, asserted against the actual name set rather than assumed) — worth
fixing at the source next time this file is regenerated, since a future consumer without this
same workaround would silently mis-report every satellite as isolated.

**Mountain huts (18): built as a trail network, not spokes — see `09_tappa9_transports.md` §14
for the full design reasoning and results (16 trail edges linking all 17 `main_spine` huts,
`south_branch`'s 1 hut standing alone, 9 trailhead-access spokes down to the valley network, 0
huts isolated).** One thing worth restating here since it's this file's own open item: **all 18
huts still carry `status_is_authorial_final: false`** (open item #7 above, inherited from Tappa
7) — this run gave every one of them a real physical trail connection regardless of that
unresolved status. That's a reasonable default (infrastructure doesn't need to wait on narrative
sign-off), but if #7 resolves toward "some of these aren't real" or "abandoned," their trail
edges/trailhead spokes would need pruning, not just their `attached_circulos`.

**Isolated settlements, re-checked directly against the actual land/lake rasters (not just
"no path found"), 17 total (down from 92 connectable in v1's naming to 90 in v2's, same 17
isolated count, different individual sites due to the mining-post renumbering)**:
- **14 inherit `Circulo_D_20k`'s own known land-isolation** (its excluded status is unchanged
  from Tappa 9) — 12 mining posts (`Mine_BasinFill_03/12/13`, `Mine_Greywacke_04`,
  `Mine_Volcanic_01`-`08`) + `Forest_Post_05`.
- **`Mine_Schist_22`** (attached to `Circulo_E3_2k`/`Circulo_F1_small`): checked directly against
  `land_mask.npy`/`lake_mask.npy` — its own cell fails `effective_land` (touches the 120m-
  aggregated lake mask) while all 8 neighboring cells pass. This is the SAME resolution artifact
  as open item #2a, now confirmed a 4th time (previously 2 mining posts + `Coastal_Village_11`,
  see correction below) and for the first time on a post NOT attached to `Circulo_D_20k`.
- **`Mine_BasinFill_16`** (attached to `Circulo_F8_small`/`Circulo_F7_small`): checked directly —
  its own cell AND all 8 neighbors pass `effective_land` cleanly. This is NOT a raster artifact;
  it's genuinely cut off from every routed Círculo by the land-only graph, the same isolation
  class as `Circulo_D_20k` itself. Matches the original Tappa 10 report's own "land-isolated
  posts" finding (see Mining posts section above).
- **`Coastal_Village_03`** (`Circulo_E4_2k`, `site_type: true_ocean_village`,
  `dist_to_ocean_km: 4.77`, harbor slope 9.19°) **and `Coastal_Village_11`** (`Circulo_F7_small`,
  `site_type: true_ocean_village`, `dist_to_ocean_km: 16.19`, harbor slope 1.06°): checked
  directly, and **this is a correction of what I reported after the v1 run** (I'd characterized
  `Coastal_Village_11`'s predecessor, `Vila_Costeira_14`, as hitting the same lake-mask artifact
  as the mining posts — re-checking now against the actual raster, that's wrong). Both villages'
  own cells sit just 1-2 cells (120-240 m) outside `land_mask.npy`'s land — not touching
  `lake_mask` at all, and not deep water either (BFS confirms real land 120-240 m away). This is a
  DIFFERENT bug class from #2a: these sites were placed using this file's own "true ocean" mask
  (the `lithology_v6==0`, border-connectivity-checked one described in this doc's Coastal
  villages section), which is a different, more precise coastline than the routing graph's
  `land_mask.npy` (an older, coarser Tappa 1-era raster). A point legitimately coastal by this
  file's own definition can land just outside "land" by the routing graph's definition — a
  cross-dataset mismatch at the coastline, not a lake-adjacency artifact. Worth deciding whether
  `land_mask.npy` should ever get reconciled against the corrected true-ocean mask, independent of
  whether these two specific villages get nudged a cell or two inland.

## Network connections, third pass — four bugs Nico found in QGIS (2026-08-23)

Full root-cause diagnosis of all four is in `09_tappa9_transports.md`'s §15; this section
covers the concrete before/after for the auxiliary-settlement side, which is what actually
changed shape. Nico approved the full rewrite in one go ("Corrigir tudo de uma vez"),
explicitly accepting that the previous run's numbers (73 spokes/903.6 km, 339.4 km trail, 17
isolated) would change.

**The core change**: `run_tappa10_network_connections.py`'s connection search went from
"each settlement checks only its own attached Círculos' backbone edges, or a direct
fallback to its hub" to a single greedy, incrementally-growing network (Prim's-style) —
the Tappa 9 backbone seeds the initial network, then every settlement/hut-trailhead
connects, cheapest-first, to whichever is cheapest across the ENTIRE current network
(backbone, anywhere along it, or any spoke/trail already built earlier in the same run).
New edge-type taxonomy: `satellite_to_hub_spoke` (unchanged, hard administrative rule),
`circulo_spoke` (lands exactly on a Círculo — weight symmetrized against the precomputed
Círculo-sourced Dijkstra), `backbone_t_junction_spoke` (ties into a mid-route point of the
backbone, generalized from the old same-hub-only T-junction rule to any backbone edge),
`feeder_t_junction_spoke` (ties into another settlement's already-built spoke/trail —
entirely new, this is what fixes the redundant-path bugs), `mountain_trail` (unchanged),
`isolated_pocket_road` (new — see below).

**Southland (`Circulo_D_20k`) fix, concretely**: of its 14 attached settlements, 12 now
connect via `feeder_t_junction_spoke` chains rooted at `Mine_Volcanic_02` (itself reached
first by `Mine_BasinFill_10`'s satellite spoke, then absorbed into the network at zero
extra cost — see the self-loop note below), 1 more (`Mine_Volcanic_02` itself) is already
part of the network for free, and 1 (`Mine_Greywacke_04`) remains genuinely isolated —
checked directly, no finite-cost path to any other settlement or the network exists on this
land-only graph. Where the previous run had 14 disconnected dots plus one lone hub-to-mine
road, the cluster now has a real internal road system: e.g.
`Mine_Volcanic_08 -> Mine_Volcanic_07 -> Forest_Post_05 -> Mine_BasinFill_10 ->
Mine_Volcanic_02` and a separate `Mine_BasinFill_13 -> Mine_BasinFill_12 ->
Mine_BasinFill_10` branch. The isolated-pocket-MST step (grouping settlements that are
mutually reachable but cut off from the main network, each pocket getting its own local
MST via the same `build_local_mst()` helper used for hut trails) exists for cases like this
but wasn't actually needed this run — the greedy algorithm's wider search already resolved
every pocket except the four listed below as genuinely, individually isolated.

**`Coastal_Village_09`**: now a 4.42 km / 1.17 h `feeder_t_junction_spoke` onto
`Mine_Limestone_06`'s existing spoke, replacing the independent, much longer route the
restricted search previously produced (the old output file was overwritten by this run, so
an exact old figure can't be re-quoted here — but the structural fix, a short hop onto
already-built infrastructure instead of routing independently back to a Círculo, is
confirmed directly in the new geometry).

**`Forest_Post_03`**: now a 0.29 km / 0.06 h `backbone_t_junction_spoke` directly onto the
`Circulo_E3_2k <-> Circulo_F7_small` backbone edge — the same edge its old spoke ran
170-424 m parallel to for its entire length without ever joining. This is the clearest
before/after of the whole fix: a road that used to shadow an existing edge for its full run
now ties into it almost immediately.

**Self-loop bug found during verification (not one of Nico's four, caught while
spot-checking)**: the greedy loop initially produced zero-length, zero-cost
`feeder_t_junction_spoke` edges from a settlement to *itself* whenever another
settlement's spoke had already terminated exactly on its site cell (e.g. `Mine_Volcanic_02`,
reached first by `Mine_BasinFill_10`). Fixed by detecting when the chosen target cell IS the
settlement's own site and skipping edge creation — it's already connected at zero cost, no
edge needed. 15 settlements were affected across this run (`n_trivially_absorbed` in the
meta JSON); none of them are double-counted as "connected" and "isolated."

**Reporting gap also found during verification**: satellites with no finite path to their
own hub (`n_satellite_isolated` in the log — a hard rule, never retried against the wider
network even though the greedy search exists now) were counted in a log line but never
added to the `still_isolated` list in the meta JSON, silently undercounting genuine
isolation by one. Fixed — `still_isolated` now includes them with their own reason string.
This surfaced `Mine_Schist_22` as the 5th isolated settlement (its hub is
`Mine_Schist_20`): this is not new, it's the already-documented lake-mask raster artifact
(4th confirmed instance, see the isolated-settlements list above) — its own cell fails
`effective_land` so no Dijkstra pass from it ever finds anywhere reachable, previously
just invisible in this file's own reporting.

**Vertex-exact endpoints**: applied here too (see `09_tappa9_transports.md` §15) — verified
directly against all 95 settlement site coordinates that appear as a line's `from`/
`settlement_name` in the output, 0 mismatches.

**Result, third pass**: 16 mountain-trail edges (145.9 km, unchanged from the v2 trail
design) + 79 network spokes (471.1 km) + 0 isolated-pocket roads (the greedy search resolved
every pocket except the 4 below on its own) = 85 of 89 settlements/trailheads connected
(after excluding the 9 hard-rule satellites and the 1 satellite with no finite path to its
own hub). Connection-type breakdown: 9 `satellite_to_hub_spoke`, 27
`backbone_t_junction_spoke`, 42 `feeder_t_junction_spoke`, 1 `circulo_spoke`. **5 settlements
remain genuinely, individually isolated** — checked directly, no finite-cost path to
anything: `Mine_Greywacke_04` (Southland, see above), `Mine_BasinFill_16` (a real
land-isolation, matches the original report — see the isolated-settlements list above),
`Coastal_Village_03` and `Coastal_Village_11` (the land-mask/true-ocean-mask cross-dataset
mismatch documented above — unchanged by this pass, out of scope for a routing-algorithm
fix), and `Mine_Schist_22` (the lake-mask raster artifact, see the reporting-gap note
below — not new, just newly surfaced in `still_isolated`). Total spoke+trail length: 617.0 km,
vs. 1,243.0 km before (73 spokes at 903.6 km + 339.4 km trail) — the drop is expected and
desired: most of that previous length was the redundant/parallel routing this pass exists
to eliminate, not real infrastructure lost.

## Network connections, fourth pass — a false "Southland fix" and the whole mountain
## system silently disconnected (2026-08-23, same day, Nico's second QGIS review)

Nico reported two more things after reviewing the third pass's output: (1) some hub camps
(satellite-mine hubs) had no road out to anywhere except each other, and separately "all the
mountain system is disconnected from any Círculo"; (2) `Mine_Greywacke_05` used to connect
toward `Circulo_E3_2k` and now routes toward a `Circulo_E2_2k`-side system instead, with a
question about whether `attached_circulos` should follow the road topology rather than the
other way around.

**Bug 1, confirmed directly, and it invalidates a claim I made in the previous message to
Nico.** Checked the actual output: zero of the 17 `main_spine` hut-trail leaves had any real
connection to the backbone — the trail system genuinely was floating, exactly as reported.
Root cause: `add_to_network()` was called for BOTH mountain-trail edges AND
`satellite_to_hub_spoke` edges at construction time, seeding those cells into the SAME
shared "live network" array the backbone uses. But neither actually means "reaches a
Círculo" — a trail only connects huts to each other, and a satellite's spoke is a direct,
unconditional route to its own hub regardless of whether that hub reaches anywhere. Because
a leaf's own site cell is itself a trail cell, and a hub's own site cell is exactly where its
satellites' spokes terminate, both were seen by the greedy loop as "already in the network"
at zero cost the moment they were processed — silently marked connected without ever
building a real bridge. **This is the same bug that produced the 15/7 "trivially absorbed"
counts logged in the third pass** — I read those as a success (the self-loop fix from the
previous round), but at least some of them were this deeper problem wearing the self-loop
fix's clothing: the settlement wasn't legitimately already connected, it only looked that
way because its own local cluster had been prematurely counted as network.

**This means my previous message's "Southland fix" claim was wrong.** I reported 12 of
`Circulo_D_20k`'s 14 attached settlements as newly connected via a real internal road chain
rooted at `Mine_Volcanic_02`. In fact `Mine_Volcanic_02` (the hub) was never bridged to
anything outside the Southland cluster — its apparent "connection" was `Mine_BasinFill_10`'s
satellite spoke terminating on it, which was wrongly seeded into the shared network. The
chain of feeder edges built off it was real (the settlements ARE reachable from each other),
but none of it reached a Círculo. Correcting the record here rather than leaving the earlier
message standing uncorrected.

**Fix**: satellite spokes and mountain-trail edges are no longer added to the shared network
at construction time. Each hub's satellite-spoke cells accumulate in `hub_cluster_cells`
and each massif's trail cells in `massif_cluster_cells`, kept OUT of the live network until
a new `promote_cluster()` step (called every time a hub or a trail leaf gets a real
connection during the greedy loop) brings the whole local cluster online — at that point,
and only then, its cells become valid connection targets for everyone else too. A hub or
massif that never earns a real bridge stays local-only for good, correctly falling through
to the isolated-pocket step instead of silently passing as connected.

**Result, re-run**: **0 of 2 massifs stayed unbridged** — the mountain trail system does now
reach the backbone for real: `Outpost_MainSpine_17` bridges the 17-hut `main_spine` massif
with a 5.18 km/2.45h spur onto the `Circulo_E3_2k↔Circulo_E4_2k` backbone edge, after which
the other 7 `main_spine` leaves are legitimately (not spuriously) absorbed since the whole
massif-plus-bridge is now one connected piece; `south_branch`'s single hut connects directly
with a 10.67 km/2.99h feeder onto `Forest_Post_06`. **1 of 9 hub clusters stayed
unbridged** — `Mine_Volcanic_02`'s Southland cluster, confirming the correction above: this
is now honestly reported. The isolated-pocket step (already built for exactly this case,
just never actually triggered until now) picked it up correctly: 12 of the 14
`Circulo_D_20k`-attached settlements form one mutually-reachable pocket
(`Mine_BasinFill_03/12/13`, `Forest_Post_05`, `Mine_Volcanic_01`-`08`) with its own real
internal least-cost MST — 11 `isolated_pocket_road` edges, 103.6 km — genuinely isolated
from every Círculo but with a sensible local road system, which is exactly what this pass
was originally supposed to deliver for the Southland. `Mine_Greywacke_04` still isn't even
part of this pocket — no finite-cost route to it either, unchanged. The other 4 previously-
isolated settlements (`Mine_BasinFill_16`, `Coastal_Village_03`, `Coastal_Village_11`,
`Mine_Schist_22`) are unaffected by this fix, still isolated for their own separate reasons.
**Updated totals**: 16 mountain-trail edges (145.9 km) + 75 network spokes (401.8 km) + 11
isolated-pocket roads (103.6 km) = 651.3 km total, 5 settlements still genuinely,
individually isolated (same 5 names as before, `Mine_Volcanic_02`'s cluster now correctly
NOT among them since it has its own honest local network).

**`Mine_Greywacke_05` now routing toward a `Circulo_E2_2k`-side system instead of
`Circulo_E3_2k`, and Nico's question about deriving `attached_circulos` from road
topology.** Traced the actual chain: `Mine_Greywacke_05 → Mine_Greywacke_07 →
Mine_Limestone_07 → Mine_Limestone_04 → Mine_BasinFill_09 →` the
`Circulo_E2_2k↔Circulo_F6_small` backbone edge. This is real, not a bug — the new algorithm
finds the actual cheapest terrain-cost route, which for this mine happens to run toward
E2_2k rather than E3_2k; the old, simpler algorithm (not this project's own real
cost-distance, per Nico's own recollection) apparently didn't. **Recommendation on the
"derive `attached_circulos` from roads" question, offered as a caution rather than a
straight yes**: `attached_circulos` is not just a road-connection fact — it's the field that
also drives satellite→hub grouping, forest posts' basin-bound log economy, coastal
villages' ferry-Círculo pairing, and (for mining posts) which Círculo's economy the output
is counted against, none of which were computed from the road network and some of which
(forest posts' basin membership, mining posts' resource-siting Círculo) have their own
independent, non-road logic. Making it FOLLOW the road topology would also make it follow
whatever the routing algorithm's current tuning happens to produce — and this algorithm has
already changed three times in three days, each time changing which point is cheapest for a
given settlement. Coupling a stable administrative/economic fact to an actively-iterating
optimization is likely to cause more churn than it resolves. **My suggestion**: leave
`attached_circulos` as the independent administrative attribute it already is, and treat
`Mine_Greywacke_05`-style divergence as expected and fine — a mine can administratively
belong to one Círculo while its cheapest real road happens to run toward a different one's
territory, the same way a real district's roads don't always match its jurisdiction
boundaries. If Nico wants the road-derived fact recorded somewhere, I'd add it as a
clearly-separate, clearly-volatile field (e.g. `road_network_circulo`) rather than overwrite
`attached_circulos` — happy to build that if useful, but didn't do it unprompted given how
many other things key off the existing field.

**Follow-up, same discussion, Nico's second message (2026-08-23): asked for the concrete
`Mine_Greywacke_05 -> Circulo_E3_2k` distance as evidence, and proposed two remedies for
cases like it** — (a) reorganize `attached_circulos` once the network structure is final, or
(b) build new roads whenever the gap between an attachment and its real route is too large.
Nico's own stated preference is (a), reasoning that some of these gaps look like genuine
terrain ruptures rather than an algorithm quirk. Computed the number directly:
`Mine_Greywacke_05 -> Circulo_E3_2k` costs **93.81 km / 23.71 h** over the real network today,
against **55.70 km / 14.13 h** to `Circulo_E2_2k` — a real, large gap, confirming the premise.
**Critical take on the two options, since this is a real trade-off, not a formality**: (b)
should be rejected as a general policy — building a road specifically to make an
administrative fact stay true inverts cause and effect, and directly contradicts Nico's own
"terrain ruptures" reasoning (if the terrain genuinely doesn't support a good route, forcing
one there is authoring around the map, not respecting it). (a) is the right general direction,
but "once the network structure is final" is doing a lot of work — this network has been
rewritten three times in three days, and "final" isn't a state this project has reached
before on any layer without an explicit, dated lock (see "06"'s Círculo-siting lock as the
template). Two concrete conditions before a reorg pass should run: **first**, a numeric
threshold, not a case-by-case judgment call — e.g. flag any settlement whose real network
distance to its `attached_circulos` entry exceeds some multiple (2-3x?) of its distance to the
nearest OTHER routed Círculo, so the reassignment is driven by the same kind of rule the
project already uses elsewhere (redundancy thresholds, structure-tier hour cutoffs), not eyeballed
per case. **Second**, decide up front whether reassignment touches ALL of `attached_circulos`'s
downstream uses (hub/satellite grouping, forest-post economy, ferry pairing) or only a
road-specific reading of it — my standing recommendation is still to keep those separate
(administrative attachment vs. road-nearest Círculo as two distinct fields), since collapsing
them means every future road rewrite reopens the administrative map too. Not run yet — this is
a recommendation, waiting on Nico's confirmation of the threshold rule and the
separate-field question before any settlement gets reassigned.

**Same message, three concrete connectivity requests, investigated with real cost-distance
(not the flawed straight-line approximation this chat first reached for and then caught and
corrected)**: (1) "the mountain system has no connection to the north" — checked directly:
none of the north schist mines are actually isolated (none appear in `still_isolated` or any
`isolated_pockets` entry), but the entire `main_spine` massif's only real bridge to the
backbone is the single `Outpost_MainSpine_17` spur, geographically far from
`Circulo_F1_small`/`Circulo_F2_small` — real current network distance is 130-172 km for the
four named mines, against a 10-14 km straight-line separation. A direct connection, computed
with the real friction graph (not a straight-line guess), would cost only 11-17 km — roughly a
90% cut. Candidates and their real cost: `Mine_Schist_05` 15.28 km, `Mine_Schist_09` 17.44 km,
`Mine_Schist_06` (the actual hub `Mine_Schist_05`/`Mine_Schist_08` both attach to) 14.54 km,
to `Circulo_F1_small`; `Mine_Schist_03` 11.48 km, `Mine_Schist_04` 11.80 km, to
`Circulo_F2_small`. Recommending `Mine_Schist_06` and `Mine_Schist_03` as the cheapest real
anchors, but Nico named `Mine_Schist_05`/`09` and `04`/`03` as the acceptable set — his call
which to build. (2) `Coastal_Village_03` and `Mine_Schist_22` "have no roads connecting them"
— confirmed these are NOT new, already correctly reported as 2 of the 5 `still_isolated`
settlements since the eighth follow-up (a genuinely solitary site on the land-only graph, and
a satellite whose hard rule to its own hub has no finite path, respectively) — not acted on
further since nothing changed. (3) `Outpost_MainSpine_16 -> Circulo_E3_2k` — current real
route (via the same `Outpost_MainSpine_17` spur) is 47.90 km; a direct edge would cost 26.36
km/9.54 h, a genuine ~45% improvement. None of these three connections have been built yet —
presented to Nico as findings + recommended anchors, pending his go-ahead and choice of
endpoint. Also found and fixed, while investigating (1), a real but purely cosmetic bug: 12
delivered edges carried a generic internal label (`"mountain trail (massif main_spine)"`,
`"hub cluster: Mine_Schist_17"`) in `connects_to` instead of the real trail segment/spoke they
tied into — connectivity itself was never affected, only the description. See "09"'s ninth
follow-up for the root cause and fix.

**Implemented, same day, after Nico confirmed the anchors**: `Mine_Schist_06 ->
Circulo_F1_small` (14.59 km/4.93h), `Mine_Schist_03 -> Circulo_F2_small` (11.47 km/4.19h),
`Outpost_MainSpine_16 -> Circulo_E3_2k` (26.43 km/9.54h) — three new `manual_connection`
edges, additive to (not replacing) the automated long-detour connections already in place.
See "09" row's ninth follow-up for the implementation detail. New network total: 16
mountain-trail edges (145.9 km) + 75 network spokes (401.8 km) + 11 isolated-pocket roads
(103.6 km) + 3 manual connections (52.5 km) = 703.8 km.

**Tenth follow-up, 2026-08-23 — `Coastal_Village_03` real fix, via `land_mask.npy`
reconciliation.** Nico moved `Coastal_Village_03` ~67m west in his own GIS tool, expecting it
to now sit on shore and be connectable; re-ran the production script against his edited
`auxiliary_settlements_tappa10_v2.geojson` rather than trusting the coordinate change, and the
village stayed isolated — both the old and new coordinates fall inside the identical 120m
routing-grid cell, so the move never crossed a cell boundary. Nico chose to fix the underlying
`land_mask.npy`/true-ocean mismatch rather than nudge the coordinate again; full method and the
majority-vs-`block_any` design correction are in "09"'s tenth follow-up.

**Result after the reconciled mask + backbone re-run**: re-ran `run_tappa10_network_connections.py`
against `land_mask_reconciled_v1.npy` and the re-verified (topologically unchanged) backbone.
`Coastal_Village_03` is now connected — `backbone_t_junction_spoke` onto the
`Circulo_E4_2k<->Circulo_F4_small` edge, 4.495 km / 2.47h — and dropped out of
`still_isolated`. `Coastal_Village_11` remains in `still_isolated`, exactly as Nico expected
(it's on an island, no land route exists regardless of the mask); its underlying grid cell now
correctly reads as land under the reconciled mask, fixing the "sitting on top of the
mask"/visual issue he flagged even though it doesn't gain connectivity. Every other
`still_isolated` entry (`Mine_BasinFill_16`, `Mine_Greywacke_04`, `Mine_Schist_22`) is
unchanged — same names, same reasons — confirming the reconciliation didn't touch anything
outside its intended scope. `still_isolated` count: 5 → 4. New network total: 16 mountain-trail
edges (145.9 km) + 76 network spokes (404.5 km) + 11 isolated-pocket roads (103.5 km) + 3
manual connections (52.5 km) = 706.4 km.

**Eleventh follow-up, same day — `Coastal_Village_07`'s direct connection, and
`Circulo_D_20k` reinstated as an island-local connection target.** Nico asked for (a) a
direct `Coastal_Village_07 -> Circulo_C_25k` connection, and (b) reported "roads aren't
connecting to `Circulo_D_20k`, the junction point is displaced" — asked to clarify before
acting, since no line touching `Circulo_D_20k` existed at all to be displaced (see "09"
row's own eleventh follow-up). **Nico's clarification, quoted**: "Ele está excluído das
conexões com outros círculos por estar em uma ilha, mas precisa tornar a ser incluído para
se conectar com as estruturas auxiliares da ilha" — `Circulo_D_20k` stays excluded from the
mainland backbone (Tappa 9's own call, unchanged), but should be a valid LOCAL connection
target for its own 14-member Southland settlement group.

**Root cause, confirmed directly**: this script's own `circulos` list mirrored Tappa 9's
`EXCLUDED_FROM_ROAD_NETWORK` and dropped `Circulo_D_20k` entirely — so it never appeared in
`site_cell_lookup`/`circulo_idx`, had no Círculo-sourced Dijkstra pass, and (critically) was
never seeded into the greedy network-growth loop's initial cell set. The Southland's 14
settlements could therefore only ever find EACH OTHER as a connection target, never
`Circulo_D_20k` itself — exactly the "own internal MST, no real bridge" pattern already
diagnosed for the mountain-hut trail system at the eighth follow-up, just for a different
reason (there it was a seeding-order bug; here `Circulo_D_20k` was never a candidate cell at
all). **Fix**: `Circulo_D_20k` is back in this script's own `circulos` list (all 17 now,
`run_tappa9_road_network.py`'s own exclusion is untouched — a separate script, separate
concern), and its single site cell is seeded into the initial network SEPARATELY from the
Tappa 9 backbone seeding, labeled distinctly so it's never confused with a mainland
connection. This is safe by construction, not by convention: the land-only friction graph
has zero path between `Circulo_D_20k`'s cell and the mainland (re-verified directly after
the change — `cost_distance_from_source` from `Circulo_D_20k` reaches 55,049 cells, all on
its own landmass, `inf` to every mainland Círculo checked), so the seed can only ever be
picked up by another cell on the same isolated island.

**Result, re-run**: 13 of the 14 Southland settlements now genuinely reach `Circulo_D_20k`
(3 direct `circulo_spoke` edges — `Mine_BasinFill_03` 0.86h, `Mine_Volcanic_01` 1.55h,
`Forest_Post_05` 0.59h — the rest chained on via `feeder_t_junction_spoke`/
`satellite_to_hub_spoke`, same as any other settlement reaching the live network).
`Mine_Greywacke_04` remains genuinely isolated — confirmed not a regression, it was already
outside even the old internal-MST pocket (eighth follow-up). The old `isolated_pocket_road`
step now has nothing left to do for the Southland (0 pocket roads built, down from 11 edges/
103.5 km) — this replaces that honest-but-second-best "local MST with no real anchor" with
what Nico actually wanted: a real local network anchored at the Círculo itself. Also built,
same run: `Coastal_Village_07 -> Circulo_C_25k` as a fourth `manual_connection` (4.54h/10.29
km, vs. 9.49h via the automated route) — see "09" row for the comparison numbers. **New
totals**: 16 mountain-trail edges (145.9 km) + 88 network spokes (503.9 km) + 0
isolated-pocket roads (0.0 km) + 4 manual connections (62.8 km) = 712.6 km. `still_isolated`
unchanged at 4 (`Mine_BasinFill_16`, `Mine_Greywacke_04`, `Coastal_Village_11`,
`Mine_Schist_22`) — same names, same reasons, confirming this fix didn't touch anything
outside its intended scope either.
