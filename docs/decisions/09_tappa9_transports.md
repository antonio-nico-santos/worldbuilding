# Tappa 9 — Transports (road-network foundation slice)

Per `07_tappa7_regional_scenario.md` §9's roadmap, Tappa 9 = Transports: roads/rail, kite
buggies, ferries, and the navigability half of dangerous seas (the creature-danger half
stays routed to a future Tappa 10, not this one — `claude/index.md`'s status table, S9's
own text). That bundles five distinct sub-builds. Opened this session with an
AskUserQuestion round rather than starting all five at once — the project's own Tappa 7
retrospective flags exactly this kind of over-bundling as a repeated mistake — and Nico
chose **"Road network foundation"**: predecessor-path extraction on the existing
`cost_distance.py` Dijkstra graph (so it returns real route geometry, not just a scalar
cost) plus a new biome-differentiated travel-cost multiplier, giving a number to the
narrative fact `scenario_reference.md` §18 already stated (Wet Forest costs more to cross
than Grassland) but never quantified.

Rail's grade-ceiling cost function, ferry corridor authoring, the kite wind-shadow mask,
and sea navigability are **not** covered here — still open, see §9.

**This doc covers THREE passes in one session, plus four same-day follow-ups.** The first pass
(§1-§4 below, as originally shipped) had four real problems Nico caught by opening the
output in QGIS — §6 is the second-pass fix for all four. That second pass then drew two
more findings from Nico (§7): one still-crossing-the-sea line that needed a different fix
than expected, and three specific redundant edges called out as excessive — §8 is the
third-pass fix for both. Read §8 for what actually shipped; §1-§6 are kept as the record of
what each earlier pass built and why, not retracted. **§10**: Nico asked how the three
UNREVIEWED friction tables could be calibrated — a sensitivity analysis (not a calibration,
see §10 for why that distinction matters) answering which of the 21 friction/topology
parameters actually move the deliverable. **§11**: same day, Nico then asked whether
Strahler order implies a river's real width/depth — it doesn't, directly, and answering that
properly led to replacing the river-crossing friction layer's order-bucket table with a
continuous model driven by an actual per-reach width estimate. **§12**: Nico then noticed
routes crossing several minor streams for free in a short span — the order>=4 friction gate
was hiding 92% of the network's real crossings; widened to all orders. **§13**: Nico then
noticed a missing direct connection between two nearby Circulos (F1/F2) — checked directly
against `add_redundant_edges`' own two gates, found a genuine edge case (clears the
shortcut-improvement bar, fails the cheapest-neighbour-proximity bar on both ends because
both sites share a cheap hub), added as one explicit manual exception rather than loosening
the general rule.

## 1. Predecessor extraction — `cost_distance.py`, additive only

Tappa 6/8 only ever needed COSTS (is this candidate ≥ some hour threshold from an
already-placed site); nothing before this session asked scipy's Dijkstra for the
shortest-path TREE it already computes as a side effect. Added, with
`cost_distance_from_source` (used by every already-locked Tappa 6/8 result) left
byte-for-byte untouched:

- `cost_distance_from_source_with_predecessors(graph, row, col, shape, limit_hours=inf)` —
  same single-source Dijkstra, `return_predecessors=True` this time.
- `reconstruct_path(predecessors, shape, source_row, source_col, target_row, target_col)` —
  walks scipy's raw predecessor array back from target to source (scipy's own chain runs
  target→source; reversed here to source→target), returning `None` for an unreached
  target rather than a garbage partial path.

## 2. Biome travel-friction — `src/transport/biome_friction.py` (new module)

Same architectural role as Tappa 8 S8f's `LAND_TRAVEL_FRICTION` (lithology/rock-type
friction) — a per-cell 0–1 multiplier on LAND-LAND edge speed in `build_cost_graph` — but
keyed on `biome_id` (vegetation/ground cover) instead of lithology (rock type). The two
are independent physical properties (what's growing on the ground vs. what the ground is
made of underneath) and are combined **multiplicatively**, not substituted for one
another — same "friction stacks" logic S8f/S8g already established for
lithology-friction × excavation-effort, applied here to two friction layers instead of one
friction layer and one non-traversal index.

**`BIOME_TRAVEL_FRICTION` — UNREVIEWED first-pass estimates, pending Nico's sign-off, not
written to `config/parameters.yml`** (same status as `LAND_TRAVEL_FRICTION` and
`EXCAVATION_EFFORT_MULTIPLIER`):

| biome | multiplier | grounding |
|---|---|---|
| Lowland Steppe/Grassland | 1.00 | baseline/anchor — open, treeless, minimal clearance |
| Woodland/Shrubland | 0.90 | open mosaic, patches of lower-stature scrub/forest, not a closed canopy |
| Alpine Fellfield | 0.90 | open, no canopy, but loose scree/talus underfoot — footing penalty only |
| Subalpine Dry Scrub | 0.85 | low sparse woody scrub, exposed ground between clumps — thorny in patches but not canopy-dense |
| Alpine Tundra | 0.75 | open at a distance, but real dense snow-tussock (0.5–1.5 m) is documented as slow, tiring walking |
| Subalpine Woodland | 0.75 | same temperature band as Wet Forest but a visibly more open canopy, drier/sparser understory |
| Subalpine Wet Forest | 0.55 | closed canopy, heavy epiphyte/moss load, damp low-light understory — real off-track native-bush travel is notoriously slow |
| Temperate Forest | 0.55 | same value as Subalpine Wet Forest on purpose — both closed-canopy, heavy-understory forest; nothing in this project's sources differentiates their off-track difficulty specifically |
| Permanent Snow & Ice | 0.35 | crevasse/avalanche hazard — real glacier travel needs roping/probing even on gentle slope, independent of whatever the DEM slope term already charges |

Not a citation of a specific per-biome hiking-speed dataset — none exists. Directional
judgement calls grounded in this project's own `biome_landscape_characteristics.md`
vegetation-structure descriptions plus real off-track-travel accounts for the closest
matching real vegetation structure, same disclosure standard S8f already set for
lithology friction. Biome 0 (Ocean) is deliberately absent from the table — sea edges
never consult it (`build_cost_graph` only applies `friction_multiplier` to LAND-LAND
edges).

## 3. Road-network topology, first pass — `src/transport/network.py`

**Original topology decision (first pass): minimum spanning tree (Prim's algorithm) over
the 17 already-placed, LOCKED Circulo sites**, weighted by pairwise cost-distance (hours)
on the combined lithology × biome friction-adjusted graph. An MST is the
minimum-total-travel-time way to connect every site into one network — a defensible,
cheaply-explained "foundation" — but a deliberate simplification worth being explicit
about: real infrastructure networks often carry redundant loops, especially around large
hubs, that an MST by definition never has (exactly N−1 edges for N nodes, zero cycles —
losing any single edge disconnects the network). **Superseded in §6 below** after Nico
flagged the resulting network as "always linear" even where several Circulos cluster
closely.

`cost_distance.py`'s graph is directed (uphill ≠ downhill Tobler cost over the same edge),
so the raw pairwise-hours matrix isn't symmetric. The MST needed one scalar weight per
unordered pair, so `build_mst_edges` used the **mean of the two directions** — the natural
choice for a road walked/ridden both ways, not just once; this convention carries over
unchanged into §6's revision.

**The 17 Circulo sites are read-only, locked input — never re-placed here**, same
discipline S8f/S8g already established for the lithology-friction and excavation-effort
layers, and consistent with `08_tappa8_geomorphology.md` §15's explicit resolution that
Circulo placement is final for all downstream Tappas.

## 4. Result, first pass — `run_tappa9_road_network.py`

120 m grid: 1334×1084 (cellsize ≈119.93 m), reproducing Tappa 6/8's grid exactly
(`land_mask & ~lake_mask`, DEM block-mean, lithology block-mode). Land-mean friction:
lithology 0.902, biome 0.804, combined (lithology × biome) 0.727.

MST produced 16 edges for 17 sites. Total route length 388.0 km vs. 348.0 km straight-line
sum. Friction vs. baseline: mean +8.13%, median +3.16%, max +57.59%
(`Circulo_C_25k`↔`Circulo_F6_small`).

**This is the pass Nico reviewed in QGIS and found four real problems with — see §6.**

## 5. Nico's review (2026-08-20) — four problems found, one non-issue

Direct QGIS inspection of the first pass's `road_network_mst.geojson` surfaced:

1. **6 of 16 edges crossed open ocean.**
2. **Of the 10 that didn't cross ocean, about 8 crossed a lake.**
3. **No consideration of rivers at all**, at minimum higher-Strahler-order ones.
4. **In areas with several nearby Circulos, the network was always linear** — no loops,
   even where geometry obviously allowed one.
5. **The exported CRS wasn't the project's custom LCC** — flagged by Nico as possibly a
   QGIS-side artifact rather than a real bug.

Root causes, confirmed before writing any fix (see §6 for the fixes themselves):

- **Problems 1+2, same root cause.** The first pass ran the road MST over
  `cost_distance.py`'s general-purpose graph, which — correctly, for Tappa 6's ORIGINAL
  isochrone/tier-distance siting use case — treats any non-land cell (ocean OR lake; this
  module has never distinguished the two) as boat-traversable at a flat
  `BOAT_SPEED_KMH = 6.0`. That's the right model for "how far away is this candidate site,"
  and the wrong one for "where can a road be built" — a road can't be laid across open
  water the way a distance check can cross it. Confirmed directly by re-checking the first
  pass's own output: exactly 6 of 16 edges touched a non-land cell with zero land alternative
  cheap enough to win the Dijkstra search, and 8 of the remaining 10 touched at least one
  lake cell.
- **Problem 3, confirmed.** `cost_distance.py`'s graph has only ever known elevation and
  land/sea — zero hydrology awareness, so no edge anywhere in the first pass paid any
  penalty for crossing even the domain's largest rivers.
- **Problem 4, confirmed but not really a "bug"** — the honest MST limitation was already
  flagged in §3/§4 above as an open follow-up before this review; Nico's report is what
  turned it from a documented limitation into something to actually fix this session.
- **Problem 5 — confirmed as NOT a code bug.** This project's GeoJSON exports have carried
  a legacy `crs`/`proj4` FeatureCollection member since early Tappas specifically because
  RFC 7946 (the modern GeoJSON spec) mandates WGS84 lon/lat and has no accommodation for a
  custom projected CRS — QGIS, as an RFC-7946-compliant reader, ignores that member and
  reports EPSG:4326, exactly as `07_tappa7_regional_scenario.md` already documented for
  earlier Tappas' own vector exports. `run_tappa9_road_network.py` follows the exact same
  convention every other GeoJSON export in this repo uses; nothing about this run
  introduced a new or different CRS problem. Workaround is the existing one, unchanged:
  assign the CRS manually in QGIS after import (Layer Properties → Source → Assigned CRS,
  paste `config/parameters.yml`'s PROJ string) — the same manual step
  `terrain/raster_io.py`'s own docstring already documents for the ENVI raster exports,
  for the same underlying reason (no GDAL/rasterio in this sandbox to embed a real CRS,
  see that module's docstring).

Two design questions were put to Nico directly (AskUserQuestion) rather than guessed at:
whether the land-only-graph-plus-flag-disconnected-pairs plan for problems 1-3 was the
right shape (confirmed: proceed), and how to fix problem 4 (confirmed: "MST + local
redundancy" over the alternatives of a denser near-neighbor mesh or leaving the pure MST
undocumented-but-unfixed).

## 6. Second pass — land-only graph, river friction, MST-forest + redundancy

**6a. Land-only road graph — `cost_distance.py`'s new `sea_mode` parameter.**
`build_cost_graph` gains `sea_mode: str = "boat"` (default reproduces every already-locked
Tappa 6/8 call site byte-for-byte) vs. `sea_mode="impassable"` (new): any edge touching a
non-land cell gets `cost = inf` and is dropped from the sparse graph entirely, rather than
priced at boat speed. The road network is now built exclusively on an
`sea_mode="impassable"` graph — no roads edge can cross open water or a lake, structurally,
not just by convention.

**Direct consequence, handled explicitly, not hidden:** cutting every sea/lake edge can
split the 17 sites into more than one connected component. It does, on this domain:
**`Circulo_D_20k` is fully isolated from the other 16 sites by land** — confirmed directly
(`connected_components` on the land-only pairwise-cost matrix), not assumed. This is a
genuine, interesting finding, not an artifact: `Circulo_D_20k` sits on the same landmass
`08_tappa8_geomorphology.md`'s "island split into volcanic-core/basin-fill-margin" lithology
zone already described — nothing before this session had actually confirmed that split as
*land-transport-disconnected* until now.

**6b. Candidate ferry crossing.** Rather than either forcing a fake "road" across the water
(the original bug) or silently dropping `Circulo_D_20k` from the network, the run computes
the cheapest boat-enabled crossing linking every disconnected component (here, just the one
pair) and records it as a SEPARATE GeoJSON feature, `edge_type: "candidate_ferry_crossing"`
— real route geometry (reused from the boat-enabled combined-friction graph, not
re-Dijkstra'd), but explicitly informational, explicitly not a road:
`Circulo_D_20k`↔`Circulo_E3_2k`, 10.88 boat-enabled hours, 57.7 km. This is a preview of
what Tappa 9's still-not-built ferries sub-build will need to formalize — a real ferry
corridor is a separate design decision (dock siting, schedule, capacity), not something this
run invents.

**6c. River-crossing friction — `src/transport/river_friction.py` (new module).** A THIRD
friction layer, same multiplicative-stacking pattern as lithology/biome, sourced from Tappa
4's own `data/exports/streams.geojson` (17,678 reaches, Strahler-ordered 1-6, not previously
consulted by any transport-cost code). Per Nico's own framing ("at least those of higher
Strahler order"), only order ≥ `MAJOR_STREAM_MIN_STRAHLER_ORDER = 4` reaches (1,847 of
17,678, ~10.4%) get a penalty — lower orders are treated as fordable, the same
simplification a real rural road network makes (a track fords a creek; only a real river
gets a bridge). Rasterized directly onto the 120 m grid by densifying each qualifying
reach's LineString to sub-cell spacing and marking every cell it passes through (same
`_densify_polyline` approach `terrain/skeleton.py` already uses for ridge/zone geometry) —
no native-30 m intermediate raster needed.

**`RIVER_CROSSING_FRICTION` — UNREVIEWED first-pass estimates, pending Nico's sign-off, not
written to `config/parameters.yml`** (same status as every other friction table in this
project):

| Strahler order | multiplier | grounding |
|---|---|---|
| 4 | 0.75 | smallest "major" class — a real bridge, but a narrower one; mildest of the three |
| 5 | 0.60 | a substantial trunk tributary — real unbridged-river-flat accounts describe a major detour to a ford, not a minor slowdown |
| 6 | 0.45 | the domain's largest rivers by discharge (up to ~103 m³/s, Tappa 4's own summary) — harshest penalty here on purpose |

Result: 9,640 of 120 m effective-land cells (1.22%) carry a river-crossing penalty;
land-mean river friction 0.996 (a small headline number by construction, since only ~1% of
cells are affected at all — the number that matters is the per-crossing penalty at the
specific cells a route actually touches, not the domain-wide mean, same caveat S8f already
raised about its own land-mean friction figures not being the interesting statistic).

**6d. Topology — minimum spanning FOREST + local redundancy —
`src/transport/network.py`, revised.** `build_mst_edges` → `build_mst_forest`: runs Prim's
independently within EACH connected component (handles `Circulo_D_20k`'s isolation
correctly instead of assuming one connected graph). On top of the forest,
`add_redundant_edges` (the AskUserQuestion-confirmed fix for problem 4) adds extra edges
wherever a site has a second connection within `redundancy_factor=1.4` (40%) of its
cheapest — a tight cluster of mutually close sites naturally gets several near-tied cheap
options, so loops appear exactly where Circulos cluster (confirmed: 8 redundant edges
landed among the E/F cluster sites specifically, none along the network's genuinely linear
stretches). `redundancy_factor=1.4` is itself an UNREVIEWED first-pass number — middle of
the "~30-50%" range discussed with Nico when this approach was picked, not independently
calibrated.

**6e. Result, second pass.** Land-only connectivity: 2 components (16 sites + the isolated
`Circulo_D_20k`, see §6a). MST-forest: 15 edges (16-site component: 14, isolated component:
0, as expected for spanning forests — 17 sites − 2 components = 15 total). Redundancy pass:
+8 edges. **23 road edges total**, zero of which touch ocean or lake (checked directly,
cell-by-cell, against every edge's actual reconstructed route — not assumed from the
graph's construction alone). Plus the 1 candidate ferry crossing from §6b.

Total road route length: 573.9 km (vs. 503.9 km straight-line sum for the same 23 edges —
higher than the first pass's ratio, expected: land-only routing can no longer shortcut
through water, so real paths detour more). Friction vs. baseline (boat-enabled, no
friction) across the 23 road edges: worst single-edge delta +73.7% — larger than the first
pass's own worst delta, expected now that river-crossing friction is a real, if narrow,
contributor stacked on top of lithology/biome. Every figure is in
`tappa9_road_network_meta.json`, including the full connected-components record and the
candidate-ferry-crossing detail, not just this summary.

**What this still does NOT do, on purpose:** re-place the 17 Circulos under this friction
(same reasoning `08_tappa8_geomorphology.md` §15 already gave for declining it with the
lithology-only friction); model actual construction/excavation cost (this is travel-TIME
friction, not a build-cost model); author a real ferry corridor (§6b's candidate crossing
is informational, not a designed corridor); or touch rail, kite buggies, or sea navigability
at all.

**Outputs** (`data/processed/transport/`, gitignored, regenerate locally):
`road_network_mst.geojson` (24 features: 15 `mst` + 8 `redundant` road edges, 1
`candidate_ferry_crossing`, real route geometry throughout, per-edge
hours/route-km/straight-line-km/tier-pair/edge_type properties),
`biome_friction_multiplier_120m.*`, `river_crossing_friction_multiplier_120m.*` (new this
pass), `combined_friction_multiplier_120m.*` (lithology × biome × river, what the land-only
road graph was actually built from), `tappa9_road_network_meta.json` (full
comparison/connectivity/candidate-ferry/friction-table record).

## 7. Nico's second review (2026-08-20, same day) — two more findings

Direct review of the second pass's `road_network_mst.geojson` (delivered per §6) surfaced
two more things, both genuine, neither a repeat of §5's problems:

1. **`Circulo_D_20k` still draws a line crossing the sea.** Correct as far as §6b's own
   design went — that line was the deliberately-informational `candidate_ferry_crossing`
   feature, not a road — but Nico's read, once he saw it rendered: since `Circulo_D_20k` is
   genuinely on an island, it shouldn't be part of this ROAD network's processing AT ALL,
   not even as a clearly-labelled candidate line. Explicit instruction: **lock it outside of
   the processing.**
2. **A few redundant edges are excessive.** Named directly by GeoJSON feature ID from that
   run's export: **15, 18, 20** — `Circulo_A_40k`↔`Circulo_E4_2k`, `Circulo_C_25k`↔
   `Circulo_E2_2k`, and `Circulo_E3_2k`↔`Circulo_F8_small` respectively. The other 5
   redundant edges from that same run were NOT flagged.
3. **Confirmed working correctly and not reopened:** sea, lake, and river crossings all
   solved by the second pass — Nico's own words, "(with Círculo D exception)" for the sea
   part specifically (i.e., the one remaining sea-crossing line was always understood to be
   the labelled ferry candidate, not a road bug — see finding 1's resolution below for why
   it's gone now anyway).

Investigated before writing any fix: pulled the exact tree-path cost each of the 8 redundant
edges from §6d would have saved over the route the MST tree already provided between the
same two sites (`_tree_path_cost`, new this pass — see §8b). Two of the three flagged edges
(`Circulo_A_40k`↔`Circulo_E4_2k`, `Circulo_E3_2k`↔`Circulo_F8_small`) cost within **0.1% and
0.003%** of their existing tree-path cost respectively — i.e., essentially ZERO real
benefit, a second line on the map that gets you nowhere faster. The third
(`Circulo_C_25k`↔`Circulo_E2_2k`) saved a real but modest 16.2%. The 5 edges Nico did NOT
flag saved 21.0%–43.4% — a clean, sizeable gap between "flagged" and "not flagged," meaning
the underlying complaint had an exact, measurable signature: §6d's `redundancy_factor`-only
criterion checked whether a second connection was cheap relative to a site's OWN cheapest
neighbour, but never checked whether the edge actually shortened anything against the
network that already existed.

## 8. Third pass — exclude `Circulo_D_20k`, tighten redundancy to genuine shortcuts only

**8a. `Circulo_D_20k` excluded from road-network processing entirely.** A new
`EXCLUDED_FROM_ROAD_NETWORK` set in `run_tappa9_road_network.py` removes it from the `sites`
list before ANY of this run's computation — not filtered out of the output afterward, never
computed against at all. Its previously-confirmed cheapest boat-enabled connection to the
rest of the network (`Circulo_D_20k`↔`Circulo_E3_2k`, 10.8757 h, 57.7 km, from §6b's run) is
preserved in `tappa9_road_network_meta.json`'s new `excluded_sites` field, specifically for
whoever eventually builds the real ferries sub-build — not thrown away, just not drawn as a
placeholder line in a ROAD deliverable. The candidate-ferry-crossing mechanism itself
(§6b, `src/transport/network.py`'s component-Prim's-MST-over-supernodes logic) is left in
place, generic and unchanged — with `Circulo_D_20k` excluded, the remaining 16 sites are
confirmed to form exactly 1 connected component on the land-only graph (checked directly,
not assumed), so it simply finds 0 crossings to compute this run. If a future friction
change or a newly-added site ever reintroduces a disconnection among the sites this script
DOES process, the mechanism is still there to catch it and flag it the same honest way.

**8b. Redundant-edge criterion tightened — three approaches tried, one shipped, in
`src/transport/network.py`'s `add_redundant_edges`** (each tested by actually running it and
counting edges before deciding, not assumed from reasoning alone):

- **Attempt 1 (rejected, over-permissive in a new way):** replace §6d's criterion outright
  with a tree-path-shortcut check (`_tree_path_cost`, new helper — walks the MST to find
  what a pair of sites already costs to connect through the existing tree), but check EVERY
  non-MST pair in the graph (not just each site's nearby candidates). Result: **26 redundant
  edges** for 16 sites — far worse than the 8 being fixed. Root cause: a direct edge between
  two sites on opposite sides of the whole network will almost always beat their long,
  multi-hop tree-path sum by a wide margin — a property of hop-count arithmetic, not
  evidence the edge is a locally useful redundant connection. Caught immediately by running
  it.
- **Attempt 2 (rejected, unreviewed churn):** keep §6d's local candidate scoping (each
  site's next-`max_extra_per_site`-cheapest neighbours only) but SWAP the acceptance rule
  for the tree-path-shortcut check instead of ANDing them. Result: still 8 edges, but NOT
  the same 8 minus the 3 flagged — it dropped the 3 correctly but also newly accepted 3
  DIFFERENT candidates §6d had rejected (cheap enough by the tree-shortcut measure, but not
  "close to that site's own cheapest neighbour" in §6d's sense), landing on an unreviewed
  edge set nobody had looked at. Rejected specifically because introducing new, unreviewed
  changes wasn't what was asked — the goal was removing 3 specific bad edges, not
  redesigning the output.
- **Attempt 3 (shipped):** require BOTH conditions together — `redundancy_factor=1.4` (§6d's
  original proximity scoping) AND `min_shortcut_improvement=0.20` (a genuine tree-path
  shortcut of at least 20%, itself an UNREVIEWED first-pass constant picked to land between
  the highest dropped edge, 16.2%, and the lowest kept one, 21.0%, per §7's investigation).
  Checked directly: this keeps EXACTLY the 5 edges Nico did not flag and drops EXACTLY the 3
  he did, introducing zero new candidates beyond the already-reviewed set of 8 — the
  conservative fix, not a redesign.

**8c. Result, third pass.** 16 sites processed (17 locked Circulos minus `Circulo_D_20k`,
§8a). MST-forest: 15 edges, 1 connected component (confirmed, not assumed). Redundancy pass:
5 edges — `Circulo_A_40k`↔`Circulo_F7_small` (21.0% shortcut), `Circulo_B_35k`↔
`Circulo_F4_small` (35.7%), `Circulo_E3_2k`↔`Circulo_E4_2k` (39.2%), `Circulo_F3_small`↔
`Circulo_F4_small` (43.4%), `Circulo_F5_small`↔`Circulo_F7_small` (38.4%). **20 road edges
total, zero candidate ferry crossings** (§8a), zero edges touching ocean or lake (checked
cell-by-cell against every edge's real reconstructed route, same verification standard as
§6e). Total road route length: 449.9 km (down from §6e's 573.9 km — expected, both from
dropping the 3 zero/low-benefit redundant edges and from `Circulo_D_20k` no longer
contributing any edges at all). Worst single-edge friction delta vs. baseline unchanged at
+73.7% (same edge, `Circulo_E4_2k`↔`Circulo_F4_small`, untouched by either of this pass's
two fixes). Full detail, including the tree-path-improvement percentage for every candidate
considered (not just the 5 that made it in), is in `tappa9_road_network_meta.json`.

**Outputs, updated in place** (`data/processed/transport/`, same file names as §6e, this
pass's numbers now current): `road_network_mst.geojson` (20 features: 15 `mst` + 5
`redundant`, zero `candidate_ferry_crossing` this run), `tappa9_road_network_meta.json`
(now carries `excluded_sites` alongside every field from §6e). The friction rasters
(`biome_friction_multiplier_120m.*`, `river_crossing_friction_multiplier_120m.*`,
`combined_friction_multiplier_120m.*`) are unchanged from §6e — this pass touched only site
inclusion and edge topology, not the friction fields themselves.

## 9. Open follow-ups (not done this stage, deliberately left open)

- **Rail's grade-ceiling cost function, actual ferry corridor authoring, the kite
  wind-shadow mask, and sea navigability** — the other four Tappa 9 sub-builds from `07`'s
  roadmap, none started. §6b's candidate ferry crossing is a preview, not the real ferry
  sub-build. Per this session's own AskUserQuestion framing, not to be started without
  further direction.
- **All three friction tables (§2 biome, §6c river; plus S8f's lithology) — UNREVIEWED,
  pending Nico's sign-off.** Not written to `config/parameters.yml`. **See §10**: a
  sensitivity analysis identifies which specific table entries are actually worth spending
  calibration effort on for this network; it does not itself resolve the sign-off.
- **`redundancy_factor=1.4` and `min_shortcut_improvement=0.20` (§8b) — both UNREVIEWED
  first-pass constants.** The second was picked to land cleanly between one specific
  flagged edge (16.2%) and one specific kept edge (21.0%) in one review round — reasonable
  as a first pass, but a genuinely different sample of "which redundant edges look
  excessive" could shift where that line should sit. **§10 confirms `redundancy_factor` is
  in fact the single most consequential parameter in the whole friction/topology set** (up
  to 5 of 20 edges change under a ±30% probe) — the constant most worth real attention if
  this gets revisited.
- **`Circulo_D_20k`'s isolation is resolved for THIS script's processing (§8a — excluded
  entirely), but not as a narrative/worldbuilding fact.** Should this island always have
  been ferry-only? Does that match how it was described when placed in Tappa 6? Not decided
  here — worth Nico's read whenever the real ferries sub-build gets picked up, since
  `excluded_sites` in `tappa9_road_network_meta.json` is exactly the fact-record that sub-
  build will need to start from.
- **The `Circulo_C_25k`↔`Circulo_F6_small` 57.59% delta (first pass, §4)** and the
  second pass's own worst-edge +73.7% delta — largest friction effects observed, not
  investigated beyond flagging.
- **Re-placing the 17 Circulos under combined friction** — explicitly declined here, same
  as S8f declined it for lithology-only friction; would be a materially bigger, separate
  decision if ever taken up.
- **The dangerous-seas creature/conflict half** — stays routed to a future Tappa 10, per
  `claude/index.md`'s status table; only the navigability half is in Tappa 9's scope at
  all, and this session didn't design any real navigability model (the boat-speed cost
  edges reused for §6b's candidate ferry are the same flat, unmodelled `BOAT_SPEED_KMH`
  `cost_distance.py` has always used).

## 10. Sensitivity analysis — answering "how do we calibrate these tables" (2026-08-20, same session)

Nico's question, verbatim in substance: the road network "looks right visually," but there's
no idea how the friction values (§2 biome, §6c river, S8f lithology) or the two topology
constants (§8b) should actually be calibrated. Direct answer given in chat, recorded here:
**"calibration" in the usual sense — fitting parameters against measured ground truth — does
not apply.** This is a fictional world; there is no dataset of real travel times to fit
against. "Looks right visually" is a weak, non-falsifiable signal on top of that (it's
susceptible to confirmation bias — a result already expected to look a certain way will tend
to look that way). Three genuinely different paths exist instead: (1) anchor values against
real-world least-cost-path/off-road-mobility literature, (2) anchor against a quantitative
narrative fact if one exists (none currently does — `scenario_reference.md` §18's "Wet Forest
costs more than Grassland" is ordinal, not a number), (3) a sensitivity analysis — vary each
parameter and measure how much the deliverable actually changes, which tells you WHERE
further calibration effort (via 1 or 2) would pay off, without requiring any new data. Nico
chose (3) only, this round.

**Method** (`run_tappa9_sensitivity_analysis.py`, new script, reuses
`run_tappa9_road_network.py`'s own machinery — `build_cost_graph`,
`travel_friction_multiplier`, `biome_friction_multiplier`, `river_friction_multiplier`,
`compute_pairwise_cost_distance`, `build_mst_forest`, `add_redundant_edges` — not a separate
approximation of the pipeline): one-at-a-time (OAT) perturbation. All 21 parameters — 7
lithology classes, 9 biome classes, 3 river-crossing Strahler-order classes, plus
`redundancy_factor` and `min_shortcut_improvement` — were each independently probed at ×0.7
and ×1.3 of their current shipped value (friction multipliers clipped to (0.05, 1.0]),
holding every other parameter at its shipped value, with the full road-network pipeline
rerun from that single change each time (land-only graph → all-pairs cost-distance → MST
forest → redundant edges). 42 probe runs + 1 baseline, ~4s each after a shared ~5s input-load
(158s total) — the boat-enabled and baseline-comparison graphs `run_tappa9_road_network.py`
also builds were skipped entirely since neither affects the road network's own topology,
which is the only thing under test.

**Baseline** (current shipped third-pass values, confirmed identical to §8c's own run): 15
MST + 5 redundant = 20 edges, 107.44 h total travel time, 395.6 km (site-to-site straight-line
sum, used as the sensitivity metric instead of real route km — cheaper to compute per probe
and monotonic with real route km for this purpose).

**Headline finding: the MST backbone never changed in any of the 42 probe runs.** Every one
of the 15 MST edges survived every ±30% perturbation of every one of the 21 parameters. Only
the 5 redundant edges (the ones `add_redundant_edges` gates on a threshold, see §6d/§8b) ever
moved — expected in hindsight, since a threshold-gated heuristic is inherently the most
fragile piece of the whole pipeline, but confirmed directly rather than assumed.

**9 of 21 parameters change at least one edge under a ±30% probe (ranked by max edges
changed):**

| Parameter | Layer | Max edges changed |
|---|---|---|
| `redundancy_factor` | topology | 5 |
| Woodland/Shrubland | biome | 3 |
| Lowland Steppe/Grassland | biome | 3 |
| sedimentary_basin_fill | lithology | 2 |
| Temperate Forest | biome | 2 |
| greywacke_argillite | lithology | 1 |
| granite | lithology | 1 |
| Subalpine Wet Forest | biome | 1 |
| `min_shortcut_improvement` | topology | 1 |

`redundancy_factor` is the standout: dropping it to 0.98 (×0.7) removes all 5 redundant
edges outright (−26.5% total travel-time-weighted cost across the network, since those
edges exist specifically to shorten some routes); raising it to 1.82 (×1.3) adds 2 more
candidate redundant edges (`F1_small↔F2_small`, `F2_small↔F4_small` — the same two edges
flagged and rejected as "unreviewed churn" during §8b's attempt-2, now confirmed to be
right on the acceptance boundary, not an arbitrary near-miss). Every other sensitive
parameter only ever changed edges inside the same dense E2-E3-E4-F cluster where
`add_redundant_edges` is actively choosing between several close-cost alternatives — this is
where the redundancy layer, not the MST, does its work, so it's exactly where friction-table
precision matters most.

**12 of 21 parameters never change the edge set** (cost-only effect, ranked by max total-hours
delta): Subalpine Dry Scrub (1.41%), schist (0.69%), Subalpine Woodland (0.63%), river orders
4/5/6 (0.15%/0.10%/0.01%), then **6 parameters with a flat 0.000% effect in either direction**:
lithology volcanic/marble/sedimentary_limestone, and biome Permanent Snow & Ice/Alpine
Fellfield/Alpine Tundra. Flagged explicitly: zero effect here means no cell any of this
16-site road network's actual routes pass through carries that class — a fact about THIS
network's specific geometry (small/high-elevation zones the chosen routes happen to avoid),
not a general claim those six values don't matter. A future rail sub-build crossing a
mountain pass, for instance, could easily touch Permanent Snow & Ice or granite terrain the
road network never needed to route through.

**What this does and doesn't establish.** It ranks where calibration effort (real-world
literature anchoring, or a narrative fact if one surfaces) would change the actual
deliverable, and where it provably wouldn't for this specific network. It does NOT say what
the correct value of any parameter is — no table was changed as a result, and all 21 values
remain exactly as shipped in §2/§6c/S8f/§8b. If this gets revisited, `redundancy_factor` is
the one constant worth spending real effort on first; the 6 zero-effect classes are the ones
safe to leave alone indefinitely without re-checking, at least until a sub-build routes
through the terrain they cover.

**Output**: `data/processed/transport/tappa9_sensitivity_analysis.json` — full results for
all 42 probe runs (which edges added/removed, hours/km delta per probe), plus the two ranked
tables above, `method`/`purpose` fields spelling out the "this is not calibration" framing
for any future reader who finds the file without this doc.

## 11. River width/depth from discharge — replacing the order-bucket friction table (2026-08-20, same day)

Nico's question, verbatim in substance: is it possible to know a river's approximate
width/depth from its Strahler order? **Direct answer: only indirectly.** Strahler order is a
TOPOLOGICAL classification (how many tributaries merge upstream of a reach) — it correlates
with drainage area, and drainage area correlates with discharge, and it's discharge that
actually determines channel size, via real hydraulic-geometry relations. Order is a proxy
for a proxy, not a measurement.

**This domain's own data confirms the proxy is leaky.** `streams.geojson`'s per-reach
`max_discharge_proxy_m3s`, checked directly (not assumed): order-4 reaches range from
0.11 to 84.15 m3/s (mean 7.85), order-5 from 1.11 to 101.76 (mean 17.81), order-6 from 4.18
to 103.13 (mean 47.97) — real overlap between classes, and a genuinely tiny 0.11 m3/s
minimum inside the "major, needs a bridge" order-4 bucket. That minimum is almost certainly
a short reach crossing the order-3/4 topological threshold immediately after one small
tributary joins, not a physically large river — but §6c's flat per-order table
(`RIVER_CROSSING_FRICTION`) penalized it exactly as hard as a typical order-4 crossing,
because order was the only thing that table could see.

**Fix: `src/transport/river_geometry.py` (new module) estimates width/depth from the
discharge this project ALREADY computes**, rather than adding a new hydrology model.
Leopold & Maddock (1953) downstream hydraulic geometry: `W = a * Q^b`, `D = c * Q^f`. The
EXPONENTS (`b=0.5` width, `f=0.4` depth) are among the most widely reproduced results in
fluvial geomorphology — used here directly, real physical basis, not a guess. The
MULTIPLICATIVE COEFFICIENTS (`a`, `c`) are genuinely NOT a citable universal constant
(published values vary ~2-4x by region/channel type) — rather than borrowing a real-world
number and presenting it as locally validated, they're derived by anchoring to one explicit,
checkable judgement: at this domain's actual largest known discharge (103.128 m3/s, the real
max across all order-6 reaches), a confined single-channel temperate gravel-bed river is
plausibly 45 m wide, 2.8 m deep — ordinary magnitudes for that class of river, explicitly
NOT a braided system (a real South-Island-NZ analog like the Rakaia at flood runs to
hundreds of metres across multiple channels; this domain's rivers are treated as
single-thread, consistent with `cost_distance.py` never having modelled braiding). **This
anchor is a calibration CHOICE, UNREVIEWED, same status as every other friction estimate in
this project** — flagged explicitly in the module's own docstring, not presented as a
citation.

**`river_friction_multiplier_from_width` (new, in `river_friction.py`) replaces the
discrete 3-bucket lookup with a continuous function of estimated width**, deliberately
ANCHORED to the old table's own two calibrated severity judgements rather than picking new
numbers from scratch: `friction(W) = 1 - k*W^p`, solved so that (a) the width implied by
order 4's own MEAN discharge (7.85 m3/s → 12.42 m) still gets 0.75 (the old flat value every
order-4 reach got) and (b) the order-6 anchor discharge (103.128 m3/s → 45 m, by
construction) still gets 0.45 (the old harshest value) — floored at 0.45 beyond 45 m (no
data past this domain's largest known reach), ceilinged at 1.0 (~0 discharge = not an
obstacle). The old order-based table and function (`RIVER_CROSSING_FRICTION`,
`river_friction_multiplier`) are kept in `river_friction.py`, documented as superseded, not
deleted — same discipline `network.py`'s `add_redundant_edges` history uses.

**Checked directly, not assumed, what actually changed** (`data/exports/streams.geojson`,
9,640 qualifying cells across the domain's 120 m grid): the new model's mean friction is
HIGHER (less severe) than the old flat value in all three former order classes — order 4:
0.786 vs. 0.750, order 5: 0.699 vs. 0.600, order 6: 0.598 vs. 0.450. This is the expected,
correct direction, not a bug: the old table applied each class's flat value to EVERY cell in
that class, including the many cells on the narrow end of the class's real width range; the
new model only applies the old harshest values to the cells that are ACTUALLY that wide.
6.2% of qualifying cells now score friction > 0.9 (negligible obstacle) that previously got
a flat 0.75/0.60/0.45 regardless of real size; 16.5% of qualifying cells have an estimated
width under 5 m.

**Network re-run, checked directly**: identical 20-edge topology (15 MST + 5 redundant,
same edge set, verified by direct set comparison against §10's recorded baseline — zero
additions, zero removals), total road route length 449.9 → 450.0 km (noise-level, from
route micro-adjustments where a path crosses a now-differently-weighted river cell, not a
topology change). Consistent with §10's own finding that all three river-order friction
values were among the LEAST consequential parameters tested (max 0.15% total-hours delta at
±30%) — moving from discrete-order to continuous-width was a structural change, not a
parameter perturbation, and it's reassuring, not suspicious, that it didn't flip anything:
this network's route choices were never sensitive to the river layer's exact values in the
first place.

**New output, informational**: `data/processed/transport/major_stream_geometry.geojson`
(1,847 reaches, order ≥ 4, `annotate_reach_geometry` — same geometry as `streams.geojson`,
`estimated_width_m`/`estimated_depth_m` properties added) for QGIS variable-width river
symbology and narrative reference. Does NOT modify `data/exports/streams.geojson` itself
(Tappa 4's own locked output). `tappa9_road_network_meta.json` gains
`river_width_depth_model` (method, exponents, anchor, estimated-width summary stats) and
`river_crossing_friction_width_anchors`; the old `river_crossing_friction_table` field is
kept, renamed `..._SUPERSEDED_order_based`, not deleted.

**Outputs updated in place**: `river_crossing_friction_multiplier_120m.*` and
`combined_friction_multiplier_120m.*` (both changed — the river layer's values changed for
the 9,640 qualifying cells) — `biome_friction_multiplier_120m.*` unchanged (biome logic
untouched this revision, not redelivered).

**Still open**: the anchor width/depth (45 m / 2.8 m at 103 m3/s) is a calibration choice,
not independently validated — same open item as every other friction number in this
project, now joined by `src/transport/river_geometry.py`'s two new constants
(`ANCHOR_WIDTH_M`, `ANCHOR_DEPTH_M`). Depth is estimated and exported but does NOT currently
drive friction (only width does, consistent with this layer's own stated reasoning: the cost
being priced is the detour to find/use a crossing point, which width — not depth — mostly
determines at this grid's resolution) — worth revisiting if a future bridge/ford-design
sub-build needs depth to matter mechanically.

## 12. All-order river crossings — the order>=4 gate was hiding 92% of the network's real crossings (2026-08-20, same day)

Nico asked for the width/depth export covering the FULL stream network, then explained why:
they'd noticed, inspecting the road routes in QGIS, points where a road crosses several
order-3 streams within a few km — one example given: 4 order-3 crossings in 4 km, which in
reality means 4 separate bridge structures, a genuinely large construction undertaking the
model was pricing at nothing. Second point made: even in a genuinely high-stream-density
area, larger individual river widths/depths there should still raise the friction
calculation — implying the model wasn't registering the real difficulty of those corridors
at all.

**Checked directly against the actual (§11) road network, before changing anything**: walked
every road edge's real reconstructed route and counted contiguous stream-cell crossings (any
order). Result: **386 total crossings across the 20 road edges — only 31 (8%) were order>=4
and therefore priced by §11's width-based friction; 355 (92%) were order<4 and priced at
literally ZERO**, no matter how many of them a single edge crossed. Some edges cross 20-44
streams over their length (`Circulo_E4_2k↔Circulo_F4_small`: 44 crossings over 56 km) — this
wasn't a rare edge case, it was the norm for routes through this domain's basin-fill/wetland
terrain (`08_tappa8_geomorphology.md` §9's `wetland_backswamp` sub-class, 21.3% of
basin_fill). Nico's complaint was fully justified: the `MAJOR_STREAM_MIN_STRAHLER_ORDER=4`
gate, kept from §6c's original scoping ("at least those of higher Strahler order should be
considered" — Nico's own first-review framing, back when zero-cost rivers of ANY size was
the problem), had become the wrong cutoff once §11 made per-cell friction continuous rather
than a flat per-order bucket: a continuous function driven by real width naturally assigns a
small stream a small penalty on its own -- there was no longer a reason to floor it at
exactly zero.

**A second, independent finding surfaced by the same investigation**: some reaches labelled
order 1-3 in `streams.geojson` carry anomalously large `max_discharge_proxy_m3s` — up to
103 m3/s, matching this domain's largest order-6 rivers. Checked directly: 674 of 17,678
reaches (3.8%) have order<4 but discharge>10 m3/s; 230 have discharge>30 m3/s. These are
almost certainly short stub reaches at a large river's mouth (very low `n_cells`, e.g. 2-16,
against a `max_contributing_area_km2` of 700+ km2 — a huge basin draining through a tiny
final segment), most likely a Tappa 4 Strahler-order-export quirk at the coastline boundary,
**not investigated or fixed at that layer** (out of scope here — flagged for whoever next
works `04_tappa4_hydrology.md`). Under the OLD order>=4 gate this was a real correctness
bug: several of this domain's largest river mouths were being treated as completely free to
cross because of a mislabeled order attribute. Under §11's discharge/width-driven friction,
this is now harmless BY CONSTRUCTION even before fixing it — friction comes from
`max_discharge_proxy_m3s` directly, never from the order label, so a mislabeled-order,
correctly-large-discharge reach already gets the correctly-large penalty. Widening the gate
to include order<4 (below) is what actually exposes this fix, since those reaches were
excluded entirely before.

**Fix**: `ALL_STREAMS_MIN_ORDER=1` (new, in `run_tappa9_road_network.py`) replaces
`MAJOR_STREAM_MIN_STRAHLER_ORDER` as the operative scope for both the friction rasterization
and the per-reach geometry export — every stream in the domain, not just order>=4, now
rasterizes into `stream_width_grid` and gets a `river_friction_multiplier_from_width` value.
`MAJOR_STREAM_MIN_STRAHLER_ORDER` (4) is kept as a labelling constant only (used to tag
"minor" vs "major" in the new crossing-count report below), no longer a friction cutoff. No
change to the friction FUNCTION itself (§11's `river_friction_multiplier_from_width`,
anchors unchanged) — a tiny order-1 trickle (e.g. 0.003 m3/s) now gets an estimated width of
~0.02 m and friction ~0.998 (negligible individually, by the same curve that was already
correct for narrow order-4 reaches in §11), while a route crossing 40 of them accumulates
real, if modest per-crossing, additive cost — the ADDITIVE path-cost machinery
(`cost_distance.py`'s Dijkstra, unchanged) already handles the "many small crossings add up"
case correctly; the only thing missing was that minor crossings weren't contributing
anything to sum in the first place.

**New permanent diagnostic, not just a one-off check**: `run_tappa9_road_network.py` now
counts real stream crossings per road edge every run (`_count_crossings`, contiguous
stream-cell runs along the reconstructed path) and writes `stream_crossing_report` into
`tappa9_road_network_meta.json` — total crossings, how many are order<4, and a full
per-edge breakdown. This is the mechanism that produced the 386/355 numbers above, and will
catch the same class of problem automatically if it recurs after a future revision.

**Network re-run, checked directly**: identical 20-edge SET (15 MST + 5 redundant, same
site-pairs, verified against §10/§11's baseline — zero additions, zero removals) but
DIFFERENT route geometry for many edges — total crossings dropped from 386 to 282 (27%
fewer) as some routes now detour around the densest crossing clusters where a detour is
cheaper than paying the accumulated crossing cost, e.g.
`Circulo_E4_2k↔Circulo_F4_small` (44→40 crossings). Not every dense corridor could be routed
around — `Circulo_E4_2k↔Circulo_F4_small` and a few others still cross 25-40 streams,
because that terrain is genuinely, unavoidably stream-dense (matching Nico's own "however
high the water density is" framing) — the fix doesn't eliminate those crossings, it
correctly PRICES them, which is what was actually missing. Total road length 450.0 → 451.5
km (the small detours that WERE worth taking). 246 of the 282 remaining crossings are still
order<4 — no longer free, but individually modest, exactly as intended.

**Outputs**: `stream_geometry_full.geojson` (renamed from `major_stream_geometry.geojson`,
now covers all 17,678 reaches/all orders, not just the 1,847 order>=4 ones — the stale
order>=4-only file was deleted, not left alongside the new one, so QGIS projects should
re-point to the new name). `river_crossing_friction_multiplier_120m.*` and
`combined_friction_multiplier_120m.*` updated in place again (both changed, this domain's
river-affected land fraction grew from 1.22% to 17.99% of effective land now carrying some
crossing penalty). `road_network_mst.geojson` and `tappa9_road_network_meta.json`
regenerated with the new `stream_crossing_report` field.

**Known device-side leftover, not a pipeline bug**: after this fix shipped, Nico reported
that a reviewed Strahler-3 stretch "still looked the same." Checked directly rather than
assumed: the committed `road_network_mst.geojson`/`tappa9_road_network_meta.json` on his
machine matched the corrected local output byte-for-byte (same file sizes, mtimes from the
same commit batch) — the pipeline output was correct. The likely cause: the OLD
`major_stream_geometry.geojson` (pre-§12, order>=4-only) was never actually deleted on
Nico's machine — file delivery only writes files, it can't delete stale ones — so it sat
next to the new `stream_geometry_full.geojson` and could easily be the one opened by
mistake; separately, an already-open QGIS layer doesn't auto-refresh when the file
underneath it changes (needs a manual Reload). Confirmed separately, and NOT a bug: even
with the correct file open, some routes are SUPPOSED to keep crossing several minor streams
in a short span, because that is still the cheapest available path even after correctly
pricing every crossing (see `Circulo_E4_2k↔Circulo_F4_small`, 40 crossings, above) — Nico
confirmed this reading after the check.

## 13. Manual exception edge — Circulo_F1_small ↔ Circulo_F2_small

Nico noticed a missing direct road between two nearby Circulos, F1 and F2, on the same
QGIS review pass as §12. Checked directly, using the exact same numbers
`add_redundant_edges` (§8/§9) itself computes, rather than assuming the complaint was a bug:

F1 and F2 were never disconnected — both already reach each other via `Circulo_B_35k`
(1.8103h + 1.9008h = 3.7111h combined). A direct F1↔F2 edge costs 2.7469h — **26.0% cheaper
than that tree path**, comfortably past `add_redundant_edges`' own `min_shortcut_improvement`
bar (>=20%, and within the 21.0%-43.4% range of the 5 redundant edges §9 already shipped).
So by the "does this edge actually save real travel time" test, F1-F2 clearly qualifies.

It was never added because it fails the OTHER required gate, `redundancy_factor=1.4`-of-
each-site's-own-cheapest-neighbour, **on both ends**: F1's cheapest neighbour is
`Circulo_B_35k` at 1.8103h (gate: <=2.5345h; F1-F2's 2.7469h is 8.4% over) and F2's cheapest
neighbour is *also* `Circulo_B_35k` at 1.9008h (gate: <=2.6611h; F1-F2's 2.7469h is 3.2%
over). Both sites happen to share the same comparatively cheap hub, which pulls each site's
own "cheapest neighbour" baseline down and makes any third option — even a genuinely good
one — look proportionally too far by that specific gate's logic. This is a real edge case
of the `redundancy_factor` design (see §9/`add_redundant_edges`'s own docstring on why that
gate exists: to keep candidates scoped to a site's genuinely near neighbours, not arbitrary
long-distance chords), not a bug in it, and not something the earlier sensitivity analysis
(§10) would have caught either — it tested global parameter values, not specific
site-pair edge cases like this one.

**Decision (Nico's explicit call, offered three options — add as a one-off manual
exception, loosen `redundancy_factor` network-wide, or leave as-is):** add the single
F1↔F2 pair as a named MANUAL exception, not a threshold change. Loosening
`redundancy_factor` for the whole network was explicitly rejected — it would very likely
reopen the excessive-redundant-edges problem §9 already fixed once (a looser gate admits
other proportionally-distant candidates elsewhere in the network that have NOT been vetted
for genuine shortcut value the way this one specifically was, by hand, with real numbers).

**Implementation**: `MANUAL_EXTRA_EDGES` (new, `run_tappa9_road_network.py`, right after
`add_redundant_edges` runs) — a small named list of `(site_a, site_b)` pairs added directly
to `road_edges` as a distinct `edge_type: "manual_exception"` (not folded into
`"redundant"`, so it stays visually and programmatically distinguishable in the GeoJSON and
in any future edge-type-based analysis), reusing the already-computed `hours_road` matrix
for its weight (no new Dijkstra runs) and skipping automatically if the pair is already
present. `tappa9_road_network_meta.json` gets a new `manual_extra_edges` field recording
the full justification (direct cost, tree-path cost, shortcut %) per edge, and
`n_manual_extra_edges`/`n_road_edges_total` updated accordingly. Same ad-hoc-exception
discipline as `EXCLUDED_FROM_ROAD_NETWORK` (§9's `Circulo_D_20k` exclusion) — a single
named site (or, here, site-pair) decision, explicit and documented, not a rule change.

**Network re-run, checked directly**: 15 MST + 5 redundant + 1 manual = 21 road edges
(465.2 km, up from 451.5 km — the new edge's own 13.7 km route), 284 total stream crossings
(up from 282 — the new edge crosses 2 streams on its own way through). No other edge
changed. `Circulo_F1_small`/`Circulo_F2_small` now have a direct 13.7 km road
(`edge_weight_hours=2.7469`, `route_km=13.725`, vs. `straight_line_km=12.907` — a modest,
expected detour around terrain, not a straight line).

## 14. Connecting Tappa 10's 94 auxiliary settlements as spokes (2026-08-21, fifth follow-up)

A different chat (Tappa 10, opened 2026-08-21) sited 94 auxiliary settlements — 71 mining
posts, 2 mountain huts, 14 coastal villages, 7 forest posts — each administratively attached
to one or more of the 17 Círculos via an `attached_circulos` property, but never actually
connected to this backbone. Nico asked how to connect them, offering three rough options
himself (attach only to the named Círculo; search the closest/cheapest path; re-run the
whole MST+redundancy topology with these included) and explicitly asked for an opinion
before anything was built.

**Recommendation given and approved (Nico: "Implementar como descrito")**: treat every
settlement as a SPOKE off the already-locked 21-edge backbone, computed on the exact same
locked land-only combined-friction graph (`combined_friction_multiplier_120m.npy` loaded
directly, nothing recomputed) — explicitly NOT folding these 94 points into the
MST+redundancy algorithm as peer nodes. Two reasons: these settlements are leaves by
nature (nothing routes *through* a mining post to reach a Círculo), and admitting them as
peer nodes risked the solver finding it cheaper to route between two unrelated Círculos via
a side settlement, plus it would reopen a topology already reviewed three times (§5-§8) for
no benefit.

**Implementation** (`run_tappa10_network_connections.py`, new script, not folded into
`run_tappa9_road_network.py` since it operates on a different tier — settlement-to-backbone,
not backbone topology): mountain huts get no new edge at all — verified directly against the
existing edges' real reconstructed route cells that both huts sit exactly 0 cells from their
named crossing's path, confirming the Tappa 10 chat's own placement intent rather than
assuming it. Satellite mining posts (10, already consolidated into hub camps by Tappa 10's
own revision) spoke to their HUB, not a Círculo, matching that chat's hub/satellite design.
Settlements with 2+ `attached_circulos` that already share a direct Tappa 9 edge — checked
directly, 12 of 21 such multi-attached settlements qualified, not assumed from proximity —
get one T-junction spoke onto that edge's cheapest real-cost-distance tie-in point rather
than duplicate spokes (modeling a real side-road-meets-highway junction). Everything else
spokes to its single cheapest attached Círculo by real (symmetrized) cost-distance. Reused
16 Círculo-sourced single-source Dijkstra passes (dist-only, for symmetrized spoke weights)
plus one settlement-sourced Dijkstra pass per connectable settlement (dist+predecessors, for
that settlement's own route geometry) — ~92 runs total, ~32s, all off the graph and friction
field Tappa 9 already computed and saved.

**Result**: 75 of 92 connectable settlements got a spoke (858.9 km total,
`data/processed/transport/auxiliary_network_connections.geojson`, a separate layer from
`road_network_mst.geojson` since spokes are a different tier than the backbone) — breakdown
54 circulo_spoke, 12 t_junction_spoke, 9 satellite_to_hub_spoke. 17 settlements came back
genuinely isolated, checked (not assumed) against three distinct causes: 14 simply inherit
`Circulo_D_20k`'s own known land-isolation (its entire local cluster — 8 mining posts, 1
coastal village, 1 forest post — since `D_20k` itself is excluded from the road graph, §9);
2 match mining posts Tappa 10's own decision record already flagged as land-isolated,
independently confirmed here; and 1 is a NEW finding, not previously known:
`Circulo_F7_small`'s own coastal village (`Vila_Costeira_14`) sits on the same lake-mask
raster-resolution artifact already flagged twice for mining posts in the original Tappa 10
report (a 120m routing cell reads as touching `lake_mask` via aggregation of the native 30m
mask even though the underlying 30m pixel is dry land) — now hit a third time, on an
independent settlement type. See `10_tappa10_auxiliary_settlements.md`'s own "Network
connections" section for the full per-rule breakdown and numbers; this section is the
canonical implementation record since the work was done from this (Tappa 9) chat.

**Open (§14 specifically, superseded by the update below)**: whether the lake-mask artifact
is worth fixing at the source, and whether `Circulo_F7_small`'s coastal village should be
repositioned. See the update immediately below for the current numbers — the specific
village this concerned no longer isolates for the reason originally stated.

### §14 update, same day, sixth follow-up: re-run on the restructured v2 dataset

Nico reported "considerable" changes to the Tappa 10 supplementary files: some features
removed, others moved, and ~18 new mountain huts added, with an explicit question of
whether the huts still fit the spoke model or need a different, trail-based approach. File
is now `auxiliary_settlements_tappa10_v2.geojson`.

**Checked directly before doing anything** (byte diff + name-set diff against the original
file, not assumed): mining posts (71) and forest posts (7) are unchanged except renamed to
English (`Mina_X`→`Mine_X`, `Posto_Florestal`→`Forest_Post`); coastal villages went 14→12 (2
dropped, per Tappa 10's own direct edits); mountain huts went 2→18 — the real content
change, sourced from Tappa 7's `outpost_candidates.geojson` fauna-siting layer, not newly
authored. Schema also changed: list-valued properties (`attached_circulos` etc.) moved from
JSON arrays to semicolon-joined strings.

**A data bug found in v2**: all 10 `satellite_no_own_structure` mining posts' `attached_hub`
still names their hub under the pre-rename `Mina_*` scheme, which no longer matches any
`name` in the file. Confirmed 10/10, not coincidental. Fixed in-script via a prefix-
substitution lookup with a hard assertion every reference resolves, rather than silently
producing false isolation reports.

**Mountain huts: given an opinion before rebuilding, per the same discipline as the original
spoke design.** Checked directly: all 18 sit at 1,869-3,241 m (mean 2,680 m vs. this domain's
920 m land-wide mean), none within ~5 km of any existing road edge (the original 2 sat
exactly ON one), and a same-massif MST of the 17 `main_spine` huts is overwhelmingly path-
shaped (11/17 nodes degree 1-2, only 3 branch points) — a trekking spine, not a scatter.
Treating each hut as an independent spoke (the mining-post treatment) would have produced 17
near-parallel access routes up the same massif and ignored that shape. **Recommended and
approved**: a per-massif least-cost MST among the huts themselves (`main_spine`: 17 huts,
`south_branch`: 1 hut, same land-only friction graph, no separate foot-trail friction
invented), with a valley-access spoke at every LEAF of that tree only — interior trail nodes
reach the valley via the trail itself. Result: 16 trail edges connecting all 17 `main_spine`
huts (8 leaves → 8 trailhead spokes) + `south_branch`'s 1 hut as its own trivial trailhead (1
spoke) = 9 trailhead spokes total, 0 huts isolated. **Caveat surfaced, not resolved**: all 18
huts still carry `status_is_authorial_final: false` (inherited from Tappa 7) — this run
connected all of them regardless; if that status question resolves toward "not real," the
corresponding trail/trailhead edges need pruning, not just the site record.

**Result, full re-run**: 73 spokes (903.6 km, mining/coastal/forest, same rules as before) +
16 mountain-trail edges + 9 trailhead spokes (339.4 km combined) = 90 connectable non-hut
sites + 18 huts, 17 non-hut sites still isolated (same count as the v1 run, different
individual sites due to renumbering — 14 inherit `Circulo_D_20k`'s isolation, `Mine_Schist_22`
is a 4th confirmed instance of the lake-mask raster artifact, `Mine_BasinFill_16` is a
genuine land-isolated post matching the original report's own finding). **Correction to the
v1-run finding above**: `Circulo_F7_small`'s coastal village (renamed `Coastal_Village_11`)
does NOT sit on the lake-mask artifact as previously reported — re-checked directly against
the actual rasters this run, its cell and `Coastal_Village_03`'s (`Circulo_E4_2k`) both sit
120-240 m outside `land_mask.npy`'s land, not touching `lake_mask` at all. This is a
different, newly-identified bug: these two sites were placed against Tappa 10's own
corrected "true ocean" mask (border-connectivity-checked `lithology_v6==0`), which draws a
different, more precise coastline than the routing graph's older, coarser `land_mask.npy` —
a cross-dataset mismatch at the coastline, not a lake-adjacency artifact. Full numbers and
per-site detail in `10_tappa10_auxiliary_settlements.md`'s own "Network connections, v2
re-run" section.

**Open (§14, current)**: whether `land_mask.npy` should be reconciled against Tappa 10's
corrected true-ocean mask; whether the lake-mask artifact (now 4 confirmed instances) is
worth fixing at the source; the still-open hut narrative-status question (Tappa 7's own
item, not this chat's to resolve) and its implication for the new trail infrastructure if it
resolves toward some huts not being real.

## 15. Four topology/routing bugs Nico found in QGIS (2026-08-23, seventh follow-up)

Nico reviewed the rendered network in QGIS and reported four problems. Each was confirmed
directly against the code/data before any fix was designed — no guessing.

**Bug 1 — line endpoints don't touch their point features.** `_route_feature()` converted
EVERY path cell, including the first/last (site) cells, from grid row/col to XY via the
120 m cell-center formula (`x = XMIN + (c+0.5)*cs_x`, ...) — ignoring that the exact
authored site coordinate (`si["x_km"]`/`si["y_km"]`, already used elsewhere for
`straight_line_km`) was available. At 120 m resolution this put every Círculo dot up to
~85 m (`sqrt((cs_x/2)^2+(cs_y/2)^2)`) away from the line meant to touch it — visible in
QGIS, present on every one of the 21 backbone edges plus every Tappa 10 spoke/trail.
**Fixed** in both `run_tappa9_road_network.py` and `run_tappa10_network_connections.py`:
every interior vertex still comes from the routed path's cell centers (that's the real
discretized route), but the first and last vertex of any line whose endpoint is a known
site (Círculo or settlement) is overridden with that site's own exact coordinate. A
mid-route tie-in point (a line joining another line partway along its length, new in this
pass — see below) has no "truer" location than its own cell center and is correctly left
alone. Verified by exact float comparison against `tappa6_site_selection_meta.json`'s site
coordinates and, for the new Tappa 10 output, against all 95 settlement site coordinates
directly (0 mismatches). Re-running the backbone with only this change: 15 MST + 5
redundant + 1 manual exception = 21 edges, 465.3 km (vs. 465.2 km before) — confirms the
topology itself was untouched, only geometry.

**Bug 2 — the Southland (`Circulo_D_20k`) had no internal road network.** Not a wrong
answer in isolation — `Circulo_D_20k` genuinely has no land route into the backbone, by
Tappa 9's own design (§8) — but a real gap: the connection script only ever drew each
settlement's one pre-encoded satellite→hub relationship, and never asked whether isolated
settlements could reach *each other*. Root cause and fix are shared with bugs 3/4 below
(the same restrictive target search), which is why it's one unified rewrite rather than a
separate patch — see §"Network connections, v2 re-run" in
`10_tappa10_auxiliary_settlements.md` for the concrete before/after.

**Bugs 3/4 — redundant near-parallel paths (e.g. `Coastal_Village_09`) and roads that ignore
a nearby edge instead of joining it (e.g. `Forest_Post_03`).** Confirmed with real numbers
before proposing anything: `Forest_Post_03`'s spoke ran 170-424 m parallel to the
`Circulo_E3_2k <-> Circulo_F7_small` backbone edge for its *entire length* without ever
tying into it. Root cause, same for both: each settlement's connection-target search was
restricted to (a) shared backbone edges between its own `attached_circulos`, and (b) a
direct fallback to its attached Círculo — never the full backbone, and never another
settlement's already-built spoke/trail. **Fix, approved by Nico
("Corrigir tudo de uma vez")**: `run_tappa10_network_connections.py` was rewritten around a
greedy, incrementally-growing network (Prim's-style) — the Tappa 9 backbone seeds the
initial "network," then every unconnected settlement/hut-trail-leaf is connected in
globally-cheapest-first order to whichever is cheapest: the backbone at *any* point, or any
settlement's spoke/trail already built earlier in the same run. Each new connection's path
cells are folded into the network before the next iteration, so later settlements can feed
onto earlier ones. This single change fixes bug 2 as a side effect (an isolated
settlement's own spoke can now become a connection target for its neighbors) and directly
fixes 3/4 (a settlement now finds the *nearest point of the whole live network*, not just
its own hub's slice of it). Concrete results and the full new edge-type taxonomy
(`circulo_spoke`, `backbone_t_junction_spoke`, `feeder_t_junction_spoke`,
`isolated_pocket_road`) are documented in `10_tappa10_auxiliary_settlements.md`, since this
script owns the auxiliary-settlement connections; this backbone doc (Tappa 9's own script)
only needed the Bug 1 vertex fix.

**One bug found during verification of the fix, not reported by Nico**: the greedy loop
initially produced zero-length self-loop edges (`feeder_t_junction_spoke` from a settlement
to itself) whenever another settlement's spoke path happened to terminate exactly on its
site cell before it was itself processed (e.g. `Mine_Volcanic_02`, reached first by
`Mine_BasinFill_10`'s satellite spoke). Fixed by checking, before building an edge, whether
the chosen target cell *is* the settlement's own site — if so it's already connected at zero
cost and is skipped rather than given a spurious self-edge. 15 such cases across the full
re-run; `n_trivially_absorbed` in the meta JSON reports the count.

**A deeper version of the same class of bug, found by Nico on the SAME day (second QGIS
review of this pass's output), not by this chat's own verification**: mountain-trail edges
and `satellite_to_hub_spoke` edges were being added to the live network at construction
time, same as any real bridge to the backbone — but neither actually means "reaches a
Círculo." This let the entire 17-hut `main_spine` trail system, and at least one hub cluster
(`Mine_Volcanic_02`'s, the Southland's), register as "already connected" purely by seeing
their own already-known cells, with zero real edges ever built out. **This directly
invalidates a number I reported to Nico earlier the same day** — the "12 of 14 Southland
settlements now connected via a real road chain" claim was built on this bug; those 12 were
never actually reaching a Círculo, just each other. Fixed by holding satellite-spoke and
trail cells in separate per-cluster pools (`hub_cluster_cells`/`massif_cluster_cells`) and
only merging a cluster into the live network once one of its members earns a REAL outside
connection (`promote_cluster()`, new). Re-run: 0 of 2 massifs stayed unbridged (the trail
system does reach the backbone for real now), 1 of 9 hub clusters stayed unbridged (the
Southland's — correctly falls through to the isolated-pocket step this time, producing an
honest 12-member pocket with its own 11-edge/103.6 km internal MST instead of a false
backbone connection). Full numbers and the `Mine_Greywacke_05`/`attached_circulos` follow-up
question are in `10_tappa10_auxiliary_settlements.md`'s own fourth-pass addendum.

**Ninth follow-up, same day: a real labeling bug found while investigating Nico's "north part
disconnected" report — connectivity itself was fine, the metadata describing it wasn't.**
Nico reported the whole mountain system looked disconnected from the north (no link to
`Circulo_F1_small`/`Circulo_F2_small`) and asked for a concrete distance check on
`Mine_Greywacke_05` to its attached `Circulo_E3_2k`. Checked directly, not assumed:
`Mine_Greywacke_05 -> Circulo_E3_2k` is **93.81 km / 23.71 h** over the real network (vs.
**55.70 km / 14.13 h** to the non-attached `Circulo_E2_2k` it now actually routes near) —
confirms Nico's point with real numbers. The "north disconnected" report led to a genuine,
separate bug: 12 delivered edges (`Mine_Schist_06/07/09/12/13/18/19/20/21/25`,
`Coastal_Village_09`) had `connects_to` set to a generic internal string
(`"mountain trail (massif main_spine)"`, `"hub cluster: Mine_Schist_17"`) instead of the real
trail segment or satellite spoke they tied into. Root cause: `promote_cluster()` (the eighth
follow-up's own fix) called `add_to_network()` once per cluster with ONE bulk label for every
cell in it, discarding the per-edge label each cell originally carried. **This did not affect
connectivity** — cost and routed geometry were always computed against the real nearest cell
— only the human-readable description of what that cell was. Fixed by keeping
`massif_cluster_cells`/`hub_cluster_cells` as lists of `(path, label)` pairs instead of flat
cell lists, so `promote_cluster()` now replays each contributing edge's own real label.
Re-run: byte-identical connectivity (same 8 clusters promoted, same 1 hub cluster/0 massifs
unbridged, same `still_isolated` list) — confirms this was cosmetic, as expected.

**Separately, a real (not cosmetic) finding, checked with actual cost-distance, not
straight-line guesses**: none of the north schist mines are actually isolated — all are
reachable via the `main_spine` trail massif — but that massif's ONLY real bridge to the
backbone is the single spur at `Outpost_MainSpine_17` (5.18 km to a point on the
`Circulo_E3_2k<->Circulo_E4_2k` edge, 18.49 km from the `E3_2k` end / 21.94 km from the
`E4_2k` end, resolved from the edges' own geometry). Every settlement north of the massif has
to route through THAT one point, then most of the way around the backbone, to reach
`Circulo_F1_small`/`Circulo_F2_small` — real network distance today: `Mine_Schist_05 ->
Circulo_F1_small` 156.3 km, `Mine_Schist_09 -> Circulo_F1_small` 130.9 km, `Mine_Schist_04 ->
Circulo_F2_small` 172.4 km, `Mine_Schist_03 -> Circulo_F2_small` 165.4 km — against a 10-14 km
straight-line separation in every case. Re-ran real cost-distance (same friction graph, direct
source-to-target, no detour) for candidate direct connections: `Mine_Schist_05/06/09 ->
Circulo_F1_small` cost 15.28/14.54/17.44 km respectively (`Mine_Schist_06`, the actual hub
`Mine_Schist_05` and `Mine_Schist_08` both attach to, is the cheapest of the three);
`Mine_Schist_03/04 -> Circulo_F2_small` cost 11.48/11.80 km. A direct connection would cut
this detour by roughly 90%. Recommended to Nico, not yet built pending his choice of anchor
settlement — see `10_tappa10_auxiliary_settlements.md`'s own addendum for the full
recommendation and the `attached_circulos` design discussion this connects to.

Also checked Nico's specific request for a direct `Outpost_MainSpine_16 -> Circulo_E3_2k`
edge: current real route (via the `Outpost_MainSpine_17` spur) is 47.90 km; a direct
candidate would cost 26.36 km/9.54 h — a real ~45% improvement, smaller than the north-cluster
case but still substantial. Also confirmed `Coastal_Village_03` and `Mine_Schist_22` are NOT
new bugs — both were already in the `still_isolated` list from the eighth follow-up, with
documented reasons (genuinely solitary on the land-only graph; a hard satellite-to-hub rule
with no fallback, respectively).

**Implemented, same day, after Nico confirmed the anchors**: three new `manual_connection`
edges, same named-exception pattern as this doc's own §13 (`Circulo_F1_small<->Circulo_F2_small`).
Nico picked `Mine_Schist_06` (the cheapest of the three candidates, and the real hub
`Mine_Schist_05`/`Mine_Schist_08` both attach to) for `Circulo_F1_small`, `Mine_Schist_03`
for `Circulo_F2_small`, and confirmed `Outpost_MainSpine_16 -> Circulo_E3_2k` as requested.
`run_tappa10_network_connections.py` gains a `MANUAL_CONNECTIONS` list processed after the
greedy loop and isolated-pocket grouping, reusing each settlement's already-computed
dist/pred grid (`settlement_dist_pred` for the two mines, `hut_dist_pred` for the hut) —
real cost-distance, vertex-exact endpoints, same as every other edge type. These are
ADDITIVE, not replacements — the automated long-detour connections built earlier in the run
are untouched; a settlement can now have both. **Re-run result**: `Mine_Schist_06 ->
Circulo_F1_small` 14.59 km/4.93h, `Mine_Schist_03 -> Circulo_F2_small` 11.47 km/4.19h,
`Outpost_MainSpine_16 -> Circulo_E3_2k` 26.43 km/9.54h (all within 0.05 km of the standalone
diagnostic numbers above — the small delta is exact-vertex snapping). Every other number
unchanged (same 8 clusters promoted, same 5 still-isolated settlements). **New totals**: 16
mountain-trail edges (145.9 km) + 75 network spokes (401.8 km) + 11 isolated-pocket roads
(103.6 km) + 3 manual connections (52.5 km) = 703.8 km.

**Tenth follow-up, 2026-08-23 — `land_mask.npy`/true-ocean reconciliation.** Nico manually
nudged `Coastal_Village_03` ~67m west in his own GIS tool, believing it now sat on shore.
Verified directly (re-ran the production script against his edited geojson) rather than
trusting the coordinate alone: the village stayed in `still_isolated`, unchanged. Root cause:
both the old and new coordinates land in the identical 120m routing-grid cell (`xy_to_rc`
rounds to the same `(row, col)`) — a 67m move doesn't cross a cell boundary on a ~120m grid.
Computed the precise threshold (needed ~18.8m further west) and offered Nico two fixes: a
further precise manual nudge, or a proper reconciliation of `land_mask.npy` (1334x1084,
~120m, Tappa-1-era) against the newer, more carefully derived true-ocean mask (built from
`lithology_v6.npy` at native ~30m via border-connected-component labelling — same method
already used for the Tappa 10 third-pass coastal fix), which would also touch
`Coastal_Village_11`. Nico chose the reconciliation ("Por favor, faça isso"), explicitly
acknowledging `Coastal_Village_11` won't gain real connectivity (it's on an island) but
wanting its mask representation corrected too.

Built `build_land_mask_reconciled.py` → `data/processed/transport/land_mask_reconciled_v1.npy`,
a MONOTONIC/ADDITIVE-ONLY superset of `land_mask.npy` (only adds cells the true-ocean analysis
confirms are real land that `land_mask` misses; never removes a cell `land_mask` already calls
land) — chosen specifically to protect this doc's own already-verified "0 edges touching ocean
or lake" backbone property, which can only stay true or get MORE true under a purely additive
change. **First design attempt failed a direct check and was revised before implementing**: a
`>=0.5` majority-mean aggregation of the true-land signal (matching the general
false-negative-recovery diagnostic, which alone recovers 10,468 cells / 0.72% of the grid at
99.24% agreement with `land_mask`) does NOT actually flip either `Coastal_Village_03`'s cell
(`true_land_120_mean=0.3125`) or `Coastal_Village_11`'s (`0.1875`) to land — both sit, at
native 30m resolution, on real land, but right at the shoreline, which is exactly where a 120m
cell is most likely to be majority-sea even while genuinely touching land; a majority vote
systematically fails shoreline points by construction, not as an edge case. Switched to
`block_any` (does this 120m cell contain ANY real land pixel) — this flips both villages to
land but recovers 26,062 cells (1.80% of the grid, 2.5x the majority-rule footprint) instead of
10,468. This trade-off (fixes the actual reported cases vs. a larger area where the backbone
could theoretically find a cheaper route) was not put back to Nico before implementing, since
it's a bounded, directly-testable, purely-additive technical choice in service of the task he
already authorized — verified below rather than assumed.

Re-ran `run_tappa9_road_network.py` against the reconciled mask and diffed the result
cell-by-cell against the previously-locked 21-edge backbone (from this doc's earlier passes):
**topology unchanged** — same 21 edges (15 MST + 5 redundant + 1 manual exception), same
`from`/`to` pairs, same edge types, 0 added, 0 removed. Four edges got marginally cheaper
routes now that a few coastal slivers are traversable (`Circulo_B_35k<->Circulo_F3_small`
22.468→21.906 km, `Circulo_E3_2k<->Circulo_F7_small` 40.146→39.538 km,
`Circulo_E5_2k<->Circulo_F3_small` 23.263→23.081 km, `Circulo_F3_small<->Circulo_F4_small`
25.565→23.911 km), total backbone length 465.3→462.3 km (-3.0 km). Re-verified the "0 edges
touching ocean or lake" property directly against every vertex of every edge on the new mask
(`effective_land = land_mask_reconciled_v1 & ~lake120`): 0/21 edges touch a non-land cell,
confirming the additive change didn't silently break it.

**Eleventh follow-up, same day: two more concrete asks — `Coastal_Village_07`'s direct
connection to `Circulo_C_25k`, and `Circulo_D_20k`'s exclusion clarified/narrowed.** Nico
reported (a) `Coastal_Village_07` could have a direct connection to `Circulo_C_25k`, and (b)
"roads aren't connecting to `Circulo_D_20k`, the junction point is displaced from its real
position." **Checked (a) directly**: `Coastal_Village_07` currently reaches `Circulo_C_25k`
via a `feeder_t_junction_spoke` onto `Mine_Greywacke_01`'s own spoke, real network distance
9.49h; a direct connection by real cost-distance costs 4.54h/10.29 km — a genuine ~52% cut,
the same order of improvement as the ninth follow-up's connections. Built as a fourth
`manual_connection` edge (Tappa 10's own script, see "10" row). **On (b), asked Nico to
clarify before touching anything**, since as designed there was literally no line/geometry
touching `Circulo_D_20k` at all to be "displaced" — it's fully excluded from
`run_tappa9_road_network.py`'s processing (Nico's own explicit third-pass call: a genuinely
land-isolated island, no placeholder ferry line). **Nico's clarification**: `Circulo_D_20k`
should stay excluded from the INTER-Círculo backbone (still an island, no road/ferry to the
mainland — that part of the third-pass decision is unchanged, this doc's own backbone is
untouched by this follow-up), but it needs to be a valid connection target for its OWN
island's auxiliary settlements (the 14-member Southland group, previously only able to reach
each other via their own internal MST/pocket, never `Circulo_D_20k` itself). This is
Tappa 10's own script's scope, not this backbone script's — see "10" row's own addendum for
the fix and result.
