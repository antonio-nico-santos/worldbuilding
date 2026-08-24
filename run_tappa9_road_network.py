"""
Tappa 9 -- road network foundation: predecessor-path extraction on Tappa 6's
cost-distance graph (src/suitability/cost_distance.py), THREE friction
layers (lithology x biome x river-crossing) combined multiplicatively, and
a first-pass road-network TOPOLOGY (minimum spanning FOREST + local
redundancy) over the 17 already-placed, LOCKED Circulo sites.

**REVISED 2026-08-20, second pass, after Nico's direct review of the first
pass caught four real problems.** See docs/decisions/09_tappa9_transports.md
S6 for the full writeup; summary of what changed and why:

1. **Sea/lake crossings (6 of 16 edges ran across open ocean, ~8 of the
   remaining 10 crossed a lake).** Root cause: the first pass ran the road
   MST over `cost_distance.py`'s general-purpose graph, which -- by design,
   for Tappa 6's ORIGINAL isochrone/tier-distance use case -- lets any
   non-land cell (ocean or lake, this module has never told them apart) be
   crossed at a flat boat speed. That's correct for "how far away is this
   candidate site," wrong for "where can a road actually be built." Fixed
   by building the road network's own graph with `sea_mode="impassable"`
   (new in `cost_distance.py`) -- no non-land edge exists in this graph at
   all. See point 2 for the direct consequence of that choice.
2. **Rivers of any size, including this domain's largest trunk rivers,
   cost nothing to cross.** `cost_distance.py` never knew anything about
   hydrology. Fixed with a THIRD friction layer, `src/transport/
   river_friction.py`: Strahler order >= 4 reaches (from Tappa 4's own
   `data/exports/streams.geojson`) get a per-order friction multiplier,
   stacked multiplicatively with lithology and biome friction, same
   pattern as the other two.
3. **The network topology was always a single unbranching tree, even where
   several Circulos cluster closely.** An MST by definition has zero
   redundant loops. Confirmed via AskUserQuestion (2026-08-20) that the
   right fix is "MST + local redundancy": `src/transport/network.py`'s
   `add_redundant_edges` adds extra connections wherever a site has a
   second option nearly as cheap as its MST-assigned one -- see that
   module's docstring.
4. **(Not a code fix)** The GeoJSON's CRS not showing as the custom LCC in
   QGIS is this project's own known, pre-existing quirk (a `crs`/`proj4`
   FeatureCollection member that RFC 7946-compliant readers, QGIS included,
   silently ignore -- see `07_tappa7_regional_scenario.md`'s own
   documentation of the same behaviour on earlier Tappas' exports) -- not
   something this script introduced or can fix without a GDAL dependency
   this sandbox doesn't have (see `terrain/raster_io.py`'s own docstring on
   why GeoTIFF/rasterio isn't available here either). Workaround unchanged
   from every other export in this project: assign the CRS manually in
   QGIS after import (Layer Properties -> Source -> Assigned CRS, paste
   `config/parameters.yml`'s PROJ string).

Point 1's consequence, handled explicitly rather than silently: cutting
sea/lake edges can (and does, on this domain) split the 17 sites into more
than one connected component -- a genuine "these sites need a FERRY to
reach the rest of the network, not a road" finding, not a bug. Any such
split is surfaced in the output metadata and, where geometry allows, drawn
as a separate `edge_type: "candidate_ferry_crossing"` feature in the
GeoJSON (using the ORIGINAL boat-enabled graph's route -- informational
only, explicitly NOT part of the road network itself, a preview of what
Tappa 9's still-not-built ferries sub-build will need to formalize).

Per `07_tappa7_regional_scenario.md` S9's roadmap, Tappa 9 = Transports:
roads/rail, kite buggies, ferries, and the navigability half of dangerous
seas. This script covers ONLY the road-network foundation slice (topology +
route geometry + the friction layers under it); rail's grade-ceiling cost
function, actual ferry corridor authoring, and the kite wind-shadow mask
are separate, not-yet-built pieces of the same Tappa -- see this run's own
decision doc (docs/decisions/09_tappa9_transports.md) for the full scope
map.

THE 17 CIRCULO SITES ARE READ-ONLY, LOCKED INPUT -- never re-placed here,
same discipline S8f/S8g already established for the lithology friction and
excavation-effort layers.

Reads:
  data/processed/climate/land_mask.npy
  data/processed/hydrology/lake_mask.npy               (native 30m)
  data/processed/dem_v3_final_30m_eroded.npy            (native 30m)
  data/processed/geomorphology/lithology_v6.npy         (native 30m)
  data/processed/suitability/biome_id_smoothed_120m.npy (native 120m)
  data/exports/streams.geojson                          (Strahler-ordered)
  data/processed/suitability/circulo_candidate_sites.geojson
  data/processed/suitability/tappa6_site_selection_meta.json
  config/parameters.yml

Writes to data/processed/transport/ (gitignored, regenerate locally):
  road_network_mst.geojson         (road edges -- MST + redundant, WITH
                                    real land-only route geometry -- plus
                                    any candidate_ferry_crossing features)
  biome_friction_multiplier_120m.*
  river_crossing_friction_multiplier_120m.*  (WIDTH-BASED as of this revision,
                                    see transport/river_geometry.py)
  combined_friction_multiplier_120m.*  (lithology x biome x river, what the
                                    road graph was actually built from)
  stream_geometry_full.geojson     (per-reach estimated width_m/depth_m for the
                                    FULL stream network, ALL orders, informational)
  tappa9_road_network_meta.json
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from geomorphology.lithology import LAND_TRAVEL_FRICTION, LITHOLOGY_CLASSES, travel_friction_multiplier
from biomes.world_biomes import BIOME_NAMES
from suitability.cost_distance import build_cost_graph, cost_distance_from_source
from suitability.terrain_metrics import block_any, block_mean, block_mode
from transport.biome_friction import BIOME_TRAVEL_FRICTION, biome_friction_multiplier
from transport.river_friction import (
    MAJOR_STREAM_MIN_STRAHLER_ORDER,
    RIVER_CROSSING_FRICTION,
    RIVER_FRICTION_WIDTH_ANCHORS,
    rasterize_major_streams,
    river_friction_multiplier,
    river_friction_multiplier_from_width,
)
from transport.river_geometry import (
    ANCHOR_DEPTH_M,
    ANCHOR_DISCHARGE_M3S,
    ANCHOR_WIDTH_M,
    DEPTH_EXPONENT,
    WIDTH_EXPONENT,
    annotate_reach_geometry,
    estimate_depth_m,
    estimate_width_m,
    rasterize_major_stream_discharge,
)
from transport.network import (
    add_redundant_edges,
    build_mst_forest,
    compute_pairwise_cost_distance,
    connected_components,
    edge_path_cells,
    _tree_path_cost,
)
from terrain.raster_io import write_envi_raw, write_prj

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 9 -- road network foundation (2nd pass: land-only + river friction + MST-forest+redundancy) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
XMIN, XMAX, YMIN, YMAX = domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"]

OUT = "data/processed/transport"
os.makedirs(OUT, exist_ok=True)

# --- reproduce Tappa 6/8's exact 120m grid -------------------------------
# **REVISED 2026-08-23, land_mask/true-ocean reconciliation (Nico's request,
# following the Coastal_Village_03/11 coastline investigation).** Loads
# land_mask_reconciled_v1.npy (build_land_mask_reconciled.py) instead of the
# raw land_mask.npy -- a MONOTONIC/ADDITIVE-ONLY superset (26,062 cells /
# 1.80% of the grid added, 0 removed) built from the newer lithology_v6-
# derived true-ocean mask. See that script's docstring for the full method
# and for why block_any (not a >=0.5 majority) was used. Because this is
# purely additive over the mask this backbone was originally verified
# against, the "0 edges touching ocean or lake" property can only stay true
# or get MORE true -- see docs/decisions/09_tappa9_transports.md for the
# direct re-run diff against the previously-locked 21-edge topology.
log("loading inputs (reproducing Tappa 6/8's grid exactly)...")
land = np.load("data/processed/transport/land_mask_reconciled_v1.npy").astype(bool)
lake30 = np.load("data/processed/hydrology/lake_mask.npy")
dem30 = np.load("data/processed/dem_v3_final_30m_eroded.npy")
lithology30 = np.load("data/processed/geomorphology/lithology_v6.npy")
biome_id = np.load("data/processed/suitability/biome_id_smoothed_120m.npy")

ny, nx = land.shape
cs_x = (XMAX - XMIN) / nx
cs_y = (YMAX - YMIN) / ny
cellsize_km = (cs_x + cs_y) / 2 / 1000.0
log(f"  120m grid: {ny}x{nx}, cellsize ({cs_x:.3f}, {cs_y:.3f}) m")

lake120 = block_any(lake30, 4)[:ny, :nx]
effective_land = land & ~lake120
dem120 = block_mean(dem30, 4)[:ny, :nx]

n_classes = len(LITHOLOGY_CLASSES)
lithology120 = block_mode(lithology30, 4, n_classes)[:ny, :nx]

assert biome_id.shape == (ny, nx), f"biome_id shape {biome_id.shape} != grid {(ny, nx)}"

# --- river-crossing friction, WIDTH-BASED, ALL ORDERS (revised AGAIN same day
# -- Nico noticed several road routes cross many order<4 streams in a short
# span (one edge crosses 44 streams over 56 km) that the order>=4-gated model
# priced at literally zero, and separately that some order<1-3-labelled
# reaches in this domain carry anomalously large discharge (a Tappa 4 export
# quirk, not fixed there -- see below). ALL_STREAMS_MIN_ORDER=1 removes the
# order gate entirely; MAJOR_STREAM_MIN_STRAHLER_ORDER (4) is kept only as a
# labelling threshold for the crossing-count breakdown below, no longer as a
# friction cutoff -- see transport/river_geometry.py + river_friction.py's
# docstrings for the full reasoning) -----------------------------------------
ALL_STREAMS_MIN_ORDER = 1
log(f"rasterizing ALL streams (Strahler order >= {ALL_STREAMS_MIN_ORDER}, no longer gated at "
    f"{MAJOR_STREAM_MIN_STRAHLER_ORDER}) onto the 120m grid (order + discharge, for width-based friction)...")
major_stream_mask, stream_order_grid = rasterize_major_streams(
    "data/exports/streams.geojson", XMIN, YMAX, cs_x, cs_y, (ny, nx),
    min_order=ALL_STREAMS_MIN_ORDER,
)
stream_discharge_grid = rasterize_major_stream_discharge(
    "data/exports/streams.geojson", XMIN, YMAX, cs_x, cs_y, (ny, nx),
    min_order=ALL_STREAMS_MIN_ORDER,
)
stream_width_grid = estimate_width_m(stream_discharge_grid)
stream_depth_grid = estimate_depth_m(stream_discharge_grid)
river_friction = river_friction_multiplier_from_width(stream_width_grid)
log(f"  {int(major_stream_mask.sum())} cells ({100 * major_stream_mask[effective_land].sum() / effective_land.sum():.2f}% "
    f"of effective land) now carry SOME river-crossing penalty (was {MAJOR_STREAM_MIN_STRAHLER_ORDER}+-only before)")
_qualifying = stream_width_grid[stream_width_grid > 0]
log(f"  estimated width across qualifying cells: min {_qualifying.min():.3f}m, "
    f"mean {_qualifying.mean():.2f}m, max {_qualifying.max():.2f}m "
    f"(Leopold-Maddock W=Q^{WIDTH_EXPONENT}, anchored {ANCHOR_WIDTH_M}m @ {ANCHOR_DISCHARGE_M3S}m3/s)")
_low_order_big = int(((stream_order_grid > 0) & (stream_order_grid < MAJOR_STREAM_MIN_STRAHLER_ORDER) & (stream_width_grid > ANCHOR_WIDTH_M * 0.3)).sum())
log(f"  {_low_order_big} cells labelled order<{MAJOR_STREAM_MIN_STRAHLER_ORDER} but with estimated width "
    f">{ANCHOR_WIDTH_M * 0.3:.0f}m -- a real Tappa 4 streams.geojson export quirk (some low-order-labelled "
    f"reaches carry large discharge, not fixed at that layer), now correctly penalized by discharge/width "
    f"regardless of the (mislabelled) order.")

log("exporting per-reach width/depth GeoJSON for the FULL stream network (informational -- QGIS "
    "symbology + narrative reference; Nico asked for the complete network, not just order>=4, after "
    "noticing routes crossing many minor streams) -- does not feed anything else in this script...")
reach_geometry = annotate_reach_geometry("data/exports/streams.geojson", ALL_STREAMS_MIN_ORDER)
with open(f"{OUT}/stream_geometry_full.geojson", "w") as f:
    json.dump(reach_geometry, f, indent=2)
log(f"  {len(reach_geometry['features'])} reaches (ALL orders) written to stream_geometry_full.geojson")

# --- combined friction: lithology (Tappa 8) x biome x river (new) ---------
log("building combined lithology x biome x river-crossing friction field...")
lith_friction = travel_friction_multiplier(lithology120)
biome_friction = biome_friction_multiplier(biome_id)
combined_friction = (lith_friction * biome_friction * river_friction).astype(np.float32)
log(f"  land-mean lithology friction: {lith_friction[effective_land].mean():.3f}")
log(f"  land-mean biome friction:     {biome_friction[effective_land].mean():.3f}")
log(f"  land-mean river friction:     {river_friction[effective_land].mean():.3f}")
log(f"  land-mean combined friction:  {combined_friction[effective_land].mean():.3f}")

# --- cost graphs: BASELINE/BOAT-ENABLED (for comparison + candidate-ferry --
# geometry) and the new LAND-ONLY ROAD graph (sea_mode="impassable") -------
log("building BASELINE cost graph (no friction, boat-enabled -- Tappa 6 reproduction)...")
t0 = time.time()
graph_baseline = build_cost_graph(dem120, effective_land, cellsize_km)
log(f"  built in {time.time() - t0:.1f}s")

log("building BOAT-ENABLED combined-friction graph (reference only, for candidate-ferry "
    "geometry -- NOT what the road network is built from)...")
t0 = time.time()
graph_boat_combined = build_cost_graph(
    dem120, effective_land, cellsize_km, friction_multiplier=combined_friction
)
log(f"  built in {time.time() - t0:.1f}s")

log("building LAND-ONLY (sea_mode='impassable') combined-friction graph -- THIS is what the "
    "road network is actually built from...")
t0 = time.time()
graph_road = build_cost_graph(
    dem120, effective_land, cellsize_km, friction_multiplier=combined_friction, sea_mode="impassable"
)
graph_build_s = time.time() - t0
log(f"  built in {graph_build_s:.1f}s")

# --- load the 17 LOCKED Circulo sites (read-only) -------------------------
log("loading the 17 already-placed Circulo sites (read-only)...")
with open("data/processed/suitability/tappa6_site_selection_meta.json") as f:
    tappa6_meta = json.load(f)

# EXCLUDED_FROM_ROAD_NETWORK -- Nico's direct call, second review (2026-08-20):
# the first (second-pass) run confirmed Circulo_D_20k is fully land-isolated
# from the other 16 sites (see docs/decisions/09_tappa9_transports.md S6a) and
# drew a "candidate_ferry_crossing" line across open water to represent that --
# Nico's read: since it's genuinely an island, LOCK it out of this run's
# processing entirely rather than draw anything for it here at all. A real
# ferry connection is Tappa 9's own not-yet-built ferries sub-build's job, not
# something this road-network script should sketch a placeholder for.
EXCLUDED_FROM_ROAD_NETWORK = {"Circulo_D_20k"}

sites = []
excluded_sites = []
for s in tappa6_meta["sites"]:
    if not s.get("placed"):
        continue
    col = int(round((s["x_km"] - XMIN / 1000.0) / cellsize_km - 0.5))
    row = int(round((YMAX / 1000.0 - s["y_km"]) / cellsize_km - 0.5))
    site = {
        "name": s["name"], "tier": s["tier"], "population": s["population"],
        "row": row, "col": col, "x_km": s["x_km"], "y_km": s["y_km"],
    }
    if s["name"] in EXCLUDED_FROM_ROAD_NETWORK:
        excluded_sites.append(site)
        continue
    sites.append(site)
log(f"  {len(sites)} placed sites loaded ({len(excluded_sites)} excluded from road-network "
    f"processing: {', '.join(s['name'] for s in excluded_sites)} -- island, land-isolated, "
    f"deferred to the ferries sub-build, see EXCLUDED_FROM_ROAD_NETWORK above)")

# --- all-pairs cost-distance + predecessors on the LAND-ONLY road graph ---
log(f"computing all-pairs cost-distance + predecessors ({len(sites)} Dijkstra runs, "
    f"LAND-ONLY road graph)...")
t0 = time.time()
hours_road, predecessors_road = compute_pairwise_cost_distance(sites, graph_road, (ny, nx))
log(f"  done in {time.time() - t0:.1f}s")

log("computing all-pairs cost-distance + predecessors on the BOAT-ENABLED combined-friction "
    "graph too (reference for candidate-ferry crossings + baseline comparison)...")
t0 = time.time()
hours_boat_combined, predecessors_boat = compute_pairwise_cost_distance(sites, graph_boat_combined, (ny, nx))
log(f"  done in {time.time() - t0:.1f}s")

log("computing all-pairs cost-distance on the plain BASELINE graph too, for comparison only "
    "(no predecessors needed -- costs only)...")
t0 = time.time()
hours_baseline = np.full((len(sites), len(sites)), np.inf)
for i, s in enumerate(sites):
    d = cost_distance_from_source(graph_baseline, s["row"], s["col"], (ny, nx))
    for j, t in enumerate(sites):
        hours_baseline[i, j] = d[t["row"], t["col"]]
log(f"  done in {time.time() - t0:.1f}s")

# --- connected components on the LAND-ONLY graph (point 1's consequence) --
components = connected_components(hours_road)
log(f"  land-only connectivity: {len(components)} connected component(s) among {len(sites)} sites "
    f"({'fully connected by land' if len(components) == 1 else 'NOT fully connected by land -- see candidate-ferry crossings below'})")
for comp in components:
    names = [sites[i]["name"] for i in comp]
    log(f"    component ({len(comp)}): {', '.join(names)}")

# --- MST FOREST + local redundancy -----------------------------------------
log("building minimum spanning FOREST over symmetrized land-only combined-friction cost-distance...")
mst_edges, _ = build_mst_forest(hours_road)
log(f"  {len(mst_edges)} MST edge(s) across {len(components)} component(s) "
    f"(expected {len(sites) - len(components)} for a full spanning forest)")

log("adding local-redundancy edges (redundancy_factor=1.4 AND min_shortcut_improvement=0.20, "
    "both required -- revised after Nico's 2nd review flagged 3 of the previous heuristic's 8 "
    "edges as excessive; see src/transport/network.py's add_redundant_edges docstring for the "
    "two intermediate attempts that were tested and rejected first)...")
redundant_edges = add_redundant_edges(
    hours_road, mst_edges, components, redundancy_factor=1.4, min_shortcut_improvement=0.20
)
log(f"  {len(redundant_edges)} redundant edge(s) added")

# --- MANUAL exception edge: Circulo_F1_small <-> Circulo_F2_small --------
# Nico noticed this connection missing on direct QGIS review of this same
# run (2026-08-20). Checked directly rather than assumed, using the exact
# numbers add_redundant_edges itself uses: F1-F2 is NOT unreachable (both
# already connect via Circulo_B_35k, 3.7111h combined) and a direct edge
# clears add_redundant_edges' own tree-path-shortcut bar with room to
# spare (2.7469h direct vs 3.7111h via the tree = 26.0% cheaper, well
# inside the 21.0%-43.4% range of the 5 redundant edges already shipped).
# It fails the OTHER required check, redundancy_factor=1.4-of-cheapest-
# neighbour, on BOTH ends: F1's cheapest neighbour is Circulo_B_35k
# (1.8103h; F1-F2's 2.7469h is 8.4% over the 1.4x=2.5345h gate) and F2's
# cheapest neighbour is ALSO Circulo_B_35k (1.9008h; 3.2% over its own
# 1.4x=2.6611h gate). Both sites route through the same comparatively
# cheap hub, which pulls each site's own "cheapest neighbour" baseline
# down and makes any third option -- even a good one -- look
# proportionally too far by that gate's specific logic, a real edge case
# of the redundancy_factor design, not a bug in it.
#
# Nico's explicit call: add this ONE pair as a manual exception rather
# than loosen redundancy_factor network-wide, which would very likely
# reopen the excessive-edges problem from the 2nd-pass review (see
# fixes_from_second_pass.excessive_redundant_edges above) -- a looser gate
# would also admit other proportionally-distant candidates elsewhere in
# the network that have NOT been checked for genuine shortcut value the
# way this one specifically was. Not a config/threshold change: a single
# named-pair addition, same ad-hoc-exception discipline as
# EXCLUDED_FROM_ROAD_NETWORK above, just in the opposite direction (adding
# a connection instead of removing a site).
MANUAL_EXTRA_EDGES = [("Circulo_F1_small", "Circulo_F2_small")]
name_to_idx = {s["name"]: k for k, s in enumerate(sites)}
sym_hours_for_manual = 0.5 * (hours_road + hours_road.T)
existing_edge_pairs = {(min(i, j), max(i, j)) for i, j, _ in mst_edges + redundant_edges}
manual_edges = []
manual_edges_meta = []
for name_a, name_b in MANUAL_EXTRA_EDGES:
    ia, ib = name_to_idx[name_a], name_to_idx[name_b]
    pair = (min(ia, ib), max(ia, ib))
    tree_cost_manual = _tree_path_cost(mst_edges, ia, ib)
    direct_cost_manual = float(sym_hours_for_manual[ia, ib])
    if pair in existing_edge_pairs:
        log(f"  MANUAL_EXTRA_EDGES: {name_a}<->{name_b} already present, skipping")
        continue
    manual_edges.append((pair[0], pair[1], direct_cost_manual))
    manual_edges_meta.append({
        "from": name_a, "to": name_b,
        "direct_cost_hours": round(direct_cost_manual, 4),
        "tree_path_cost_hours": round(float(tree_cost_manual), 4) if tree_cost_manual else None,
        "tree_path_shortcut_improvement_pct": round(
            100 * (tree_cost_manual - direct_cost_manual) / tree_cost_manual, 1
        ) if tree_cost_manual else None,
        "reason": "Nico noticed this connection missing on direct QGIS review. Clears "
        "add_redundant_edges' own >=20% tree-path-shortcut bar (26.0%) but fails "
        "redundancy_factor=1.4 on both endpoints (8.4% / 3.2% over each site's own "
        "cheapest-neighbour gate) because both sites route cheaply through the same hub "
        "(Circulo_B_35k). Added as a manual exception rather than loosening "
        "redundancy_factor network-wide -- see this script's own inline comment at "
        "MANUAL_EXTRA_EDGES for the full reasoning.",
    })
    log(f"  MANUAL_EXTRA_EDGES: added {name_a}<->{name_b} ({direct_cost_manual:.4f}h, "
        f"{manual_edges_meta[-1]['tree_path_shortcut_improvement_pct']}% shortcut over the "
        f"tree path) -- Nico's explicit call")

road_edges = (
    [(i, j, w, "mst") for i, j, w in mst_edges]
    + [(i, j, w, "redundant") for i, j, w in redundant_edges]
    + [(i, j, w, "manual_exception") for i, j, w in manual_edges]
)

log("reconstructing real land-only route geometry for every road edge (reusing the already-"
    "computed predecessor arrays -- no new Dijkstra runs)...")
road_paths = edge_path_cells([(i, j, w) for i, j, w, _ in road_edges], sites, predecessors_road, (ny, nx))

# --- crossing-count check (Nico's own diagnostic, made permanent): how many -
# DISTINCT stream crossings does each road edge's real route make, and how --
# many of those were previously priced at zero (order < MAJOR_STREAM_MIN_..)?
log("counting distinct stream crossings per road edge (checked directly against the real route, "
    "not estimated) -- a contiguous run of stream cells along the path counts as ONE crossing...")


def _count_crossings(path, order_grid):
    crossings = []
    run_orders = []
    in_stream = False
    for r, c in path:
        o = int(order_grid[r, c])
        if o > 0 and not in_stream:
            run_orders = [o]
            in_stream = True
        elif o > 0:
            run_orders.append(o)
        elif in_stream:
            crossings.append(max(run_orders))
            in_stream = False
    if in_stream:
        crossings.append(max(run_orders))
    return crossings


total_crossings = 0
total_minor_crossings = 0  # order < MAJOR_STREAM_MIN_STRAHLER_ORDER -- was zero-cost before this revision
edge_crossing_report = []
for (i, j, w_hours, etype), path in zip(road_edges, road_paths):
    orders = _count_crossings(path, stream_order_grid)
    n_minor = sum(1 for o in orders if o < MAJOR_STREAM_MIN_STRAHLER_ORDER)
    total_crossings += len(orders)
    total_minor_crossings += n_minor
    edge_crossing_report.append({
        "from": sites[i]["name"], "to": sites[j]["name"], "edge_type": etype,
        "n_crossings": len(orders), "n_minor_crossings_order_lt_4": n_minor,
        "crossing_orders": orders,
    })
log(f"  {total_crossings} total stream crossings across all {len(road_edges)} road edges "
    f"({total_minor_crossings} of them order<{MAJOR_STREAM_MIN_STRAHLER_ORDER}, priced at ZERO before "
    f"this revision -- now priced via width, same as every other crossing).")

# --- candidate ferry crossings: cheapest inter-component connections, -----
# using the BOAT-ENABLED graph's real route geometry (informational only) --
ferry_edges = []
ferry_paths = []
if len(components) > 1:
    log(f"{len(components)} land-disconnected component(s) found -- computing candidate ferry "
        f"crossings (boat-enabled graph, informational only, NOT part of the road network)...")
    comp_of_site = {}
    for ci, comp in enumerate(components):
        for idx in comp:
            comp_of_site[idx] = ci
    sym_boat = 0.5 * (hours_boat_combined + hours_boat_combined.T)
    # cheapest single link between every pair of components -> MST over the
    # component-level graph, so every component ends up linked into one
    # candidate-ferry-augmented network, same "minimum total connection
    # cost" logic build_mst_forest already uses at the site level
    n_comp = len(components)
    comp_best = {}  # (ca, cb) -> (cost, i, j)
    for i in range(len(sites)):
        for j in range(len(sites)):
            if i == j:
                continue
            ci, cj = comp_of_site[i], comp_of_site[j]
            if ci == cj:
                continue
            key = (min(ci, cj), max(ci, cj))
            cost = sym_boat[i, j]
            if key not in comp_best or cost < comp_best[key][0]:
                comp_best[key] = (cost, i, j)
    # Prim's over the n_comp component-supernodes using comp_best as edge weights
    in_tree = {0}
    remaining = set(range(1, n_comp))
    while remaining:
        best = None
        for ca in in_tree:
            for cb in remaining:
                key = (min(ca, cb), max(ca, cb))
                if key not in comp_best:
                    continue
                cost, i, j = comp_best[key]
                if best is None or cost < best[0]:
                    best = (cost, i, j, cb)
        cost, i, j, cb = best
        ferry_edges.append((i, j, float(cost)))
        in_tree.add(cb)
        remaining.discard(cb)
    ferry_paths = edge_path_cells(ferry_edges, sites, predecessors_boat, (ny, nx))
    log(f"  {len(ferry_edges)} candidate ferry crossing(s) linking all {n_comp} components")
else:
    log("  all sites in one land-connected component -- no candidate ferry crossings needed")

# --- assemble GeoJSON -------------------------------------------------------


def _route_feature(i, j, w_hours, path, edge_type, hours_matrix):
    si, sj = sites[i], sites[j]
    # Nico found (2026-08-21, network-connections review): every road edge's
    # LineString endpoint snapped to its 120m CELL CENTER instead of the
    # site's own exact authored coordinate, so the point layer (Circulo
    # dots) and the line layer visibly didn't touch in QGIS -- up to ~85m
    # apart at this grid resolution. Fixed here: every INTERIOR vertex
    # still comes from the routed path's cell centers (that's the real
    # route geometry, correctly discretized), but the first and last
    # vertices are overridden with the two sites' own exact x_km/y_km, so
    # the line's endpoints exactly match the point features they connect.
    coords = []
    for r, c in path:
        x = XMIN + (c + 0.5) * cs_x
        y = YMAX - (r + 0.5) * cs_y
        coords.append([x, y])
    if coords:
        coords[0] = [si["x_km"] * 1000.0, si["y_km"] * 1000.0]
        coords[-1] = [sj["x_km"] * 1000.0, sj["y_km"] * 1000.0]
    seg_km = sum(
        float(np.hypot(coords[k + 1][0] - coords[k][0], coords[k + 1][1] - coords[k][1]))
        for k in range(len(coords) - 1)
    ) / 1000.0
    tier_pair = sorted({t for t in (si["tier"], sj["tier"]) if t})
    return seg_km, {
        "type": "Feature",
        "properties": {
            "edge_type": edge_type,
            "from": si["name"], "to": sj["name"],
            "tier_pair": tier_pair,
            "hours_i_to_j": round(float(hours_matrix[i, j]), 4),
            "hours_j_to_i": round(float(hours_matrix[j, i]), 4),
            "edge_weight_hours": round(float(w_hours), 4),
            "route_km": round(seg_km, 3),
            "n_cells": len(path),
            "straight_line_km": round(
                float(np.hypot(si["x_km"] - sj["x_km"], si["y_km"] - sj["y_km"])), 3
            ),
        },
        "geometry": {"type": "LineString", "coordinates": coords},
    }


features = []
total_road_km = 0.0
for (i, j, w_hours, etype), path in zip(road_edges, road_paths):
    seg_km, feat = _route_feature(i, j, w_hours, path, etype, hours_road)
    total_road_km += seg_km
    features.append(feat)

total_ferry_km = 0.0
for (i, j, w_hours), path in zip(ferry_edges, ferry_paths):
    seg_km, feat = _route_feature(i, j, w_hours, path, "candidate_ferry_crossing", hours_boat_combined)
    total_ferry_km += seg_km
    features.append(feat)

geojson = {
    "type": "FeatureCollection",
    "name": "road_network_mst",
    "crs": {"type": "proj4", "properties": {"proj4": CRS_PROJ4}},
    "features": features,
}
with open(f"{OUT}/road_network_mst.geojson", "w") as f:
    json.dump(geojson, f, indent=2)

# --- export friction rasters (float32, matching S8f's convention) ---------
write_envi_raw(
    f"{OUT}/biome_friction_multiplier_120m", biome_friction, XMIN, YMIN, (cs_x + cs_y) / 2,
    "Tappa 9 biome travel-friction multiplier (BIOME_TRAVEL_FRICTION, "
    "src/transport/biome_friction.py) -- UNREVIEWED first-pass estimates, not locked.",
    dtype="f4",
)
np.save(f"{OUT}/biome_friction_multiplier_120m.npy", biome_friction)
write_prj(f"{OUT}/biome_friction_multiplier_120m.prj", CRS_PROJ4)

write_envi_raw(
    f"{OUT}/river_crossing_friction_multiplier_120m", river_friction, XMIN, YMIN, (cs_x + cs_y) / 2,
    "Tappa 9 river-crossing travel-friction multiplier -- WIDTH-BASED, ALL Strahler orders "
    "(river_friction_multiplier_from_width, src/transport/river_friction.py, driven by "
    f"estimated channel width from src/transport/river_geometry.py), Strahler order >= "
    f"{ALL_STREAMS_MIN_ORDER} -- UNREVIEWED first-pass estimates, not locked.",
    dtype="f4",
)
np.save(f"{OUT}/river_crossing_friction_multiplier_120m.npy", river_friction)
write_prj(f"{OUT}/river_crossing_friction_multiplier_120m.prj", CRS_PROJ4)

write_envi_raw(
    f"{OUT}/combined_friction_multiplier_120m", combined_friction, XMIN, YMIN, (cs_x + cs_y) / 2,
    "Tappa 9 combined (lithology x biome x river-crossing) travel-friction multiplier -- the "
    "field the LAND-ONLY road-network graph was actually built from.",
    dtype="f4",
)
np.save(f"{OUT}/combined_friction_multiplier_120m.npy", combined_friction)
write_prj(f"{OUT}/combined_friction_multiplier_120m.prj", CRS_PROJ4)

# --- baseline-vs-combined comparison across the road edges (transparency, -
# same spirit as S8f's baseline-vs-friction comparison) --------------------
edge_comparisons = []
worst_delta_pct = None
for i, j, w_hours, etype in road_edges:
    b = float(hours_baseline[i, j])
    c = float(hours_road[i, j])
    delta_pct = round(100 * (c - b) / b, 2) if b > 0 and np.isfinite(c) else None
    if delta_pct is not None and (worst_delta_pct is None or delta_pct > worst_delta_pct):
        worst_delta_pct = delta_pct
    edge_comparisons.append({
        "from": sites[i]["name"], "to": sites[j]["name"], "edge_type": etype,
        "baseline_hours_boat_enabled": round(b, 4), "road_hours_land_only": round(c, 4),
        "delta_pct": delta_pct,
    })

meta = {
    "scope_note": (
        "Road-network FOUNDATION only, THIRD PASS (2026-08-20) after two rounds of Nico's "
        "review. Predecessor-path extraction on cost_distance.py's Tobler-hiking-function "
        "graph; THREE multiplicative friction layers (lithology S8f, biome S2, river-crossing "
        "S6); a LAND-ONLY graph (sea_mode='impassable') so the road network itself never "
        "crosses open water or a lake; a minimum spanning FOREST + tree-path-shortcut "
        "local-redundancy topology over 16 of the 17 LOCKED Circulo sites (Circulo_D_20k "
        "excluded, see excluded_sites below); and candidate ferry crossings (informational, "
        "separate edge_type) for any INCLUDED site pair the land-only graph still can't "
        "connect (none, this run). Rail's grade-ceiling cost function, actual ferry corridor "
        "authoring, and the kite wind-shadow mask are separate, not-yet-built pieces of "
        "Tappa 9 -- see docs/decisions/09_tappa9_transports.md."
    ),
    "excluded_sites": [
        {
            "name": s["name"], "tier": s["tier"], "population": s["population"],
            "reason": "Confirmed fully land-transport-isolated from the other 16 sites "
            "(2nd-pass run, land-only connected-components check). Nico's explicit call "
            "(3rd-pass review): a genuinely isolated island Circulo doesn't belong in a ROAD "
            "network's processing at all -- excluded here entirely rather than represented as "
            "a placeholder candidate-ferry-crossing line. Cheapest known boat-enabled "
            "connection to the rest of the network (2nd-pass run, boat-enabled combined-"
            "friction graph, for the future ferries sub-build's reference): "
            "Circulo_D_20k<->Circulo_E3_2k, 10.8757 h, 57.7 km.",
        }
        for s in excluded_sites
    ],
    "fixes_from_first_pass": {
        "sea_lake_crossings": "FIXED -- road graph now built with sea_mode='impassable' "
        "(cost_distance.py); no road edge can cross open water or a lake. First pass had 6/16 "
        "edges crossing open ocean and ~8/10 of the rest crossing a lake.",
        "river_crossings": "FIXED -- new river_crossing_friction_multiplier_120m layer "
        f"(Strahler order >= {MAJOR_STREAM_MIN_STRAHLER_ORDER}), stacked multiplicatively with "
        "lithology and biome friction.",
        "linear_topology": "FIXED -- MST + local redundancy (add_redundant_edges, revised "
        "again in this 3rd pass) adds loops where a nearby second connection ALSO provides a "
        "genuine, non-trivial shortcut over the existing tree route -- see "
        "src/transport/network.py's docstring for the exact criterion and two intermediate "
        "attempts that were tested and rejected first.",
        "crs_in_qgis": "NOT a code bug -- pre-existing project-wide GeoJSON export quirk "
        "(RFC 7946 readers ignore the crs/proj4 member), same as every other GeoJSON export "
        "in this project; manual CRS assignment in QGIS is the existing workaround.",
    },
    "fixes_from_second_pass": {
        "circulo_d_isolation": "FIXED per Nico's explicit call -- Circulo_D_20k excluded "
        "entirely from this run's site list (see excluded_sites above) rather than drawn as a "
        "candidate_ferry_crossing feature. The 2nd pass's own candidate-ferry-crossing "
        "mechanism (src/transport/network.py + the ferry-crossing block in this script) is "
        "left in place, generic and unchanged, for any FUTURE included-site disconnection --  "
        "it simply finds 0 components to bridge this run, since Circulo_D_20k (the only "
        "disconnected site) is no longer in the input at all.",
        "excessive_redundant_edges": "FIXED -- 3 of the 2nd pass's 8 redundant edges "
        "(GeoJSON feature IDs 15, 18, 20 in that run's export) were flagged directly by Nico "
        "as excessive. Root cause confirmed: the old redundancy criterion (within "
        "redundancy_factor=1.4 of a site's own cheapest neighbour ALONE) never checked whether "
        "the edge actually shortened anything versus the route the MST tree already provided "
        "-- two of the three flagged edges cost within 0.1% of the existing tree path (zero "
        "real benefit). `add_redundant_edges` now requires BOTH the original proximity check "
        "AND a genuine tree-path-shortcut improvement (>=20%) -- tested and confirmed this "
        "keeps exactly the 5 edges Nico did not flag and drops exactly the 3 he did, without "
        "introducing any new, unreviewed edges (two other approaches -- an O(n^2) candidate "
        "search, and swapping the criterion outright instead of ANDing it -- were tried first "
        "and rejected for over- or under-shooting; see src/transport/network.py's "
        "add_redundant_edges docstring for the full history and numbers).",
    },
    "topology": "minimum spanning FOREST (Prim's per connected component) + local-redundancy "
    "edges (BOTH redundancy_factor=1.4-of-cheapest-neighbour AND a >=20% tree-path-shortcut "
    "improvement required -- only added if the direct edge beats the existing MST-tree-path "
    "cost between the same two sites by at "
    "least that fraction), over symmetrized (mean of both directions) combined-friction "
    "cost-distance (hours) on the LAND-ONLY graph -- see src/transport/network.py's module "
    "docstring. PLUS a small set of MANUAL exception edges added by explicit, individually "
    "justified decision where the automatic rule's two gates (proximity-to-cheapest-neighbour "
    "AND tree-path-shortcut) disagreed on a specific pair -- see manual_extra_edges below; "
    "these are NOT a change to the automatic rule itself.",
    "connected_components": [
        {"members": [sites[i]["name"] for i in comp], "size": len(comp)} for comp in components
    ],
    "n_connected_components": len(components),
    "candidate_ferry_crossings": [
        {
            "from": sites[i]["name"], "to": sites[j]["name"],
            "boat_enabled_hours": round(float(w), 4),
            "note": "informational only -- NOT part of the road network; links two land-"
            "disconnected components (among the INCLUDED sites only -- see excluded_sites "
            "above) at their cheapest boat-enabled crossing point. Real ferry corridor "
            "authoring is a separate, not-yet-built Tappa 9 sub-build.",
        }
        for i, j, w in ferry_edges
    ],
    "river_crossing_friction_model": "WIDTH-BASED (revised same day as the sensitivity "
    "analysis -- supersedes the order-bucket table this field used to hold). "
    "friction(W) = 1 - k*W^p, continuous in estimated channel width W (metres), anchored to "
    "this model's own predecessor's two severity judgements so the 'typical' order-4/order-6 "
    "case is unchanged: see river_crossing_friction_width_anchors below. Floored at "
    f"{RIVER_FRICTION_WIDTH_ANCHORS['friction_at_order6_anchor']} (this domain's largest known "
    "river), ceilinged at 1.0 (~0 discharge -> not an obstacle). See "
    "src/transport/river_friction.py's river_friction_multiplier_from_width and "
    "src/transport/river_geometry.py for the full derivation.",
    "river_crossing_friction_width_anchors": RIVER_FRICTION_WIDTH_ANCHORS,
    "river_crossing_friction_table_SUPERSEDED_order_based": {
        f"strahler_order_{k}": v for k, v in RIVER_CROSSING_FRICTION.items()
    },
    "river_crossing_min_order": ALL_STREAMS_MIN_ORDER,
    "river_crossing_min_order_history": f"Was {MAJOR_STREAM_MIN_STRAHLER_ORDER} (order>=4 only) through "
    "the first width-based revision. Widened to ALL orders (1+) same day, after Nico noticed road "
    "routes crossing many order<4 streams in a short span with zero cost -- see "
    "stream_crossing_report below for the actual count this caught (was previously invisible).",
    "river_width_depth_model": {
        "method": "Leopold & Maddock (1953) downstream hydraulic geometry: "
        "W = ANCHOR_WIDTH_M * (Q/ANCHOR_DISCHARGE_M3S)^WIDTH_EXPONENT, same form for depth.",
        "width_exponent": WIDTH_EXPONENT, "depth_exponent": DEPTH_EXPONENT,
        "anchor_discharge_m3s": ANCHOR_DISCHARGE_M3S, "anchor_width_m": ANCHOR_WIDTH_M,
        "anchor_depth_m": ANCHOR_DEPTH_M,
        "status": "Exponents are widely-reproduced textbook values (real physical basis). "
        "The anchor width/depth at anchor_discharge_m3s is a CALIBRATION CHOICE (a plausible "
        "confined-channel magnitude for this domain's largest known discharge), not a cited "
        "coefficient -- UNREVIEWED, same status as every other friction estimate in this "
        "project. See src/transport/river_geometry.py's module docstring.",
        "estimated_width_m": {
            "min": round(float(_qualifying.min()), 3), "mean": round(float(_qualifying.mean()), 3),
            "max": round(float(_qualifying.max()), 3),
        },
        "low_order_high_discharge_cells": _low_order_big,
        "low_order_high_discharge_note": f"cells labelled order<{MAJOR_STREAM_MIN_STRAHLER_ORDER} in "
        "streams.geojson but with estimated width above 30% of this domain's largest known river -- a "
        "real Tappa 4 export quirk (some low-order-labelled reaches, almost certainly short stub "
        "reaches at a large river's mouth, carry anomalously large max_discharge_proxy_m3s), NOT fixed "
        "at that layer -- but no longer a road-network correctness problem, since friction now comes "
        "from discharge/width directly, not from the (sometimes wrong) order label.",
        "per_reach_export": f"data/processed/transport/stream_geometry_full.geojson "
        f"({len(reach_geometry['features'])} reaches, ALL orders (>= {ALL_STREAMS_MIN_ORDER}), "
        "estimated_width_m/estimated_depth_m properties added, geometry unmodified).",
    },
    "stream_crossing_report": {
        "method": "For every road edge's REAL reconstructed route, count contiguous runs of "
        "stream-covered cells (any order) as one crossing each -- checked directly against the "
        "route geometry, not estimated. Nico's own diagnostic (noticed several routes crossing many "
        "order<4 streams within a few km, e.g. 4 order-3 crossings in 4km), made permanent here.",
        "total_crossings_all_road_edges": total_crossings,
        "total_minor_crossings_order_lt_4": total_minor_crossings,
        "note": f"{total_minor_crossings} of {total_crossings} total crossings ({round(100*total_minor_crossings/total_crossings, 1) if total_crossings else 0}%) "
        f"are order<{MAJOR_STREAM_MIN_STRAHLER_ORDER} -- these were priced at literally ZERO before "
        "this revision (river_crossing_min_order was 4), regardless of how many of them a single road "
        "edge crossed. Now priced via the same continuous width-based function as every other "
        "crossing -- individually small for a minor stream, but a route crossing many of them in a "
        "short span now correctly accumulates real (if modest per-crossing) extra cost, addressing "
        "Nico's specific complaint that a road could cross several minor streams in a few km 'for "
        "free' where a detour might have been cheaper overall.",
        "per_edge": edge_crossing_report,
    },
    "biome_travel_friction_table": {
        BIOME_NAMES[code]: mult for code, mult in BIOME_TRAVEL_FRICTION.items()
    },
    "friction_status": "UNREVIEWED first-pass estimates for all three layers (lithology, "
    "biome, river-crossing) -- directional judgement calls, not a cited hiking-speed dataset. "
    "River-crossing friction is now width-based (see above), grounded in a real physical "
    "relation (Leopold-Maddock exponents) applied to this project's own computed discharge, "
    "one derivational step better-grounded than lithology/biome's still-direct judgement "
    "calls -- but the width model's own anchor coefficient remains a calibration choice, not "
    "independently measured. Pending Nico's sign-off; not written to config/parameters.yml -- "
    "same status as Tappa 8's LAND_TRAVEL_FRICTION and EXCAVATION_EFFORT_MULTIPLIER.",
    "friction_combination": "combined = lithology_friction * biome_friction * "
    "river_crossing_friction (all <=1.0, independent physical properties, combined "
    "multiplicatively, same 'friction stacks' logic as Tappa 8 S8f/S8g).",
    "resolution_m": [cs_x, cs_y],
    "manual_extra_edges": manual_edges_meta,
    "n_sites": len(sites),
    "n_mst_edges": len(mst_edges),
    "n_redundant_edges": len(redundant_edges),
    "n_manual_extra_edges": len(manual_edges),
    "n_road_edges_total": len(road_edges),
    "n_candidate_ferry_crossings": len(ferry_edges),
    "total_road_route_km": round(total_road_km, 2),
    "total_ferry_route_km": round(total_ferry_km, 2),
    "total_road_straight_line_km": round(
        sum(f["properties"]["straight_line_km"] for f in features if f["properties"]["edge_type"] != "candidate_ferry_crossing"), 2
    ),
    "graph_build_seconds": graph_build_s,
    "worst_edge_delta_pct_road_vs_baseline": worst_delta_pct,
    "road_edge_baseline_vs_combined_comparison": edge_comparisons,
    "connectivity_check": f"{len(components)} connected component(s) among {len(sites)} sites "
    "on the LAND-ONLY graph, checked directly (not assumed) before MST-forest construction. "
    + ("Fully connected by land." if len(components) == 1 else
       f"NOT fully connected by land -- {len(ferry_edges)} candidate ferry crossing(s) computed "
       "to link the remaining components (informational, see candidate_ferry_crossings above)."),
    "sites": [{k: v for k, v in s.items() if k not in ("row", "col")} for s in sites],
}
with open(f"{OUT}/tappa9_road_network_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log(f"=== DONE in {time.time() - t_start:.1f}s -- {len(mst_edges)} MST + {len(redundant_edges)} "
    f"redundant + {len(manual_edges)} manual exception = {len(road_edges)} road edges "
    f"({total_road_km:.1f} km), {len(ferry_edges)} candidate ferry crossing(s) "
    f"({total_ferry_km:.1f} km) ===")
