"""
Tappa 10 -- connecting the auxiliary settlements (mining posts, mountain
huts, coastal villages, forest posts -- `auxiliary_settlements_tappa10_v2.geojson`)
to the road network Tappa 9 already built and locked (`road_network_mst.geojson`,
21 edges: 15 MST + 5 redundant + 1 manual exception).

**THIRD implementation, 2026-08-22, after Nico reviewed the second run in QGIS
and reported four concrete problems.** All four checked directly against the
actual geometry/data before fixing anything -- summary of what changed and
why:

1. **Line endpoints didn't touch the Circulo/settlement point features --
   "a little bit dislocated," on EVERY road, including the Tappa 9 backbone.**
   Root cause: every LineString's first/last vertex came from its path's
   source/target CELL CENTER (120m grid), never the site's own exact
   authored coordinate -- up to ~85m apart at this resolution. Fixed in both
   this script and `run_tappa9_road_network.py`: every line's two endpoints
   now use the exact site coordinate when the endpoint IS a site (a Circulo
   or a settlement); interior vertices still come from the real routed path
   (correct, and inherently grid-resolution). A tie-in point in the MIDDLE
   of another line has no more precise "true" location than its own cell
   center -- that's not a bug, so those points are left as-is.

2. **The Southland pocket around the excluded `Circulo_D_20k` had almost no
   internal road network** -- 14 settlements attached to it are correctly
   unreachable from the main 16-Circulo backbone (D_20k itself is a real
   land-isolated island, per Tappa 9), but the SECOND run only ever drew one
   edge among them (a satellite's spoke to its own hub) because that was the
   only relationship the source data already encoded. Nothing ever checked
   whether these leftover-isolated sites could at least reach EACH OTHER.
   Fixed: after the main connection pass (item 3 below), any settlement
   still unreachable from the main network is grouped with every OTHER
   still-unreachable settlement it CAN reach (mutual real cost-distance,
   same land-only graph) into a connected pocket, and each pocket of 2+
   members gets its own internal least-cost MST -- same machinery as the
   mountain-hut trail network, reused rather than reinvented. A genuinely
   solitary site (no other unreachable neighbor within reach either) stays
   isolated, honestly.

3 & 4. **Redundant, near-parallel paths, and a spoke that ran alongside an
   existing road for its ENTIRE length without ever joining it** (Nico's
   examples: `Coastal_Village_09` traced its own 30.2 km route to
   `Circulo_E2_2k` instead of a ~1.6 km hop onto `Mine_Limestone_07`'s
   already-built 21.5 km spoke to the SAME Circulo; `Forest_Post_03`'s spoke
   ran 170-420 m from the `Circulo_E3_2k<->Circulo_F7_small` road for its
   whole 6.1 km length because that road isn't between any of its OWN
   `attached_circulos`, so the old T-junction rule could never even
   consider it). **Same root cause for both, confirmed directly**: every
   settlement was connected independently against the fixed Tappa 9
   backbone only, filtered to its own `attached_circulos` -- never checking
   the Tappa 9 backbone MORE broadly, and never checking any OTHER
   settlement's already-built spoke as a possible shortcut. Nico's own
   suggestion (item 4) was exactly the fix: **the road network itself
   (Tappa 9's backbone) is now always a valid connection target for every
   settlement, not just edges between its own attached Círculos** -- and
   this script now also lets a settlement connect to another settlement's
   ALREADY-BUILT spoke/trail when that's cheaper than reaching the backbone
   directly.

**New algorithm, replacing the old independent-per-settlement + attached-
Círculos-filtered design**: a GREEDY, incrementally-growing network. Start
with the Tappa 9 backbone's own cells as the "current network." Repeatedly
find the (still-unconnected settlement, target cell) pair with the globally
CHEAPEST real cost-distance to the current network -- the target can be any
backbone cell, any Círculo, or any cell already claimed by a PREVIOUSLY
connected settlement's own spoke/trail in this same run -- connect it, add
its new path's cells to the network, repeat. This is the same idea as
Prim's algorithm for a minimum spanning tree, generalized to "grow a network
out from an already-fixed backbone," and it naturally produces exactly the
behavior Nico asked for: a settlement near an already-built spoke hops onto
it instead of duplicating the route, and a settlement near ANY part of the
Tappa 9 backbone (not just its own named Círculos) can tie in there. Mountain
huts keep their own separate per-massif trail network (built first, since it
has nothing to do with the valley backbone), then each trail's LEAF nodes
join this same greedy pool for their valley-access connection.

**Symmetrization note**: the old design could symmetrize every spoke's
weight (mean of both directions) because every target was always a Círculo,
and Círculo-sourced dist-only Dijkstra passes were already computed for
every settlement (Step A below). A generalized target (a point mid-way
along another settlement's own spoke) doesn't have that reverse distance
available without a full Dijkstra pass FROM that arbitrary point, which
would defeat the point of the greedy algorithm's speed. So: weight is
symmetrized (same convention as before) whenever the connection happens to
land exactly on a Círculo; otherwise the reported weight is one-way real
cost-distance, flagged as such via `weight_is_symmetrized`.

**Graph reused, not rebuilt**: the exact LAND-ONLY (`sea_mode="impassable"`)
combined-friction graph Tappa 9 already built and locked --
`combined_friction_multiplier_120m.npy` loaded directly from disk, no
friction recomputed.

Writes: `data/processed/transport/auxiliary_network_connections.geojson`
(spoke + trail + local-pocket-road LineStrings -- a separate layer from
`road_network_mst.geojson`) and `tappa10_network_connections_meta.json`
(full accounting: connection type per settlement, still-isolated list,
per-massif trail topology, isolated-pocket topology).
"""
from __future__ import annotations

import json
import time

import numpy as np

from params import load_params
from suitability.cost_distance import (
    build_cost_graph,
    cost_distance_from_source,
    cost_distance_from_source_with_predecessors,
    reconstruct_path,
)
from suitability.terrain_metrics import block_any, block_mean

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 10 -- connecting auxiliary settlements to the Tappa 9 road network "
    "(greedy network-growth, third implementation) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
XMIN, XMAX, YMIN, YMAX = domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"]

OUT = "data/processed/transport"

# --- reproduce the exact grid + LAND-ONLY combined-friction graph Tappa 9 --
# already built and locked -- load the friction field from disk rather than
# recompute it (lithology x biome x river-crossing, all already decided).
log("loading inputs (reproducing Tappa 6/8/9's grid exactly)...")
# **REVISED 2026-08-23** -- see run_tappa9_road_network.py's matching comment
# and build_land_mask_reconciled.py's docstring: loads the monotonic/
# additive-only land_mask_reconciled_v1.npy instead of the raw land_mask.npy,
# so Tappa 10's own connection-finding uses the identical grid the re-run
# backbone was built from.
land = np.load("data/processed/transport/land_mask_reconciled_v1.npy").astype(bool)
lake30 = np.load("data/processed/hydrology/lake_mask.npy")
dem30 = np.load("data/processed/dem_v3_final_30m_eroded.npy")

ny, nx = land.shape
cs_x = (XMAX - XMIN) / nx
cs_y = (YMAX - YMIN) / ny
cellsize_km = (cs_x + cs_y) / 2 / 1000.0

lake120 = block_any(lake30, 4)[:ny, :nx]
effective_land = land & ~lake120
dem120 = block_mean(dem30, 4)[:ny, :nx]

combined_friction = np.load(f"{OUT}/combined_friction_multiplier_120m.npy")
assert combined_friction.shape == (ny, nx), "combined_friction shape mismatch vs. this grid"

log("building the LAND-ONLY combined-friction graph (identical to Tappa 9's -- no friction "
    "recomputed, loaded straight from combined_friction_multiplier_120m.npy)...")
t0 = time.time()
graph_road = build_cost_graph(
    dem120, effective_land, cellsize_km, friction_multiplier=combined_friction, sea_mode="impassable"
)
log(f"  built in {time.time() - t0:.1f}s")


def xy_to_rc(x, y):
    col = int(round((x - XMIN) / cs_x - 0.5))
    row = int(round((YMAX - y) / cs_y - 0.5))
    return max(0, min(ny - 1, row)), max(0, min(nx - 1, col))


def rc_to_xy(r, c):
    x = XMIN + (c + 0.5) * cs_x
    y = YMAX - (r + 0.5) * cs_y
    return x, y


def path_to_coords(path, source_xy, target_xy=None):
    """Real routed path -> LineString coords, with exact endpoints (see
    module docstring, item 1). `source_xy` is always an exact site
    coordinate. `target_xy` is an exact site coordinate too when the target
    IS a known site (a Círculo, or another settlement reached exactly);
    pass None when the target is a mid-route tie-in point on another line,
    which has no truer location than its own cell center."""
    coords = [list(rc_to_xy(r, c)) for r, c in path]
    if coords:
        coords[0] = [float(source_xy[0]), float(source_xy[1])]
        if target_xy is not None:
            coords[-1] = [float(target_xy[0]), float(target_xy[1])]
    return coords


# --- the 17 Círculos, ALL of them (2026-08-23, TENTH-follow-up addendum:
# `Circulo_D_20k` reinstated into THIS script's own `circulos` list) --------
# `Circulo_D_20k` stays excluded from Tappa 9's own backbone
# (`run_tappa9_road_network.py`'s `EXCLUDED_FROM_ROAD_NETWORK`, UNCHANGED --
# it is still a genuine land-isolated island, no mainland road/ferry line is
# drawn for it there) -- that exclusion is about the INTER-Círculo backbone
# and is not touched here. But per Nico's explicit clarification, D_20k's
# OWN 14-settlement Southland pocket should be able to connect TO Círculo D
# itself, not just to each other -- the previous full exclusion (mirroring
# Tappa 9's) meant D_20k never appeared in `site_cell_lookup`/`circulo_idx`
# at all, so no settlement's spoke could ever land ON it, and the greedy
# network-growth loop had no seed cell there for anything to connect to in
# the first place. Fixed by including D_20k in `circulos` (below) AND
# separately seeding its own single site cell into the initial network
# (see the `net_r`/`net_c` seeding block) -- but NOT adding any edge from it
# to `road_network_mst.geojson`'s backbone cells. Since the land-only
# friction graph has zero path from D_20k's cell to the mainland (Tappa 9's
# own "fully land-transport-isolated" finding, re-verified below after this
# change), seeding just this one cell is safe by construction: it can only
# ever be reached by other cells on the SAME isolated landmass -- i.e.
# exactly the Southland settlements -- never by a mainland settlement.
with open("data/processed/suitability/tappa6_site_selection_meta.json") as f:
    tappa6_meta = json.load(f)

circulos = []
for s in tappa6_meta["sites"]:
    if not s.get("placed"):
        continue
    row = int(round((YMAX / 1000.0 - s["y_km"]) / cellsize_km - 0.5))
    col = int(round((s["x_km"] - XMIN / 1000.0) / cellsize_km - 0.5))
    circulos.append({"name": s["name"], "row": row, "col": col,
                      "x": s["x_km"] * 1000.0, "y": s["y_km"] * 1000.0})
circulo_idx = {c["name"]: i for i, c in enumerate(circulos)}
circulo_names = {c["name"] for c in circulos}
log(f"  {len(circulos)} Círculos loaded (all 17 -- Circulo_D_20k is back in this script's own "
    f"list as a LOCAL/island connection target only; it still has no mainland backbone edge, "
    f"see Tappa 9's own EXCLUDED_FROM_ROAD_NETWORK)")

# --- the v2 auxiliary settlements -------------------------------------------
with open("data/processed/suitability/auxiliary_settlements_tappa10_v2.geojson") as f:
    aux_fc = json.load(f)


def _split_list_field(v):
    """v2's schema flattens list-valued properties to semicolon-joined
    strings -- normalize back to a list here, once."""
    if v is None or v == "":
        return []
    if isinstance(v, list):
        return v
    return [part.strip() for part in str(v).split(";") if part.strip()]


settlements = []
for feat in aux_fc["features"]:
    p = dict(feat["properties"])
    p["attached_circulos"] = _split_list_field(p.get("attached_circulos"))
    p["resources"] = _split_list_field(p.get("resources"))
    x, y = feat["geometry"]["coordinates"]
    row, col = xy_to_rc(x, y)
    settlements.append({**p, "row": row, "col": col, "x": x, "y": y})
name_to_settlement = {s["name"]: s for s in settlements}
log(f"  {len(settlements)} auxiliary settlements loaded (v2)")

# --- fix stale attached_hub references (v2 data-quality issue) -------------
HUB_NAME_PREFIX_FIXES = [
    ("Mina_Calcario", "Mine_Limestone"),
    ("Mina_Granito", "Mine_Granite"),
    ("Mina_", "Mine_"),
]


def resolve_stale_hub_name(raw_name):
    if not raw_name:
        return None
    if raw_name in name_to_settlement:
        return raw_name
    for old_prefix, new_prefix in HUB_NAME_PREFIX_FIXES:
        if raw_name.startswith(old_prefix):
            candidate = new_prefix + raw_name[len(old_prefix):]
            if candidate in name_to_settlement:
                return candidate
    return None


n_hub_fixed = 0
n_hub_unresolved = 0
for s in settlements:
    raw_hub = s.get("attached_hub")
    if not raw_hub:
        continue
    resolved = resolve_stale_hub_name(raw_hub)
    if resolved is None:
        n_hub_unresolved += 1
        log(f"  WARNING: {s['name']}'s attached_hub {raw_hub!r} does not resolve -- leaving "
            f"as-is, will show up as isolated below")
    elif resolved != raw_hub:
        n_hub_fixed += 1
        s["attached_hub"] = resolved
log(f"  attached_hub stale-name fix: {n_hub_fixed} references corrected, {n_hub_unresolved} unresolved")

# --- site lookup by exact cell: every Círculo + every settlement, for
# exact-vertex snapping and for recognizing "this connection lands exactly
# on a known site" (as opposed to a mid-route tie-in point) -----------------
site_cell_lookup = {}
for c in circulos:
    site_cell_lookup[(c["row"], c["col"])] = ("circulo", c["name"], c["x"], c["y"])
for s in settlements:
    site_cell_lookup[(s["row"], s["col"])] = ("settlement", s["name"], s["x"], s["y"])

# --- Tappa 9 backbone -> seed the growing "current network" cell set ------
with open(f"{OUT}/road_network_mst.geojson") as f:
    road_fc = json.load(f)

net_r, net_c, net_label = [], [], []
for feat in road_fc["features"]:
    p = feat["properties"]
    if p["edge_type"] not in ("mst", "redundant", "manual_exception"):
        continue
    label = f"Circulo backbone: {p['from']} <-> {p['to']}"
    for x, y in feat["geometry"]["coordinates"]:
        r, c = xy_to_rc(x, y)
        net_r.append(r)
        net_c.append(c)
        net_label.append(label)
log(f"  {len(net_r)} backbone cells seeded as the initial connectable network "
    f"({sum(1 for pf in road_fc['features'] if pf['properties']['edge_type'] in ('mst', 'redundant', 'manual_exception'))} edges)")

# --- Circulo_D_20k's own site cell -> seeded SEPARATELY, NOT as part of the
# mainland backbone above (2026-08-23, tenth-follow-up addendum). This is
# the actual fix for "Southland settlements can't reach Círculo D itself":
# a single-cell seed, its own island-local anchor. It is provably safe to
# add unconditionally (no adjacency/reachability check needed here) because
# the land-only friction graph has already been confirmed (Tappa 9's own
# 2nd-pass land-only connected-components check, re-verified below) to have
# NO path at all between D_20k's cell and the mainland network's cells --
# so this seed can only ever be picked as "nearest" by another cell on the
# same isolated landmass. `site_cell_lookup` already maps D_20k's
# (row, col) to a "circulo" entry (built above, now that it's back in
# `circulos`), so a settlement landing exactly here is correctly typed
# `circulo_spoke`, not folded into the "Circulo backbone: ..." label.
circulo_d = next((c for c in circulos if c["name"] == "Circulo_D_20k"), None)
if circulo_d is not None:
    net_r.append(circulo_d["row"])
    net_c.append(circulo_d["col"])
    net_label.append("Circulo (island-local only, no mainland road): Circulo_D_20k")
    log("  +1 cell seeded for Circulo_D_20k itself (island-local anchor only, "
        "no connection to the mainland backbone)")


def add_to_network(path, label):
    for r, c in path:
        net_r.append(r)
        net_c.append(c)
        net_label.append(label)


def nearest_network_cell(dist_arr):
    """Cheapest (row, col, cost, label) among all cells currently in the
    growing network, for one settlement's precomputed dist array."""
    costs = dist_arr[np.array(net_r), np.array(net_c)]
    idx = int(np.argmin(costs))
    return net_r[idx], net_c[idx], float(costs[idx]), net_label[idx]


# --- Step A: one Círculo-sourced Dijkstra run per entry in `circulos` (17,
# since the tenth-follow-up addendum above -- Circulo_D_20k's own run will
# correctly return `inf` to every mainland settlement, real hours only to
# its own island's), dist-only, sampled at every settlement's cell -- gives
# circulo->settlement hours for symmetrized spoke weights whenever a
# connection happens to land exactly on a Círculo.
log(f"running {len(circulos)} Círculo-sourced Dijkstra passes (dist-only, for symmetrized "
    f"weights when a connection lands exactly on a Círculo)...")
t0 = time.time()
circulo_to_settlement_hours = np.full((len(circulos), len(settlements)), np.inf)
for ci, c in enumerate(circulos):
    dist = cost_distance_from_source(graph_road, c["row"], c["col"], (ny, nx))
    for si, s in enumerate(settlements):
        circulo_to_settlement_hours[ci, si] = dist[s["row"], s["col"]]
log(f"  done in {time.time() - t0:.1f}s")
settlement_index = {s["name"]: i for i, s in enumerate(settlements)}


def build_local_mst(members, dist_pred_by_name):
    """Least-cost MST among `members` (list of site dicts with name/row/col/
    x/y), using each member's own precomputed (dist, pred) full-grid arrays.
    Real symmetrized cost-distance as edge weight (Prim's, trivial at this
    size). Returns (edges, degree) where edges is a list of
    (i, j, weight_hours, path_cells)."""
    n = len(members)
    edges = []
    if n <= 1:
        return edges, {0: 0} if n == 1 else {}
    in_tree = [0]
    remaining = set(range(1, n))
    while remaining:
        best = None
        for i in in_tree:
            dist_i, _ = dist_pred_by_name[members[i]["name"]]
            for j in remaining:
                dist_j, _ = dist_pred_by_name[members[j]["name"]]
                r_j, c_j = members[j]["row"], members[j]["col"]
                r_i, c_i = members[i]["row"], members[i]["col"]
                w = 0.5 * (float(dist_i[r_j, c_j]) + float(dist_j[r_i, c_i]))
                if best is None or w < best[0]:
                    best = (w, i, j)
        w, i, j = best
        edges.append((i, j, w))
        in_tree.append(j)
        remaining.discard(j)
    degree = {i: 0 for i in range(n)}
    for i, j, w in edges:
        degree[i] += 1
        degree[j] += 1
    return edges, degree


# --- mountain huts: per-massif trail MST (unrelated to the valley backbone,
# built first), then each tree's LEAVES join the general connection pool
# below for their own valley-access spoke -----------------------------------
log("mountain huts: building per-massif trail networks...")
huts = [s for s in settlements if s["settlement_type"] == "mountain_hut"]
massifs = {}
for h in huts:
    massifs.setdefault(h.get("massif") or "unassigned", []).append(h)
log(f"  {len(huts)} huts across {len(massifs)} massif(s): "
    f"{ {m: len(v) for m, v in massifs.items()} }")

t0 = time.time()
hut_dist_pred = {}
for h in huts:
    dist_h, pred_h = cost_distance_from_source_with_predecessors(graph_road, h["row"], h["col"], (ny, nx))
    hut_dist_pred[h["name"]] = (dist_h, pred_h)
log(f"  {len(huts)} hut-sourced Dijkstra passes done in {time.time() - t0:.1f}s")

trail_features = []
trail_topology_report = {}
hut_leaves = []  # site dicts needing a valley-access connection

# Trail cells are NOT added to the shared `net_r`/`net_c` network immediately --
# a massif's trail only connects its own huts to EACH OTHER, not to any Círculo,
# until one of its leaves finds a real valley-access bridge below. Seeding trail
# cells into the shared network at construction time was a real bug (found by
# Nico, 2026-08-23 second report): every leaf's own site is itself a trail cell,
# so it would immediately find itself "already in the network" at zero cost and
# get silently skipped -- leaving the ENTIRE mountain system floating,
# disconnected from every Círculo, with no valley spoke ever built. Fix: keep
# each massif's cells in `massif_cluster_cells` until `promote_cluster` (below,
# in the greedy loop) brings the whole massif online the moment ANY of its
# leaves builds a real connection to the outside network. Each massif/hub entry
# is a list of (path, label) pairs -- one per contributing trail edge / satellite
# spoke -- rather than one flat cell list, so that once promoted, every cell
# keeps ITS OWN real origin label instead of collapsing to one generic
# "mountain trail (massif X)" / "hub cluster: Y" string. That collapse was a
# real, separate bug found while investigating Nico's 2026-08-23 "north part
# looks disconnected" report: it didn't affect connectivity (the cost/route were
# still computed against the correct nearest real cell), but it made
# `connects_to` on 12 delivered edges an uninformative internal label instead of
# naming the real trail segment or satellite spoke the settlement actually ties
# into -- see the "ninth follow-up" addendum in the decision doc.
massif_cluster_cells = {}
hut_leaf_massif = {}

for massif_name, members in massifs.items():
    edges, degree = build_local_mst(members, hut_dist_pred)
    for i, j, w in edges:
        hi, hj = members[i], members[j]
        dist_i, pred_i = hut_dist_pred[hi["name"]]
        path = reconstruct_path(pred_i, (ny, nx), hi["row"], hi["col"], hj["row"], hj["col"])
        coords = path_to_coords(path, (hi["x"], hi["y"]), (hj["x"], hj["y"]))
        route_km = sum(
            float(np.hypot(coords[k + 1][0] - coords[k][0], coords[k + 1][1] - coords[k][1]))
            for k in range(len(coords) - 1)
        ) / 1000.0
        trail_features.append({
            "type": "Feature",
            "properties": {
                "edge_type": "mountain_trail",
                "massif": massif_name, "from": hi["name"], "to": hj["name"],
                "edge_weight_hours": round(w, 4), "weight_is_symmetrized": True,
                "route_km": round(route_km, 3), "n_cells": len(path),
            },
            "geometry": {"type": "LineString", "coordinates": coords},
        })
        massif_cluster_cells.setdefault(massif_name, []).append(
            (path, f"mountain trail: {hi['name']} <-> {hj['name']}"))
    n = len(members)
    leaves = [i for i in range(n) if degree.get(i, 0) <= 1]
    hut_leaves.extend(members[i] for i in leaves)
    for i in leaves:
        hut_leaf_massif[members[i]["name"]] = massif_name
    trail_topology_report[massif_name] = {
        "n_huts": n, "n_trail_edges": len(edges), "n_leaves": len(leaves),
        "hut_names": [m["name"] for m in members], "leaf_names": [members[i]["name"] for i in leaves],
    }
    log(f"  massif '{massif_name}': {n} huts, {len(edges)} trail edge(s), {len(leaves)} leaf/leaves")

log(f"  {len(trail_features)} mountain-trail edges built; {len(hut_leaves)} leaves need a "
    f"valley-access connection")

# --- precompute dist+pred for every mining/coastal/forest settlement ------
SPOKE_TYPES = {"mining_post", "coastal_village", "forest_post"}
log(f"running one Dijkstra pass per mining/coastal/forest settlement "
    f"({sum(1 for s in settlements if s['settlement_type'] in SPOKE_TYPES)} sites)...")
t0 = time.time()
settlement_dist_pred = {}
n_done = 0
for s in settlements:
    if s["settlement_type"] not in SPOKE_TYPES:
        continue
    n_done += 1
    settlement_dist_pred[s["name"]] = cost_distance_from_source_with_predecessors(
        graph_road, s["row"], s["col"], (ny, nx))
    if n_done % 20 == 0:
        log(f"  {n_done} settlements processed...")
log(f"  done in {time.time() - t0:.1f}s")

# --- satellites: hard administrative rule, connect to their HUB first,
# not part of the greedy network-growth pool below --------------------------
connected_features = []
connection_summary = {}


def bump(counter, key):
    connection_summary[key] = connection_summary.get(key, 0) + 1


satellites = [s for s in settlements
              if s["settlement_type"] in SPOKE_TYPES
              and s.get("structure_tier") == "satellite_no_own_structure" and s.get("attached_hub")]
n_satellite_isolated = 0
satellite_isolated_report = []
# Same fix as the mountain-trail one, generalized: a satellite's spoke is a
# DIRECT satellite->hub route (hard rule, computed regardless of whether the
# hub itself reaches anywhere) -- it says nothing about whether the hub (and
# so the whole little satellite cluster) is actually connected to a Círculo.
# Seeding it into the shared network immediately let a hub with no real
# outside connection still register as "already in the network" the moment
# it was itself processed below (found by Nico, 2026-08-23: "hubs of two
# mines that only connect them and don't have road to other places"). Fix:
# keep each hub's satellite-spoke cells in `hub_cluster_cells` until
# `promote_cluster` (in the greedy loop) brings the whole cluster online the
# moment the HUB itself builds a real connection to the outside network.
hub_cluster_cells = {}
for s in satellites:
    hub = name_to_settlement.get(s["attached_hub"])
    if hub is None:
        n_satellite_isolated += 1
        satellite_isolated_report.append({
            "name": s["name"], "settlement_type": s["settlement_type"],
            "attached_circulos": s.get("attached_circulos"),
            "reason": f"attached_hub '{s.get('attached_hub')}' does not resolve to any "
            "known settlement -- unrelated to the routing graph, a data-reference problem.",
        })
        continue
    dist_s, pred_s = settlement_dist_pred[s["name"]]
    one_way_hours = float(dist_s[hub["row"], hub["col"]])
    if not np.isfinite(one_way_hours):
        n_satellite_isolated += 1
        satellite_isolated_report.append({
            "name": s["name"], "settlement_type": s["settlement_type"],
            "attached_circulos": s.get("attached_circulos"),
            "reason": f"no finite-cost land route to its own hub ('{s['attached_hub']}') on "
            "this land-only graph -- a hard administrative rule (satellite -> own hub only) "
            "with no fallback, so this is NOT retried against the wider network.",
        })
        continue
    path = reconstruct_path(pred_s, (ny, nx), s["row"], s["col"], hub["row"], hub["col"])
    coords = path_to_coords(path, (s["x"], s["y"]), (hub["x"], hub["y"]))
    route_km = sum(
        float(np.hypot(coords[k + 1][0] - coords[k][0], coords[k + 1][1] - coords[k][1]))
        for k in range(len(coords) - 1)
    ) / 1000.0
    bump(connection_summary, "satellite_to_hub_spoke")
    connected_features.append({
        "type": "Feature",
        "properties": {
            "edge_type": "satellite_to_hub_spoke", "settlement_name": s["name"],
            "settlement_type": s["settlement_type"], "structure_tier": s.get("structure_tier"),
            "connects_to": hub["name"], "one_way_hours": round(one_way_hours, 4),
            "edge_weight_hours": round(one_way_hours, 4), "weight_is_symmetrized": False,
            "route_km": round(route_km, 3), "n_cells": len(path),
        },
        "geometry": {"type": "LineString", "coordinates": coords},
    })
    hub_cluster_cells.setdefault(hub["name"], []).append((path, f"spoke: {s['name']}"))
log(f"  {len(satellites) - n_satellite_isolated} satellite-to-hub spokes built "
    f"({n_satellite_isolated} satellites isolated -- no finite path to their own hub)")

# --- greedy network-growth: every non-satellite mining/coastal/forest
# settlement + every hut trail leaf, connecting to whichever is cheapest --
# the Tappa 9 backbone directly, OR any spoke/trail already built THIS RUN.
satellite_names = {s["name"] for s in satellites}
pool = [s for s in settlements
        if s["settlement_type"] in SPOKE_TYPES and s["name"] not in satellite_names]
pool_dist_pred = dict(settlement_dist_pred)
pool.extend(hut_leaves)
pool_dist_pred.update(hut_dist_pred)
unconnected = {s["name"]: s for s in pool}

promoted_hub_clusters = set()
promoted_massifs = set()
n_clusters_promoted = 0


def promote_cluster(name):
    """Once `name` has a real (or trivially-already-there) connection to the
    live network, bring its whole local cluster online too: if `name` is a
    hub, every satellite spoke that already terminates on it is now REALLY
    reaching the outside world, not just each other; if `name` is a hut trail
    leaf, its whole massif (reachable via the tree itself) is now REALLY
    reaching the outside world too. See the two comments above (satellite
    loop, trail loop) for why these were kept OUT of the network until now."""
    global n_clusters_promoted
    if name in hub_cluster_cells and name not in promoted_hub_clusters:
        for path, label in hub_cluster_cells[name]:
            add_to_network(path, label)
        promoted_hub_clusters.add(name)
        n_clusters_promoted += 1
    massif = hut_leaf_massif.get(name)
    if massif is not None and massif not in promoted_massifs:
        for path, label in massif_cluster_cells.get(massif, []):
            add_to_network(path, label)
        promoted_massifs.add(massif)
        n_clusters_promoted += 1


log(f"greedy network-growth over {len(unconnected)} settlements/trailheads...")
t0 = time.time()
n_connected_this_phase = 0
n_trivially_absorbed = 0
while unconnected:
    best = None  # (cost, name, target_r, target_c, label)
    for name, s in unconnected.items():
        dist_s, _ = pool_dist_pred[name]
        r, c, cost, label = nearest_network_cell(dist_s)
        if best is None or cost < best[0]:
            best = (cost, name, r, c, label)
    cost, name, target_r, target_c, label = best
    if not np.isfinite(cost):
        break  # every remaining settlement is now unreachable from the current network
    s = unconnected.pop(name)
    site_info = site_cell_lookup.get((target_r, target_c))
    if site_info is not None and site_info[1] == name:
        # This settlement's OWN site cell was already absorbed into the growing
        # network -- typically because another settlement's spoke path (e.g. a
        # satellite_to_hub_spoke) terminates exactly here, on this settlement's
        # own site. It is therefore already connected at zero cost; emitting a
        # "feeder_t_junction_spoke" edge from it to itself would be a spurious
        # zero-length self-loop. Nothing further to route -- just mark it done.
        n_trivially_absorbed += 1
        promote_cluster(name)
        continue
    dist_s, pred_s = pool_dist_pred[name]
    if site_info is not None and site_info[0] == "circulo":
        edge_type = "circulo_spoke"
        target_xy = (site_info[2], site_info[3])
        connects_to = site_info[1]
        ci = circulo_idx[site_info[1]]
        symmetrized_hours = 0.5 * (cost + circulo_to_settlement_hours[ci, settlement_index.get(name, -1)]) \
            if name in settlement_index else None
    elif site_info is not None:  # lands exactly on another settlement's own site cell
        edge_type = "feeder_t_junction_spoke"
        target_xy = (site_info[2], site_info[3])
        connects_to = site_info[1]
        symmetrized_hours = None
    elif label.startswith("Circulo backbone:"):
        edge_type = "backbone_t_junction_spoke"
        target_xy = None
        connects_to = label.replace("Circulo backbone: ", "")
        symmetrized_hours = None
    else:
        edge_type = "feeder_t_junction_spoke"
        target_xy = None
        connects_to = label.replace("mountain trail: ", "").replace("spoke: ", "")
        symmetrized_hours = None

    path = reconstruct_path(pred_s, (ny, nx), s["row"], s["col"], target_r, target_c)
    coords = path_to_coords(path, (s["x"], s["y"]), target_xy)
    route_km = sum(
        float(np.hypot(coords[k + 1][0] - coords[k][0], coords[k + 1][1] - coords[k][1]))
        for k in range(len(coords) - 1)
    ) / 1000.0
    weight_hours = symmetrized_hours if symmetrized_hours is not None else cost
    bump(connection_summary, edge_type)
    connected_features.append({
        "type": "Feature",
        "properties": {
            "edge_type": edge_type, "settlement_name": s["name"],
            "settlement_type": s["settlement_type"], "structure_tier": s.get("structure_tier"),
            "connects_to": connects_to, "one_way_hours": round(cost, 4),
            "edge_weight_hours": round(weight_hours, 4),
            "weight_is_symmetrized": symmetrized_hours is not None,
            "route_km": round(route_km, 3), "n_cells": len(path),
        },
        "geometry": {"type": "LineString", "coordinates": coords},
    })
    add_to_network(path, f"spoke: {s['name']}")
    promote_cluster(name)
    n_connected_this_phase += 1
    if n_connected_this_phase % 20 == 0:
        log(f"  {n_connected_this_phase} connected, {len(unconnected)} remaining...")

log(f"  done in {time.time() - t0:.1f}s -- {n_connected_this_phase} settlements/trailheads "
    f"connected to the growing network, {n_trivially_absorbed} already absorbed at zero cost "
    f"(another spoke terminated exactly on their site), {len(unconnected)} unreachable from it entirely")
log(f"  {n_clusters_promoted} local cluster(s) (hub-satellite groups / hut massifs) promoted "
    f"into the live network -- {len(hub_cluster_cells) - len(promoted_hub_clusters)} hub "
    f"cluster(s) and {len(massifs) - len(promoted_massifs)} massif(s) never got a real outside "
    f"connection and stayed local-only")
log(f"  connection type breakdown: {connection_summary}")

# --- isolated-pocket local roads: among what's left, group by MUTUAL real
# reachability (still on this same land-only graph, just cut off from the
# main network) and build one local MST per pocket of 2+ members ----------
leftover = list(unconnected.values())
pocket_features = []
pocket_report = []
still_isolated = []
if leftover:
    log(f"grouping {len(leftover)} leftover-unreachable settlements into local pockets...")
    n = len(leftover)
    adj = [[] for _ in range(n)]
    for i in range(n):
        dist_i, _ = pool_dist_pred[leftover[i]["name"]]
        for j in range(n):
            if i == j:
                continue
            r_j, c_j = leftover[j]["row"], leftover[j]["col"]
            if np.isfinite(dist_i[r_j, c_j]):
                adj[i].append(j)
    visited = [False] * n
    components = []
    for i in range(n):
        if visited[i]:
            continue
        stack, comp = [i], []
        visited[i] = True
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
        components.append(comp)
    log(f"  {len(components)} pocket(s) found ({sum(1 for c in components if len(c) > 1)} "
        f"connectable internally, {sum(1 for c in components if len(c) == 1)} genuinely solitary)")

    for comp in components:
        members = [leftover[i] for i in comp]
        if len(members) == 1:
            m = members[0]
            still_isolated.append({
                "name": m["name"], "settlement_type": m["settlement_type"],
                "attached_circulos": m.get("attached_circulos"),
                "reason": "no finite-cost land route to the main network OR to any other "
                "still-unreachable settlement -- genuinely solitary on this land-only graph.",
            })
            continue
        edges, degree = build_local_mst(members, pool_dist_pred)
        for i, j, w in edges:
            mi, mj = members[i], members[j]
            dist_i, pred_i = pool_dist_pred[mi["name"]]
            path = reconstruct_path(pred_i, (ny, nx), mi["row"], mi["col"], mj["row"], mj["col"])
            coords = path_to_coords(path, (mi["x"], mi["y"]), (mj["x"], mj["y"]))
            route_km = sum(
                float(np.hypot(coords[k + 1][0] - coords[k][0], coords[k + 1][1] - coords[k][1]))
                for k in range(len(coords) - 1)
            ) / 1000.0
            bump(connection_summary, "isolated_pocket_road")
            pocket_features.append({
                "type": "Feature",
                "properties": {
                    "edge_type": "isolated_pocket_road", "from": mi["name"], "to": mj["name"],
                    "edge_weight_hours": round(w, 4), "weight_is_symmetrized": True,
                    "route_km": round(route_km, 3), "n_cells": len(path),
                },
                "geometry": {"type": "LineString", "coordinates": coords},
            })
        pocket_report.append({
            "n_members": len(members), "n_edges": len(edges),
            "member_names": [m["name"] for m in members],
        })
        log(f"  pocket ({len(members)} members): {[m['name'] for m in members]} -- "
            f"{len(edges)} local road(s) built")
else:
    log("  nothing left unreachable -- every settlement connected to the main network")

still_isolated.extend(satellite_isolated_report)
if satellite_isolated_report:
    log(f"  +{len(satellite_isolated_report)} satellite(s) isolated from their own hub "
        f"(hard rule, not retried against the network) added to still_isolated")

# --- Ninth follow-up (2026-08-23): three manual connections Nico requested
# after real cost-distance confirmed the north schist cluster's only bridge
# to the backbone is a single, distant spur (130-172 km real route against a
# 10-14 km straight-line separation). Same "named manual exception" pattern
# Tappa 9 already established for its own Circulo_F1_small<->Circulo_F2_small
# edge (see run_tappa9_road_network.py's MANUAL_EXTRA_EDGES): added directly
# rather than tuning the greedy algorithm's cluster-promotion logic, which by
# design only ever finds ONE bridge per cluster and has no notion of "also
# build a second, closer one." These do NOT replace the automated connections
# already built above (still present, unchanged) -- they're additional,
# deliberately shorter routes to a specific named Círculo, each confirmed
# with real cost-distance (not straight-line) before being added. Anchors
# chosen with Nico: Mine_Schist_06 (the real hub Mine_Schist_05/08 both
# attach to, cheapest of the three candidates offered) for Circulo_F1_small,
# Mine_Schist_03 for Circulo_F2_small, and Outpost_MainSpine_16 (Nico's own
# explicit pick, no alternative offered) for Circulo_E3_2k.
log("building 4 manual connections (Nico's explicit picks, ninth + eleventh follow-up)...")
MANUAL_CONNECTIONS = [
    ("Mine_Schist_06", "Circulo_F1_small", "Nico's explicit request, 2026-08-23 ninth "
     "follow-up -- the settlement's automated connection (still present, unchanged) reaches "
     "the backbone only via a single distant massif bridge; this direct edge is real "
     "cost-distance to its named Círculo, roughly 90% cheaper than that route."),
    ("Mine_Schist_03", "Circulo_F2_small", "Nico's explicit request, 2026-08-23 ninth "
     "follow-up -- same north-cluster long-detour finding as Mine_Schist_06 above, ~90% "
     "cheaper than the automated route."),
    ("Outpost_MainSpine_16", "Circulo_E3_2k", "Nico's explicit request, 2026-08-23 ninth "
     "follow-up -- the automated connection reaches the backbone only via a single distant "
     "hub-cluster bridge; this direct edge is real cost-distance to its named Círculo, "
     "roughly 45% cheaper than that route."),
    # Eleventh follow-up (2026-08-23): Nico asked for a direct
    # Coastal_Village_07 -> Circulo_C_25k connection. Currently 13.58 km
    # away as a feeder_t_junction_spoke onto Mine_Greywacke_01's own spoke,
    # a real network distance of 9.49h to Circulo_C_25k; a direct connection
    # costs 4.54h/10.29 km by real cost-distance -- a genuine ~52% cut.
    ("Coastal_Village_07", "Circulo_C_25k", "Nico's explicit request, 2026-08-23 eleventh "
     "follow-up -- the automated connection (still present, unchanged) reaches Circulo_C_25k "
     "only via a 13.58 km feeder hop onto Mine_Greywacke_01's own spoke, a real network "
     "distance of 9.49h; this direct edge is real cost-distance to Circulo_C_25k, 4.54h -- "
     "roughly 52% cheaper than that route."),
]
manual_features = []
manual_report = []
for settlement_name, circulo_name, manual_reason in MANUAL_CONNECTIONS:
    s = name_to_settlement[settlement_name]
    c = circulos[circulo_idx[circulo_name]]
    if s["settlement_type"] == "mountain_hut":
        dist_s, pred_s = hut_dist_pred[settlement_name]
    else:
        dist_s, pred_s = settlement_dist_pred[settlement_name]
    cost_hr = float(dist_s[c["row"], c["col"]])
    path = reconstruct_path(pred_s, (ny, nx), s["row"], s["col"], c["row"], c["col"])
    coords = path_to_coords(path, (s["x"], s["y"]), (c["x"], c["y"]))
    route_km = sum(
        float(np.hypot(coords[k + 1][0] - coords[k][0], coords[k + 1][1] - coords[k][1]))
        for k in range(len(coords) - 1)
    ) / 1000.0
    bump(connection_summary, "manual_connection")
    manual_features.append({
        "type": "Feature",
        "properties": {
            "edge_type": "manual_connection", "settlement_name": settlement_name,
            "settlement_type": s["settlement_type"], "structure_tier": s.get("structure_tier"),
            "connects_to": circulo_name, "one_way_hours": round(cost_hr, 4),
            "edge_weight_hours": round(cost_hr, 4), "weight_is_symmetrized": False,
            "route_km": round(route_km, 3), "n_cells": len(path),
            "manual_reason": manual_reason,
        },
        "geometry": {"type": "LineString", "coordinates": coords},
    })
    manual_report.append({
        "settlement_name": settlement_name, "connects_to": circulo_name,
        "route_km": round(route_km, 3), "edge_weight_hours": round(cost_hr, 4),
    })
    log(f"  {settlement_name} -> {circulo_name}: {route_km:.2f} km / {cost_hr:.2f} h")

# --- assemble output ---------------------------------------------------
all_features = trail_features + connected_features + pocket_features + manual_features
out_fc = {
    "type": "FeatureCollection",
    "name": "auxiliary_network_connections",
    "crs": {"type": "proj4", "properties": {"proj4": CRS_PROJ4}},
    "features": all_features,
}
with open(f"{OUT}/auxiliary_network_connections.geojson", "w") as f:
    json.dump(out_fc, f, indent=2)

total_trail_km = sum(f["properties"]["route_km"] for f in trail_features)
total_spoke_km = sum(f["properties"]["route_km"] for f in connected_features)
total_pocket_km = sum(f["properties"]["route_km"] for f in pocket_features)

meta = {
    "scope_note": "Connects the Tappa 10 v2 auxiliary settlements to the already-locked "
    "Tappa 9 road network via a GREEDY, incrementally-growing network (Prim's-style): each "
    "settlement connects to whichever is cheapest -- the Tappa 9 backbone at ANY point (not "
    "just edges between its own attached Círculos), or another settlement's ALREADY-BUILT "
    "spoke/trail -- fixing redundant near-parallel paths and roads that ran alongside an "
    "existing edge without joining it. Mountain huts keep a separate per-massif trail "
    "network built first; only their trail leaves join the general pool. Settlements left "
    "unreachable from the main network are grouped by mutual reachability into local "
    "pockets (e.g. the Circulo_D_20k cluster) and each pocket of 2+ gets its own internal "
    "least-cost MST, so a genuinely isolated region still gets a sensible local road system "
    "instead of disconnected dots. See this script's own module docstring for the full "
    "reasoning behind each of the four problems this run fixed.",
    "graph": "identical LAND-ONLY combined-friction graph Tappa 9 built and locked -- "
    "combined_friction_multiplier_120m.npy loaded from disk, nothing recomputed.",
    "vertex_fix": "every LineString's endpoints now use the exact site coordinate when the "
    "endpoint is a known site (Círculo or settlement) -- previously snapped to the 120m "
    "cell center, visibly offset from the point features in QGIS. Applied here AND in "
    "run_tappa9_road_network.py (the backbone itself).",
    "v2_data_fixes": {
        "attached_hub_stale_names_fixed": n_hub_fixed,
        "attached_hub_unresolved": n_hub_unresolved,
    },
    "trail_topology": trail_topology_report,
    "isolated_pockets": pocket_report,
    "n_settlements_total": len(settlements),
    "n_settlements_by_type": {
        t: sum(1 for s in settlements if s["settlement_type"] == t)
        for t in ("mining_post", "mountain_hut", "coastal_village", "forest_post")
    },
    "connection_type_breakdown": connection_summary,
    "n_trail_edges": len(trail_features),
    "n_network_spokes": len(connected_features),
    "n_trivially_absorbed": n_trivially_absorbed,
    "n_clusters_promoted": n_clusters_promoted,
    "n_hub_clusters_never_bridged": len(hub_cluster_cells) - len(promoted_hub_clusters),
    "n_massifs_never_bridged": len(massifs) - len(promoted_massifs),
    "n_pocket_roads": len(pocket_features),
    "n_still_isolated": len(still_isolated),
    "still_isolated": still_isolated,
    "manual_connections": manual_report,
    "n_manual_connections": len(manual_features),
    "total_trail_route_km": round(total_trail_km, 2),
    "total_spoke_route_km": round(total_spoke_km, 2),
    "total_pocket_route_km": round(total_pocket_km, 2),
    "total_manual_route_km": round(sum(f["properties"]["route_km"] for f in manual_features), 2),
}
with open(f"{OUT}/tappa10_network_connections_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

total_manual_km = sum(f["properties"]["route_km"] for f in manual_features)
log(f"=== DONE in {time.time() - t_start:.1f}s -- {len(trail_features)} mountain-trail edges "
    f"({total_trail_km:.1f} km) + {len(connected_features)} network spokes "
    f"({total_spoke_km:.1f} km) + {len(pocket_features)} isolated-pocket roads "
    f"({total_pocket_km:.1f} km) + {len(manual_features)} manual connections "
    f"({total_manual_km:.1f} km), {len(still_isolated)} settlement(s) genuinely solitary ===")
