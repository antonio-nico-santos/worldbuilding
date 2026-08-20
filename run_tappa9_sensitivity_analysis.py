"""
Tappa 9 -- sensitivity analysis for the three UNREVIEWED friction tables
(lithology S8f, biome S9-first-pass, river-crossing S9-second-pass) plus the
two road-topology constants (`redundancy_factor`, `min_shortcut_improvement`
in `src/transport/network.py`'s `add_redundant_edges`).

WHY THIS EXISTS (Nico, 2026-08-20, same day as the third pass): asked how
these tables could be "calibrated." Direct answer, given here rather than in
chat only: there is no empirical dataset to fit against -- this is a
fictional world, not a real place with measured travel times. "Calibration"
in the usual model-fitting sense does not apply. What DOES apply, without
requiring any new data, is asking a narrower, answerable question: **for
each individual parameter, does changing it -- within a plausible range --
change anything the deliverable actually shows?** A parameter whose
plausible range never changes a single edge in the road network doesn't need
more calibration effort. A parameter that flips the topology on a small
perturbation is exactly where effort (real-world literature anchoring, or
narrative-fact anchoring, see the decision doc) would actually pay off.

METHOD: one-at-a-time (OAT) perturbation. Holding every other parameter at
its current shipped value, each individual table entry (a single lithology
class's multiplier, a single biome's multiplier, a single Strahler-order
river multiplier) or topology constant is perturbed to two alternates
(x0.7 and x1.3 of its current value, clipped to a friction multiplier's
valid range (0.05, 1.0]) and the FULL road-network pipeline -- combined
friction -> land-only graph -> all-pairs cost-distance -> MST forest ->
redundant edges -- is rerun from that single change. This is NOT a claim
that +-30% brackets the "true" uncertainty (there is no true value); it's a
fixed, consistent probe width applied uniformly so results are comparable
across all 21 parameters.

Deliberately reuses this project's already-shipped machinery
(`build_cost_graph`, `travel_friction_multiplier`, `biome_friction_multiplier`,
`river_friction_multiplier`, `compute_pairwise_cost_distance`,
`build_mst_forest`, `add_redundant_edges`) rather than reimplementing
anything -- the road-network graph this script perturbs IS the one
`run_tappa9_road_network.py` ships, not an approximation of it. Skips that
script's BOAT-ENABLED and BASELINE graphs entirely (candidate-ferry-crossing
geometry and the baseline-vs-friction comparison) since neither affects the
land-only road network's own topology, which is the only thing under test
here -- this keeps each perturbation run to ~4s instead of ~28s.

Reads the same inputs as run_tappa9_road_network.py. Writes:
  data/processed/transport/tappa9_sensitivity_analysis.json
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from geomorphology.lithology import LAND_TRAVEL_FRICTION, travel_friction_multiplier
from biomes.world_biomes import BIOME_NAMES
from suitability.cost_distance import build_cost_graph
from suitability.terrain_metrics import block_any, block_mean, block_mode
from transport.biome_friction import BIOME_TRAVEL_FRICTION, biome_friction_multiplier
from transport.river_friction import (
    MAJOR_STREAM_MIN_STRAHLER_ORDER,
    RIVER_CROSSING_FRICTION,
    rasterize_major_streams,
    river_friction_multiplier,
)
from transport.network import add_redundant_edges, build_mst_forest, compute_pairwise_cost_distance

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 9 -- friction/topology sensitivity analysis (OAT, x0.7 / x1.3 probes) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
XMIN, XMAX, YMIN, YMAX = domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"]

OUT = "data/processed/transport"
os.makedirs(OUT, exist_ok=True)

# --- reproduce the exact grid run_tappa9_road_network.py uses --------------
log("loading inputs...")
land = np.load("data/processed/climate/land_mask.npy").astype(bool)
lake30 = np.load("data/processed/hydrology/lake_mask.npy")
dem30 = np.load("data/processed/dem_v3_final_30m_eroded.npy")
lithology30 = np.load("data/processed/geomorphology/lithology_v6.npy")
biome_id = np.load("data/processed/suitability/biome_id_smoothed_120m.npy")

ny, nx = land.shape
cs_x = (XMAX - XMIN) / nx
cs_y = (YMAX - YMIN) / ny
cellsize_km = (cs_x + cs_y) / 2 / 1000.0

lake120 = block_any(lake30, 4)[:ny, :nx]
effective_land = land & ~lake120
dem120 = block_mean(dem30, 4)[:ny, :nx]

from geomorphology.lithology import LITHOLOGY_CLASSES
n_classes = max(LAND_TRAVEL_FRICTION.keys()) + 1
lithology120 = block_mode(lithology30, 4, n_classes)[:ny, :nx]

major_stream_mask, stream_order_grid = rasterize_major_streams(
    "data/exports/streams.geojson", XMIN, YMAX, cs_x, cs_y, (ny, nx),
    min_order=MAJOR_STREAM_MIN_STRAHLER_ORDER,
)

with open("data/processed/suitability/tappa6_site_selection_meta.json") as f:
    tappa6_meta = json.load(f)

EXCLUDED_FROM_ROAD_NETWORK = {"Circulo_D_20k"}
sites = []
for s in tappa6_meta["sites"]:
    if not s.get("placed") or s["name"] in EXCLUDED_FROM_ROAD_NETWORK:
        continue
    col = int(round((s["x_km"] - XMIN / 1000.0) / cellsize_km - 0.5))
    row = int(round((YMAX / 1000.0 - s["y_km"]) / cellsize_km - 0.5))
    sites.append({
        "name": s["name"], "tier": s["tier"], "population": s["population"],
        "row": row, "col": col, "x_km": s["x_km"], "y_km": s["y_km"],
    })
log(f"  {len(sites)} sites loaded (Circulo_D_20k excluded, same as the road-network script)")

BASE_REDUNDANCY_FACTOR = 1.4
BASE_MIN_SHORTCUT_IMPROVEMENT = 0.20


def run_pipeline(lith_table, biome_table, river_table, redundancy_factor, min_shortcut_improvement):
    """One full (fast) pipeline pass: combined friction -> land-only graph ->
    all-pairs cost-distance -> MST forest -> redundant edges. Returns
    (edge_set, n_mst, n_redundant, total_hours, total_straight_line_km)."""
    lith_friction = travel_friction_multiplier(lithology120, multipliers=lith_table)
    biome_fric = biome_friction_multiplier(biome_id, multipliers=biome_table)
    river_fric = river_friction_multiplier(stream_order_grid, multipliers=river_table)
    combined = (lith_friction * biome_fric * river_fric).astype(np.float32)

    graph_road = build_cost_graph(
        dem120, effective_land, cellsize_km, friction_multiplier=combined, sea_mode="impassable"
    )
    hours_road, _ = compute_pairwise_cost_distance(sites, graph_road, (ny, nx))
    mst_edges, components = build_mst_forest(hours_road)
    redundant_edges = add_redundant_edges(
        hours_road, mst_edges, components,
        redundancy_factor=redundancy_factor, min_shortcut_improvement=min_shortcut_improvement,
    )
    all_edges = [(i, j) for i, j, _ in mst_edges] + [(i, j) for i, j, _ in redundant_edges]
    edge_set = frozenset(tuple(sorted(e)) for e in all_edges)
    total_hours = sum(w for _, _, w in mst_edges) + sum(w for _, _, w in redundant_edges)
    total_straight_km = sum(
        float(np.hypot(sites[i]["x_km"] - sites[j]["x_km"], sites[i]["y_km"] - sites[j]["y_km"]))
        for i, j in edge_set
    )
    return {
        "edge_set": edge_set,
        "n_mst": len(mst_edges),
        "n_redundant": len(redundant_edges),
        "n_components": len(components),
        "total_hours": round(float(total_hours), 4),
        "total_straight_line_km": round(total_straight_km, 2),
    }


def edge_name_set(edge_set):
    return sorted(f"{sites[i]['name']}<->{sites[j]['name']}" for i, j in edge_set)


log("running BASELINE pipeline (current shipped values)...")
t0 = time.time()
baseline = run_pipeline(
    LAND_TRAVEL_FRICTION, BIOME_TRAVEL_FRICTION, RIVER_CROSSING_FRICTION,
    BASE_REDUNDANCY_FACTOR, BASE_MIN_SHORTCUT_IMPROVEMENT,
)
per_run_s = time.time() - t0
log(f"  done in {per_run_s:.1f}s -- {baseline['n_mst']} MST + {baseline['n_redundant']} redundant "
    f"= {len(baseline['edge_set'])} edges, {baseline['total_hours']:.2f}h total, "
    f"{baseline['total_straight_line_km']:.1f} km (straight-line sum)")

CLASS_LABELS = {
    1: "lithology:sedimentary_basin_fill", 2: "lithology:greywacke_argillite",
    3: "lithology:schist", 4: "lithology:volcanic", 5: "lithology:marble",
    6: "lithology:sedimentary_limestone", 7: "lithology:granite",
}

param_specs = []
for code, val in LAND_TRAVEL_FRICTION.items():
    param_specs.append(("lithology", CLASS_LABELS[code], code, val))
for code, val in BIOME_TRAVEL_FRICTION.items():
    param_specs.append(("biome", f"biome:{BIOME_NAMES[code]}", code, val))
for code, val in RIVER_CROSSING_FRICTION.items():
    param_specs.append(("river", f"river:strahler_order_{code}", code, val))

log(f"running {len(param_specs)} friction parameters x 2 probes (x0.7 / x1.3) "
    f"+ 2 topology constants x 2 probes -- estimated {(len(param_specs) * 2 + 4) * per_run_s:.0f}s total...")

results = []
for layer, label, code, base_val in param_specs:
    for probe_name, factor in (("x0.7", 0.7), ("x1.3", 1.3)):
        probe_val = min(1.0, max(0.05, round(base_val * factor, 4)))
        if layer == "lithology":
            table = dict(LAND_TRAVEL_FRICTION); table[code] = probe_val
            r = run_pipeline(table, BIOME_TRAVEL_FRICTION, RIVER_CROSSING_FRICTION,
                              BASE_REDUNDANCY_FACTOR, BASE_MIN_SHORTCUT_IMPROVEMENT)
        elif layer == "biome":
            table = dict(BIOME_TRAVEL_FRICTION); table[code] = probe_val
            r = run_pipeline(LAND_TRAVEL_FRICTION, table, RIVER_CROSSING_FRICTION,
                              BASE_REDUNDANCY_FACTOR, BASE_MIN_SHORTCUT_IMPROVEMENT)
        else:  # river
            table = dict(RIVER_CROSSING_FRICTION); table[code] = probe_val
            r = run_pipeline(LAND_TRAVEL_FRICTION, BIOME_TRAVEL_FRICTION, table,
                              BASE_REDUNDANCY_FACTOR, BASE_MIN_SHORTCUT_IMPROVEMENT)
        added = r["edge_set"] - baseline["edge_set"]
        removed = baseline["edge_set"] - r["edge_set"]
        results.append({
            "parameter": label, "layer": layer, "base_value": base_val,
            "probe": probe_name, "probe_value": probe_val,
            "edges_added": edge_name_set(added), "edges_removed": edge_name_set(removed),
            "n_edge_changes": len(added) + len(removed),
            "total_hours_delta_pct": round(100 * (r["total_hours"] - baseline["total_hours"]) / baseline["total_hours"], 3),
            "total_km_delta_pct": round(100 * (r["total_straight_line_km"] - baseline["total_straight_line_km"]) / baseline["total_straight_line_km"], 3),
        })
    log(f"  {label}: done")

for const_name, base_val, kind in (
    ("redundancy_factor", BASE_REDUNDANCY_FACTOR, "redundancy_factor"),
    ("min_shortcut_improvement", BASE_MIN_SHORTCUT_IMPROVEMENT, "min_shortcut_improvement"),
):
    for probe_name, factor in (("x0.7", 0.7), ("x1.3", 1.3)):
        probe_val = round(base_val * factor, 4)
        kwargs = {"redundancy_factor": BASE_REDUNDANCY_FACTOR, "min_shortcut_improvement": BASE_MIN_SHORTCUT_IMPROVEMENT}
        kwargs[kind] = probe_val
        r = run_pipeline(LAND_TRAVEL_FRICTION, BIOME_TRAVEL_FRICTION, RIVER_CROSSING_FRICTION, **kwargs)
        added = r["edge_set"] - baseline["edge_set"]
        removed = baseline["edge_set"] - r["edge_set"]
        results.append({
            "parameter": f"topology:{const_name}", "layer": "topology", "base_value": base_val,
            "probe": probe_name, "probe_value": probe_val,
            "edges_added": edge_name_set(added), "edges_removed": edge_name_set(removed),
            "n_edge_changes": len(added) + len(removed),
            "total_hours_delta_pct": round(100 * (r["total_hours"] - baseline["total_hours"]) / baseline["total_hours"], 3),
            "total_km_delta_pct": round(100 * (r["total_straight_line_km"] - baseline["total_straight_line_km"]) / baseline["total_straight_line_km"], 3),
        })
    log(f"  topology:{const_name}: done")

# --- rank parameters by impact ----------------------------------------------
by_param = {}
for r in results:
    p = r["parameter"]
    by_param.setdefault(p, {"layer": r["layer"], "max_edge_changes": 0, "max_abs_hours_delta_pct": 0.0})
    by_param[p]["max_edge_changes"] = max(by_param[p]["max_edge_changes"], r["n_edge_changes"])
    by_param[p]["max_abs_hours_delta_pct"] = max(by_param[p]["max_abs_hours_delta_pct"], abs(r["total_hours_delta_pct"]))

topology_sensitive = sorted(
    [p for p, v in by_param.items() if v["max_edge_changes"] > 0],
    key=lambda p: -by_param[p]["max_edge_changes"],
)
topology_robust = sorted(
    [p for p, v in by_param.items() if v["max_edge_changes"] == 0],
    key=lambda p: -by_param[p]["max_abs_hours_delta_pct"],
)

log(f"=== {len(topology_sensitive)}/{len(by_param)} parameters change at least one edge under "
    f"a x0.7/x1.3 probe; {len(topology_robust)} never change the edge set (cost-only effect) ===")
for p in topology_sensitive:
    log(f"  SENSITIVE: {p} (up to {by_param[p]['max_edge_changes']} edge change(s))")

out = {
    "method": "One-at-a-time (OAT) perturbation, x0.7 and x1.3 of each parameter's current "
    "shipped value (friction multipliers clipped to (0.05, 1.0]), full road-network pipeline "
    "rerun per probe (combined friction -> land-only graph -> all-pairs cost-distance -> MST "
    "forest -> redundant edges), compared against the baseline (currently shipped) edge set.",
    "purpose": "Answers 'which of these UNREVIEWED friction/topology values actually matter to "
    "the deliverable' -- NOT a claim about what the 'correct' values are (no ground truth exists "
    "for a fictional world). A parameter that never changes the edge set across this probe width "
    "is low-priority for further calibration effort; a parameter that does is where real-world "
    "literature anchoring or a narrative-stated fact would be worth spending effort on.",
    "n_sites": len(sites), "n_parameters_tested": len(by_param),
    "baseline": {
        "n_mst_edges": baseline["n_mst"], "n_redundant_edges": baseline["n_redundant"],
        "n_edges_total": len(baseline["edge_set"]), "total_hours": baseline["total_hours"],
        "total_straight_line_km": baseline["total_straight_line_km"],
        "edges": edge_name_set(baseline["edge_set"]),
    },
    "topology_sensitive_parameters": [
        {"parameter": p, "layer": by_param[p]["layer"], "max_edge_changes": by_param[p]["max_edge_changes"]}
        for p in topology_sensitive
    ],
    "topology_robust_parameters_by_cost_impact": [
        {"parameter": p, "layer": by_param[p]["layer"],
         "max_abs_total_hours_delta_pct": round(by_param[p]["max_abs_hours_delta_pct"], 3)}
        for p in topology_robust
    ],
    "all_probe_results": results,
}
with open(f"{OUT}/tappa9_sensitivity_analysis.json", "w") as f:
    json.dump(out, f, indent=2)

log(f"=== DONE in {time.time() - t_start:.1f}s -- wrote {OUT}/tappa9_sensitivity_analysis.json ===")
