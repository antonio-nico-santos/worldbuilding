"""
Tappa 8 -- transport lithology-cost multiplier (Nico's post-Tappa-8 follow-up
list, item 1: "Create a weight cost to calculate how this lithology could
help/difficult the construction of the transport system").

SCOPE, confirmed with Nico via AskUserQuestion before building this: a
rock-type FRICTION WEIGHT layered onto Tappa 6's existing cost-distance graph
(src/suitability/cost_distance.py -- Tobler's hiking function on land, a flat
boat speed on any sea-touching edge), NOT full road/rail network design
(predecessor-path extraction, biome-differentiated costs, a rail grade
ceiling -- Nico explicitly declined that larger scope, calling it "starting
Tappa 9"). See geomorphology/lithology.py's LAND_TRAVEL_FRICTION /
travel_friction_multiplier() docstring for the actual per-class values,
their (non-quantitative) real-world grounding, and an important reframe:
Nico's own phrasing was about "construction" (excavation/foundation/tunnel
cost), which this does NOT model -- there is no citable unit-cost data for
that anywhere in this project. What's built here is the nearest defensible
proxy the EXISTING cost graph can represent: a travel-TIME friction factor.

IMPORTANT -- this does NOT touch anything already locked. Tappa 6's actual
Circulo placement (circulo_candidate_sites.geojson,
tappa6_site_selection_meta.json) is read-only input here, never
overwritten. This script builds a SEPARATE friction-adjusted cost graph
alongside the original baseline graph and directly compares pairwise
cost-distance (hours) between the 17 ALREADY-PLACED sites under both, so
Nico can see exactly how much this new layer would have changed things --
re-running site_selection.place_circulos() WITH this friction (i.e.
actually re-placing the 17 Circulos under it) is a distinct, separate
follow-up decision, not made here.

Reads:
  data/processed/climate/land_mask.npy
  data/processed/hydrology/lake_mask.npy       (native 30m, same block-ANY
                                                convention as Tappa 6)
  data/processed/dem_v3_final_30m_eroded.npy   (native 30m, block-MEAN to
                                                120m, identical to Tappa 6)
  data/processed/geomorphology/lithology_v6.npy (native 30m, categorical --
                                                block-MODE (majority vote)
                                                to 120m, NEW downsample kind
                                                for this pipeline)
  data/processed/suitability/tappa6_site_selection_meta.json (the 17
                                                ALREADY-PLACED sites' x_km/
                                                y_km/tier -- read-only)

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  transport_friction_multiplier_120m.*   (float32 raster, the multiplier
                                          field itself, for QGIS/visual QA)
  tappa8_transport_friction_meta.json    (multiplier table + citations,
                                          class areas, and the full
                                          baseline-vs-friction pairwise
                                          comparison table)
"""
import sys, time, json
from functools import partial

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import LAND_TRAVEL_FRICTION, travel_friction_multiplier, LITHOLOGY_CLASSES
from suitability.cost_distance import BOAT_SPEED_KMH, build_cost_graph, cost_distance_from_source
from suitability.terrain_metrics import block_mean, block_any, block_mode

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- transport lithology-cost multiplier (Tappa 6 cost_distance follow-up) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
XMIN, XMAX, YMIN, YMAX = domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"]

OUT = "data/processed/suitability"

# --- reproduce Tappa 6's exact 120m grid -------------------------------------------
log("loading Tappa 6 inputs (reproducing run_tappa6_site_selection.py's grid exactly)...")
land = np.load("data/processed/climate/land_mask.npy").astype(bool)
lake30 = np.load("data/processed/hydrology/lake_mask.npy")
dem30 = np.load("data/processed/dem_v3_final_30m_eroded.npy")
lithology30 = np.load("data/processed/geomorphology/lithology_v6.npy")

ny, nx = land.shape
cs_x = (XMAX - XMIN) / nx
cs_y = (YMAX - YMIN) / ny
cellsize_km = (cs_x + cs_y) / 2 / 1000.0
log(f"  120m grid: {ny}x{nx}, cellsize ({cs_x:.3f}, {cs_y:.3f}) m")

lake120 = block_any(lake30, 4)[:ny, :nx]
effective_land = land & ~lake120
dem120 = block_mean(dem30, 4)[:ny, :nx]
n_effective_land = int(effective_land.sum())

log("downsampling lithology_v6 (30m, categorical, 8 classes) to 120m via block-mode "
    "majority vote -- NEW block_mode() in terrain_metrics.py, first categorical consumer "
    "of that module...")
n_classes = len(LITHOLOGY_CLASSES)
lithology120 = block_mode(lithology30, 4, n_classes)[:ny, :nx]

friction = travel_friction_multiplier(lithology120)

mapped_codes = set(LAND_TRAVEL_FRICTION.keys())
unmapped_on_land = effective_land & ~np.isin(lithology120, list(mapped_codes))
n_unmapped_on_land = int(unmapped_on_land.sum())
log(f"  {n_unmapped_on_land}/{n_effective_land} effective-land cells "
    f"({100 * n_unmapped_on_land / max(1, n_effective_land):.3f}%) had no LAND_TRAVEL_FRICTION "
    f"class code (land_mask/lithology_v6 coastline don't perfectly agree) -- fell back to "
    f"friction=1.0 (neutral), not an error.")

class_areas_km2 = {
    name: float((lithology120 == code).sum() * cs_x * cs_y / 1e6)
    for code, name in LITHOLOGY_CLASSES.items()
}
log(f"  120m lithology class areas (km2): {class_areas_km2}")

# --- baseline vs friction-adjusted cost graphs --------------------------------------
log("building BASELINE cost graph (no friction -- exact Tappa 6 reproduction)...")
t0 = time.time()
graph_baseline = build_cost_graph(dem120, effective_land, cellsize_km)
graph_baseline_s = time.time() - t0
log(f"  built in {graph_baseline_s:.1f}s")

log("building FRICTION-ADJUSTED cost graph (LAND_TRAVEL_FRICTION applied)...")
t0 = time.time()
graph_friction = build_cost_graph(dem120, effective_land, cellsize_km, friction_multiplier=friction)
graph_friction_s = time.time() - t0
log(f"  built in {graph_friction_s:.1f}s")

cd_baseline = partial(cost_distance_from_source, graph_baseline, shape=(ny, nx))
cd_friction = partial(cost_distance_from_source, graph_friction, shape=(ny, nx))

# --- compare against the 17 ALREADY-PLACED (locked) Circulo sites -------------------
log("loading the 17 already-placed Circulo sites (read-only -- not re-placing them)...")
with open(f"{OUT}/tappa6_site_selection_meta.json") as f:
    tappa6_meta = json.load(f)

sites = []
for s in tappa6_meta["sites"]:
    if not s.get("placed"):
        continue
    col = int(round((s["x_km"] - XMIN / 1000.0) / cellsize_km - 0.5))
    row = int(round((YMAX / 1000.0 - s["y_km"]) / cellsize_km - 0.5))
    sites.append({"name": s["name"], "tier": s["tier"], "row": row, "col": col})
log(f"  {len(sites)} placed sites loaded")

TIER_MIN_HOURS = {
    frozenset({"large"}): 12.0,
    frozenset({"medium"}): 8.0,
    frozenset({"small"}): 2.5,
    frozenset({"large", "medium"}): 6.0,
}

log("recomputing pairwise cost-distance (hours) for every tier-constrained pair, "
    "baseline vs friction-adjusted...")
comparisons = []
worst_margin_change_hours = None
n_would_now_violate = 0
recomputed_baseline_min_margin = None

for i in range(len(sites)):
    for j in range(i + 1, len(sites)):
        a, b = sites[i], sites[j]
        if not (a["tier"] and b["tier"]):
            continue
        hr_req = TIER_MIN_HOURS.get(frozenset({a["tier"], b["tier"]}))
        if hr_req is None:
            continue

        hrs_baseline = float(cd_baseline(a["row"], a["col"])[b["row"], b["col"]])
        hrs_friction = float(cd_friction(a["row"], a["col"])[b["row"], b["col"]])
        margin_baseline = hrs_baseline - hr_req
        margin_friction = hrs_friction - hr_req
        delta_hours = hrs_friction - hrs_baseline
        now_violates = hrs_friction < hr_req

        if recomputed_baseline_min_margin is None or margin_baseline < recomputed_baseline_min_margin:
            recomputed_baseline_min_margin = margin_baseline
        margin_change = margin_friction - margin_baseline
        if worst_margin_change_hours is None or margin_change < worst_margin_change_hours:
            worst_margin_change_hours = margin_change
        if now_violates:
            n_would_now_violate += 1

        comparisons.append({
            "a": a["name"], "b": b["name"],
            "tier_pair": sorted({a["tier"], b["tier"]}),
            "required_hours": hr_req,
            "baseline_hours": round(hrs_baseline, 4),
            "friction_hours": round(hrs_friction, 4),
            "delta_hours": round(delta_hours, 4),
            "delta_pct": round(100 * delta_hours / hrs_baseline, 2) if hrs_baseline > 0 else None,
            "baseline_margin_hours": round(margin_baseline, 4),
            "friction_margin_hours": round(margin_friction, 4),
            "would_now_violate_tier_rule": now_violates,
        })

log(f"  {len(comparisons)} tier-constrained pairs re-checked")
log(f"  sanity check -- recomputed baseline smallest margin: {recomputed_baseline_min_margin:.4f}h "
    f"(Tappa 6's own already-committed record: "
    f"{tappa6_meta['min_distance_verification']['smallest_tier_margin_hours']:.4f}h)")
log(f"  worst margin CHANGE from adding friction: {worst_margin_change_hours:+.4f}h")
log(f"  pairs that would now VIOLATE their tier rule under friction: {n_would_now_violate} "
    f"(0 in the original, already-locked Tappa 6 run)")

# --- export ---------------------------------------------------------------------
log("exporting friction multiplier raster (120m, float32) + comparison meta...")
np.save(f"{OUT}/transport_friction_multiplier_120m.npy", friction.astype(np.float32))
write_envi_raw(
    f"{OUT}/transport_friction_multiplier_120m", friction.astype(np.float32),
    xmin=XMIN, ymin=YMIN, cellsize=(cs_x + cs_y) / 2,
    description=(
        "Tappa 8 transport lithology-cost multiplier: per-cell Tobler-speed friction factor "
        "derived from lithology_v6 (block-mode downsampled to 120m). 1.0 = basin_fill baseline "
        "(no penalty), lower = slower overland travel. See LAND_TRAVEL_FRICTION in "
        "geomorphology/lithology.py for per-class values, citations, and caveats -- UNREVIEWED "
        "first-pass estimates, not locked."
    ),
    dtype="f4",
)
write_prj(f"{OUT}/transport_friction_multiplier_120m.prj", CRS_PROJ4)

meta = {
    "scope_note": (
        "A travel-TIME friction multiplier layered onto Tappa 6's existing Tobler-hiking-"
        "function cost graph (src/suitability/cost_distance.py), NOT a road/rail construction-"
        "cost model (no excavation/foundation/tunnel cost data exists anywhere in this project "
        "to build one from) and NOT a re-run of Tappa 6's actual site placement -- "
        "circulo_candidate_sites.geojson and tappa6_site_selection_meta.json are read-only "
        "inputs here, never modified. See run_tappa8_transport_friction.py's module docstring "
        "and lithology.py's LAND_TRAVEL_FRICTION docstring for the full reasoning."
    ),
    "land_travel_friction_table": {
        LITHOLOGY_CLASSES[code]: mult for code, mult in LAND_TRAVEL_FRICTION.items()
    },
    "friction_status": "UNREVIEWED first-pass estimates -- directional judgement calls grounded "
    "in real-world travel-difficulty accounts per rock type, not a cited per-lithology hiking-"
    "speed dataset (none was found or is claimed to exist). Pending Nico's sign-off; not written "
    "to config/parameters.yml.",
    "boat_speed_kmh_unaffected_by_friction": BOAT_SPEED_KMH,
    "resolution_m": [cs_x, cs_y],
    "class_areas_km2_120m": class_areas_km2,
    "unmapped_on_land_cells": n_unmapped_on_land,
    "unmapped_on_land_of_effective_land_pct": round(100 * n_unmapped_on_land / max(1, n_effective_land), 4),
    "graph_build_seconds": {"baseline": graph_baseline_s, "friction_adjusted": graph_friction_s},
    "sanity_check": {
        "recomputed_baseline_smallest_tier_margin_hours": recomputed_baseline_min_margin,
        "tappa6_committed_smallest_tier_margin_hours":
            tappa6_meta["min_distance_verification"]["smallest_tier_margin_hours"],
        "note": "these two should match closely -- confirms this script's baseline graph is a "
        "faithful reproduction of Tappa 6's own, before friction is layered on.",
    },
    "tier_pair_comparison": comparisons,
    "worst_margin_change_hours": worst_margin_change_hours,
    "n_pairs_that_would_now_violate_tier_rule": n_would_now_violate,
    "final_note": "Nothing about Tappa 6's actual, already-locked Circulo placement changes as a "
    "result of this script. If Nico wants to see where the 17 sites WOULD land if "
    "place_circulos() were re-run WITH this friction from scratch (a materially bigger, "
    "separate decision -- it could move sites, not just change measured margins), that is a "
    "distinct follow-up, not implied by running this script.",
}
with open(f"{OUT}/tappa8_transport_friction_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
