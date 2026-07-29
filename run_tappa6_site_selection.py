"""
Tappa 6 (final deliverable) -- Circulo candidate site selection, greedy by
population, on top of suitability_circulo_120m. See
src/suitability/site_selection.py's module docstring for the full method,
the density assumption, and honest limitations (coarse resolution for the
smallest Circulos, greedy-not-jointly-optimal placement, square-window
approximation of a circular footprint pending Tappa 7's actual layout).

COST-DISTANCE UPGRADE (this version): the earlier straight-line-km tier
minimums (60km large-large / 25km medium-medium) are replaced by
COST-DISTANCE (walking/boat TIME, hours) minimums, applied now to all three
tiers including a new "small" tier -- see src/suitability/cost_distance.py
for the friction model (Tobler's hiking function on land, a flat boat speed
crossing any non-land cell) and site_selection.py's docstring for why two
separate mechanisms (straight-line footprint non-overlap + cost-distance
tier rules) coexist rather than one replacing the other.

LAKE EXCLUSION + LARGE-MEDIUM RULE (this version, round 2): Nico caught two
more real problems by inspecting the previous run's output directly. (1) 7
of the 17 sites (not the 4 he'd spotted by eye -- checked precisely here)
landed ON a lake, 2 of them large Circulos with 53-95% of their own site
window covered by lake water -- lake_mask.npy (Tappa 4 hydrology, native
30m) existed in the repo already but was never wired into site_selection,
since land_mask only excludes ocean, not inland lakes. Fixed by building an
EFFECTIVE land mask (land_mask & ~lake_mask_120m) and using THAT, not the
raw land_mask, for every site-selection computation (window screening, the
cost graph, everything) -- lakes now behave like any other non-land cell
(crossable by boat, not buildable on). (2) The nearest large-medium pair
was only ~3-5km apart in the previous run (cross-tier pairs had NO minimum
at all, falling back to the tiny footprint-buffer spacing) -- added
frozenset({"large","medium"}): 6.0h to TIER_MIN_HOURS, feasibility-tested
the same way as every other threshold in this file (see below).

Reads:
  data/processed/suitability/suitability_circulo_120m.npy
  data/processed/climate/land_mask.npy
  data/processed/hydrology/lake_mask.npy     (native 30m, Tappa 4 -- block-
                                              ANY downsampled to 120m here,
                                              same hazard-conservative
                                              convention as stream_mask)
  data/processed/dem_v3_final_30m_eroded.npy  (native 30m, block-mean
                                               downsampled to 120m here for
                                               the cost graph's elevation
                                               term -- NOT the same
                                               downsample as slope_pct_120m,
                                               which uses block-MAX; see
                                               cost_distance.py's docstring)
  data/processed/biomes/biome_id.npy          (informational only -- which
                                               biome each site lands in,
                                               a preview of the Tappa 7
                                               architecture-style lookup,
                                               NOT a formal Tappa 6 output)

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  circulo_candidate_sites.geojson   (Point features, one per Circulo)
  circulo_claimed_footprint_120m    (bool raster -- union of (a) every
                                     site's straight-line footprint-buffer
                                     circle and (b) for tiered sites, the
                                     actual cost-distance ISOCHRONE below
                                     its tier's hour threshold, which is an
                                     irregular terrain/coast-following
                                     shape, not a circle -- a sanity-check
                                     visual, not a formal layer)
  tappa6_site_selection_meta.json
"""

from __future__ import annotations

import json
import os
import time
from functools import partial

import numpy as np

from src.biomes.world_biomes import BIOME_NAMES
from src.suitability.cost_distance import (
    BOAT_SPEED_KMH,
    build_cost_graph,
    cost_distance_from_source,
)
from src.suitability.site_selection import (
    BUFFER_FACTOR,
    DENSITY_PPL_KM2,
    MIN_LAND_FRACTION,
    place_circulos,
)
from src.suitability.terrain_metrics import block_any, block_mean
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
PROJ4 = (
    "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

# Populations from the Tappa 6 planning chat. The 8 smallest only had a
# combined total (5,000) specified, not individual sizes -- split evenly
# here (625 each), a simplifying assumption, see site_selection.py docstring.
#
# tier (3rd element): ALL THREE tiers now carry a cost-distance minimum
# (see TIER_MIN_HOURS below) -- large/medium were the original "greek
# city-state" idea, small is new, added after Nico pointed out some small
# Circulos were landing under 2km apart ("praticamente integrado a ele").
CIRCULOS = [
    ("Circulo_A_40k", 40000, "large"),
    ("Circulo_B_35k", 35000, "large"),
    ("Circulo_C_25k", 25000, "large"),
    ("Circulo_D_20k", 20000, "large"),
] + [(f"Circulo_E{i+1}_2k", 2000, "medium") for i in range(5)] + [
    (f"Circulo_F{i+1}_small", 625, "small") for i in range(8)
]

# Cost-distance (HOURS) minimum separation, same-tier pairs only -- see
# site_selection.py's docstring for why cross-tier pairs fall back to
# ordinary footprint spacing instead. Feasibility checked directly (same
# method as the earlier km version): each tier tested ALONE (ignoring the
# other tiers) across a range of thresholds using this exact cost graph --
#   large (4 sites):  4/4 placed up to 15h; 16h+ drops to 3/4 -- picked
#                     12h, 3h of margin below the ceiling.
#   medium (5 sites): 5/5 placed up to 12h; 15h+ drops to 4/5 -- picked
#                     8h, 4h+ of margin below the ceiling.
#   small (8 sites):  8/8 placed at every threshold tested up to 8h (never
#                     found a ceiling in the tested range) -- picked 2.5h,
#                     Nico's own recommended midpoint of "2-3 horas de
#                     caminhada", comfortably inside the feasible range.
# All three together (interleaved, the real placement order): still
# 17/17 placed, 0 violations, smallest margin ~0.014h (~49s) -- see
# min_distance_verification in the meta this run writes.
#
# These hour values are a fresh judgement call, not a literal unit
# conversion of the old 60km/25km -- at Tobler's flat-ground pace (~5km/h)
# 12h would be ~60km and 8h ~40km AS IF the whole route were flat land,
# which is roughly why these numbers were chosen (keeps the original
# "day-trip-ish for a large city, less for a medium one" feel) but actual
# achieved cost-distance on this terrain is NOT simply hours*5km/h --
# slopes slow land routes down and BOAT_SPEED_KMH (6km/h, faster than
# Tobler's flat pace) can make a sea-crossing route cheaper per km than an
# equivalent land route, so the real relationship between "hours" and "km
# apart on the map" varies a lot by where the two sites actually are.
#
# large-medium (NEW, round 2): Nico noticed a medium Circulo only ~4km
# from a large one -- cross-tier pairs previously had NO minimum at all
# beyond ordinary footprint spacing. Feasibility-tested the same way,
# together with the other 3 tier rules and the lake exclusion below (not
# in isolation, since by this point the other rules already constrain
# the map a lot): 17/17 still placed up to 10h; 11h+ drops to 13/17.
# Picked 6h -- comfortably below the ceiling, and intentionally between
# small-small (2.5h) and medium-medium (8h): a medium town should feel
# more independent from a nearby metropolis than from another village,
# but doesn't need the same full separation two mediums need from each
# other. Result: nearest large-medium pair went from ~3-5km / <1h to
# 29km / 6.19h -- see min_distance_verification for the exact figures
# this run measured.
TIER_MIN_HOURS = {
    frozenset({"large"}): 12.0,
    frozenset({"medium"}): 8.0,
    frozenset({"small"}): 2.5,
    frozenset({"large", "medium"}): 6.0,
}


def main():
    t0 = time.time()
    suit = np.load("data/processed/suitability/suitability_circulo_120m.npy").astype(np.float64)
    land = np.load("data/processed/climate/land_mask.npy").astype(bool)
    lake30 = np.load("data/processed/hydrology/lake_mask.npy")
    biome_id = np.load("data/processed/biomes/biome_id.npy")
    dem30 = np.load("data/processed/dem_v3_final_30m_eroded.npy")

    ny, nx = land.shape
    cs_x = (XMAX - XMIN) / nx
    cs_y = (YMAX - YMIN) / ny
    cellsize_km = (cs_x + cs_y) / 2 / 1000.0

    # lake_mask: block-ANY, not block-mean -- same hazard-conservative
    # reasoning terrain_metrics.py already uses for stream_mask/slope (a
    # thin arm of a lake shouldn't disappear inside an otherwise-dry 120m
    # cell). effective_land is what EVERY site-selection computation below
    # uses in place of the raw land_mask -- lakes are on-land elevation-wise
    # but not buildable, so they must not count as valid site area, valid
    # window land-fraction, OR valid on-foot terrain in the cost graph
    # (a lake still gets crossed at BOAT_SPEED_KMH, it just can't host a
    # Circulo).
    lake120 = block_any(lake30, 4)[:ny, :nx]
    effective_land = land & ~lake120

    # 120m DEM for the cost graph's elevation/slope term -- block-MEAN, not
    # block-max (see cost_distance.py's docstring: we want a realistic
    # average edge slope between neighbouring 120m cells here, not a
    # hazard-conservative worst case within a cell).
    dem120 = block_mean(dem30, 4)[:ny, :nx]

    t_graph = time.time()
    cost_graph = build_cost_graph(dem120, effective_land, cellsize_km)
    cost_distance_fn = partial(cost_distance_from_source, cost_graph, shape=(ny, nx))
    graph_build_s = time.time() - t_graph

    results = place_circulos(
        suit, effective_land, cellsize_km, CIRCULOS, TIER_MIN_HOURS,
        cost_graph=cost_graph, cost_distance_fn=cost_distance_fn,
        xmin_km=XMIN / 1000.0, ymax_km=YMAX / 1000.0,
    )

    features = []
    claimed_viz = np.zeros((ny, nx), dtype=bool)
    for r in results:
        if not r["placed"]:
            continue
        row, col = r["row"], r["col"]
        x, y = r["x_km"] * 1000.0, r["y_km"] * 1000.0  # GeoJSON in project meters
        biome_name = BIOME_NAMES[int(biome_id[row, col])]
        features.append({
            "type": "Feature",
            "properties": {
                "name": r["name"],
                "population": r["population"],
                "tier": r["tier"],
                "radius_km": round(r["radius_km"], 3),
                "window_cells": r["window_cells"],
                "mean_suitability": round(r["mean_suitability"], 4),
                "land_fraction": round(r["land_fraction"], 4),
                "biome_at_site": biome_name,
            },
            "geometry": {"type": "Point", "coordinates": [x, y]},
        })
        # base visualization: every site's ordinary footprint-buffer circle
        # (straight-line, mechanism 1 -- physical non-overlap)
        viz_radius_km = r["radius_km"] * BUFFER_FACTOR
        rc = int(np.ceil(viz_radius_km / cellsize_km))
        r0, r1 = max(0, row - rc), min(ny, row + rc + 1)
        c0, c1 = max(0, col - rc), min(nx, col + rc + 1)
        yy, xx = np.ogrid[r0:r1, c0:c1]
        dist2 = (yy - row) ** 2 + (xx - col) ** 2
        claimed_viz[r0:r1, c0:c1] |= dist2 <= (viz_radius_km / cellsize_km) ** 2

        # NOTE, tried and reverted: shading each tiered site's actual
        # cost-distance isochrone (hours_from_site < tier_thresh) into this
        # raster too, on top of the footprint circle. Checked directly: a
        # SINGLE large Circulo's 12h isochrone alone already covers ~24% of
        # all land, and with 17 sites of overlapping isochrones (4 large @
        # 12h + 5 medium @ 8h + 8 small @ 2.5h) the UNION covers 99.998% of
        # the entire grid -- i.e. this visualization would be almost solid,
        # not a useful "what's still open" picture. This does NOT mean the
        # tier rule is toothless: enforcement is a per-PAIR lookup (this
        # candidate vs. that specific already-placed site), verified
        # correct below (0 violations) -- it's the "shade every isochrone
        # ever used" idea for the visual that doesn't work at this hour
        # scale on this map, not the underlying constraint.

    geojson = {
        "type": "FeatureCollection",
        "name": "circulo_candidate_sites",
        "crs": {"type": "proj4", "properties": {"proj4": PROJ4}},
        "features": features,
    }

    out = "data/processed/suitability"
    os.makedirs(out, exist_ok=True)
    with open(f"{out}/circulo_candidate_sites.geojson", "w") as f:
        json.dump(geojson, f, indent=2)

    claimed_i2 = claimed_viz.astype(np.int16)
    np.save(f"{out}/circulo_claimed_footprint_120m.npy", claimed_i2)
    write_envi_raw(
        f"{out}/circulo_claimed_footprint_120m", claimed_i2, XMIN, YMIN, cs_x,
        "Tappa6 circulo_claimed_footprint_120m", dtype="i2",
    )
    write_prj(f"{out}/circulo_claimed_footprint_120m.prj", PROJ4)

    n_placed = sum(1 for r in results if r["placed"])
    n_failed = len(results) - n_placed

    # verify every pairwise requirement actually holds in the OUTPUT, not
    # just trust the placement logic: straight-line footprint non-overlap
    # for every pair, PLUS cost-distance (hours) for ANY pair whose tier
    # combination (same-tier OR the large-medium cross rule) is a key in
    # TIER_MIN_HOURS, PLUS a direct lake check on every placed site.
    placed = [r for r in results if r["placed"]]
    min_gap_km_found = None
    min_margin_hours_found = None
    violations = []
    for i in range(len(placed)):
        for j in range(i + 1, len(placed)):
            a, b = placed[i], placed[j]
            d_km = float(np.hypot(a["x_km"] - b["x_km"], a["y_km"] - b["y_km"]))
            req_km = max(a["radius_km"] * BUFFER_FACTOR, b["radius_km"] * BUFFER_FACTOR)
            if min_gap_km_found is None or (d_km - req_km) < min_gap_km_found:
                min_gap_km_found = d_km - req_km
            if d_km < req_km - 1e-6:
                violations.append({
                    "type": "footprint_km", "a": a["name"], "b": b["name"],
                    "dist_km": d_km, "required_km": req_km,
                })

            if a["tier"] and b["tier"]:
                hr_req = TIER_MIN_HOURS.get(frozenset({a["tier"], b["tier"]}))
                if hr_req is not None:
                    dist_hrs = float(
                        cost_distance_fn(a["row"], a["col"])[b["row"], b["col"]]
                    )
                    margin = dist_hrs - hr_req
                    if min_margin_hours_found is None or margin < min_margin_hours_found:
                        min_margin_hours_found = margin
                    if dist_hrs < hr_req - 1e-6:
                        violations.append({
                            "type": "tier_hours", "a": a["name"], "b": b["name"],
                            "dist_hours": dist_hrs, "required_hours": hr_req,
                        })

    n_on_lake = 0
    for r in placed:
        if bool(lake120[r["row"], r["col"]]):
            n_on_lake += 1
            violations.append({"type": "on_lake", "name": r["name"]})

    meta = {
        "assumptions": {
            "density_ppl_km2": DENSITY_PPL_KM2,
            "buffer_factor": BUFFER_FACTOR,
            "min_land_fraction": MIN_LAND_FRACTION,
            "tier_min_hours": {
                "large-large": 12.0, "medium-medium": 8.0, "small-small": 2.5,
                "large-medium": 6.0,
            },
            "boat_speed_kmh": BOAT_SPEED_KMH,
            "cost_model": "Tobler's hiking function on land (signed slope from a "
            "120m block-MEAN DEM), flat BOAT_SPEED_KMH on any edge touching a "
            "non-land cell -- see cost_distance.py.",
            "lake_exclusion": "land_mask & ~lake_mask_120m (lake_mask.npy, Tappa 4 "
            "hydrology, block-ANY downsampled) used for ALL site-selection "
            "computations in place of the raw land_mask -- added after Nico "
            "found sites landing inside lakes in the previous run (checked "
            "directly here: 7 of 17, not the 4 he'd spotted visually, 2 of "
            "them large Circulos with 53-95% lake coverage in their own "
            "window). 0 of 17 land on a lake in this run -- see "
            "min_distance_verification for the direct re-check.",
            "note": "8 smallest Circulos' individual populations were not specified "
            "(only their 5,000 combined total) -- split evenly (625 each) here. "
            "Cost-distance (hours) replaces the earlier straight-line km tier "
            "minimums, applied now to all 3 tiers (large/medium/small) plus a "
            "new large-medium cross-tier rule (6h, after Nico noticed a medium "
            "Circulo only ~4km from a large one in the previous run) -- Nico "
            "asked to test whether this changes the same-biome clustering "
            "outcome, and to give the small Circulos (previously untiered) "
            "their own minimum, expressed in walking hours since that's the "
            "natural unit for 'not practically integrated'.",
        },
        "resolution_m": [cs_x, cs_y],
        "cost_graph_build_seconds": graph_build_s,
        "n_circulos": len(CIRCULOS),
        "n_placed": n_placed,
        "n_failed_to_place": n_failed,
        "min_distance_verification": {
            "smallest_footprint_margin_km": min_gap_km_found,
            "smallest_tier_margin_hours": min_margin_hours_found,
            "n_sites_on_lake": n_on_lake,
            "n_violations": len(violations),
            "violations": violations,
        },
        "sites": [
            {k: v for k, v in r.items() if k not in ("row", "col")} for r in results
        ],
    }
    with open(f"{out}/tappa6_site_selection_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time() - t0:.1f}s -- {n_placed}/{len(CIRCULOS)} placed, {n_failed} failed")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
