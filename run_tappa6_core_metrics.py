"""
Tappa 6 (step 1 of 2) -- core suitability layers + the Povo Silencioso
exclusion distance, at the 120 m working grid.

Reads:
  data/processed/dem_v3_final_30m_eroded.npy          (Tappa 1, native 30 m)
  data/processed/hydrology/contributing_area_km2.npy  (Tappa 4, native 30 m)
  data/processed/hydrology/stream_mask.npy            (Tappa 4, native 30 m)
  data/processed/climate/land_mask.npy                (Tappa 2, 120 m)

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  slope_pct_120m / slope_suitability_120m
  dist_to_stream_km_120m / water_suitability_120m
  twi_120m / agriculture_suitability_120m
  dist_to_povo_silencioso_km_120m / povo_silencioso_exclusion_120m
  tappa6_core_metrics_meta.json

Each .npy also ships as ENVI .bin/.hdr/.prj (see src/terrain/raster_io.py),
same convention as every prior stage. Solar exposure and the biome
suitability lookup are NOT computed here -- see the Tappa 6 decision doc's
open-items list.

water_suitability_120m and povo_silencioso_exclusion_120m were added after
the initial layers 1-4 build, once it was clear the raw distance fields
(dist_to_stream_km, dist_to_povo_silencioso_km) needed a 0-1 transform
before they could enter a weighted composite alongside slope_suitability/
agriculture_suitability/solar_suitability/biome_suitability. See
water_suitability's and povo_silencioso_exclusion_factor's docstrings in
src/suitability/terrain_metrics.py for the honest caveats on both --
water_suitability in particular reads as near-flat over most of the land
(this world's stream network is dense), by design, not a bug.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
from scipy import ndimage

from src.suitability.terrain_metrics import (
    agriculture_suitability,
    block_any,
    block_max,
    block_mean,
    compute_slope_pct,
    distance_to_stream_km,
    povo_silencioso_distance_km,
    povo_silencioso_exclusion_factor,
    slope_suitability,
    topographic_wetness_index,
    water_suitability,
)
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
FACTOR = 4  # native 30 m -> working 120 m
PROJ4 = (
    "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)
# Identified once, by hand, in the Tappa 6 planning chat -- see
# povo_silencioso_distance_km's docstring for why this is a fixed list, not a
# geometric re-derivation.
POVO_SILENCIOSO_LABELS = [10, 7, 6, 47, 20]

# Placeholder, uncalibrated slope-suitability knee (see slope_suitability's
# docstring) -- revisit once the settlement-size downstream filter exists.
SLOPE_GENTLE_PCT = 5.0
SLOPE_HARD_LIMIT_PCT = 30.0

# Placeholder water-suitability knee, set from this world's own dist-to-
# stream percentiles (0.5km ~ where "trivially close" stops covering nearly
# everyone; 5km ~ p99.5) -- see water_suitability's docstring.
WATER_GENTLE_KM = 0.5
WATER_HARD_LIMIT_KM = 5.0

# Placeholder Povo Silencioso "respect the territory" exclusion buffer --
# not derived from any specified treaty/lore distance -- see
# povo_silencioso_exclusion_factor's docstring.
PS_HARD_BUFFER_KM = 5.0
PS_SOFT_BUFFER_KM = 15.0


def main():
    t0 = time.time()
    dem30 = np.load("data/processed/dem_v3_final_30m_eroded.npy")
    contrib30 = np.load("data/processed/hydrology/contributing_area_km2.npy")
    stream30 = np.load("data/processed/hydrology/stream_mask.npy")
    land120 = np.load("data/processed/climate/land_mask.npy").astype(bool)

    ny120, nx120 = land120.shape
    cs_x120 = (XMAX - XMIN) / nx120
    cs_y120 = (YMAX - YMIN) / ny120

    # --- slope ---
    slope_pct_30 = compute_slope_pct(dem30, cellsize_m=30.0)
    slope_pct_120 = block_max(slope_pct_30, FACTOR)  # hazard-safe: max, not mean
    slope_pct_120_mean_check = block_mean(slope_pct_30, FACTOR)
    slope_suit = slope_suitability(slope_pct_120, SLOPE_GENTLE_PCT, SLOPE_HARD_LIMIT_PCT)

    # --- distance to stream ---
    stream_120 = block_any(stream30, FACTOR)
    dist_stream_km = distance_to_stream_km(stream_120, (cs_y120, cs_x120), factor=1)
    water_suit = water_suitability(dist_stream_km, WATER_GENTLE_KM, WATER_HARD_LIMIT_KM)

    # --- TWI / agriculture proxy ---
    contrib_120 = block_mean(contrib30, FACTOR)
    twi = topographic_wetness_index(contrib_120, slope_pct_120, cellsize_m=cs_x120)
    twi = np.where(land120 & np.isfinite(twi), twi, np.nan)
    agri_suit = agriculture_suitability(twi, land120)

    # --- Povo Silencioso exclusion ---
    labeled, _ = ndimage.label(land120, structure=np.ones((3, 3)))
    archipelago_mask = np.isin(labeled, POVO_SILENCIOSO_LABELS)
    other_land = land120 & ~archipelago_mask
    dist_ps_km = povo_silencioso_distance_km(land120, POVO_SILENCIOSO_LABELS, (cs_y120, cs_x120))
    ps_exclusion = povo_silencioso_exclusion_factor(dist_ps_km, PS_HARD_BUFFER_KM, PS_SOFT_BUFFER_KM)

    out = "data/processed/suitability"
    os.makedirs(out, exist_ok=True)

    layers = {
        "slope_pct_120m": slope_pct_120,
        "slope_suitability_120m": slope_suit,
        "dist_to_stream_km_120m": dist_stream_km,
        "water_suitability_120m": water_suit,
        "twi_120m": twi,
        "agriculture_suitability_120m": agri_suit,
        "dist_to_povo_silencioso_km_120m": dist_ps_km,
        "povo_silencioso_exclusion_120m": ps_exclusion,
    }
    for name, arr in layers.items():
        arr32 = arr.astype(np.float32)
        np.save(f"{out}/{name}.npy", arr32)
        write_envi_raw(f"{out}/{name}", arr32, XMIN, YMIN, cs_x120, f"Tappa6 {name}", dtype="f4")
        write_prj(f"{out}/{name}.prj", PROJ4)

    land_finite_twi = land120 & np.isfinite(twi)
    meta = {
        "resolution_m": [cs_x120, cs_y120],
        "grid_shape": list(land120.shape),
        "slope": {
            "method": "np.gradient at native 30m; block-MAX downsample to 120m",
            "gentle_pct_threshold": SLOPE_GENTLE_PCT,
            "hard_limit_pct_threshold": SLOPE_HARD_LIMIT_PCT,
            "land_mean_slope_pct_blockmax": float(slope_pct_120[land120].mean()),
            "land_mean_slope_pct_blockmean_forcomparison": float(
                slope_pct_120_mean_check[land120].mean()
            ),
            "land_slope_pct_p50_p95_p99": [
                float(v) for v in np.percentile(slope_pct_120[land120], [50, 95, 99])
            ],
        },
        "distance_to_stream_km": {
            "land_mean": float(dist_stream_km[land120].mean()),
            "land_median": float(np.median(dist_stream_km[land120])),
            "land_max": float(dist_stream_km[land120].max()),
        },
        "water_suitability": {
            "gentle_km_threshold": WATER_GENTLE_KM,
            "hard_limit_km_threshold": WATER_HARD_LIMIT_KM,
            "land_frac_within_0.5km": float((dist_stream_km[land120] <= 0.5).mean()),
            "land_frac_within_1.0km": float((dist_stream_km[land120] <= 1.0).mean()),
            "land_frac_within_5.0km": float((dist_stream_km[land120] <= 5.0).mean()),
            "land_mean_score": float(water_suit[land120].mean()),
            "land_p5_p50_p95_score": [
                float(v) for v in np.percentile(water_suit[land120], [5, 50, 95])
            ],
            "caveat": "near-flat over most of the land by design -- dense stream network, see water_suitability's docstring",
        },
        "twi_agriculture_proxy": {
            "land_p5_p50_p95": [float(v) for v in np.percentile(twi[land_finite_twi], [5, 50, 95])],
            "caveat": "geomorphic proxy only -- no real soil/pedology data exists in this project",
        },
        "povo_silencioso_exclusion": {
            "archipelago_labels": POVO_SILENCIOSO_LABELS,
            "archipelago_area_km2": float(archipelago_mask.sum() * (cs_x120 / 1000) * (cs_y120 / 1000)),
            "land_dist_km_mean_excl_target": float(dist_ps_km[other_land].mean()),
            "land_dist_km_median_excl_target": float(np.median(dist_ps_km[other_land])),
            "hard_buffer_km": PS_HARD_BUFFER_KM,
            "soft_buffer_km": PS_SOFT_BUFFER_KM,
            "other_land_frac_within_10km": float((dist_ps_km[other_land] <= 10.0).mean()),
            "other_land_area_km2_within_10km": float(
                (dist_ps_km[other_land] <= 10.0).sum() * (cs_x120 / 1000) * (cs_y120 / 1000)
            ),
            "caveat": "localised corner-of-the-map effect by construction, only ~1% of non-archipelago land falls within 10km of the archipelago at all",
        },
    }
    with open(f"{out}/tappa6_core_metrics_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time() - t0:.1f}s")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
