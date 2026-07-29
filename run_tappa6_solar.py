"""
Tappa 6 (step 2 of 2) -- solar exposure / insolation proxy at the 120 m
working grid. See src/suitability/solar.py's module docstring for the full
method and honest limitations.

Reads:
  data/processed/climate/surface_elevation_m.npy  (Tappa 2, 120 m)
  data/processed/climate/land_mask.npy            (Tappa 2, 120 m)
  config/parameters.yml                            (CRS / LCC parameters)

Writes to data/processed/suitability/ (gitignored, regenerate locally --
the horizon-angle precompute is the expensive step, ~65 s for 16 directions
x 10 km reach on this grid; the 12-month integration is another ~110 s):
  slope_deg_120m / aspect_deg_120m         (energy-averaging slope/aspect --
                                             NOT the same as the hazard
                                             slope_pct_120m from step 1,
                                             which deliberately used block-
                                             MAX; see module docstring)
  horizon_mean_deg_120m                    (diagnostic: mean horizon angle
                                             across the 16 directions)
  monthly_insolation_MJm2_120m              (12, ny, nx)
  annual_insolation_MJm2_120m
  june_insolation_MJm2_120m                 (southern-winter worst-month
                                             diagnostic, relevant given
                                             this world's reliance on solar
                                             power)
  solar_suitability_annual_120m / solar_suitability_june_120m  (0-1)
  tappa6_solar_meta.json
"""

from __future__ import annotations

import json
import time

import numpy as np

from src.climate.grid import latitude_grid
from src.suitability.solar import (
    MONTH_NAMES,
    annual_insolation_MJ_m2,
    horizon_angles_deg,
    monthly_insolation_MJ_m2_day,
    normalize_suitability,
)
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
LAT_0, LAT_1, LAT_2, LON_0 = -44.0, -44.48, -43.52, 42.0
PROJ4 = (
    f"+proj=lcc +lat_1={LAT_1} +lat_2={LAT_2} +lat_0={LAT_0} +lon_0={LON_0} "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)
HORIZON_N_DIRS = 16
HORIZON_MAX_DIST_M = 10000.0


def slope_aspect_deg(elevation_m, cs_y, cs_x):
    """Ordinary (non-block-max) gradient -- appropriate here because
    insolation is an energy-AVERAGING quantity, unlike the step-1 hazard
    slope, which deliberately used block-max to avoid hiding locally steep
    ground. Two different quantities from the same DEM family, by design,
    not an inconsistency."""
    gy, gx = np.gradient(elevation_m, cs_y, cs_x)
    slope_deg = np.degrees(np.arctan(np.sqrt(gx**2 + gy**2)))
    dzdN, dzdE = -gy, gx  # row0=North convention -> d/d(north) = -d/d(row)
    aspect_deg = np.degrees(np.arctan2(-dzdE, -dzdN)) % 360  # downslope-facing azimuth, 0=N cw
    return slope_deg, aspect_deg


def main():
    t0 = time.time()
    elev = np.load("data/processed/climate/surface_elevation_m.npy").astype(np.float64)
    land = np.load("data/processed/climate/land_mask.npy").astype(bool)
    ny, nx = elev.shape
    cs_x = (XMAX - XMIN) / nx
    cs_y = (YMAX - YMIN) / ny
    cellsize = (cs_x + cs_y) / 2

    lat_grid = latitude_grid(ny, nx, YMAX, XMIN, cellsize, lat_1=LAT_1, lat_2=LAT_2, lat_0=LAT_0)
    slope_deg, aspect_deg = slope_aspect_deg(elev, cs_y, cs_x)

    print(f"[{time.time()-t0:.0f}s] computing horizon angles ({HORIZON_N_DIRS} dirs, {HORIZON_MAX_DIST_M/1000:.0f} km)...")
    horizon_deg, horizon_az = horizon_angles_deg(elev, cellsize, HORIZON_N_DIRS, HORIZON_MAX_DIST_M)

    print(f"[{time.time()-t0:.0f}s] integrating 12 representative months...")
    monthly = monthly_insolation_MJ_m2_day(lat_grid, elev, slope_deg, aspect_deg, horizon_deg, horizon_az)
    annual = annual_insolation_MJ_m2(monthly).astype(np.float32)
    june = monthly[5]

    solar_suit_annual, norm_annual = normalize_suitability(annual, land)
    solar_suit_june, norm_june = normalize_suitability(june, land)

    out = "data/processed/suitability"
    import os

    os.makedirs(out, exist_ok=True)
    layers = {
        "slope_deg_120m": slope_deg,
        "aspect_deg_120m": aspect_deg,
        "horizon_mean_deg_120m": horizon_deg.mean(axis=0),
        "annual_insolation_MJm2_120m": annual,
        "june_insolation_MJm2_120m": june,
        "solar_suitability_annual_120m": solar_suit_annual,
        "solar_suitability_june_120m": solar_suit_june,
    }
    for name, arr in layers.items():
        arr32 = arr.astype(np.float32)
        np.save(f"{out}/{name}.npy", arr32)
        write_envi_raw(f"{out}/{name}", arr32, XMIN, YMIN, cs_x, f"Tappa6 {name}", dtype="f4")
        write_prj(f"{out}/{name}.prj", PROJ4)
    np.save(f"{out}/monthly_insolation_MJm2_120m.npy", monthly)

    # aspect validation on the real terrain: north vs south facing at matched
    # moderate slope, should show a clear N > S asymmetry at this latitude
    moderate = land & (slope_deg > 8) & (slope_deg < 25)
    north_facing = moderate & (np.abs(((aspect_deg + 180) % 360) - 180) < 30)
    south_facing = moderate & (np.abs(((aspect_deg - 180 + 180) % 360) - 180) < 30)

    meta = {
        "resolution_m": [cs_x, cs_y],
        "method_summary": "FAO-56 elevation-corrected clear-sky Rso, scaled by a "
        "numerically-integrated tilted/shaded-vs-horizontal ratio (16-direction "
        "10km horizon shading, 0.5hr steps, 12 mid-month representative days); "
        "diffuse = 15% clear-sky fraction x slope+terrain sky-view-factor. See "
        "src/suitability/solar.py docstring for full method and limitations.",
        "monthly_land_mean_MJm2_day": {
            MONTH_NAMES[i]: float(monthly[i][land].mean()) for i in range(12)
        },
        "annual_land_mean_MJm2": float(annual[land].mean()),
        "annual_land_p5_p50_p95_MJm2": [float(v) for v in np.percentile(annual[land], [5, 50, 95])],
        "horizon": {
            "n_dirs": HORIZON_N_DIRS,
            "max_dist_m": HORIZON_MAX_DIST_M,
            "land_mean_horizon_deg": float(horizon_deg[:, land].mean()),
            "land_max_horizon_deg": float(horizon_deg[:, land].max()),
        },
        "aspect_validation": {
            "description": "moderate-slope (8-25 deg) cells, N-facing vs S-facing, "
            "annual insolation -- should show N > S at this southern-hemisphere "
            "latitude (opposite of northern-hemisphere intuition)",
            "north_facing_mean_MJm2": float(annual[north_facing].mean()),
            "south_facing_mean_MJm2": float(annual[south_facing].mean()),
            "ratio_N_over_S": float(annual[north_facing].mean() / annual[south_facing].mean()),
        },
        "normalization_pctl_2_98": {"annual": norm_annual, "june": norm_june},
    }
    with open(f"{out}/tappa6_solar_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time()-t0:.1f}s")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
