#!/usr/bin/env python3
"""
Tappa 4 driver — hydrology: depression-filled DEM, D8 flow direction,
(uniform and precipitation-weighted) flow accumulation, a stream network by
contributing-area threshold, a lake/depression mask, and major drainage
basins. Runs at the DEM's native 30 m resolution (5334x4334), not the 120 m
climate working grid -- see docs/decisions/04_tappa4_hydrology.md for why.

    python run_tappa4_hydrology.py [--dem PATH] [--climate-dir DIR] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from hydrology.flow import (                                   # noqa: E402
    accumulate_flow,
    direction_codes,
    label_basins,
    priority_flood_d8,
)
from hydrology.weighting import upsample_precip_to_dem          # noqa: E402
from terrain.raster_io import write_envi_raw, write_prj         # noqa: E402

XMIN, XMAX, YMIN, YMAX = -65000.0, 65000.0, -80000.0, 80000.0
RES_M = 30.0
CELL_KM2 = (RES_M / 1000.0) ** 2
PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

STREAM_THRESHOLD_KM2 = 0.3        # see decision doc S3 for the calibration
LAKE_FILL_THRESHOLD_M = 2.0       # depression-fill raise > this = mapped as a lake
MAJOR_BASIN_MIN_KM2 = 20.0        # basins smaller than this lumped as "minor/direct coastal"
SECONDS_PER_YEAR = 365.25 * 86400.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", default="data/processed/dem_v3_final_30m_eroded.npy")
    ap.add_argument("--climate-dir", default="data/processed/climate")
    ap.add_argument("--out", default="data/processed/hydrology")
    args = ap.parse_args()

    t0 = time.time()
    dem = np.load(args.dem).astype(np.float64)
    ny, nx = dem.shape
    land = dem > 0.0
    print(f"[{time.time()-t0:6.1f}s] loaded DEM {ny}x{nx}, land fraction {land.mean():.4f}")

    # --- 1. depression fill + D8 direction (single priority-flood pass) --
    seed = dem <= 0.0
    seed[0, :] = seed[-1, :] = seed[:, 0] = seed[:, -1] = True
    filled, receiver, pop_order = priority_flood_d8(dem, seed, epsilon=1e-4)
    print(f"[{time.time()-t0:6.1f}s] priority-flood fill + D8 direction done")
    assert np.all(filled >= dem - 1e-9)
    assert not np.any(np.isnan(filled))

    codes = direction_codes(receiver, (ny, nx))
    print(f"[{time.time()-t0:6.1f}s] direction codes assigned")

    # --- 2. lake / depression mask ----------------------------------------
    raise_m = (filled - dem).reshape(ny, nx)
    lake_mask = (raise_m > LAKE_FILL_THRESHOLD_M) & land
    print(f"[{time.time()-t0:6.1f}s] lake mask: {lake_mask.sum()} cells "
          f"({lake_mask.sum()*CELL_KM2:.2f} km2), max fill raise {raise_m.max():.1f} m")

    # --- 3. flow accumulation: uniform (contributing area) ----------------
    w_area = land.astype(np.float64)   # ocean contributes no land runoff
    accum_cells = accumulate_flow(receiver, pop_order, w_area).reshape(ny, nx)
    area_km2 = accum_cells * CELL_KM2
    print(f"[{time.time()-t0:6.1f}s] uniform flow accumulation done, "
          f"max contributing area {area_km2.max():.1f} km2")

    # --- 4. flow accumulation: precipitation-weighted ----------------------
    precip_monthly = np.load(Path(args.climate_dir) / "precipitation_monthly_mm.npy")
    precip_annual = precip_monthly.sum(axis=0)
    precip_30 = upsample_precip_to_dem(precip_annual, (ny, nx))
    w_precip = np.where(land, precip_30, 0.0)
    accum_precip_mm = accumulate_flow(receiver, pop_order, w_precip).reshape(ny, nx)
    # discharge proxy: treat accumulated (precip_mm * cell) as if it all
    # became streamflow with no losses -- an explicit UPPER BOUND (no
    # infiltration/ET/routing lag), not a real discharge estimate. See
    # decision doc S5 for why this is stated as a proxy, not a discharge model.
    discharge_proxy_m3s = accum_precip_mm * (RES_M ** 2 / 1000.0) / SECONDS_PER_YEAR
    print(f"[{time.time()-t0:6.1f}s] precipitation-weighted accumulation done, "
          f"max discharge proxy {discharge_proxy_m3s.max():.1f} m3/s")

    # --- 5. stream network ---------------------------------------------
    stream_mask = (area_km2 >= STREAM_THRESHOLD_KM2) & land
    stream_length_km = stream_mask.sum() * (RES_M / 1000.0)
    land_area_km2 = land.sum() * CELL_KM2
    drainage_density = stream_length_km / land_area_km2
    print(f"[{time.time()-t0:6.1f}s] stream network: {stream_mask.sum()} cells, "
          f"~{stream_length_km:.0f} km, drainage density ~{drainage_density:.2f} km/km2")

    # --- 6. drainage basins ----------------------------------------------
    basin_raw = label_basins(receiver, pop_order).reshape(ny, nx)
    basin_ids, basin_counts = np.unique(basin_raw[land], return_counts=True)
    basin_area_km2 = basin_counts * CELL_KM2
    major = basin_ids[basin_area_km2 >= MAJOR_BASIN_MIN_KM2]
    # renumber major basins 1..N by descending size for a stable, compact
    # raster; everything else (countless trivial direct-coastal micro-
    # catchments) collapses to 0 ("minor/direct coastal drainage").
    order_desc = major[np.argsort(-basin_area_km2[np.isin(basin_ids, major)])]
    remap = {int(old): i + 1 for i, old in enumerate(order_desc)}
    basin_labeled = np.zeros((ny, nx), dtype=np.int32)
    mask_major = np.isin(basin_raw, order_desc) & land
    # vectorized remap via searchsorted on the sorted major-id array
    sorter = np.argsort(order_desc)
    idx_in_major = sorter[np.searchsorted(order_desc, basin_raw[mask_major], sorter=sorter)]
    basin_labeled[mask_major] = idx_in_major + 1
    n_major = len(order_desc)
    print(f"[{time.time()-t0:6.1f}s] basins: {len(basin_ids)} total, "
          f"{n_major} major (>= {MAJOR_BASIN_MIN_KM2} km2), "
          f"{(basin_area_km2[np.isin(basin_ids, major)]).sum():.0f} km2 of "
          f"{land_area_km2:.0f} km2 land in major basins")

    # --- 7. windward/leeward river-size asymmetry check (Tappa 2 S1 wind) -
    dist_to_coast = ndimage.distance_transform_edt(land) * RES_M / 1000.0  # crude proxy, see note
    p33, p67 = np.percentile(precip_30[land], [33.33, 66.67])
    wet = land & (precip_30 >= p67)
    dry = land & (precip_30 <= p33)
    wet_stream = stream_mask & wet
    dry_stream = stream_mask & dry
    asym = {
        "wet_side_mean_discharge_proxy_m3s": float(discharge_proxy_m3s[wet_stream].mean()) if wet_stream.any() else None,
        "dry_side_mean_discharge_proxy_m3s": float(discharge_proxy_m3s[dry_stream].mean()) if dry_stream.any() else None,
        "wet_side_max_discharge_proxy_m3s": float(discharge_proxy_m3s[wet_stream].max()) if wet_stream.any() else None,
        "dry_side_max_discharge_proxy_m3s": float(discharge_proxy_m3s[dry_stream].max()) if dry_stream.any() else None,
    }
    print(f"[{time.time()-t0:6.1f}s] windward/leeward discharge asymmetry: {asym}")

    # --- 8. write outputs ----------------------------------------------
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ymin = YMAX - ny * RES_M

    np.save(out / "filled_dem_30m.npy", filled.astype(np.float32))
    np.save(out / "flow_direction_code.npy", codes)
    np.save(out / "contributing_area_km2.npy", area_km2.astype(np.float32))
    np.save(out / "discharge_proxy_m3s.npy", discharge_proxy_m3s.astype(np.float32))
    np.save(out / "stream_mask.npy", stream_mask)
    np.save(out / "lake_mask.npy", lake_mask)
    np.save(out / "basin_labeled.npy", basin_labeled)

    write_envi_raw(str(out / "contributing_area_km2"), area_km2.astype(np.float32),
                    XMIN, ymin, RES_M, "Tappa 4 D8 contributing area, float32 km2", dtype="f4")
    write_envi_raw(str(out / "discharge_proxy_m3s"), discharge_proxy_m3s.astype(np.float32),
                    XMIN, ymin, RES_M, "Tappa 4 mean annual discharge proxy (no infiltration/ET/routing), float32 m3/s", dtype="f4")
    write_envi_raw(str(out / "stream_mask"), stream_mask.astype(np.int16),
                    XMIN, ymin, RES_M, f"Tappa 4 stream network, area>={STREAM_THRESHOLD_KM2} km2, int16 (0/1)", dtype="i2")
    write_envi_raw(str(out / "lake_mask"), lake_mask.astype(np.int16),
                    XMIN, ymin, RES_M, f"Tappa 4 lake/depression mask, fill raise>{LAKE_FILL_THRESHOLD_M} m, int16 (0/1)", dtype="i2")
    write_envi_raw(str(out / "basin_labeled"), basin_labeled,
                    XMIN, ymin, RES_M, f"Tappa 4 major drainage basins (>= {MAJOR_BASIN_MIN_KM2} km2, 0=minor/coastal), int16", dtype="i2")
    write_envi_raw(str(out / "flow_direction_code"), codes,
                    XMIN, ymin, RES_M, "Tappa 4 D8 flow direction, ESRI convention (E1 SE2 S4 SW8 W16 NW32 N64 NE128, 0=outlet), int16", dtype="i2")
    for stem in ("contributing_area_km2", "discharge_proxy_m3s", "stream_mask", "lake_mask", "basin_labeled", "flow_direction_code"):
        write_prj(str(out / f"{stem}.prj"), PROJ4)

    meta = {
        "grid": {"ny": ny, "nx": nx, "res_m": RES_M, "xmin": XMIN, "ymax": YMAX},
        "params": {
            "stream_threshold_km2": STREAM_THRESHOLD_KM2,
            "lake_fill_threshold_m": LAKE_FILL_THRESHOLD_M,
            "major_basin_min_km2": MAJOR_BASIN_MIN_KM2,
            "priority_flood_epsilon_m": 1e-4,
        },
        "summary": {
            "land_area_km2": float(land_area_km2),
            "max_contributing_area_km2": float(area_km2.max()),
            "max_discharge_proxy_m3s": float(discharge_proxy_m3s.max()),
            "stream_cells": int(stream_mask.sum()),
            "stream_length_km_approx": float(stream_length_km),
            "drainage_density_km_per_km2": float(drainage_density),
            "lake_cells": int(lake_mask.sum()),
            "lake_area_km2": float(lake_mask.sum() * CELL_KM2),
            "n_basins_total": int(len(basin_ids)),
            "n_basins_major": int(n_major),
            "major_basin_area_km2_sum": float((basin_area_km2[np.isin(basin_ids, major)]).sum()),
        },
        "windward_leeward_asymmetry": asym,
        "runtime_s": time.time() - t0,
    }
    (out / "tappa4_hydrology_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{time.time()-t0:6.1f}s] wrote {out}")
    print(json.dumps(meta["summary"], indent=2))


if __name__ == "__main__":
    main()
