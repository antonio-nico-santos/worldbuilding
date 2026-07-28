#!/usr/bin/env python3
"""
Tappa 4 addendum — seasonal / intermittent flow. Monthly discharge and a
per-stream-cell "months flowing" (0-12) classification, driven by a real
monthly snowpack simulation rather than raw monthly precipitation. See
docs/decisions/04_tappa4_hydrology.md S12.

Re-runs the priority-flood fill/direction pass (~335 s) rather than
reloading `filled_dem_30m.npy` from disk: that file was saved as float32,
whose ~0.0001-0.001 m precision at this DEM's elevation range collides
with the priority-flood's 1e-4 m epsilon tie-breaker, silently corrupting
the strict-descent ordering for ~22K cells (confirmed directly, not
theoretical -- see decision doc S12). The receiver/pop_order this stage
needs have to come from the same in-memory float64 pass Tappa 4's main
run used, not a round trip through a lossy export.

    python run_tappa4_seasonal.py [--dem PATH] [--climate-dir DIR] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from climate.snow import SnowParams, monthly_water_input          # noqa: E402
from hydrology.flow import accumulate_flow, priority_flood_d8      # noqa: E402
from hydrology.weighting import upsample_precip_to_dem             # noqa: E402
from terrain.raster_io import write_envi_raw, write_prj            # noqa: E402

XMIN, YMAX, RES_M = -65000.0, 80000.0, 30.0
PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()
SECONDS_PER_YEAR = 365.25 * 86400.0
FLOWING_FRACTION = 0.10   # "flowing" if month's discharge >= this fraction
                          # of the cell's own mean annual discharge (S12)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", default="data/processed/dem_v3_final_30m_eroded.npy")
    ap.add_argument("--climate-dir", default="data/processed/climate")
    ap.add_argument("--hydro-dir", default="data/processed/hydrology")
    ap.add_argument("--out", default="data/processed/hydrology/seasonal")
    args = ap.parse_args()

    t0 = time.time()
    dem = np.load(args.dem).astype(np.float64)
    ny, nx = dem.shape
    land = dem > 0.0
    stream_mask = np.load(Path(args.hydro_dir) / "stream_mask.npy")
    print(f"[{time.time()-t0:6.1f}s] loaded DEM {ny}x{nx}, "
          f"{stream_mask.sum()} stream cells")

    # --- re-run priority-flood for a full-precision receiver/pop_order ---
    seed = dem <= 0.0
    seed[0, :] = seed[-1, :] = seed[:, 0] = seed[:, -1] = True
    filled, receiver, pop_order = priority_flood_d8(dem, seed, epsilon=1e-4)
    print(f"[{time.time()-t0:6.1f}s] priority-flood re-run done "
          f"(in-memory float64, not reloaded from disk)")

    # sanity check: reproduce the ORIGINAL annual discharge_proxy_m3s with
    # this fresh receiver/pop_order and the same weighting it was built
    # with -- must match closely, or something upstream has changed.
    precip_monthly = np.load(Path(args.climate_dir) / "precipitation_monthly_mm.npy")
    temp_monthly = np.load(Path(args.climate_dir) / "temperature_monthly_c.npy")
    precip_annual_30 = upsample_precip_to_dem(precip_monthly.sum(axis=0), (ny, nx))
    w_check = np.where(land, precip_annual_30, 0.0)
    accum_check = accumulate_flow(receiver, pop_order, w_check).reshape(ny, nx)
    disch_check = accum_check * (RES_M ** 2 / 1000.0) / SECONDS_PER_YEAR
    orig_disch = np.load(Path(args.hydro_dir) / "discharge_proxy_m3s.npy")
    max_abs_diff = np.abs(disch_check - orig_disch).max()
    corr = np.corrcoef(disch_check.ravel(), orig_disch.ravel())[0, 1]
    print(f"[{time.time()-t0:6.1f}s] sanity check vs original annual run: "
          f"max abs diff {max_abs_diff:.4f} m3/s, correlation {corr:.6f}")
    assert corr > 0.999 and max_abs_diff < 1.0, (
        "reconstructed flow graph does not reproduce the original annual "
        "discharge -- do not trust the monthly routing below")

    # --- monthly snowpack simulation (climate/snow.py) --------------------
    p = SnowParams()
    water_input_mm, melt_mm, rain_mm, snowpack_end_mm = monthly_water_input(
        precip_monthly, temp_monthly, DAYS, p, spinup_cycles=3)
    print(f"[{time.time()-t0:6.1f}s] monthly snowpack simulation done "
          f"(3 spin-up cycles), mean end-of-year snowpack "
          f"{snowpack_end_mm[np.load(Path(args.climate_dir)/'land_mask.npy')].mean():.1f} mm w.e.")

    # --- route each month's water input through the SAME flow network ----
    # NOTE: divide by THAT MONTH's seconds (days_in_month * 86400), not a
    # full year -- a month's accumulated water volume spread over a full
    # year understates its actual mean rate during the month by ~11.8-12.9x
    # (caught via a magnitude sanity check against the annual run: monthly
    # maxima came out an order of magnitude below the annual max, which is
    # not physically possible for a mean-monthly vs. mean-annual rate).
    monthly_discharge = np.zeros((12, ny, nx), dtype=np.float32)
    for m in range(12):
        water_30 = upsample_precip_to_dem(water_input_mm[m], (ny, nx))
        w = np.where(land, water_30, 0.0)
        accum = accumulate_flow(receiver, pop_order, w).reshape(ny, nx)
        seconds_this_month = DAYS[m] * 86400.0
        monthly_discharge[m] = accum * (RES_M ** 2 / 1000.0) / seconds_this_month
        print(f"[{time.time()-t0:6.1f}s]   {MONTHS[m]}: max discharge "
              f"{monthly_discharge[m].max():.1f} m3/s")

    # --- classify months flowing, stream cells only ------------------------
    mean_annual = monthly_discharge.mean(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        flowing = monthly_discharge >= (FLOWING_FRACTION * mean_annual)
    months_flowing = np.where(stream_mask, flowing.sum(axis=0), 0).astype(np.int16)
    print(f"[{time.time()-t0:6.1f}s] months_flowing classified "
          f"(threshold: {FLOWING_FRACTION*100:.0f}% of each cell's own annual mean)")

    n_stream = stream_mask.sum()
    dist = np.bincount(months_flowing[stream_mask], minlength=13)
    print("  months-flowing distribution over stream cells:",
          {i: int(dist[i]) for i in range(13) if dist[i] > 0})
    perennial = int((months_flowing == 12).sum())
    intermittent = int(((months_flowing >= 1) & (months_flowing <= 11)).sum())
    print(f"  perennial (12/12): {perennial} ({perennial/n_stream*100:.1f}%), "
          f"intermittent (1-11): {intermittent} ({intermittent/n_stream*100:.1f}%)")

    # --- write outputs -----------------------------------------------------
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ymin = YMAX - ny * RES_M

    np.save(out / "monthly_discharge_proxy_m3s.npy", monthly_discharge)
    np.save(out / "months_flowing.npy", months_flowing)

    write_envi_raw(str(out / "monthly_discharge_proxy_m3s"), monthly_discharge,
                    XMIN, ymin, RES_M,
                    "Tappa 4 seasonal addendum: monthly discharge proxy (snowmelt-timed), float32 m3/s, 12 bands",
                    dtype="f4", band_names=MONTHS)
    write_envi_raw(str(out / "months_flowing"), months_flowing,
                    XMIN, ymin, RES_M,
                    f"Tappa 4 seasonal addendum: months flowing (0-12), threshold {FLOWING_FRACTION*100:.0f}% of own annual mean, int16",
                    dtype="i2")
    for stem in ("monthly_discharge_proxy_m3s", "months_flowing"):
        write_prj(str(out / f"{stem}.prj"), PROJ4)

    meta = {
        "grid": {"ny": ny, "nx": nx, "res_m": RES_M, "xmin": XMIN, "ymax": YMAX},
        "params": {
            "flowing_fraction_of_annual_mean": FLOWING_FRACTION,
            "snow_params": vars(p),
            "spinup_cycles": 3,
        },
        "sanity_check_vs_original_annual_run": {
            "max_abs_diff_m3s": float(max_abs_diff),
            "correlation": float(corr),
        },
        "summary": {
            "n_stream_cells": int(n_stream),
            "perennial_cells": perennial,
            "perennial_fraction": perennial / n_stream,
            "intermittent_cells": intermittent,
            "intermittent_fraction": intermittent / n_stream,
            "months_flowing_distribution": {int(i): int(dist[i]) for i in range(13)},
        },
        "runtime_s": time.time() - t0,
    }
    (out / "tappa4_seasonal_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{time.time()-t0:6.1f}s] wrote {out}")
    print(json.dumps(meta["summary"], indent=2))


if __name__ == "__main__":
    main()
