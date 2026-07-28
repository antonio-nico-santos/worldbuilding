#!/usr/bin/env python3
"""
Tappa 3 driver — months-with-snow, seasonality index, and a mass-balance
permanent-snow / ELA proxy, all derived from Tappa 2's monthly temperature
and precipitation stacks. See docs/decisions/03_tappa3_snow.md for the
model and its validation; see src/climate/snow.py for the maths.

    python run_tappa3_snow.py [--climate-dir DIR] [--out DIR] [--params PATH]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from climate.grid import latitude_grid                        # noqa: E402
from climate.snow import (                                     # noqa: E402
    SnowParams,
    annual_mass_balance,
    monthly_snow_rain,
    months_with_snow,
    seasonality_index,
)
from climate.temperature import TemperatureParams, temperature_year  # noqa: E402
from terrain.raster_io import write_envi_raw, write_prj        # noqa: E402

XMIN, XMAX, YMIN, YMAX = -65000.0, 65000.0, -80000.0, 80000.0
LCC = dict(lat_1=-44.48, lat_2=-43.52, lat_0=-44.0)
PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()

_S_MAP = {
    "t50_snow_rain_c": "t50_snow_rain_c",
    "sigma_day_c": "sigma_day_c",
    "degree_day_factor_mm_per_c_day": "degree_day_factor_mm_per_c_day",
    "min_snow_month_mm": "min_snow_month_mm",
}


def load_snow_params(path: str) -> SnowParams:
    p = SnowParams()
    if not path or not Path(path).exists():
        return p
    import yaml
    cfg = (yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}).get("climate", {}).get("snow") or {}
    for key, field in _S_MAP.items():
        if key in cfg:
            setattr(p, field, float(cfg[key]))
    return p


def ela_by_band(elevation_m: np.ndarray, balance_mm: np.ndarray, mask: np.ndarray, bin_m: float = 100.0):
    """Bin cells (restricted to `mask`) by elevation, average balance per
    bin, and linearly interpolate the elevation where the binned mean
    balance crosses zero (accumulation = melt). Returns (ela_m, bin_table)
    where bin_table is a list of (elev_mid, mean_balance, n_cells) for
    plotting/inspection. Returns (None, table) if the balance never
    crosses zero within the sampled elevation range.
    """
    z = elevation_m[mask]
    b = balance_mm[mask]
    if z.size == 0:
        return None, []
    lo = np.floor(z.min() / bin_m) * bin_m
    hi = np.ceil(z.max() / bin_m) * bin_m
    edges = np.arange(lo, hi + bin_m, bin_m)
    table = []
    for i in range(len(edges) - 1):
        sel = (z >= edges[i]) & (z < edges[i + 1])
        n = int(sel.sum())
        if n == 0:
            continue
        table.append((float((edges[i] + edges[i + 1]) / 2.0), float(b[sel].mean()), n))
    ela = None
    for (z0, b0, _), (z1, b1, _) in zip(table, table[1:]):
        # balance rises with elevation (colder -> less melt, often more
        # orographic snow), so the ELA crossing goes negative -> non-negative
        if b0 < 0 and b1 >= 0:
            # linear interpolation to the zero crossing between the two bins
            ela = z0 + (0.0 - b0) * (z1 - z0) / (b1 - b0)
            break
    return ela, table


def run(temp_c, precip_mm, elevation_m, land, snow_params, tag=""):
    """One full Tappa-3 derived-metrics pass, factored out so the lapse-rate
    sensitivity check (S4 of the decision doc) can reuse it on an
    alternately-generated temperature stack without duplicating logic."""
    snow_mm, rain_mm = monthly_snow_rain(precip_mm, temp_c, snow_params)
    n_snow_months = months_with_snow(snow_mm, snow_params)
    seasonality_c = seasonality_index(temp_c)
    accum_mm, melt_mm, balance_mm = annual_mass_balance(precip_mm, temp_c, DAYS, snow_params)
    permanent_mask = (balance_mm >= 0.0) & land

    # naive Tappa-2-style comparison metric: all 12 monthly means < 0 C
    naive_permanent_mask = (temp_c.max(axis=0) < 0.0) & land

    # windward/leeward split via annual precipitation tercile, land cells only
    annual_precip = precip_mm.sum(axis=0)
    land_precip = annual_precip[land]
    p33, p67 = np.percentile(land_precip, [33.33, 66.67])
    wet_mask = land & (annual_precip >= p67)
    dry_mask = land & (annual_precip <= p33)

    ela_wet, table_wet = ela_by_band(elevation_m, balance_mm, wet_mask)
    ela_dry, table_dry = ela_by_band(elevation_m, balance_mm, dry_mask)

    result = {
        "snow_mm": snow_mm, "rain_mm": rain_mm,
        "n_snow_months": n_snow_months, "seasonality_c": seasonality_c,
        "accum_mm": accum_mm, "melt_mm": melt_mm, "balance_mm": balance_mm,
        "permanent_mask": permanent_mask, "naive_permanent_mask": naive_permanent_mask,
        "wet_mask": wet_mask, "dry_mask": dry_mask,
    }
    summary = {
        "permanent_snow_area_km2": float(permanent_mask.sum() * (0.12 ** 2)),  # 120 m cells -> km2 (0.12 km/side)
        "permanent_snow_fraction_of_land": float(permanent_mask.sum() / land.sum()),
        "naive_permanent_snow_area_km2": float(naive_permanent_mask.sum() * (0.12 ** 2)),
        "naive_permanent_snow_fraction_of_land": float(naive_permanent_mask.sum() / land.sum()),
        "ela_windward_m": ela_wet,
        "ela_leeward_m": ela_dry,
        "ela_differential_m": (ela_dry - ela_wet) if (ela_wet is not None and ela_dry is not None) else None,
        "n_snow_months_land_mean": float(n_snow_months[land].mean()),
        "n_snow_months_land_median": float(np.median(n_snow_months[land])),
        "seasonality_c_land_mean": float(seasonality_c[land].mean()),
        "seasonality_c_land_min": float(seasonality_c[land].min()),
        "seasonality_c_land_max": float(seasonality_c[land].max()),
    }
    if tag:
        print(f"  [{tag}] permanent snow: {summary['permanent_snow_area_km2']:.1f} km2 "
              f"({summary['permanent_snow_fraction_of_land']*100:.2f}% of land), "
              f"naive: {summary['naive_permanent_snow_area_km2']:.1f} km2 "
              f"({summary['naive_permanent_snow_fraction_of_land']*100:.2f}%)  "
              f"ELA windward={ela_wet} m, leeward={ela_dry} m")
    return result, summary, {"wet": table_wet, "dry": table_dry}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climate-dir", default="data/processed/climate")
    ap.add_argument("--out", default="data/processed/climate")
    ap.add_argument("--params", default="config/parameters.yml")
    args = ap.parse_args()

    t0 = time.time()
    cdir = Path(args.climate_dir)
    temp_c = np.load(cdir / "temperature_monthly_c.npy")
    precip_mm = np.load(cdir / "precipitation_monthly_mm.npy")
    elevation_m = np.load(cdir / "surface_elevation_m.npy")
    land = np.load(cdir / "land_mask.npy")
    ny, nx = elevation_m.shape
    res = 120.0
    print(f"[{time.time()-t0:6.1f}s] loaded Tappa 2 outputs: {ny}x{nx} @ {res:.0f} m, land {land.mean():.3f}")

    snow_params = load_snow_params(args.params)
    print(f"         snow params: {snow_params}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    result, summary, band_tables = run(temp_c, precip_mm, elevation_m, land, snow_params, tag="v1 (locked)")

    # --- sensitivity check: lapse_seasonal_amplitude_c_per_km, flagged in
    # 02_tappa2_climate.md S8 as needing a re-check here since it moves
    # summit seasonality and therefore the permanent-snow area. -----------
    print("[sensitivity] re-running with lapse_seasonal_amplitude_c_per_km = 0.0 (vs locked 0.3)")
    dist = np.load(cdir / "distance_to_coast_km.npy")
    lat = latitude_grid(ny, nx, YMAX, XMIN, res, **LCC)
    tp_flat = TemperatureParams(lapse_seasonal_amplitude_c_per_km=0.0)
    temp_c_flat = temperature_year(elevation_m, lat, dist, tp_flat)
    _, summary_flat, _ = run(temp_c_flat, precip_mm, elevation_m, land, snow_params, tag="sensitivity: amplitude=0.0")

    # --- save arrays -------------------------------------------------------
    np.save(out / "months_with_snow.npy", result["n_snow_months"])
    np.save(out / "seasonality_index_c.npy", result["seasonality_c"])
    np.save(out / "snow_accum_mm.npy", result["accum_mm"])
    np.save(out / "snow_melt_mm.npy", result["melt_mm"])
    np.save(out / "mass_balance_mm.npy", result["balance_mm"])
    np.save(out / "permanent_snow_mask.npy", result["permanent_mask"])

    ymin = YMAX - ny * res
    write_envi_raw(str(out / "months_with_snow"), result["n_snow_months"].astype(np.int16),
                   XMIN, ymin, res, "Tappa 3 months with modelled snow (0-12), int16", dtype="i2")
    write_envi_raw(str(out / "seasonality_index_c"), result["seasonality_c"],
                   XMIN, ymin, res, "Tappa 3 seasonality index, warmest-coldest month mean, float32 C", dtype="f4")
    write_envi_raw(str(out / "mass_balance_mm"), result["balance_mm"],
                   XMIN, ymin, res, "Tappa 3 annual snow mass balance, float32 mm w.e.", dtype="f4")
    write_envi_raw(str(out / "permanent_snow_mask"), result["permanent_mask"].astype(np.int16),
                   XMIN, ymin, res, "Tappa 3 permanent snow / perennial firn proxy mask, int16 (0/1)", dtype="i2")
    for stem in ("months_with_snow", "seasonality_index_c", "mass_balance_mm", "permanent_snow_mask"):
        write_prj(str(out / f"{stem}.prj"), PROJ4)

    meta = {
        "grid": {"ny": ny, "nx": nx, "res_m": res, "xmin": XMIN, "ymax": YMAX},
        "snow_params": vars(snow_params),
        "summary": summary,
        "sensitivity_lapse_seasonal_amplitude_0": summary_flat,
        "elevation_balance_bands": band_tables,
    }
    (out / "tappa3_snow_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{time.time()-t0:6.1f}s] wrote {out}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
