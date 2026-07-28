#!/usr/bin/env python3
"""
Tappa 2 driver — monthly temperature and precipitation for the whole domain.

Reads the eroded v3 DEM, coarsens it to the 120 m climate working grid, and
writes twelve monthly layers of each variable plus an annual summary.
See docs/decisions/02_tappa2_climate.md for the model and its calibration.

    python run_tappa2_climate.py [--dem PATH] [--out DIR] [--factor 4]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from climate import precipitation as pp                       # noqa: E402
from terrain.raster_io import write_envi_raw, write_prj       # noqa: E402
from climate.grid import coarsen, latitude_grid               # noqa: E402
from climate.temperature import (                             # noqa: E402
    TemperatureParams,
    distance_to_coast_km,
    sea_level_temperature,
    temperature_month,
)

# Tappa 0 domain, in the project's custom LCC (config/parameters.yml)
XMIN, XMAX, YMIN, YMAX = -65000.0, 65000.0, -80000.0, 80000.0
LCC = dict(lat_1=-44.48, lat_2=-43.52, lat_0=-44.0)
PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()


# Mapping from config/parameters.yml keys to the dataclass fields. Kept
# explicit rather than name-matched: the YAML uses the project's existing
# `_c_per_1000m` / signed-lapse convention (negative = cooling with height),
# the code uses a positive C/km, and silently coercing between the two is
# exactly the kind of thing that produces an upside-down mountain.
_T_MAP = {
    "t_ref_sea_level_c": "t_ref_sea_level_c",
    "lat_ref_deg": "lat_ref_deg",
    "dt_dlat_c_per_deg": "dt_dlat_c_per_deg",
    "amplitude_coastal_c": "amplitude_coastal_c",
    "amplitude_inland_c": "amplitude_inland_c",
    "continentality_scale_km": "continentality_scale_km",
    "peak_month_coastal": "peak_month_coastal",
    "peak_month_shift_inland": "peak_month_shift_inland",
    "lapse_peak_month": "lapse_peak_month",
    "lapse_seasonal_amplitude_c_per_1000m": "lapse_seasonal_amplitude_c_per_km",
}
_P_MAP = {
    "wind_direction_from_deg": "wind_from_deg",
    "wind_speed_ms": "wind_speed_ms",
    "wind_speed_seasonal_amplitude_ms": "wind_speed_seasonal_amplitude_ms",
    "wind_peak_month": "wind_peak_month",
    "nm_moist_stability_per_s": "nm_moist_stability",
    "tau_c_s": "tau_c_s",
    "tau_f_s": "tau_f_s",
    "orographic_duty_cycle": "orographic_duty_cycle",
    "background_precip_mm_per_month": "background_precip_mm_per_month",
    "lee_floor_fraction": "lee_floor_fraction",
    "fft_pad_km": "pad_km",
}


def load_climate_params(path: str):
    """Read config/parameters.yml into the two dataclasses.

    The README calls that file the single source of truth; this makes it
    actually true for Tappa 2. Missing keys fall back to the dataclass
    defaults, and passing an empty path skips the file entirely.
    """
    tp, pp_ = TemperatureParams(), pp.PrecipParams()
    factor = 4
    if not path or not Path(path).exists():
        print(f"         (no parameter file at {path!r}; using dataclass defaults)")
        return tp, pp_, factor
    import yaml

    cfg = (yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}).get("climate") or {}
    res = cfg.get("working_resolution_m")
    dem_res = 30.0
    if res:
        if res % dem_res:
            raise ValueError(f"working_resolution_m {res} is not a multiple of the DEM's {dem_res:.0f} m")
        factor = int(res // dem_res)
    for key, field in _T_MAP.items():
        if key in cfg.get("temperature", {}):
            setattr(tp, field, float(cfg["temperature"][key]))
    if "lapse_rate_c_per_1000m" in cfg.get("temperature", {}):
        # signed in the file (negative = cooling with height), positive in code
        tp.lapse_rate_c_per_km = abs(float(cfg["temperature"]["lapse_rate_c_per_1000m"]))
    for key, field in _P_MAP.items():
        if key in cfg.get("precipitation", {}):
            setattr(pp_, field, float(cfg["precipitation"][key]))
    print(f"         parameters from {path}")
    return tp, pp_, factor


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", default="data/processed/dem_v3_final_30m_eroded.npy")
    ap.add_argument("--out", default="data/processed/climate")
    ap.add_argument("--factor", type=int, default=0,
                    help="DEM coarsening factor; 0 = derive from parameters.yml")
    ap.add_argument("--params", default="config/parameters.yml",
                    help="parameter file; pass '' to use the dataclass defaults")
    args = ap.parse_args()

    tp, pparams, factor = load_climate_params(args.params)
    if args.factor:
        factor = args.factor
    args.factor = factor

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res = 30.0 * args.factor

    t0 = time.time()
    dem = np.load(args.dem)
    surface, uplift, land_fraction = coarsen(dem, args.factor)
    del dem
    ny, nx = surface.shape
    land = surface > 0.0
    print(f"[{time.time()-t0:6.1f}s] grid {ny}x{nx} @ {res:.0f} m, land {land.mean():.3f}")

    lat = latitude_grid(ny, nx, YMAX, XMIN, res, **LCC)
    dist = distance_to_coast_km(land, res)
    print(f"         lat {lat.min():.3f}..{lat.max():.3f}, max dist-to-coast {dist.max():.1f} km")

    # Domain-mean SEA-LEVEL temperature per month drives the moisture terms.
    # Sea level, not surface: Cw and Hw describe the incoming airstream over
    # the ocean upwind, not the air already sitting on a 3 km summit.
    t_sl = np.array([float(sea_level_temperature(m, lat, dist, tp).mean()) for m in range(1, 13)])
    u_month = np.array([float(pp.wind_speed(m, pparams)) for m in range(1, 13)])
    rho_s = np.array([float(pp.saturation_vapour_density(t)) for t in t_sl])
    reference_flux = float((u_month * rho_s).mean())

    temp = np.empty((12, ny, nx), np.float32)
    prec = np.empty((12, ny, nx), np.float32)
    states = []
    for i, m in enumerate(range(1, 13)):
        temp[i] = temperature_month(m, surface, lat, dist, tp)
        prec[i], st = pp.precip_month(
            uplift, res, m, t_sl[i], tp.lapse_rate_c_per_km, pparams,
            days=DAYS[i], reference_flux=reference_flux,
        )
        states.append({k: float(v) for k, v in st.items()})
        print(f"[{time.time()-t0:6.1f}s] {MONTHS[i]}  T_sl {t_sl[i]:5.2f}C  U {st['wind_speed_ms']:4.1f} m/s"
              f"  Cw_eff {st['cw_effective_kg_m3']*1e3:.3f} g/m3  P {prec[i][land].mean():6.1f} mm"
              f"  T_land {temp[i][land].mean():5.2f}C")

    np.save(out / "temperature_monthly_c.npy", temp)  # noqa: E501  (.npy vs ENVI .bin, no clash)
    np.save(out / "precipitation_monthly_mm.npy", prec)
    np.save(out / "surface_elevation_m.npy", surface.astype(np.float32))
    np.save(out / "distance_to_coast_km.npy", dist)
    np.save(out / "land_mask.npy", land)

    # --- ENVI export for QGIS ------------------------------------------
    # ymin is derived from the grid's own shape, not from the ROI constant:
    # 1334 rows x 120 m overhangs YMIN by 80 m (the same 20 m overhang the
    # 30 m DEM already had, times the coarsening factor). Anchoring the
    # header to the array keeps the north edge exactly on YMAX.
    ymin = YMAX - ny * res
    # Temperature ships as float32, in plain degrees Celsius. Tappa 1 exported
    # int16 because halving a 92 MB DEM was worth one documented unit
    # convention; that argument does not carry over here. This stack is 69 MB
    # instead of 35, it is gitignored, and it regenerates in 19 s — so the
    # 35 MB buys away a raster whose pixels read "-821" in the QGIS identify
    # tool and silently mean -8.21 C. Wrong trade. Precipitation stays int16
    # because whole millimetres need no scaling to be readable.
    write_envi_raw(str(out / "temperature_monthly_c"), temp,
                   XMIN, ymin, res, "Tappa 2 monthly mean temperature, float32 degrees C",
                   dtype="f4", band_names=MONTHS)
    write_envi_raw(str(out / "precipitation_monthly_mm"), np.round(prec),
                   XMIN, ymin, res, "Tappa 2 monthly precipitation, int16 = mm",
                   dtype="i2", band_names=MONTHS)
    ann_p = prec.sum(0)
    ann_t = temp.mean(0)
    write_envi_raw(str(out / "annual_precipitation_mm"), np.round(ann_p),
                   XMIN, ymin, res, "Tappa 2 annual precipitation total, int16 = mm", dtype="i2")
    write_envi_raw(str(out / "annual_mean_temperature_c"), ann_t,
                   XMIN, ymin, res, "Tappa 2 annual mean temperature, float32 degrees C", dtype="f4")
    for stem in ("temperature_monthly_c", "precipitation_monthly_mm",
                 "annual_precipitation_mm", "annual_mean_temperature_c"):
        write_prj(str(out / f"{stem}.prj"), PROJ4)

    meta = {
        "grid": {"ny": ny, "nx": nx, "res_m": res, "xmin": XMIN, "ymax": YMAX, "ymin": float(YMAX - ny * res),
                 "coarsen_factor": args.factor, "source_dem": args.dem},
        "parameter_source": args.params or "dataclass defaults",
        "temperature_params": {k: v for k, v in vars(tp).items() if k != "notes"},
        "precip_params": vars(pparams),
        "monthly_state": states,
        "reference_flux": reference_flux,
        "summary": {
            "annual_precip_mm": {
                "min": float(ann_p[land].min()), "p05": float(np.percentile(ann_p[land], 5)),
                "median": float(np.median(ann_p[land])), "p95": float(np.percentile(ann_p[land], 95)),
                "max": float(ann_p[land].max()), "mean": float(ann_p[land].mean()),
            },
            "annual_mean_temp_c": {
                "min": float(ann_t[land].min()), "median": float(np.median(ann_t[land])),
                "max": float(ann_t[land].max()), "mean": float(ann_t[land].mean()),
            },
            "warmest_month_mean_c": float(temp.max(0)[land].mean()),
            "coldest_month_mean_c": float(temp.min(0)[land].mean()),
        },
    }
    (out / "climate_run_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{time.time()-t0:6.1f}s] wrote {out}")
    print(json.dumps(meta["summary"], indent=2))


if __name__ == "__main__":
    main()
