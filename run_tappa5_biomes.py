#!/usr/bin/env python3
"""
Tappa 5 driver — biome classification from Tappa 2's monthly temperature and
precipitation stacks, via Holdridge Life Zones (src/biomes/holdridge.py)
fine-tuned to this world's own data (src/biomes/world_biomes.py). Permanent
snow/ice is overridden from Tappa 3's mass-balance mask. See
docs/decisions/05_tappa5_biomes.md for the full decision record, including
two rejected drafts (a 10-class quartile scheme with an unreadable 96 km^2
sliver, and Holdridge's own unmodified moisture bands, which classify
nothing on this world drier than "Subhumid").

    python run_tappa5_biomes.py [--climate-dir DIR] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from biomes.holdridge import BELT_NAMES                        # noqa: E402
from biomes.world_biomes import (                               # noqa: E402
    BIOME_NAMES, BIOME_COLORS_HEX, classify_world_biomes,
)
from climate.snow import SnowParams, annual_mass_balance        # noqa: E402
from terrain.raster_io import write_envi_raw, write_prj         # noqa: E402

XMIN, XMAX, YMIN, YMAX = -65000.0, 65000.0, -80000.0, 80000.0
PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
CELL_KM2 = 0.12 ** 2  # 120 m working grid


def hillshade(elev: np.ndarray, res_m: float = 120.0, azimuth_deg: float = 315.0, altitude_deg: float = 45.0) -> np.ndarray:
    dy, dx = np.gradient(elev, res_m, res_m)
    slope = np.pi / 2 - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(-dx, dy)
    az, alt = np.deg2rad(azimuth_deg), np.deg2rad(altitude_deg)
    return np.clip(np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect), 0, 1)


def fragmentation_stats(biome_id: np.ndarray) -> dict:
    stats = {}
    for i in range(1, len(BIOME_NAMES)):
        mask = biome_id == i
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        lbl, n = ndimage.label(mask, structure=np.ones((3, 3)))
        sizes = ndimage.sum(mask, lbl, index=np.arange(1, n + 1))
        stats[BIOME_NAMES[i]] = {
            "area_km2": round(n_cells * CELL_KM2, 1),
            "n_patches": int(n),
            "median_patch_km2": round(float(np.median(sizes)) * CELL_KM2, 4),
            "largest_patch_km2": round(float(sizes.max()) * CELL_KM2, 1),
            "fraction_area_in_patches_lt_1km2": round(float(sizes[sizes < 7].sum() / sizes.sum()), 4),
        }
    return stats


def validate_against_nz(elev, land, annual_precip, biome_id, dry_scrub_id=6) -> dict:
    """Windward/leeward split (annual-precip tercile, same method as
    Tappa 3/4) applied to the biome layer itself, plus the specific
    asymmetry the wind-direction decision (02_tappa2_climate.md S1) predicts:
    a dry-adapted band should be wider on the leeward (dry, NE) side than
    the windward (wet, SW) side."""
    land_precip = annual_precip[land]
    p33, p67 = np.percentile(land_precip, [33.33, 66.67])
    wet_mask = land & (annual_precip >= p67)
    dry_mask = land & (annual_precip <= p33)

    def biome_frac(mask):
        n = mask.sum()
        return {BIOME_NAMES[i]: round(float(np.sum(mask & (biome_id == i)) / n), 4)
                for i in range(1, len(BIOME_NAMES)) if n}

    dry_scrub_wet_side = int(np.sum(wet_mask & (biome_id == dry_scrub_id)))
    dry_scrub_dry_side = int(np.sum(dry_mask & (biome_id == dry_scrub_id)))
    return {
        "windward_wet_tercile_biome_fraction": biome_frac(wet_mask),
        "leeward_dry_tercile_biome_fraction": biome_frac(dry_mask),
        "dry_scrub_km2_windward_side": round(dry_scrub_wet_side * CELL_KM2, 1),
        "dry_scrub_km2_leeward_side": round(dry_scrub_dry_side * CELL_KM2, 1),
        "dry_scrub_leeward_to_windward_ratio": (
            round(dry_scrub_dry_side / dry_scrub_wet_side, 2) if dry_scrub_wet_side
            else ("leeward-only, zero windward cells" if dry_scrub_dry_side else None)
        ),
    }


def render_overview(elev, land, biome_id, out_png: Path, title: str):
    hs = hillshade(elev)
    n = len(BIOME_NAMES)
    cmap = ListedColormap(BIOME_COLORS_HEX)
    norm = BoundaryNorm(np.arange(-0.5, n + 0.5, 1), cmap.N)
    fig, ax = plt.subplots(figsize=(11, 13.5), dpi=170)
    ax.imshow(hs, cmap="gray", vmin=0.15, vmax=1.0)
    ax.imshow(biome_id, cmap=cmap, norm=norm, alpha=0.72)
    ax.set_axis_off()
    ax.set_title(title, fontsize=11)
    areas = [np.sum(biome_id == i) * CELL_KM2 for i in range(n)]
    handles = [Patch(facecolor=BIOME_COLORS_HEX[i], edgecolor="none",
                     label=f"{BIOME_NAMES[i]}  ({areas[i]:.0f} km2)")
               for i in range(n) if areas[i] > 0]
    ax.legend(handles=handles, loc="lower left", fontsize=8, framealpha=0.9, title="Biome")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climate-dir", default="data/processed/climate")
    ap.add_argument("--out", default="data/processed/biomes")
    ap.add_argument("--assets", default="docs/decisions/assets/05_tappa5_biomes")
    args = ap.parse_args()

    t0 = time.time()
    cdir = Path(args.climate_dir)
    temp_c = np.load(cdir / "temperature_monthly_c.npy")
    precip_mm = np.load(cdir / "precipitation_monthly_mm.npy")
    land = np.load(cdir / "land_mask.npy").astype(bool)
    elev = np.load(cdir / "surface_elevation_m.npy")
    ny, nx = elev.shape
    res = 120.0
    print(f"[{time.time()-t0:6.1f}s] loaded Tappa 2 outputs: {ny}x{nx} @ {res:.0f} m, land {land.mean():.3f}")

    # Recompute Tappa 3's permanent-snow mask directly from the snow module
    # (cheap, ~1 s) rather than depend on a possibly-stale saved copy.
    snow_params = SnowParams()
    _, _, balance_mm = annual_mass_balance(precip_mm, temp_c, DAYS, snow_params)
    permanent_snow = (balance_mm >= 0.0)
    print(f"[{time.time()-t0:6.1f}s] permanent snow (mass balance): "
          f"{(permanent_snow & land).sum()*CELL_KM2:.1f} km2")

    result = classify_world_biomes(temp_c, precip_mm, land, permanent_snow)
    print(f"[{time.time()-t0:6.1f}s] classified. moisture tercile edges "
          f"(PET ratio): {result.moisture_tercile_edges}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    assets = Path(args.assets)
    assets.mkdir(parents=True, exist_ok=True)

    # --- save arrays ---------------------------------------------------
    np.save(out / "biome_id.npy", result.biome_id)
    np.save(out / "biotemperature_c.npy", result.holdridge.biotemperature_c)
    np.save(out / "pet_ratio.npy", result.holdridge.pet_ratio)
    np.save(out / "moisture_idx.npy", result.moisture_idx)

    ymin = YMAX - ny * res
    write_envi_raw(str(out / "biome_id"), result.biome_id.astype(np.int16),
                   XMIN, ymin, res, "Tappa 5 biome classification, int16 (0=ocean, 1-9=land biomes, see meta json)", dtype="i2")
    write_envi_raw(str(out / "biotemperature_c"), result.holdridge.biotemperature_c,
                   XMIN, ymin, res, "Tappa 5 Holdridge biotemperature, float32 C", dtype="f4")
    write_envi_raw(str(out / "pet_ratio"), result.holdridge.pet_ratio,
                   XMIN, ymin, res, "Tappa 5 Holdridge PET ratio (PET/annual precip), float32", dtype="f4")
    for stem in ("biome_id", "biotemperature_c", "pet_ratio"):
        write_prj(str(out / f"{stem}.prj"), PROJ4)

    # --- diagnostics / validation ---------------------------------------
    frag = fragmentation_stats(result.biome_id)
    annual_precip = precip_mm.sum(axis=0)
    nz_check = validate_against_nz(elev, land, annual_precip, result.biome_id)

    belt_km2 = {BELT_NAMES[b]: round(float(np.sum(land & (result.holdridge.belt_idx == b))) * CELL_KM2, 1)
                for b in range(len(BELT_NAMES))}
    areas_km2 = {BIOME_NAMES[i]: round(float(np.sum(result.biome_id == i)) * CELL_KM2, 1)
                 for i in range(len(BIOME_NAMES))}

    meta = {
        "grid": {"ny": ny, "nx": nx, "res_m": res, "xmin": XMIN, "ymax": YMAX},
        "moisture_tercile_pet_ratio_edges": result.moisture_tercile_edges,
        "biotemperature_c_land": {
            "min": float(result.holdridge.biotemperature_c[land].min()),
            "max": float(result.holdridge.biotemperature_c[land].max()),
            "mean": float(result.holdridge.biotemperature_c[land].mean()),
        },
        "pet_ratio_land": {
            "min": float(result.holdridge.pet_ratio[land].min()),
            "max": float(result.holdridge.pet_ratio[land].max()),
            "mean": float(result.holdridge.pet_ratio[land].mean()),
        },
        "holdridge_belt_km2_pre_fold": belt_km2,
        "biome_area_km2": areas_km2,
        "fragmentation": frag,
        "permanent_snow_km2": round(float((permanent_snow & land).sum()) * CELL_KM2, 1),
        "holdridge_polar_belt_km2": belt_km2.get("Polar"),
        "validation_vs_nz": nz_check,
    }
    (out / "tappa5_biomes_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{time.time()-t0:6.1f}s] wrote {out}")
    print(json.dumps({"biome_area_km2": areas_km2, "validation_vs_nz": nz_check}, indent=2))

    # --- overview figure --------------------------------------------------
    render_overview(elev, land, result.biome_id, assets / "05_overview.png",
                     "Tappa 5 — biomes (Holdridge belt x world-recalibrated moisture tercile,\n"
                     "permanent snow overridden by Tappa 3 mass balance)")
    print(f"[{time.time()-t0:6.1f}s] wrote {assets / '05_overview.png'}")


if __name__ == "__main__":
    main()
