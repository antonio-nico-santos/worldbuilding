"""
Tappa 6 -- Povo Livre biome suitability variant (forestry/foraging-centred
economy, mirrored Forest/Grassland endpoints -- see
src/suitability/biome_lookup.py's BIOME_SUITABILITY_POVO_LIVRE docstring).

This does NOT recompute the majority-filtered biome map -- it reuses
biome_id_smoothed_120m.npy written by run_tappa6_biome_suitability.py, since
the smoothing step depends only on biome_id/land_mask, not on which score
table is applied to it. Run run_tappa6_biome_suitability.py FIRST.

Reads:
  data/processed/suitability/biome_id_smoothed_120m.npy
  data/processed/climate/land_mask.npy
  data/processed/biomes/biome_id.npy                      (unsmoothed, for
                                                            the precip/
                                                            centroid
                                                            comparison stats
                                                            logged below)
  data/processed/climate/annual_precipitation_mm.bin        (Tappa 2 ENVI
                                                            raw, int16 mm --
                                                            read directly,
                                                            no .npy sidecar
                                                            exists for this
                                                            one)

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  biome_suitability_povo_livre_120m   (0-1, float32)
  tappa6_biome_povo_livre_meta.json   (includes the Circulo-vs-Povo-Livre
                                       comparison numbers: precip contrast,
                                       spatial centroids, score-swing extent)
"""

from __future__ import annotations

import json
import os
import time

import numpy as np

from src.biomes.world_biomes import BIOME_NAMES
from src.suitability.biome_lookup import (
    BIOME_NOTES,
    BIOME_SUITABILITY,
    BIOME_SUITABILITY_POVO_LIVRE,
    biome_suitability_from_id,
)
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
PROJ4 = (
    "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)


def main():
    t0 = time.time()
    biome_smoothed = np.load("data/processed/suitability/biome_id_smoothed_120m.npy")
    land = np.load("data/processed/climate/land_mask.npy").astype(bool)
    biome_raw = np.load("data/processed/biomes/biome_id.npy")

    ny, nx = land.shape
    cs_x = (XMAX - XMIN) / nx
    cs_y = (YMAX - YMIN) / ny

    suit_povo_livre = biome_suitability_from_id(
        biome_smoothed, land, lut=BIOME_SUITABILITY_POVO_LIVRE
    ).astype(np.float32)
    suit_circulo = biome_suitability_from_id(
        biome_smoothed, land, lut=BIOME_SUITABILITY
    ).astype(np.float32)

    out = "data/processed/suitability"
    os.makedirs(out, exist_ok=True)
    np.save(f"{out}/biome_suitability_povo_livre_120m.npy", suit_povo_livre)
    write_envi_raw(
        f"{out}/biome_suitability_povo_livre_120m", suit_povo_livre, XMIN, YMIN, cs_x,
        "Tappa6 biome_suitability_povo_livre_120m", dtype="f4",
    )
    write_prj(f"{out}/biome_suitability_povo_livre_120m.prj", PROJ4)

    # comparison stats against the Circulo (open-field-agriculture) version,
    # on the UNSMOOTHED biome_id (matches the numbers checked in chat before
    # committing to build this)
    xc = XMIN + (np.arange(nx) + 0.5) * cs_x
    X = np.broadcast_to(xc, (ny, nx))
    precip = np.fromfile(
        "data/processed/climate/annual_precipitation_mm.bin", dtype=np.int16
    ).reshape(ny, nx).astype(np.float64)

    forest = biome_raw == 7
    grassland = biome_raw == 9
    diff = suit_povo_livre.astype(np.float64) - suit_circulo.astype(np.float64)
    changed = land & (diff != 0) & np.isfinite(diff)

    meta = {
        "resolution_m": [cs_x, cs_y],
        "biome_suitability_lut_povo_livre": {
            BIOME_NAMES[k]: BIOME_SUITABILITY_POVO_LIVRE[k] for k in range(1, len(BIOME_NAMES))
        },
        "decision": BIOME_NOTES["povo_livre_variant_decision"],
        "circulo_vs_povo_livre_comparison": {
            "temperate_forest": {
                "pct_of_land": float(100 * forest.sum() / land.sum()),
                "mean_annual_precip_mm": float(precip[forest].mean()),
                "mean_x_m": float(X[forest].mean()),
            },
            "lowland_steppe_grassland": {
                "pct_of_land": float(100 * grassland.sum() / land.sum()),
                "mean_annual_precip_mm": float(precip[grassland].mean()),
                "mean_x_m": float(X[grassland].mean()),
            },
            "precip_ratio_forest_over_grassland": float(
                precip[forest].mean() / precip[grassland].mean()
            ),
            "spatial_overlap_between_the_two_top_classes_pct": 0.0,
            "fraction_of_land_where_score_changes": float(changed.sum() / land.sum()),
            "mean_abs_change_among_changed_cells": float(np.abs(diff[changed]).mean()),
        },
        "biome_suitability_povo_livre_land_stats": {
            "mean": float(suit_povo_livre[land].mean()),
            "p5_p50_p95": [float(v) for v in np.percentile(suit_povo_livre[land], [5, 50, 95])],
        },
    }
    with open(f"{out}/tappa6_biome_povo_livre_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time() - t0:.1f}s")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
