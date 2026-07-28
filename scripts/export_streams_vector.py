#!/usr/bin/env python3
"""
Export Tappa 4's raster stream network as a smoothed vector (GeoJSON).

Reads `flow_direction_code` + `stream_mask` (+ `contributing_area_km2` /
`discharge_proxy_m3s` for per-reach attributes) from
`data/processed/hydrology/`, segments the network into reaches at
confluences, Chaikin-smooths each reach, and writes
`data/exports/streams.geojson` -- lightweight, web-ready, meant to be
committed (unlike `data/processed/`, which is gitignored).

No re-run of the hydrology pipeline: this is a pure export step over
already-computed rasters. See docs/decisions/04_tappa4_hydrology.md S11.

    python scripts/export_streams_vector.py [--hydro-dir DIR] [--out PATH]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hydrology.vectorize import (                     # noqa: E402
    cell_to_xy,
    chaikin_smooth,
    segment_reaches,
    simplify_rdp,
    strahler_order,
)

XMIN, YMAX, RES_M = -65000.0, 80000.0, 30.0
CHAIKIN_ITERATIONS = 4
SIMPLIFY_TOLERANCE_M = 15.0   # half a cell; removes Chaikin's own redundant
                               # near-collinear points, not real shape detail


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hydro-dir", default="data/processed/hydrology")
    ap.add_argument("--out", default="data/exports/streams.geojson")
    args = ap.parse_args()

    t0 = time.time()
    hdir = Path(args.hydro_dir)
    codes = np.load(hdir / "flow_direction_code.npy")
    stream = np.load(hdir / "stream_mask.npy")
    area_km2 = np.load(hdir / "contributing_area_km2.npy")
    discharge = np.load(hdir / "discharge_proxy_m3s.npy")
    ny, nx = codes.shape
    print(f"[{time.time()-t0:5.1f}s] loaded rasters {ny}x{nx}")

    reaches, node_kind = segment_reaches(codes, stream)
    print(f"[{time.time()-t0:5.1f}s] segmented into {len(reaches)} reaches "
          f"({sum(1 for v in node_kind.values() if v=='head')} heads, "
          f"{sum(1 for v in node_kind.values() if v=='confluence')} confluences)")

    orders = strahler_order(reaches, node_kind)
    print(f"[{time.time()-t0:5.1f}s] Strahler order assigned, max order "
          f"{max(orders.values())}")

    area_flat = area_km2.ravel()
    disch_flat = discharge.ravel()

    features = []
    for i, reach in enumerate(reaches):
        xy = cell_to_xy(reach, nx, XMIN, YMAX, RES_M)
        smoothed = chaikin_smooth(xy, iterations=CHAIKIN_ITERATIONS)
        smoothed = simplify_rdp(smoothed, SIMPLIFY_TOLERANCE_M)
        features.append({
            "type": "Feature",
            "properties": {
                "reach_id": i,
                "strahler_order": orders[i],
                "n_cells": int(len(reach)),
                "max_contributing_area_km2": round(float(area_flat[reach].max()), 4),
                "mean_discharge_proxy_m3s": round(float(disch_flat[reach].mean()), 4),
                "max_discharge_proxy_m3s": round(float(disch_flat[reach].max()), 4),
            },
            "geometry": {
                "type": "LineString",
                "coordinates": [[round(float(x), 1), round(float(y), 1)] for x, y in smoothed],
            },
        })

    fc = {
        "type": "FeatureCollection",
        "crs_note": (
            "Coordinates are in the project's custom 'Fictional World LCC' "
            "CRS (metres), NOT WGS84 lon/lat -- same convention as "
            "data/input/*.geojson. Assign the CRS manually in QGIS "
            "(Layer Properties -> Source -> Assigned CRS) after loading; "
            "see config/parameters.yml for the PROJ4 string."
        ),
        "generated_by": "scripts/export_streams_vector.py",
        "smoothing": f"Chaikin corner-cutting, {CHAIKIN_ITERATIONS} iterations",
        "features": features,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(fc))
    print(f"[{time.time()-t0:5.1f}s] wrote {out} ({len(features)} features, "
          f"{out.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
