"""
Tappa 9 -- per-reach river width/depth estimation from `streams.geojson`'s
already-computed discharge proxy, via downstream hydraulic geometry.

Added same day as `river_friction.py`'s first shipped version, in direct
response to Nico's own observation: Strahler order is a TOPOLOGICAL
classification (how many tributaries merge upstream), not a physical
measurement -- it correlates with channel size only indirectly, through
drainage area and therefore discharge. This domain's own data confirms the
proxy is imperfect: a handful of order-4 reaches carry as little as
0.11 m3/s (max_discharge_proxy_m3s), a tiny fraction of the 7.85 m3/s
order-4 MEAN -- almost certainly a short reach that crosses the order-3/4
topological threshold immediately after one small tributary joins, not a
reach that is actually physically large. `river_friction.py`'s original
per-order bucket (§6c) would have penalized that reach exactly as harshly
as a typical order-4 river; this module exists so the friction layer can
be driven by an estimate of the reach's actual physical size instead.

**PHYSICAL BASIS: Leopold & Maddock (1953) downstream hydraulic
geometry.** Channel width and depth scale with discharge as power laws:

    W = a * Q^b      (width, metres)
    D = c * Q^f      (depth, metres)

The EXPONENTS (b ~ 0.5 for width, f ~ 0.4 for depth) are among the most
widely reproduced results in fluvial geomorphology -- consistent across
regions and channel types because they follow from a simple physical
argument (a channel widens roughly with the square root of the flow it
must pass) more than from any one dataset. This module uses those two
textbook exponent values directly.

**The MULTIPLICATIVE COEFFICIENTS (a, c) are NOT a citable universal
constant** -- published values vary roughly 2-4x across regions, channel
substrate (bedrock vs. alluvial), and confinement (single-thread vs.
braided), because they capture local channel-forming history, not just
discharge. Rather than presenting a borrowed real-world coefficient as if
it were locally validated for THIS fictional world, this module derives
`a`/`c` by anchoring to ONE explicit, checkable magnitude judgement,
consistent with how every other "no ground truth, need a number" table in
this project has been built (see e.g. river_friction.py's own severity
picks): at this domain's largest known discharge (`ANCHOR_DISCHARGE_M3S`,
103 m3/s -- the actual max `max_discharge_proxy_m3s` across all order-6
reaches, not a round number), a confined, non-braided temperate gravel-bed
river of that size is plausibly `ANCHOR_WIDTH_M` = 45 m wide and
`ANCHOR_DEPTH_M` = 2.8 m deep -- ordinary magnitudes for a large single-
channel temperate river (NOT a braided system like a real South-Island-NZ
Rakaia/Waimakariri at flood, which run to hundreds of metres across
multiple channels; this domain's rivers are treated as single-thread,
consistent with `cost_distance.py` never having modelled channel braiding).
**This anchor is a calibration CHOICE, not a literature citation --
UNREVIEWED, same status as every other friction estimate in this
project.** Every other reach's width/depth is scaled from this one anchor
by the exponent relation, so the anchor's absolute value matters for
narrative/QGIS-symbology use, but NOT for `river_friction.py`'s relative
severity ranking between reaches (a ratio of two Q^b terms cancels `a`
out) -- see that module for how the anchor's uncertainty is kept from
leaking into the friction multiplier's operative behaviour.

**Discharge basis**: `max_discharge_proxy_m3s` (the reach's own downstream-
end value, not `mean_discharge_proxy_m3s`, which blends a reach's
upstream/downstream ends) -- same field `river_friction.py`'s own
docstring already anchors order 6's characterization on, for consistency.
`max_discharge_proxy_m3s` is itself flagged in `04_tappa4_hydrology.md` S3
as an UPPER BOUND (no infiltration/evapotranspiration loss, unrouted) --
so width/depth estimates inherit that same upper-bound character, not
independently.
"""
from __future__ import annotations

import json
import math

import numpy as np

from transport.river_friction import _densify_polyline

__all__ = [
    "WIDTH_EXPONENT",
    "DEPTH_EXPONENT",
    "ANCHOR_DISCHARGE_M3S",
    "ANCHOR_WIDTH_M",
    "ANCHOR_DEPTH_M",
    "estimate_width_m",
    "estimate_depth_m",
    "rasterize_major_stream_discharge",
    "annotate_reach_geometry",
]

# Leopold & Maddock (1953) downstream hydraulic geometry exponents --
# reproduced across fluvial geomorphology texts (e.g. Knighton, "Fluvial
# Forms and Processes"; Leopold, Wolman & Miller, "Fluvial Processes in
# Geomorphology") as representative midpoint values. See module docstring.
WIDTH_EXPONENT = 0.50   # b
DEPTH_EXPONENT = 0.40   # f

# Calibration anchor -- see module docstring. NOT a literature citation.
ANCHOR_DISCHARGE_M3S = 103.128  # this domain's actual max max_discharge_proxy_m3s (order 6)
ANCHOR_WIDTH_M = 45.0
ANCHOR_DEPTH_M = 2.8


def estimate_width_m(discharge_m3s):
    """W = ANCHOR_WIDTH_M * (Q / ANCHOR_DISCHARGE_M3S) ** WIDTH_EXPONENT.
    Works on scalars or numpy arrays. Discharge <= 0 -> width 0 (no channel)."""
    q = np.asarray(discharge_m3s, dtype=np.float64)
    q_safe = np.maximum(q, 0.0)
    w = ANCHOR_WIDTH_M * (q_safe / ANCHOR_DISCHARGE_M3S) ** WIDTH_EXPONENT
    return np.where(q > 0, w, 0.0)


def estimate_depth_m(discharge_m3s):
    """D = ANCHOR_DEPTH_M * (Q / ANCHOR_DISCHARGE_M3S) ** DEPTH_EXPONENT.
    Same conventions as `estimate_width_m`."""
    q = np.asarray(discharge_m3s, dtype=np.float64)
    q_safe = np.maximum(q, 0.0)
    d = ANCHOR_DEPTH_M * (q_safe / ANCHOR_DISCHARGE_M3S) ** DEPTH_EXPONENT
    return np.where(q > 0, d, 0.0)


def rasterize_major_stream_discharge(
    streams_geojson_path: str,
    xmin: float,
    ymax: float,
    cellsize_x: float,
    cellsize_y: float,
    shape: tuple[int, int],
    min_order: int,
) -> np.ndarray:
    """Rasterize `max_discharge_proxy_m3s` for every reach with
    `strahler_order >= min_order`, onto the same (ny, nx) grid convention
    `river_friction.py`'s `rasterize_major_streams` uses (row 0 = north,
    reuses that module's own `_densify_polyline` directly rather than a
    third copy of the same ~15-line algorithm). At a cell shared by more
    than one qualifying reach (a confluence), keeps the MAX discharge seen
    -- for discharge specifically (unlike order, see that function's own
    caveat) this is not an approximation: discharge strictly increases
    downstream in this project's tree-structured D8 drainage network, so
    the higher-discharge reach at a shared cell IS the physically correct
    value, not a stand-in.

    Returns a float32 grid, 0.0 where no qualifying reach passes through.
    """
    ny, nx = shape
    discharge_grid = np.zeros(shape, dtype=np.float32)
    spacing_m = min(cellsize_x, cellsize_y) / 2.0

    with open(streams_geojson_path) as f:
        features = json.load(f)["features"]

    for feat in features:
        p = feat["properties"]
        order = p.get("strahler_order", 0)
        if order < min_order:
            continue
        q = float(p.get("max_discharge_proxy_m3s", 0.0))
        coords = np.array(feat["geometry"]["coordinates"], dtype=np.float64)
        dense = _densify_polyline(coords, spacing_m)
        cols = np.clip(((dense[:, 0] - xmin) / cellsize_x).astype(np.int64), 0, nx - 1)
        rows = np.clip(((ymax - dense[:, 1]) / cellsize_y).astype(np.int64), 0, ny - 1)
        discharge_grid[rows, cols] = np.maximum(discharge_grid[rows, cols], q)

    return discharge_grid


def annotate_reach_geometry(streams_geojson_path: str, min_order: int) -> dict:
    """Read every reach with `strahler_order >= min_order` and return a
    GeoJSON FeatureCollection (same LineString geometry, unmodified) with
    `estimated_width_m`/`estimated_depth_m` properties added -- for QGIS
    symbology (variable-width river rendering) and narrative/reference use.
    Does NOT modify `data/exports/streams.geojson` itself (that file is
    Tappa 4's own locked output); the caller writes this to a separate
    Tappa-9-owned file.
    """
    with open(streams_geojson_path) as f:
        fc = json.load(f)

    out_features = []
    for feat in fc["features"]:
        p = feat["properties"]
        if p.get("strahler_order", 0) < min_order:
            continue
        q = float(p.get("max_discharge_proxy_m3s", 0.0))
        w = float(estimate_width_m(q))
        d = float(estimate_depth_m(q))
        new_props = dict(p)
        new_props["estimated_width_m"] = round(w, 2)
        new_props["estimated_depth_m"] = round(d, 2)
        out_features.append({
            "type": "Feature",
            "properties": new_props,
            "geometry": feat["geometry"],
        })

    return {
        "type": "FeatureCollection",
        "name": "major_stream_geometry",
        "crs": fc.get("crs"),
        "features": out_features,
    }
