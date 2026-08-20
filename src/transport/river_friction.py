"""
Tappa 9 -- river-crossing friction, a THIRD friction layer stacked
multiplicatively on top of lithology (Tappa 8 S8f) and biome (Tappa 9 S2).

Added in direct response to Nico's review of the first road-network pass:
that pass gave zero cost to crossing a river of any size, because
`cost_distance.py`'s graph has never known anything about hydrology (only
elevation and land/sea). This module fixes that gap the same way every
other friction layer in this project gets built: a per-cell multiplier,
combined multiplicatively with whatever else is already stacked, grounded
in the project's own existing data (`data/exports/streams.geojson`, Tappa 4's
Strahler-ordered stream vector, 1-6) rather than inventing a new hydrology
model.

**Only HIGHER-order streams get a multiplier, per Nico's own framing** ("at
least those of higher Strahler order should be considered"). Order
`MAJOR_STREAM_MIN_STRAHLER_ORDER` (4) and up -- the top 3 of 6 Strahler
classes on this domain, `1847 / 17678` reaches (~10.4%) -- are treated as
needing a real bridge; orders 1-3 (headwater/tributary streams, `n_cells`
in the low tens per reach in the source vector) are left unpenalized as
fordable on foot/by cart, the same simplification real rural road networks
make (a farm track fords a creek; it bridges a river).

**RASTERIZATION, not a pre-existing raster.** `data/processed/hydrology/`
has a native-30m `stream_mask.npy`, but it carries no per-cell Strahler
order (order lives only on the vector's per-REACH attributes in
`streams.geojson`), so it can't answer "which cells are >= order 4" on its
own. This module densifies each qualifying reach's LineString to
sub-cell spacing (same `_densify_polyline`-style approach
`terrain/skeleton.py` already uses for ridge/zone geometry -- nearest-cell
snapping, not exact point-to-segment distance, an acceptable approximation
at 120 m when spacing is well under a cell width) and marks every 120 m
cell it passes through, directly at working resolution -- no native-30m
intermediate needed.

**RIVER_CROSSING_FRICTION -- UNREVIEWED first-pass estimates, pending
Nico's sign-off, not written to `config/parameters.yml`** (same status as
every other friction table in this project): scaled by order within the
"major" band itself, not a single flat number, on the reasoning that a
order-4 reach and a order-6 trunk river are not the same crossing --
order 6 carries the domain's largest discharge (`max_discharge_proxy_m3s`
up to ~103 m3/s in Tappa 4's own summary) and is the least likely to have
an easy natural crossing point.

**REVISED same day, superseded as the OPERATIVE model (kept below, not
retracted).** Nico asked directly whether Strahler order can tell you a
river's actual width/depth -- answer: only indirectly (order is a
topological proxy for discharge, and it's discharge, via real hydraulic-
geometry relations, that determines channel size -- see
`transport/river_geometry.py`'s full writeup). That module estimates a
per-reach width from the discharge this project ALREADY computed
(`max_discharge_proxy_m3s`, Tappa 4), rather than binning by order. This
exposed a real weakness in the order-only table above: some order-4
reaches on this domain carry as little as 0.11 m3/s (vs. the class's own
7.85 m3/s mean) -- a reach that crosses the order-3/4 topological
threshold right after one small tributary joins, not a physically large
river -- and the flat per-order table penalized it exactly as hard as a
typical order-4 crossing. `river_friction_multiplier_from_width` (below)
replaces the discrete 3-bucket lookup with a continuous function of
estimated width, ANCHORED to this table's own two calibrated severity
judgements (0.75 at the order-4 class's mean-discharge-implied width,
0.45 at the order-6 anchor's width) so the "typical" case for each old
class is unchanged and only the WITHIN-class spread is now differentiated
-- see that function's own docstring for the exact derivation.
`river_friction_multiplier` (order-based) and `RIVER_CROSSING_FRICTION`
are kept below as the documented predecessor, not deleted -- same
discipline `network.py`'s `add_redundant_edges` history uses.
"""
from __future__ import annotations

import json
import math

import numpy as np

__all__ = [
    "MAJOR_STREAM_MIN_STRAHLER_ORDER",
    "RIVER_CROSSING_FRICTION",
    "rasterize_major_streams",
    "river_friction_multiplier",
    "RIVER_FRICTION_WIDTH_ANCHORS",
    "river_friction_multiplier_from_width",
]

MAJOR_STREAM_MIN_STRAHLER_ORDER = 4

# Keyed on strahler_order (only orders >= MAJOR_STREAM_MIN_STRAHLER_ORDER are
# ever consulted -- rasterize_major_streams never marks a lower-order cell in
# the first place, so this table doesn't need entries below 4).
RIVER_CROSSING_FRICTION: dict[int, float] = {
    4: 0.75,  # smallest "major" class -- a real bridge, but a narrower one;
              # mildest penalty of the three, same "anchor near the low end"
              # convention the biome/lithology tables use for their most
              # common/least severe class
    5: 0.60,  # a substantial trunk tributary -- real accounts of unbridged
              # river-flat travel describe this as a major deviation
              # (detour to a ford/crossing point), not a minor slowdown
    6: 0.45,  # the domain's largest rivers by discharge -- harshest
              # penalty in this table on purpose, though still nowhere near
              # as severe as e.g. biome friction's Permanent-Snow-&-Ice
              # 0.35, since a real bridge, once built, restores normal speed
              # at the crossing POINT itself; what this number actually
              # prices in is the cost of detouring to find/use it, at the
              # single-cell resolution this grid can represent
}


def _densify_polyline(coords: np.ndarray, max_spacing_m: float) -> np.ndarray:
    """Same approach as `terrain/skeleton.py`'s own `_densify_polyline` --
    linear resampling so consecutive points are at most `max_spacing_m`
    apart, an acceptable stand-in for exact point-to-segment distance when
    `max_spacing_m` is well under the target cell size (here, well under
    120 m)."""
    out = [coords[0]]
    for i in range(len(coords) - 1):
        p0, p1 = coords[i], coords[i + 1]
        seg_len = np.hypot(*(p1 - p0))
        if seg_len <= max_spacing_m:
            out.append(p1)
            continue
        n_steps = int(np.ceil(seg_len / max_spacing_m))
        for s in range(1, n_steps + 1):
            out.append(p0 + (p1 - p0) * (s / n_steps))
    return np.array(out)


def rasterize_major_streams(
    streams_geojson_path: str,
    xmin: float,
    ymax: float,
    cellsize_x: float,
    cellsize_y: float,
    shape: tuple[int, int],
    min_order: int = MAJOR_STREAM_MIN_STRAHLER_ORDER,
) -> tuple[np.ndarray, np.ndarray]:
    """Rasterize every stream reach with `strahler_order >= min_order` onto
    a (ny, nx) grid at the given origin/cellsize -- returns
    (major_stream_mask, strahler_order_grid): a bool mask (True = a
    qualifying reach passes through this cell) and an int8 grid carrying
    the HIGHEST qualifying order present in each cell (0 where
    major_stream_mask is False), so `river_friction_multiplier` can apply
    per-order values rather than one flat number.

    `xmin`/`ymax`/`cellsize_x`/`cellsize_y` must match the target grid's own
    convention exactly (row 0 = north/ymax, same as every other raster in
    this pipeline) -- callers pass the same domain constants used to build
    the road-network's own 120 m grid.
    """
    ny, nx = shape
    mask = np.zeros(shape, dtype=bool)
    order_grid = np.zeros(shape, dtype=np.int8)
    spacing_m = min(cellsize_x, cellsize_y) / 2.0

    with open(streams_geojson_path) as f:
        features = json.load(f)["features"]

    n_reaches = 0
    for feat in features:
        order = feat["properties"].get("strahler_order", 0)
        if order < min_order:
            continue
        n_reaches += 1
        coords = np.array(feat["geometry"]["coordinates"], dtype=np.float64)
        dense = _densify_polyline(coords, spacing_m)
        cols = np.clip(((dense[:, 0] - xmin) / cellsize_x).astype(np.int64), 0, nx - 1)
        rows = np.clip(((ymax - dense[:, 1]) / cellsize_y).astype(np.int64), 0, ny - 1)
        mask[rows, cols] = True
        # keep the HIGHEST order seen at any cell shared by multiple reaches
        # (confluences) -- a cell at a confluence of two order-4 reaches
        # immediately downstream IS the order-5 reach in real Strahler
        # ordering, but we don't re-derive that here; taking max across
        # whatever reaches actually rasterize onto the cell is a reasonable,
        # cheap stand-in, not a claim of exact hydrological correctness at
        # confluence cells specifically
        order_grid[rows, cols] = np.maximum(order_grid[rows, cols], order)

    return mask, order_grid


def river_friction_multiplier(
    order_grid: np.ndarray,
    multipliers: dict[int, float] = RIVER_CROSSING_FRICTION,
    default: float = 1.0,
) -> np.ndarray:
    """Map each cell's major-stream order (0 = no qualifying stream, from
    `rasterize_major_streams`) to a friction multiplier. Cells with order 0
    (the overwhelming majority of the grid) get `default` (neutral) --
    same fallback convention as `biome_friction_multiplier` and
    `travel_friction_multiplier`."""
    out = np.full(order_grid.shape, default, dtype=np.float32)
    for order, mult in multipliers.items():
        out[order_grid == order] = mult
    return out


# --- Width-based friction (supersedes the order-bucket model above) --------
#
# Two anchor points, DELIBERATELY reused from RIVER_CROSSING_FRICTION's own
# already-made severity judgements rather than picking new numbers from
# scratch:
#   - order 4's mean discharge (7.8521 m3/s, this domain's actual value)
#     implies an estimated width of ~12.42 m (river_geometry.estimate_width_m)
#     -> friction 0.75, the OLD flat value every order-4 reach got. This
#     anchors the new curve so the "typical" order-4 crossing is UNCHANGED,
#     not re-litigated.
#   - order 6's anchor discharge (103.128 m3/s, ANCHOR_DISCHARGE_M3S) is by
#     construction ANCHOR_WIDTH_M (45.0 m) wide -> friction 0.45, the OLD
#     harshest value.
# A power law friction(W) = 1 - k*W^p is the unique monotonic curve through
# (0, 1.0) and both anchor points (checked directly at each import via the
# two anchors below, not assumed to fit) -- floored at 0.45 beyond 45 m
# (this domain's largest known reach; no data to extrapolate past it) and
# ceilinged at 1.0 (a reach with ~0 discharge is not an obstacle at all).
RIVER_FRICTION_WIDTH_ANCHORS = {
    "width_m_at_typical_order4": 12.4247,  # river_geometry.estimate_width_m(7.8521)
    "friction_at_typical_order4": 0.75,
    "width_m_at_order6_anchor": 45.0,      # == river_geometry.ANCHOR_WIDTH_M
    "friction_at_order6_anchor": 0.45,
}
_w1 = RIVER_FRICTION_WIDTH_ANCHORS["width_m_at_typical_order4"]
_f1 = RIVER_FRICTION_WIDTH_ANCHORS["friction_at_typical_order4"]
_w2 = RIVER_FRICTION_WIDTH_ANCHORS["width_m_at_order6_anchor"]
_f2 = RIVER_FRICTION_WIDTH_ANCHORS["friction_at_order6_anchor"]
_WIDTH_FRICTION_P = math.log((1 - _f2) / (1 - _f1)) / math.log(_w2 / _w1)
_WIDTH_FRICTION_K = (1 - _f1) / (_w1 ** _WIDTH_FRICTION_P)
_WIDTH_FRICTION_FLOOR = _f2  # 0.45 -- no data past the order-6 anchor width


def river_friction_multiplier_from_width(
    width_grid: np.ndarray,
    default: float = 1.0,
) -> np.ndarray:
    """Continuous version of `river_friction_multiplier`: maps each cell's
    ESTIMATED WIDTH (metres, 0 = no qualifying stream -- from
    `river_geometry.rasterize_major_stream_discharge` +
    `river_geometry.estimate_width_m`, already scoped to
    `strahler_order >= MAJOR_STREAM_MIN_STRAHLER_ORDER` upstream) to a
    friction multiplier via `friction(W) = 1 - _WIDTH_FRICTION_K *
    W ** _WIDTH_FRICTION_P`, floored at `_WIDTH_FRICTION_FLOOR` (0.45,
    this domain's largest known river) and ceilinged at 1.0.

    Cells with width 0 (no qualifying stream) get `default` -- same
    fallback convention as every other friction function in this project.
    A cell WITH a qualifying stream but a genuinely small estimated width
    (this domain has order-4-by-topology reaches under 2 m wide, see this
    module's docstring) now correctly gets a friction close to 1.0 instead
    of the old flat 0.75 -- the entire point of this revision.
    """
    w = np.asarray(width_grid, dtype=np.float64)
    friction = 1.0 - _WIDTH_FRICTION_K * np.power(np.maximum(w, 0.0), _WIDTH_FRICTION_P)
    friction = np.clip(friction, _WIDTH_FRICTION_FLOOR, 1.0)
    out = np.where(w > 0, friction, default).astype(np.float32)
    return out
