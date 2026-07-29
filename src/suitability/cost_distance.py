"""
Tappa 6 -- cost-distance (walking/boat TIME, not straight-line km) between
Circulo sites. Built after Nico's own follow-up idea: barriers (mountains,
open sea) should count for MORE effective separation than the same straight-
line km over flat terrain, and flat open ground should count for LESS -- the
"greek city-state" intuition is really about travel TIME, not map distance.

Two travel modes, chosen per edge of an 8-connected grid graph spanning the
WHOLE domain (land AND sea -- see BOAT_SPEED_KMH below, Nico chose to model
sea crossing as a real (if slow) travel mode, not an impassable barrier):

- LAND-LAND edges: Tobler's hiking function, a standard, widely-used
  empirical hiking-speed model from slope:
      speed_kmh = 6 * exp(-3.5 * abs(slope_ratio + 0.05))
  where slope_ratio = rise/run (signed, tangent of the slope angle) IN THE
  DIRECTION OF TRAVEL for that edge -- this makes the cost graph genuinely
  DIRECTED (uphill costs more than downhill over the same edge, and Tobler's
  peak speed is actually a gentle ~2.86% downhill grade, not flat ground).
  At slope_ratio=0 this gives ~5.04 km/h, matching the well-known "average
  human walking pace on flat ground" calibration point Tobler's function is
  built around.
- Any edge touching a non-land (sea) cell: a flat BOAT_SPEED_KMH, regardless
  of direction or the other endpoint's elevation. HONEST PLACEHOLDER: 6.0
  km/h, deliberately modest (comparable to a brisk flat walking pace) --
  simple rowing/small-sail boats, no motors, no wind/current/tide modelling
  (none of that data exists in this project). A land<->sea edge (launching/
  landing) also just uses the boat speed for the whole edge -- launch/land
  friction is not modelled separately. Revisit this number if the scenario
  ever specifies actual watercraft technology.

  IMPORTANT SIDE EFFECT, worth knowing before reading any cost-distance
  result: 6.0 km/h is actually slightly FASTER than Tobler's flat-ground
  walking speed (~5.04 km/h), and open water has no slope penalty at all.
  That means, AS PARAMETERISED, water is not really a barrier in this
  model -- it is closer to a slight travel SHORTCUT relative to walking
  over any terrain with real relief. This is a direct, deliberate
  consequence of Nico choosing "mar com custo de barco" over "mar como
  barreira forte" -- but it does mean the original mountains-and-sea-both-
  separate-more intuition only actually holds for MOUNTAINS in this
  version; sea crossings are cheap. If sea is meant to still act as some
  kind of barrier (while remaining crossable), lowering BOAT_SPEED_KMH
  well below ~5 km/h would restore that -- not done here since Nico chose
  boat-mode specifically, not the barrier option.

Elevation for the slope term comes from a 120 m DEM, block-MEAN downsampled
from the native 30 m DEM (mean, not block-max like slope_pct_120m elsewhere
in this pipeline -- here we want a realistic average edge slope between
neighbouring 120 m cells, not a hazard-conservative worst-case within a
cell -- see terrain_metrics.py's module docstring for why block-max is
used THERE instead).

Performance note: the full graph has ~1.45M nodes (this world's 120 m grid)
and ~8 directed edges per node (~11.6M edges). A single-source Dijkstra over
that (scipy.sparse.csgraph.dijkstra) is what cost_distance_from_source runs
once per already-placed Circulo during greedy placement -- `limit_hours`
lets scipy stop expanding once accumulated cost exceeds the limit, which
matters a lot here since we only ever need to know whether a candidate is
>= some threshold, not its exact cost past that point.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

__all__ = [
    "BOAT_SPEED_KMH",
    "tobler_speed_kmh",
    "build_cost_graph",
    "cost_distance_from_source",
]

BOAT_SPEED_KMH = 6.0
_MIN_SPEED_KMH = 0.05  # floor so cost=dist/speed never overflows on near-vertical edges

_DIRECTIONS = [
    (-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
    (-1, -1, 2**0.5), (-1, 1, 2**0.5), (1, -1, 2**0.5), (1, 1, 2**0.5),
]


def tobler_speed_kmh(slope_ratio: np.ndarray) -> np.ndarray:
    """Tobler's hiking function. `slope_ratio` = signed rise/run (tangent) in
    the direction of travel. Returns km/h, floored at _MIN_SPEED_KMH."""
    speed = 6.0 * np.exp(-3.5 * np.abs(slope_ratio + 0.05))
    return np.maximum(speed, _MIN_SPEED_KMH)


def build_cost_graph(
    dem_120m: np.ndarray, land_mask: np.ndarray, cellsize_km: float
) -> csr_matrix:
    """Directed 8-connected cost graph (edge weight = travel time in HOURS)
    over the full ny*nx grid (land and sea both included as nodes -- sea
    nodes are only ever traversed at BOAT_SPEED_KMH, never scored as a
    Circulo site elsewhere in this pipeline, but they must exist as graph
    nodes for a sea crossing to be a path at all).
    """
    ny, nx = land_mask.shape
    node_id = np.arange(ny * nx, dtype=np.int64).reshape(ny, nx)
    dist_m = cellsize_km * 1000.0

    rows_all, cols_all, data_all = [], [], []
    for dy, dx, factor in _DIRECTIONS:
        edge_dist_km = cellsize_km * factor
        edge_dist_m = dist_m * factor

        r0, r1 = max(0, -dy), ny - max(0, dy)
        c0, c1 = max(0, -dx), nx - max(0, dx)
        sr0, sr1 = max(0, dy), ny - max(0, -dy)
        sc0, sc1 = max(0, dx), nx - max(0, -dx)

        src_elev = dem_120m[r0:r1, c0:c1]
        dst_elev = dem_120m[sr0:sr1, sc0:sc1]
        src_land = land_mask[r0:r1, c0:c1]
        dst_land = land_mask[sr0:sr1, sc0:sc1]

        slope_ratio = (dst_elev - src_elev) / edge_dist_m
        speed_land = tobler_speed_kmh(slope_ratio)
        cost_land = edge_dist_km / speed_land
        cost_sea = edge_dist_km / BOAT_SPEED_KMH

        both_land = src_land & dst_land
        cost = np.where(both_land, cost_land, cost_sea)

        src_id = node_id[r0:r1, c0:c1]
        dst_id = node_id[sr0:sr1, sc0:sc1]

        rows_all.append(src_id.ravel())
        cols_all.append(dst_id.ravel())
        data_all.append(cost.ravel())

    rows = np.concatenate(rows_all)
    cols = np.concatenate(cols_all)
    data = np.concatenate(data_all)
    n = ny * nx
    return csr_matrix((data, (rows, cols)), shape=(n, n))


def cost_distance_from_source(
    graph: csr_matrix, row: int, col: int, shape: tuple[int, int], limit_hours: float = np.inf
) -> np.ndarray:
    """Single-source Dijkstra (hours) from (row, col) to every cell, reshaped
    to `shape`. Cells not reached within `limit_hours` come back as np.inf --
    treat that as "definitely far enough", not a computation failure.
    """
    ny, nx = shape
    src_id = row * nx + col
    dist = dijkstra(graph, directed=True, indices=[src_id], limit=limit_hours)[0]
    return dist.reshape(ny, nx)
