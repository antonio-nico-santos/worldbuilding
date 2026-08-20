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

LITHOLOGY FRICTION (optional, added for the transport lithology-cost
multiplier follow-up -- see geomorphology/lithology.py's
`travel_friction_multiplier` docstring for the actual per-class values,
citations, and honest caveats): `build_cost_graph`'s new `friction_multiplier`
parameter is a generic per-cell 0-1-ish multiplier applied to LAND-LAND edge
speed only (sea/boat edges are untouched -- lithology is only defined on
land). Deliberately kept lithology-agnostic HERE -- this module still only
knows about elevation and land/sea, exactly like before; it has no import of
anything from geomorphology, matching how it already doesn't know about lakes
either (that exclusion happens in the calling run script via `land_mask`, not
here). The caller assembles whatever friction field it wants (lithology,
later maybe biome or vegetation) and passes the array in; passing None
(default) reproduces the EXACT previous behaviour bit-for-bit -- this keeps
every already-locked Tappa 6 result (site_selection's 17/17 placement, 0
violations) reproducible without this parameter ever being touched, since
lithology friction is a NEW, separate, not-yet-signed-off layer, not a retro-
active correction to what's already committed.

An edge's friction is the arithmetic mean of its two endpoint cells' own
multipliers (a property of the ground being crossed, not of travel
direction -- unlike slope, which already IS directional and stays that way).

Performance note: the full graph has ~1.45M nodes (this world's 120 m grid)
and ~8 directed edges per node (~11.6M edges). A single-source Dijkstra over
that (scipy.sparse.csgraph.dijkstra) is what cost_distance_from_source runs
once per already-placed Circulo during greedy placement -- `limit_hours`
lets scipy stop expanding once accumulated cost exceeds the limit, which
matters a lot here since we only ever need to know whether a candidate is
>= some threshold, not its exact cost past that point.

PREDECESSOR EXTRACTION / PATH RECONSTRUCTION (Tappa 9, additive only --
`cost_distance_from_source` above is untouched, byte-for-byte, so every
already-locked Tappa 6/8 result stays reproducible). Tappa 6/8 only ever
needed COSTS (is this candidate >= some hour threshold from an already-
placed site) -- Tappa 9's road-network work needs the actual PATH a route
would follow, not just its length. scipy's dijkstra already computes a
shortest-path tree as a side effect of finding the distances; the only
thing missing was asking for it (`return_predecessors=True`) and having a
function to walk it back into a cell sequence. See
`cost_distance_from_source_with_predecessors` and `reconstruct_path` below.
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
    "cost_distance_from_source_with_predecessors",
    "reconstruct_path",
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
    dem_120m: np.ndarray,
    land_mask: np.ndarray,
    cellsize_km: float,
    friction_multiplier: np.ndarray | None = None,
    sea_mode: str = "boat",
) -> csr_matrix:
    """Directed 8-connected cost graph (edge weight = travel time in HOURS)
    over the full ny*nx grid (land and sea both included as nodes -- sea
    nodes are only ever traversed at BOAT_SPEED_KMH, never scored as a
    Circulo site elsewhere in this pipeline, but they must exist as graph
    nodes for a sea crossing to be a path at all).

    `friction_multiplier`: optional, same shape as `land_mask` -- a per-cell
    multiplier (<=1.0 slows travel, >1.0 would speed it up, though nothing in
    this project currently produces a >1.0 value) applied to LAND-LAND edge
    speed only, via the mean of the edge's two endpoint multipliers. None
    (default) reproduces prior behaviour exactly -- see module docstring's
    "LITHOLOGY FRICTION" section.

    `sea_mode` (Tappa 9 road-network fix, additive -- default UNCHANGED from
    every already-locked Tappa 6/8 call site): "boat" (default) reproduces
    the original behaviour byte-for-byte -- any edge touching a non-land
    cell (ocean OR lake; this module has never distinguished the two, both
    are simply "not land_mask") costs `edge_dist_km / BOAT_SPEED_KMH`,
    appropriate for Tappa 6's isochrone/tier-distance siting use case, where
    a real (if slow) boat crossing is a legitimate way to reach a site.
    "impassable" is NEW: any edge touching a non-land cell gets cost=inf
    (the edge is effectively removed from the graph) -- for a *road*
    network specifically, a route can't be "built" across open water or a
    lake the way it can be walked/boated for a distance check; those
    crossings belong to a real ferry connection (a separate, not-yet-built
    Tappa 9 sub-build), not a road. Using "impassable" can leave some
    node-pairs at infinite cost (no all-land path exists between them) --
    that is not a bug, it is the correct signal that those two points need
    a ferry, not a road; the caller must handle disconnected components
    (see `src/transport/network.py`'s minimum-spanning-FOREST, not -tree).
    """
    if sea_mode not in ("boat", "impassable"):
        raise ValueError(f"sea_mode must be 'boat' or 'impassable', got {sea_mode!r}")
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
        if friction_multiplier is not None:
            src_fric = friction_multiplier[r0:r1, c0:c1]
            dst_fric = friction_multiplier[sr0:sr1, sc0:sc1]
            edge_friction = 0.5 * (src_fric + dst_fric)
            speed_land = np.maximum(speed_land * edge_friction, _MIN_SPEED_KMH)
        cost_land = edge_dist_km / speed_land
        cost_sea = (
            edge_dist_km / BOAT_SPEED_KMH
            if sea_mode == "boat"
            else np.inf  # "impassable" -- broadcasts fine, scipy's csr_matrix
                         # just stores a finite-looking inf entry; dijkstra
                         # treats it as an edge that's never worth taking
        )

        both_land = src_land & dst_land
        cost = np.where(both_land, cost_land, cost_sea)

        src_id = node_id[r0:r1, c0:c1]
        dst_id = node_id[sr0:sr1, sc0:sc1]

        # "impassable" mode: drop inf-cost entries from the sparse matrix
        # entirely rather than storing them -- an explicit stored inf would
        # still occupy space and, worse, would make the node LOOK connected
        # to csr_matrix/scipy machinery that treats "has a stored entry" and
        # "reachable" as different questions in some code paths (e.g.
        # `.nnz`-based checks); dropping the edge is the same "no edge here"
        # a plain non-land cell that's never adjacent to anything already
        # represents.
        if sea_mode == "impassable":
            keep = np.isfinite(cost.ravel())
            rows_all.append(src_id.ravel()[keep])
            cols_all.append(dst_id.ravel()[keep])
            data_all.append(cost.ravel()[keep])
        else:
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


def cost_distance_from_source_with_predecessors(
    graph: csr_matrix, row: int, col: int, shape: tuple[int, int], limit_hours: float = np.inf
) -> tuple[np.ndarray, np.ndarray]:
    """Same single-source Dijkstra as `cost_distance_from_source`, but also
    returns the predecessor array scipy computes as a side effect of
    building the shortest-path tree -- needed to reconstruct an actual
    route, not just read off its cost. Returns (dist_grid, predecessors),
    where `predecessors` is scipy's raw flat-node-id array (shape (ny*nx,),
    -9999 = "never reached this node from the source" INCLUDING the source
    cell itself, which is scipy's own convention, not a bug) -- pass it to
    `reconstruct_path` rather than indexing it directly.
    """
    ny, nx = shape
    src_id = row * nx + col
    dist, predecessors = dijkstra(
        graph, directed=True, indices=[src_id], limit=limit_hours, return_predecessors=True
    )
    return dist[0].reshape(ny, nx), predecessors[0]


def reconstruct_path(
    predecessors: np.ndarray,
    shape: tuple[int, int],
    source_row: int,
    source_col: int,
    target_row: int,
    target_col: int,
) -> list[tuple[int, int]] | None:
    """Walk `predecessors` (from `cost_distance_from_source_with_predecessors`,
    same source cell it was computed from) back from the target to the
    source, returning an ordered list of (row, col) cells SOURCE -> TARGET
    (scipy's own chain runs target -> source, reversed here to the more
    useful direction for building a route polyline).

    Returns None if the target was never reached within whatever
    `limit_hours` the predecessor array was computed with -- checked by
    confirming the walked-back chain actually terminates AT the source
    (scipy's -9999 sentinel means "no predecessor recorded", which is true
    both for a genuinely unreached node and, correctly, for the source cell
    itself; the two are told apart here by checking which node the walk
    stopped at, not just whether it stopped).
    """
    ny, nx = shape
    source_id = source_row * nx + source_col
    node = target_row * nx + target_col
    if node == source_id:
        return [(source_row, source_col)]

    path_ids = [node]
    seen = {node}
    while predecessors[node] != -9999:
        node = int(predecessors[node])
        if node in seen:
            # defensive only -- a valid shortest-path tree is acyclic by
            # construction, this should never trigger
            return None
        seen.add(node)
        path_ids.append(node)

    if node != source_id:
        return None  # target unreachable (or outside limit_hours)

    path_ids.reverse()
    return [(pid // nx, pid % nx) for pid in path_ids]
