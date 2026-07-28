"""
D8 flow routing: combined depression-filling + flow-direction assignment via
the Priority-Flood+Epsilon algorithm (Barnes, Lehman & Mulla 2014, "Priority-
flood: An optimal depression-filling and watershed-labeling algorithm for
digital elevation models"), plus a single reverse-topological pass for
(optionally weighted) flow accumulation.

Why this algorithm specifically, and why one pass does both jobs:

r.watershed's job is really two things GRASS keeps separate internally
(depression filling, then D8/flow accumulation) but which the priority-flood
algorithm naturally computes together. It is a Dijkstra-style flood fill: seed
a min-heap with every drainage sink (here: the ocean, plus the domain
boundary, since this DEM's land reaches all four edges -- see
02_tappa2_climate.md S4f), then repeatedly pop the lowest open cell and
"flood" into its unclosed neighbors. Two things fall out for free:

1. **Depression filling**: any neighbor at or below the popping cell's
   elevation is raised to `popped_elevation + epsilon` before being closed.
   This guarantees strictly-decreasing elevation along every flow path back
   to a sink -- no pits, and (because of the epsilon, not just filling to
   dead-flat) no true flats either, so there is no separate flat-resolution
   pass needed (the classic problem with naive priority-flood + a bare D8
   pass afterwards).
2. **Flow direction**: whichever cell was popped when a neighbor gets closed
   IS that neighbor's single receiver (the neighbor was only reachable by
   flooding in from there) -- D8 direction is a byproduct of the fill, not a
   second pass over the (now-filled) DEM.

The pop order itself is a valid reverse-topological order of the flow DAG:
every cell is popped strictly after its receiver (the cell floods in FROM its
receiver). So flow accumulation is a single O(n) pass over that order,
reversed, summing each cell's accumulated weight into its receiver's --
no sorting, no iteration to convergence.

No GRASS/numba/rasterio in this sandbox (same constraint as Tappa 1-3, see
noise.py); everything here is numpy plus a plain-Python heapq driving the
per-cell traversal (the one part that cannot be vectorized -- it is an
inherently sequential graph traversal). Validated at small scale first (see
`scripts/test_flow_small.py`) before committing to the full 5334x4334 domain.
"""

from __future__ import annotations

import heapq

import numpy as np

# 8-connected offsets: (drow, dcol), and the true cell-to-cell distance for
# each (used nowhere in the fill/direction itself -- D8 direction only needs
# ordering, not slope -- but returned per-cell for callers that want it, e.g.
# a future stream-power / slope-dependent step).
_NEIGHBORS = [
    (-1, -1), (-1, 0), (-1, 1),
    (0, -1),           (0, 1),
    (1, -1),  (1, 0),  (1, 1),
]


def priority_flood_d8(
    dem: np.ndarray,
    seed_mask: np.ndarray,
    epsilon: float = 1e-4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Depression-fill `dem` and assign D8 flow directions in one pass.

    `seed_mask` marks the cells flow is allowed to drain OUT of (ocean cells,
    domain-boundary land cells) -- these are pushed onto the queue first and
    never modified.

    Returns `(filled, receiver, pop_order)`:
      - `filled`: float64, same shape as `dem`, depression-filled (strictly
        decreasing along every flow path to a seed).
      - `receiver`: int64, flat index (row*nx+col) of the single downstream
        neighbor each cell drains to. Seed cells receive themselves (their
        own flat index) -- callers must check `receiver[i] == i` to detect
        an outlet.
      - `pop_order`: int64, flat indices in the order cells were closed --
        a valid reverse-topological order of the flow DAG (every cell
        appears strictly after its receiver). Consumed by
        `accumulate_flow`.
    """
    ny, nx = dem.shape
    n = ny * nx
    dem_flat = dem.astype(np.float64).ravel()
    filled = dem_flat.copy()
    receiver = np.full(n, -1, dtype=np.int64)
    closed = np.zeros(n, dtype=bool)

    seed_idx = np.flatnonzero(seed_mask.ravel())
    if seed_idx.size == 0:
        raise ValueError("seed_mask has no True cells -- nothing for flow to drain to")

    heap: list[tuple[float, int]] = []
    for i in seed_idx:
        closed[i] = True
        receiver[i] = i
        heap.append((filled[i], i))
    heapq.heapify(heap)

    pop_order = np.empty(n, dtype=np.int64)
    n_popped = 0

    while heap:
        elev, idx = heapq.heappop(heap)
        # stale entries: a cell can only be pushed once in this algorithm
        # (it's closed at push time, see below), so no staleness to filter --
        # heap holds exactly one entry per closed cell.
        pop_order[n_popped] = idx
        n_popped += 1
        row, col = divmod(idx, nx)
        for dr, dc in _NEIGHBORS:
            r, c = row + dr, col + dc
            if r < 0 or r >= ny or c < 0 or c >= nx:
                continue
            j = r * nx + c
            if closed[j]:
                continue
            if filled[j] <= elev:
                filled[j] = elev + epsilon
            closed[j] = True
            receiver[j] = idx
            heapq.heappush(heap, (filled[j], j))

    assert n_popped == n, f"priority flood only reached {n_popped}/{n} cells -- disconnected region?"
    return filled.reshape(ny, nx), receiver, pop_order


def label_basins(receiver: np.ndarray, pop_order: np.ndarray) -> np.ndarray:
    """Label every cell with the flat index of the outlet (seed) it
    ultimately drains to -- one forward pass over `pop_order` (seeds/roots
    first, by construction), the mirror image of `accumulate_flow`'s reverse
    pass: a cell's basin is only known once its receiver's basin is, and
    receivers are always popped (and therefore visited here) before their
    contributors.
    """
    n = receiver.size
    basin = np.empty(n, dtype=np.int64)
    for idx in pop_order:
        r = receiver[idx]
        basin[idx] = idx if r == idx else basin[r]
    return basin


# ESRI/ArcGIS-convention D8 direction codes (power-of-two, the most widely
# recognised encoding for a GIS audience): E=1, SE=2, S=4, SW=8, W=16,
# NW=32, N=64, NE=128. 0 marks an outlet (a cell that is its own receiver).
_DIRECTION_CODE = {
    (0, 1): 1, (1, 1): 2, (1, 0): 4, (1, -1): 8,
    (0, -1): 16, (-1, -1): 32, (-1, 0): 64, (-1, 1): 128,
}


def direction_codes(receiver: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """D8 flow direction as ESRI-convention codes, from the flat `receiver`
    index array. Outlets (self-receiving cells) get 0."""
    ny, nx = shape
    n = ny * nx
    idx = np.arange(n)
    row, col = np.divmod(idx, nx)
    rrow, rcol = np.divmod(receiver, nx)
    dr = np.clip(rrow - row, -1, 1)
    dc = np.clip(rcol - col, -1, 1)
    code = np.zeros(n, dtype=np.uint8)
    for (ddr, ddc), c in _DIRECTION_CODE.items():
        mask = (dr == ddr) & (dc == ddc) & (receiver != idx)
        code[mask] = c
    return code.reshape(ny, nx)


def accumulate_flow(receiver: np.ndarray, pop_order: np.ndarray, weight: np.ndarray) -> np.ndarray:
    """Flow accumulation: `weight` summed over each cell's full upstream
    contributing area, via one reverse pass over `pop_order` (see module
    docstring for why this order is exactly right, no sorting needed).
    """
    accum = weight.astype(np.float64).ravel().copy()
    n = accum.size
    recv = receiver
    # reverse pop order: most-upstream-processed-last cells first, so every
    # contributor has already added itself into `accum` before its receiver
    # is visited and, in turn, pushes its (now-complete) total downstream.
    for idx in pop_order[::-1]:
        r = recv[idx]
        if r != idx:
            accum[r] += accum[idx]
    return accum
