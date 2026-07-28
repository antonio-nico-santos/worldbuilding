"""
Raster stream network -> smoothed vector polylines.

`stream_mask` (from run_tappa4_hydrology.py) only says "is this cell part of
a channel" -- it has no notion of individual reaches, so it can't be
smoothed or drawn as a set of lines directly: every tributary's cells would
need to be walked into one tangled, self-overlapping path. This module
segments the network into reaches first (splitting at confluences, exactly
the same node/edge structure a real hydrographic vector layer -- e.g.
NHDPlus -- uses), THEN smooths each reach independently with Chaikin's
corner-cutting algorithm.

No re-run of the priority-flood needed: `flow_direction_code` (D8, ESRI
convention) already encodes the full receiver graph losslessly, so the
receiver array is reconstructed from it directly.
"""

from __future__ import annotations

import numpy as np

_DELTA = {
    1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
    16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1),
}


def receiver_from_codes(codes: np.ndarray) -> np.ndarray:
    """Reconstruct the flat receiver-index array from ESRI-convention D8
    direction codes. Cells with code 0 (outlets) or whose target falls
    outside the grid receive themselves (self-receiving, same convention
    `flow.py` uses for seeds)."""
    ny, nx = codes.shape
    n = ny * nx
    idx = np.arange(n)
    row, col = np.divmod(idx, nx)
    receiver = idx.copy()
    for code, (dr, dc) in _DELTA.items():
        m = (codes.ravel() == code)
        r, c = row[m] + dr, col[m] + dc
        valid = (r >= 0) & (r < ny) & (c >= 0) & (c < nx)
        target = np.where(valid, r * nx + c, idx[m])
        receiver[m] = np.where(valid, target, idx[m])
    return receiver


def segment_reaches(codes: np.ndarray, stream_mask: np.ndarray):
    """Split the stream network into reaches: maximal chains of stream
    cells between two "nodes" (a channel head, a confluence, or an
    outlet), matching a standard hydrographic network data model.

    Returns `(reaches, node_kind)`:
      - `reaches`: list of arrays of flat cell indices, head-to-outlet
        order, each including both of its endpoint nodes (so adjacent
        reaches share exactly one vertex at a confluence -- topologically
        connected, not just visually coincident).
      - `node_kind`: dict flat-index -> "head" | "confluence" | "outlet",
        for callers that want to label endpoints (e.g. for Strahler order).
    """
    ny, nx = codes.shape
    n = ny * nx
    stream_flat = stream_mask.ravel()
    receiver = receiver_from_codes(codes)

    # in-degree within the stream network only: how many OTHER stream
    # cells drain directly into this one.
    in_degree = np.zeros(n, dtype=np.int32)
    stream_idx = np.flatnonzero(stream_flat)
    for i in stream_idx:
        r = receiver[i]
        if r != i and stream_flat[r]:
            in_degree[r] += 1

    is_head = stream_flat & (in_degree == 0)
    is_confluence = stream_flat & (in_degree >= 2)
    starts = np.flatnonzero(is_head | is_confluence)

    node_kind = {}
    for i in np.flatnonzero(is_head):
        node_kind[int(i)] = "head"
    for i in np.flatnonzero(is_confluence):
        node_kind[int(i)] = "confluence"

    reaches = []
    for start in starts:
        # a confluence cell is itself the first vertex of every downstream
        # reach that begins there (one per cell, by construction, since
        # each cell has exactly one receiver) -- walk forward from it.
        chain = [int(start)]
        cur = start
        while True:
            r = receiver[cur]
            if r == cur:
                node_kind[int(r)] = "outlet"
                break
            chain.append(int(r))
            if not stream_flat[r]:
                node_kind[int(r)] = "outlet"
                break
            if is_confluence[r] or is_head[r]:
                # reached another node -- reach ends here (inclusive),
                # the next reach starting at `r` is a separate list entry
                break
            cur = r
        if len(chain) >= 2:
            reaches.append(np.array(chain, dtype=np.int64))
    return reaches, node_kind


def strahler_order(reaches: list[np.ndarray], node_kind: dict) -> dict:
    """Strahler order per reach, computed over the (small) reach graph --
    tractable here specifically because reach segmentation already
    collapsed millions of cells down to thousands of reaches. Standard
    rule: a head reach is order 1; at a confluence, two reaches of equal
    order N produce order N+1 downstream, otherwise the max order carries
    through unchanged.
    """
    # map: starting node (head/confluence) -> reach; ending node -> list of
    # incoming reaches (for confluences, exactly the reaches that end there)
    by_start = {int(r[0]): i for i, r in enumerate(reaches)}
    incoming: dict[int, list[int]] = {}
    for i, r in enumerate(reaches):
        incoming.setdefault(int(r[-1]), []).append(i)

    order = [None] * len(reaches)

    def resolve(i: int) -> int:
        if order[i] is not None:
            return order[i]
        start = int(reaches[i][0])
        if node_kind.get(start) == "head":
            order[i] = 1
        else:
            ups = [resolve(j) for j in incoming.get(start, [])]
            if not ups:
                order[i] = 1
            elif len(ups) == 1:
                order[i] = ups[0]
            else:
                m = max(ups)
                order[i] = m + 1 if ups.count(m) >= 2 else m
        return order[i]

    for i in range(len(reaches)):
        resolve(i)
    return {i: order[i] for i in range(len(reaches))}


def chaikin_smooth(points: np.ndarray, iterations: int = 4) -> np.ndarray:
    """Chaikin corner-cutting: each iteration replaces every edge with two
    points 1/4 and 3/4 along it, converging to a smooth curve. Endpoints
    are preserved exactly (both left un-cut) so adjacent reaches still
    share a coincident vertex at confluences after smoothing."""
    pts = points.astype(np.float64)
    for _ in range(iterations):
        if len(pts) < 3:
            break
        new_pts = [pts[0]]
        for i in range(len(pts) - 1):
            p0, p1 = pts[i], pts[i + 1]
            new_pts.append(0.75 * p0 + 0.25 * p1)
            new_pts.append(0.25 * p0 + 0.75 * p1)
        new_pts.append(pts[-1])
        pts = np.array(new_pts)
    return pts


def simplify_rdp(points: np.ndarray, tolerance: float) -> np.ndarray:
    """Ramer-Douglas-Peucker simplification. Applied AFTER Chaikin
    smoothing, not instead of it: Chaikin needs the original jagged
    vertices to know where the curve should bend, but its own output
    packs in far more near-collinear points than the resulting smooth
    curve needs to render correctly -- this strips those back out.
    `tolerance` in the same units as `points` (metres here); endpoints are
    always preserved, so reaches still meet exactly at confluences.

    Iterative (explicit stack), not recursive: a near-straight reach --
    and this DEM has some genuinely long straight stretches, see
    01_tappa1_terrain.md's ridge-shelf geometry -- is exactly the
    worst case for naive recursive RDP (O(n) recursion depth instead of
    the usual O(log n)), which would blow past Python's default recursion
    limit on this project's longer post-Chaikin reaches (up to ~4,200
    points before simplification)."""
    pts = points.astype(np.float64)
    n = len(pts)
    if n < 3:
        return pts
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        i0, i1 = stack.pop()
        if i1 - i0 < 2:
            continue
        start, end = pts[i0], pts[i1]
        line = end - start
        line_len2 = np.dot(line, line)
        seg = pts[i0 + 1:i1]
        if line_len2 == 0:
            d = np.linalg.norm(seg - start, axis=1)
        else:
            t = np.clip(((seg - start) @ line) / line_len2, 0.0, 1.0)
            proj = start + t[:, None] * line
            d = np.linalg.norm(seg - proj, axis=1)
        idx_local = int(np.argmax(d))
        if d[idx_local] > tolerance:
            idx = i0 + 1 + idx_local
            keep[idx] = True
            stack.append((i0, idx))
            stack.append((idx, i1))
    return pts[keep]


def cell_to_xy(flat_idx: np.ndarray, nx: int, xmin: float, ymax: float, res_m: float) -> np.ndarray:
    """Flat cell indices -> (x, y) cell-center world coordinates, in the
    project's CRS (metres) -- same convention as raster_io's map info."""
    row, col = np.divmod(flat_idx, nx)
    x = xmin + (col + 0.5) * res_m
    y = ymax - (row + 0.5) * res_m
    return np.stack([x, y], axis=1)
