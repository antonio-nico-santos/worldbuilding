"""
Tappa 8 -- real (DEM-grounded) ridge crest extraction, replacing the
authored-polyline buffer that made lithology look artificially smooth
(Nico's catch).

v1 of this module thresholded `ridge_accum_cells` (inverted-DEM flow
accumulation, see run_tappa8_ridge_extraction.py -- the same trick as
Tappa 4's stream extraction, just fed -DEM) by a single percentile over
each ridge's ENTIRE shelf-reach corridor. That was wrong in a way that
only showed up at full-domain scale, not in the small validation crop:
a corridor sized to the shelf reach (up to ~24 km either side of the
authored line, millions of cells) is mountainous terrain generally, full
of real minor ridges and spurs that have nothing to do with the named
massif -- even a strict top-1% cut over that many cells still pulled in a
wide, diffuse scatter of those secondary ridges, not a single coherent
line. Schist area exploded to 80% of land as a direct result (every
scattered point got its own 15 km buffer).

FIX: arclength-bin the authored line (same technique Tappa 7 already used
for outpost placement along the Spine/South Branch), and within each
bin's own LOCAL search window pick the single highest-accumulation land
cell. This keeps the search anchored to where the named ridge actually
is at every point along its length, instead of searching the whole
corridor at once -- one real point per bin, forming a coherent sequence
that follows the true terrain while staying tied to the ridge's own
identity.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from terrain.skeleton import _densify_polyline


def extract_real_crest(
    ridge,
    ridge_coords: np.ndarray,
    dem: np.ndarray,
    ridge_accum_cells: np.ndarray,
    land_mask: np.ndarray,
    xmin: float,
    ymax: float,
    cellsize_m: float,
    bin_spacing_m: float = 1000.0,
    search_radius_km: float | None = None,
):
    """`ridge_coords`: the ridge's own raw authored LineString coordinates
    (not densified -- this function does its own, coarser, arclength
    binning). `search_radius_km` defaults to the ridge's own `falloff_km`
    (the real crest shouldn't need to wander farther from the authored
    line than the schist radius itself, or the corridor concept stops
    meaning anything).

    Returns (tree, crest_elev, info) -- a cKDTree over one real point per
    bin (skipping bins where the local window has no land), the matching
    elevations, and a diagnostic dict.
    """
    radius_m = (search_radius_km or ridge.falloff_km) * 1000.0
    bins = _densify_polyline(ridge_coords, bin_spacing_m)

    ny, nx = dem.shape
    pts_xy = []
    pts_elev = []
    n_empty_bins = 0

    for bx, by in bins:
        col_c = int((bx - xmin) / cellsize_m)
        row_c = int((ymax - by) / cellsize_m)
        r_cells = int(np.ceil(radius_m / cellsize_m))
        r0, r1 = max(0, row_c - r_cells), min(ny, row_c + r_cells + 1)
        c0, c1 = max(0, col_c - r_cells), min(nx, col_c + r_cells + 1)
        if r0 >= r1 or c0 >= c1:
            n_empty_bins += 1
            continue

        sub_accum = ridge_accum_cells[r0:r1, c0:c1]
        sub_land = land_mask[r0:r1, c0:c1]
        # also enforce the circular radius within the bounding box, not
        # just a square window
        rr, cc = np.meshgrid(np.arange(r0, r1), np.arange(c0, c1), indexing="ij")
        yy_sub = ymax - rr * cellsize_m
        xx_sub = xmin + cc * cellsize_m
        within_r = np.hypot(xx_sub - bx, yy_sub - by) <= radius_m
        valid = sub_land & within_r
        if not valid.any():
            n_empty_bins += 1
            continue

        masked_accum = np.where(valid, sub_accum, -1.0)
        local_argmax = np.unravel_index(np.argmax(masked_accum), masked_accum.shape)
        real_row, real_col = r0 + local_argmax[0], c0 + local_argmax[1]
        pts_xy.append((xmin + real_col * cellsize_m, ymax - real_row * cellsize_m))
        pts_elev.append(float(dem[real_row, real_col]))

    info = {
        "n_bins": len(bins),
        "n_points": len(pts_xy),
        "n_empty_bins": n_empty_bins,
        "search_radius_km": radius_m / 1000.0,
        "bin_spacing_m": bin_spacing_m,
    }
    if len(pts_xy) < 2:
        return None, None, {**info, "status": "too few real-crest points"}

    tree = cKDTree(np.array(pts_xy))
    return tree, np.array(pts_elev, dtype=np.float32), {**info, "status": "ok"}


def extract_real_crest_cross_section(
    ridge,
    ridge_coords: np.ndarray,
    dem: np.ndarray,
    land_mask: np.ndarray,
    xmin: float,
    ymax: float,
    cellsize_m: float,
    bin_spacing_m: float = 1000.0,
    cross_section_halfwidth_m: float = 2000.0,
    sample_spacing_m: float | None = None,
):
    """Option A -- "cross-section snapping" (Nico's second catch: Option B's
    real-crest points, even radius-limited to 2.5 km, still read as
    artificial once zoomed in -- straight 1 km chords between points found
    by a 2D window search, a "low-poly" look, not organic).

    At each arclength station along the AUTHORED line, this samples the DEM
    along a short line segment PERPENDICULAR to the local line direction
    (a literal cross-section) and picks the single highest-elevation land
    cell on that segment -- i.e. "the summit of the spine surface" the way
    Nico originally described it, read directly off elevation rather than
    inferred from flow-accumulation topology.

    This is a different, tighter constraint than Option B's radius-limited
    2D window: the search is confined to a 1D perpendicular strip anchored
    to each station's own position and local direction, so a station's
    point is geometrically tied to its neighbors by construction (it can't
    jump sideways along the ridge the way a 2D window search could) --
    ordering/coherence falls out of the method itself rather than needing a
    radius cap to prevent it.

    `cross_section_halfwidth_m=2000.0` is a first-pass placeholder (same
    order of magnitude as Option B's radius fix, chosen to keep the crest
    from wandering implausibly far off the authored line -- not
    independently calibrated). `sample_spacing_m` defaults to `cellsize_m`
    (full native-grid resolution along the cross-section).

    Returns (tree, crest_elev, info) -- same contract as
    `extract_real_crest`.
    """
    sample_spacing_m = sample_spacing_m or cellsize_m
    bins = np.asarray(_densify_polyline(ridge_coords, bin_spacing_m), dtype=np.float64)
    n = len(bins)
    ny, nx = dem.shape

    n_side = max(1, int(round(cross_section_halfwidth_m / sample_spacing_m)))
    offsets = np.arange(-n_side, n_side + 1) * sample_spacing_m  # (n_offsets,)

    pts_xy = []
    pts_elev = []
    n_empty_bins = 0

    for i in range(n):
        if n == 1:
            tangent = np.array([1.0, 0.0])
        elif i == 0:
            tangent = bins[1] - bins[0]
        elif i == n - 1:
            tangent = bins[-1] - bins[-2]
        else:
            tangent = bins[i + 1] - bins[i - 1]
        tnorm = np.hypot(tangent[0], tangent[1])
        tangent = tangent / tnorm if tnorm > 1e-6 else np.array([1.0, 0.0])
        perp = np.array([-tangent[1], tangent[0]])

        sample_xy = bins[i][None, :] + offsets[:, None] * perp[None, :]
        cols = np.round((sample_xy[:, 0] - xmin) / cellsize_m).astype(int)
        rows = np.round((ymax - sample_xy[:, 1]) / cellsize_m).astype(int)
        in_bounds = (rows >= 0) & (rows < ny) & (cols >= 0) & (cols < nx)
        if not in_bounds.any():
            n_empty_bins += 1
            continue
        rows_v, cols_v = rows[in_bounds], cols[in_bounds]
        land_v = land_mask[rows_v, cols_v]
        if not land_v.any():
            n_empty_bins += 1
            continue

        elev_v = dem[rows_v, cols_v]
        elev_masked = np.where(land_v, elev_v, -np.inf)
        best = int(np.argmax(elev_masked))
        real_row, real_col = int(rows_v[best]), int(cols_v[best])
        pts_xy.append((xmin + real_col * cellsize_m, ymax - real_row * cellsize_m))
        pts_elev.append(float(dem[real_row, real_col]))

    info = {
        "n_bins": n,
        "n_points": len(pts_xy),
        "n_empty_bins": n_empty_bins,
        "cross_section_halfwidth_m": cross_section_halfwidth_m,
        "sample_spacing_m": sample_spacing_m,
        "bin_spacing_m": bin_spacing_m,
        "method": "cross-section snapping (Option A)",
    }
    if len(pts_xy) < 2:
        return None, None, {**info, "status": "too few real-crest points"}

    tree = cKDTree(np.array(pts_xy))
    return tree, np.array(pts_elev, dtype=np.float32), {**info, "status": "ok"}


def extract_real_crest_network(
    ridge,
    ridge_coords: np.ndarray,
    dem: np.ndarray,
    ridge_accum_cells: np.ndarray,
    land_mask: np.ndarray,
    xmin: float,
    ymax: float,
    cellsize_m: float,
    corridor_halfwidth_km: float = 2.0,
    accum_percentile_in_corridor: float = 90.0,
    dense_spacing_m: float = 200.0,
):
    """Nico's third catch: even Option A's cross-section snapping (v3) is
    still fundamentally anchored to the AUTHORED line's own station
    structure (one point every bin_spacing_m, geometry effectively
    "applying DEM values onto the authored line"). Requested instead: use
    the DEM's OWN ridge line as the actual source geometry, and only use
    the authored line to say WHICH part of the terrain belongs to this
    named massif.

    Method: `ridge_accum_cells` (inverted-DEM flow accumulation -- see
    run_tappa8_ridge_extraction.py; high values there are real ridge/divide
    cells, the exact mirror of how high accumulation marks real stream
    cells in Tappa 4) is thresholded by its OWN in-corridor percentile,
    inside a narrow corridor around the authored line
    (`corridor_halfwidth_km`, NOT the full falloff/shelf reach -- that
    wide-corridor version is what failed the first time this stage tried a
    percentile cut, see this module's own top-of-file history). Every land
    cell in the corridor at or above that percentile is kept as a real
    crest cell -- the union is the DEM's actual (branching, irregular)
    ridge network inside this named massif's own footprint, not a sparse
    point-per-bin sequence.

    Both `corridor_halfwidth_km=2.0` and `accum_percentile_in_corridor=90.0`
    are first-pass placeholders (2 km matches the scale that already killed
    cross-ridge jump artifacts in v2/v3; 90th percentile is an unreviewed
    starting guess, not calibrated against how thick a real crest network
    should read at this map's scale).

    Returns (tree, crest_elev, info) -- tree built over EVERY matching
    cell (could be hundreds to thousands of points, not one per bin), plus
    each point's own DEM elevation and a diagnostic dict. Restricted to a
    local bounding box around the ridge's own extent (+corridor) for
    memory/time, not the full domain.
    """
    dense = _densify_polyline(ridge_coords, dense_spacing_m)
    corridor_tree = cKDTree(dense)

    corridor_m = corridor_halfwidth_km * 1000.0
    x0 = float(dense[:, 0].min() - corridor_m)
    x1 = float(dense[:, 0].max() + corridor_m)
    y0 = float(dense[:, 1].min() - corridor_m)
    y1 = float(dense[:, 1].max() + corridor_m)

    ny, nx = dem.shape
    col0 = max(0, int((x0 - xmin) / cellsize_m))
    col1 = min(nx, int(np.ceil((x1 - xmin) / cellsize_m)))
    row0 = max(0, int((ymax - y1) / cellsize_m))
    row1 = min(ny, int(np.ceil((ymax - y0) / cellsize_m)))

    info_base = {
        "corridor_halfwidth_km": corridor_halfwidth_km,
        "accum_percentile_in_corridor": accum_percentile_in_corridor,
        "bbox_rows": [row0, row1], "bbox_cols": [col0, col1],
    }
    if row0 >= row1 or col0 >= col1:
        return None, None, {**info_base, "status": "empty bounding box"}

    sub_dem = dem[row0:row1, col0:col1]
    sub_land = land_mask[row0:row1, col0:col1]
    sub_accum = ridge_accum_cells[row0:row1, col0:col1]

    rr, cc = np.meshgrid(np.arange(row0, row1), np.arange(col0, col1), indexing="ij")
    xx_sub = xmin + cc * cellsize_m
    yy_sub = ymax - rr * cellsize_m
    sub_xy = np.column_stack([xx_sub.ravel(), yy_sub.ravel()])

    dist_m, _ = corridor_tree.query(sub_xy)
    dist_m = dist_m.reshape(sub_dem.shape)
    in_corridor = sub_land & (dist_m <= corridor_m)

    if not in_corridor.any():
        return None, None, {**info_base, "status": "no land in corridor"}

    threshold = np.percentile(sub_accum[in_corridor], accum_percentile_in_corridor)
    crest_mask = in_corridor & (sub_accum >= threshold)

    n_crest = int(crest_mask.sum())
    info = {
        **info_base,
        "n_corridor_land_cells": int(in_corridor.sum()),
        "accum_threshold": float(threshold),
        "n_crest_cells": n_crest,
    }
    if n_crest < 2:
        return None, None, {**info, "status": "too few real-crest cells"}

    crest_rows, crest_cols = np.where(crest_mask)
    pts_xy = np.column_stack([
        xmin + (col0 + crest_cols) * cellsize_m,
        ymax - (row0 + crest_rows) * cellsize_m,
    ])
    pts_elev = sub_dem[crest_rows, crest_cols].astype(np.float32)

    tree = cKDTree(pts_xy)
    return tree, pts_elev, {**info, "status": "ok"}


def extract_real_crest_network_local(
    ridge,
    ridge_coords: np.ndarray,
    dem: np.ndarray,
    ridge_accum_cells: np.ndarray,
    land_mask: np.ndarray,
    xmin: float,
    ymax: float,
    cellsize_m: float,
    corridor_halfwidth_km: float = 2.0,
    window_spacing_m: float = 1000.0,
    window_halflength_m: float = 1500.0,
    accum_percentile_in_window: float = 75.0,
):
    """Fixes a real bug found in `extract_real_crest_network`'s GLOBAL
    corridor-wide percentile: checked coverage along the Spine's arclength
    after the fact and found real-crest cells only in the first 10.2 km of
    a 64.4 km line -- 84% of the ridge had ZERO qualifying cells, not a
    plausible geological saddle. Root cause: `ridge_accum_cells` is an
    UPSTREAM-CONTRIBUTING-AREA proxy (inverted-DEM flow accumulation), and
    contributing area along a ridge trunk is NOT uniform the way it is
    along a real stream -- it's dominated by wherever the inverted network
    happens to converge most, not "how ridge-like this point is" locally.
    A single percentile over the WHOLE corridor is entirely captured by
    that one high-convergence stretch, starving everywhere else. This is
    the exact same category of mistake as this module's very first failed
    attempt (global percentile over the whole shelf-reach corridor) --
    just one level less wide this time, so it wasn't as visually obvious
    until arclength coverage was actually checked.

    FIX: same medicine as before -- make the percentile LOCAL. At each
    arclength station (`window_spacing_m` apart), take the top
    `accum_percentile_in_window` of cells within a small along-line window
    (`window_halflength_m` either side of the station) intersected with the
    lateral corridor, and OR the result into a running crest mask. Every
    station gets its own locally-relevant cutoff, so coverage no longer
    depends on how the whole ridge's accumulation happens to be
    distributed -- guaranteed continuous coverage the same way the earlier
    per-bin methods guaranteed it, but keeping this function's actual
    improvement over those (a multi-cell, non-single-point swath per
    station, preserving real texture instead of one argmax pixel).

    All four numeric knobs are first-pass placeholders, not calibrated.
    """
    dense = _densify_polyline(ridge_coords, 200.0)
    seg_vec = np.diff(dense, axis=0)
    seg_len = np.hypot(seg_vec[:, 0], seg_vec[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])
    total_len = cum[-1]

    corridor_m = corridor_halfwidth_km * 1000.0
    x0 = float(dense[:, 0].min() - corridor_m)
    x1 = float(dense[:, 0].max() + corridor_m)
    y0 = float(dense[:, 1].min() - corridor_m)
    y1 = float(dense[:, 1].max() + corridor_m)

    ny, nx = dem.shape
    col0 = max(0, int((x0 - xmin) / cellsize_m))
    col1 = min(nx, int(np.ceil((x1 - xmin) / cellsize_m)))
    row0 = max(0, int((ymax - y1) / cellsize_m))
    row1 = min(ny, int(np.ceil((ymax - y0) / cellsize_m)))

    info_base = {
        "corridor_halfwidth_km": corridor_halfwidth_km,
        "window_spacing_m": window_spacing_m,
        "window_halflength_m": window_halflength_m,
        "accum_percentile_in_window": accum_percentile_in_window,
        "total_line_length_km": total_len / 1000.0,
    }
    if row0 >= row1 or col0 >= col1:
        return None, None, {**info_base, "status": "empty bounding box"}

    sub_dem = dem[row0:row1, col0:col1]
    sub_land = land_mask[row0:row1, col0:col1]
    sub_accum = ridge_accum_cells[row0:row1, col0:col1]

    rr, cc = np.meshgrid(np.arange(row0, row1), np.arange(col0, col1), indexing="ij")
    xx_sub = xmin + cc * cellsize_m
    yy_sub = ymax - rr * cellsize_m

    crest_mask = np.zeros(sub_dem.shape, dtype=bool)
    n_stations = max(1, int(np.ceil(total_len / window_spacing_m)) + 1)
    n_empty_windows = 0
    n_ok_windows = 0

    for s in range(n_stations):
        arclen = min(s * window_spacing_m, total_len)
        station_xy = dense[np.searchsorted(cum, arclen)]
        along_lo, along_hi = arclen - window_halflength_m, arclen + window_halflength_m

        # local along-line sub-window of the dense authored points, then
        # the corridor test is still done against the FULL dense line (not
        # just this window's points) so lateral distance stays exact
        dist_to_station = np.hypot(xx_sub - station_xy[0], yy_sub - station_xy[1])
        # along-line window approximated by straight-line distance to the
        # station point being within window_halflength_m + corridor_m --
        # a generous circular stand-in for the along-line band, cheap to
        # compute per station without re-projecting every cell onto the line
        local = sub_land & (dist_to_station <= (window_halflength_m + corridor_m))
        if not local.any():
            n_empty_windows += 1
            continue

        threshold = np.percentile(sub_accum[local], accum_percentile_in_window)
        window_crest = local & (sub_accum >= threshold)
        crest_mask |= window_crest
        n_ok_windows += 1

    # still enforce the true lateral corridor (the circular per-station
    # window above is a superset near corridor edges/ends)
    dense_tree = cKDTree(dense)
    sub_xy = np.column_stack([xx_sub.ravel(), yy_sub.ravel()])
    dist_to_line, _ = dense_tree.query(sub_xy)
    dist_to_line = dist_to_line.reshape(sub_dem.shape)
    crest_mask &= (dist_to_line <= corridor_m)

    n_crest = int(crest_mask.sum())
    info = {
        **info_base,
        "n_stations": n_stations, "n_ok_windows": n_ok_windows, "n_empty_windows": n_empty_windows,
        "n_crest_cells": n_crest,
    }
    if n_crest < 2:
        return None, None, {**info, "status": "too few real-crest cells"}

    crest_rows, crest_cols = np.where(crest_mask)
    pts_xy = np.column_stack([
        xmin + (col0 + crest_cols) * cellsize_m,
        ymax - (row0 + crest_rows) * cellsize_m,
    ])
    pts_elev = sub_dem[crest_rows, crest_cols].astype(np.float32)

    tree = cKDTree(pts_xy)
    return tree, pts_elev, {**info, "status": "ok"}


def height_normalized_falloff_m(crest_elev: np.ndarray, base_falloff_km: float,
                                 low_multiplier: float = 0.6, high_multiplier: float = 1.0):
    """Per-crest-point falloff, scaled by that point's own elevation
    relative to the REST OF THIS RIDGE's real crest points only (not a
    world-wide or cross-ridge normalization) -- the highest point on the
    ridge gets `high_multiplier * base_falloff_km` (schist band at full
    width), the lowest point gets `low_multiplier * base_falloff_km`
    (narrower band), linearly interpolated in between. Nico's explicit
    request: "1x base value on the highest point, 0.6x at the lowest point
    of the ridge." Both multipliers are exactly what was asked for, not
    independently calibrated beyond that.
    """
    elev = np.asarray(crest_elev, dtype=np.float64)
    lo, hi = float(elev.min()), float(elev.max())
    if hi - lo < 1e-6:
        t = np.ones_like(elev)
    else:
        t = (elev - lo) / (hi - lo)
    multiplier = low_multiplier + (high_multiplier - low_multiplier) * t
    return (base_falloff_km * 1000.0 * multiplier).astype(np.float32)
