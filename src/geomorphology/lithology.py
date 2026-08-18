"""
Tappa 8 -- Geomorphology: lithology classification.

Direction was locked in Tappa 7 (docs/decisions/07_tappa7_regional_scenario.md
S1) as decision-only; this module is the first actual implementation. Four
rock classes, reusing existing Tappa 1 skeleton geometry rather than
authoring anything new, combined by a straightforward priority order (a
categorical paint order, NOT Tappa 1's `max`-compositing of continuous
elevation contributions -- see the precedence note below):

    volcanic > schist > greywacke/argillite > sedimentary basin fill

- schist (metamorphic core): within each ridge's own `falloff_km` of its
  axis -- the exact distance field Tappa 1 uses for `structure()`,
  recomputed here (not reused from the DEM pipeline directly, since Tappa 1
  never persisted a standalone distance-to-ridge raster, only the elevation
  it drives). Real analog: metamorphic grade increases toward the Alpine
  Fault in the Haast Schist / Southern Alps.
- greywacke/argillite (flank): beyond the schist radius, out to each
  ridge's own SHELF reach (`falloff_km * shelf_multiplier`). The
  multipliers are the ones `run_tappa1_terrain.py` actually passed to
  `generate_dem` for the v3 run on disk -- reproduced here as
  `SHELF_MULTIPLIERS`, NOT re-derived. "Island" is deliberately absent from
  that dict (it was never added when the ridge was authored), so it falls
  back to `RidgeField`'s own class default of 3.0 -- a real fact about the
  DEM this world's coastline was actually generated from, not a choice
  made in this stage.
- volcanic: the SW Island landmass, identified by connected-component
  labeling of the native-resolution land mask -- the same method Tappa 7
  already used and validated at 120 m (794.6 km2, 16.44 km minimum water
  gap to the mainland). Recomputed here at native 30 m per this stage's own
  resolution decision, not reused from the 120 m biome grid.
- sedimentary basin fill: the default label for anything the other three
  don't claim. A cell inside one of the six authored plateau/plains zone
  polygons is flagged `basin_fill_grounded` -- those six zones were
  checked directly against `stream_mask` in Tappa 7 S1 (mean distance to
  the nearest stream 0.25-0.72 km across all six) and are a real,
  citable Canterbury Plains analog. Everything else that falls to this
  class by elimination (background-noise coastline slivers, remote
  offshore islets, ridge terrain beyond every shelf reach) paints the
  same class value but has NO such grounding -- see the decision doc for
  why conflating the two would overstate how well-supported this class
  is everywhere it appears. This grounding split is this module's own
  addition, not something S1 asked for explicitly.

Polygon membership uses a HARD boundary (`Path.contains_points`), not the
`edge_transition_km`-blended smootherstep weight `ZoneField.blend_weight`
computes for elevation continuity in Tappa 1 -- a rock-type boundary is a
structural fact, not something that should fade smoothly the way an
elevation target does. This is a deliberate divergence from how Tappa 1
uses the same zone polygons, not an oversight.
"""

from __future__ import annotations

import numpy as np

from terrain.skeleton import build_ridge_fields, build_zone_fields, load_geojson

# Reproduced verbatim from run_tappa1_terrain.py -- NOT re-derived here.
# "Island" is absent on purpose (see module docstring): it inherits
# RidgeField's own default (3.0), a genuine fact about the DEM already on
# disk, not a Tappa 8 decision.
SHELF_MULTIPLIERS = {
    "Spine": 1.6,
    "North branch (Big Brother)": 1.3,
    "West Branch (Little Brother)": 1.3,
    "South Branch": 1.3,
}
DEFAULT_SHELF_MULTIPLIER = 3.0

LITHOLOGY_CLASSES = {
    0: "ocean",
    1: "sedimentary_basin_fill",
    2: "greywacke_argillite",
    3: "schist",
    4: "volcanic",
}
CLASS_OCEAN = 0
CLASS_BASIN_FILL = 1
CLASS_GREYWACKE = 2
CLASS_SCHIST = 3
CLASS_VOLCANIC = 4


def load_ridges(path: str):
    return build_ridge_fields(
        load_geojson(path),
        shelf_multipliers=SHELF_MULTIPLIERS,
        default_shelf_multiplier=DEFAULT_SHELF_MULTIPLIER,
    )


def load_zones(path: str):
    return build_zone_fields(load_geojson(path))


def _grid_xy(xmin, xmax, ymin, ymax, ny, nx):
    """Cell-center coordinates, row 0 = north (same convention as the DEM
    array and every other raster in this pipeline). Uses the domain
    xmin/ymax directly rather than the array's own slightly-overhanging
    true top edge (Tappa 2 S2's ~20 m overhang) -- negligible next to the
    km-scale falloff/shelf distances this module compares against."""
    y = ymax - (np.arange(ny) + 0.5) * ((ymax - ymin) / ny)
    x = xmin + (np.arange(nx) + 0.5) * ((xmax - xmin) / nx)
    xx, yy = np.meshgrid(x, y)
    return xx, yy


def identify_landmasses(land_mask: np.ndarray, cellsize_m: float):
    """Connected-component label the land mask (8-connectivity) and return
    (labeled, mainland_label, sw_island_label, component_areas_km2).

    Same method Tappa 7 already validated at 120 m (mainland 8,879 km2,
    second component 794.6 km2 = the SW Island, 16.44 km min water gap,
    then 207 much smaller islets totalling 230.5 km2) -- recomputed here at
    native 30 m. Picks the two largest components by area as mainland/SW
    island; asserts the second component's area is in the right ballpark
    rather than assuming the ranking transfers cleanly across a 4x
    resolution change.
    """
    structure = np.ones((3, 3), dtype=int)
    labeled, n = ndimage_label(land_mask, structure)
    areas_cells = np.bincount(labeled.ravel())
    areas_cells[0] = 0  # background/ocean
    order = np.argsort(areas_cells)[::-1]
    cell_km2 = (cellsize_m / 1000.0) ** 2
    areas_km2 = areas_cells * cell_km2
    mainland_label = int(order[0])
    sw_island_label = int(order[1])
    return labeled, mainland_label, sw_island_label, areas_km2


def ndimage_label(mask, structure):
    from scipy import ndimage

    return ndimage.label(mask, structure=structure)


def classify(
    dem: np.ndarray,
    ridges,
    zones,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    row_chunk_size: int = 500,
    real_crest_trees: dict | None = None,
    real_crest_falloff_m: dict | None = None,
    warp_dx: np.ndarray | None = None,
    warp_dy: np.ndarray | None = None,
):
    """`warp_dx`/`warp_dy`: optional (ny, nx) meter displacement fields,
    applied ONLY to the query coordinates used for ridge distance/falloff
    (schist/greywacke boundary noise), NOT to the plateau-zone membership
    test below -- Nico's "apply noise at the end" request, a domain-warp of
    the query points rather than the source geometry itself, so it composes
    with whichever real-crest source (authored line, Option A, or the DEM
    ridge-network) is in play. See `boundary_noise.py` for how these are
    generated -- an independent value-noise implementation for this stage,
    NOT a reuse of Tappa 1's actual domain_warp (that code lives in
    terrain.generate/terrain.erosion, which were never staged into this
    session).

    `real_crest_trees`: optional {ridge.name: cKDTree} override -- when
    given for a ridge, its REAL DEM-grounded crest tree (from
    `real_crest.extract_real_crest`/`extract_real_crest_cross_section`) is
    queried instead of the ridge's own authored-polyline `tree`. See
    real_crest.py for why (Nico's catch: a constant-width buffer around a
    hand-authored line reads as an artificial "buffer," not a real
    metamorphic belt).

    `real_crest_falloff_m`: optional {ridge.name: np.ndarray} of PER-POINT
    falloff values (meters), aligned index-for-index with that ridge's
    `real_crest_trees[ridge.name]` tree points (e.g. from
    `real_crest.height_normalized_falloff_m`). When given for a ridge, each
    grid cell's falloff/shelf reach is taken from whichever real-crest
    point is nearest to it (the same point the distance itself is measured
    against), NOT the ridge's flat `falloff_km` constant -- this is Nico's
    height-normalized falloff request (highest crest point = full base
    falloff, lowest = a fraction of it). Falls back to the ridge's own
    scalar `falloff_km` when omitted, exactly as before."""
    """Returns a dict of (ny, nx) arrays: lithology (uint8, class codes
    above), basin_fill_grounded (bool), schist_grade (float32, 0-1, NaN
    outside the schist class), min_ridge_dist_km (float32, nearest-ridge
    distance regardless of class -- kept for QA/debug, not itself an
    output layer), plus landmass diagnostics.

    Chunked over rows (same reasoning as Tappa 1 S6: several ridges x
    several zones x a 23.1M-cell grid means several float64 temporaries
    that don't all need to coexist for the whole grid at once).
    """
    ny, nx = dem.shape
    land_mask = dem > 0

    labeled, mainland_label, sw_island_label, landmass_areas_km2 = identify_landmasses(
        land_mask, cellsize_m=(xmax - xmin) / nx
    )
    volcanic_mask = labeled == sw_island_label

    lithology = np.zeros((ny, nx), dtype=np.uint8)
    basin_fill_grounded = np.zeros((ny, nx), dtype=bool)
    schist_grade = np.full((ny, nx), np.nan, dtype=np.float32)
    min_ridge_dist_km = np.full((ny, nx), np.nan, dtype=np.float32)

    ln2 = np.log(2.0)

    for row0 in range(0, ny, row_chunk_size):
        row1 = min(row0 + row_chunk_size, ny)
        xx, yy = _grid_xy(xmin, xmax, ymin, ymax, ny, nx)
        xx = xx[row0:row1]
        yy = yy[row0:row1]
        xy = np.column_stack([xx.ravel(), yy.ravel()])

        if warp_dx is not None and warp_dy is not None:
            wdx = warp_dx[row0:row1].ravel()
            wdy = warp_dy[row0:row1].ravel()
            xy_ridge = xy + np.column_stack([wdx, wdy])
        else:
            xy_ridge = xy

        chunk_land = land_mask[row0:row1]
        chunk_volcanic = volcanic_mask[row0:row1]

        schist_any = np.zeros(xy.shape[0], dtype=bool)
        greywacke_any = np.zeros(xy.shape[0], dtype=bool)
        grade_max = np.zeros(xy.shape[0], dtype=np.float32)
        dist_min_km = np.full(xy.shape[0], np.inf, dtype=np.float32)

        for ridge in ridges:
            tree = (real_crest_trees or {}).get(ridge.name, ridge.tree)
            dist_m, nearest_idx = tree.query(xy_ridge, k=1)
            per_point_falloff = (real_crest_falloff_m or {}).get(ridge.name)
            if per_point_falloff is not None:
                falloff_m = per_point_falloff[nearest_idx]
            else:
                falloff_m = ridge.falloff_km * 1000.0
            shelf_m = falloff_m * ridge.shelf_multiplier
            schist_any |= dist_m <= falloff_m
            greywacke_any |= dist_m <= shelf_m
            grade = np.exp(-ln2 * (dist_m / falloff_m) ** 2).astype(np.float32)
            grade_max = np.maximum(grade_max, np.where(dist_m <= falloff_m, grade, 0.0))
            dist_min_km = np.minimum(dist_min_km, (dist_m / 1000.0).astype(np.float32))

        zone_any = np.zeros(xy.shape[0], dtype=bool)
        for zone in zones:
            if zone.feature_type != "plateau":
                continue
            zone_any |= zone.path.contains_points(xy)

        schist_any = schist_any.reshape(row1 - row0, nx)
        greywacke_any = greywacke_any.reshape(row1 - row0, nx)
        zone_any = zone_any.reshape(row1 - row0, nx)
        grade_max = grade_max.reshape(row1 - row0, nx)
        dist_min_km = dist_min_km.reshape(row1 - row0, nx)

        cls = np.where(
            chunk_volcanic,
            CLASS_VOLCANIC,
            np.where(schist_any, CLASS_SCHIST, np.where(greywacke_any, CLASS_GREYWACKE, CLASS_BASIN_FILL)),
        )
        cls = np.where(chunk_land, cls, CLASS_OCEAN).astype(np.uint8)

        lithology[row0:row1] = cls
        basin_fill_grounded[row0:row1] = zone_any & (cls == CLASS_BASIN_FILL)
        grade_out = np.where(cls == CLASS_SCHIST, grade_max, np.nan).astype(np.float32)
        schist_grade[row0:row1] = grade_out
        min_ridge_dist_km[row0:row1] = np.where(chunk_land, dist_min_km, np.nan)

    return {
        "lithology": lithology,
        "basin_fill_grounded": basin_fill_grounded,
        "schist_grade": schist_grade,
        "min_ridge_dist_km": min_ridge_dist_km,
        "labeled_landmasses": labeled,
        "mainland_label": mainland_label,
        "sw_island_label": sw_island_label,
        "landmass_areas_km2": landmass_areas_km2,
    }


def jade_eligible_mask(lithology: np.ndarray, schist_grade: np.ndarray, grade_percentile: float = 80.0):
    """Jade/pounamu (and the co-located gold-quartz-mica veins, see
    resources.py) sit in the highest metamorphic grades near the axis, not
    uniformly across the whole schist class. Threshold is a PERCENTILE of
    this world's own schist-grade distribution (top 20% by default), not
    an absolute grade value -- same "calibrate against this world's own
    data" move already used for the biome moisture terciles (Tappa 5) and
    the windward/leeward precipitation tercile (Tappa 3/4). 80th percentile
    is a first-pass placeholder, not visually reviewed yet -- flagged in
    the decision doc for QA the same way Tappa 6's slope thresholds were.
    """
    schist_mask = lithology == CLASS_SCHIST
    if not schist_mask.any():
        return np.zeros_like(schist_mask)
    threshold = np.nanpercentile(schist_grade[schist_mask], grade_percentile)
    return schist_mask & (schist_grade >= threshold)


def place_jade_pods(
    lithology: np.ndarray,
    schist_grade: np.ndarray,
    xx: np.ndarray,
    yy: np.ndarray,
    cellsize_m: float,
    n_pods: int = 10,
    min_separation_km: float = 5.0,
    radius_range_m: tuple[float, float] = (300.0, 800.0),
    grade_percentile: float = 80.0,
    candidate_pool_size: int = 500,
    seed: int = 13,
):
    """Real jade/pounamu deposits are NOT a smooth function of distance to
    a ridge axis -- if they were, prospecting would be trivial (Nico's
    catch). The source citation itself says "isolated pods", plural and
    discrete, not a painted band. This function keeps the percentile
    threshold as a SUITABILITY zone (geologically plausible ground, per
    `jade_eligible_mask`) and then stochastically places a small number of
    discrete pod footprints within it, weighted toward higher grade but
    NOT deterministic -- two runs with different seeds would put pods in
    different specific spots within the same plausible zone, which is the
    honest way to represent "we know roughly where to look, not exactly
    where the ore is."

    `n_pods=10`, the radius range, and `min_separation_km=5.0` are
    placeholders, not derived from a named citation -- flagged the same
    way as every other first-pass, undated numeric choice in this project
    (e.g. Tappa 3's sigma_day_c). `min_separation_km` reuses the same 5 km
    scale already used elsewhere in this project for "these should read as
    genuinely separate" buffers (Povo Silencioso hard buffer, Tappa 7's
    cove minimum spacing) rather than inventing a new one.

    Returns (pod_mask, suitable_mask, pod_centers_xy, pod_radii_m).
    """
    suitable_mask = jade_eligible_mask(lithology, schist_grade, grade_percentile)
    ys, xs = np.where(suitable_mask)
    if len(ys) == 0:
        return np.zeros_like(suitable_mask), suitable_mask, [], []

    rng = np.random.default_rng(seed)
    weights = schist_grade[ys, xs].astype(np.float64)
    weights = weights / weights.sum()

    pool_size = min(candidate_pool_size, len(ys))
    pool_idx = rng.choice(len(ys), size=pool_size, replace=False, p=weights)

    chosen_xy = []
    chosen_rc = []
    for i in pool_idx:
        cx, cy = xx[ys[i], xs[i]], yy[ys[i], xs[i]]
        if all(np.hypot(cx - px, cy - py) >= min_separation_km * 1000.0 for px, py in chosen_xy):
            chosen_xy.append((cx, cy))
            chosen_rc.append((ys[i], xs[i]))
        if len(chosen_xy) >= n_pods:
            break

    radii_m = rng.uniform(radius_range_m[0], radius_range_m[1], size=len(chosen_xy))

    pod_mask = np.zeros_like(suitable_mask)
    ny, nx = suitable_mask.shape
    for (cx, cy), r_m in zip(chosen_xy, radii_m):
        r_cells = int(np.ceil(r_m / cellsize_m)) + 1
        # locate the pod center's own row/col via the matching chosen_rc entry
        ry, rx = chosen_rc[chosen_xy.index((cx, cy))]
        r0, r1 = max(0, ry - r_cells), min(ny, ry + r_cells + 1)
        c0, c1 = max(0, rx - r_cells), min(nx, rx + r_cells + 1)
        sub_xx, sub_yy = xx[r0:r1, c0:c1], yy[r0:r1, c0:c1]
        disk = (np.hypot(sub_xx - cx, sub_yy - cy) <= r_m) & suitable_mask[r0:r1, c0:c1]
        pod_mask[r0:r1, c0:c1] |= disk

    return pod_mask, suitable_mask, chosen_xy, radii_m
