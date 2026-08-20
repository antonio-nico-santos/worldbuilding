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

import json

import numpy as np
from matplotlib.path import Path

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
CLASS_MARBLE = 5
CLASS_SEDIMENTARY_LIMESTONE = 6
CLASS_GRANITE = 7

LITHOLOGY_CLASSES.update({
    CLASS_MARBLE: "marble",
    CLASS_SEDIMENTARY_LIMESTONE: "sedimentary_limestone",
    CLASS_GRANITE: "granite",
})

# --- Authored lithology zones (v6) --------------------------------------------------
# marble / sedimentary_limestone / granite, added for the Urban Scale materials-economy
# follow-up (decision doc S8). Unlike schist/greywacke/basin_fill/volcanic, none of
# these three are derivable from elevation+relief (terrain_relief.classify_from_terrain)
# -- their presence is a fact of geological age/depositional history (Oligocene
# limestone cover, Ordovician Takaka-terrane marble, Miocene Coromandel-style granite
# intrusion), orthogonal to current terrain shape. So they're hand-authored polygons
# (data/input/lithology_authoral.geojson, drawn by Nico), composited ON TOP of the v5
# DEM-native result by apply_authoral_zones() below -- NOT part of classify_from_terrain
# itself, and not consumed by the legacy classify()/build_zone_fields path above either.

AUTHORAL_FEATURE_TYPE_TO_CLASS = {
    "marble": CLASS_MARBLE,
    "sedimentary_limestone": CLASS_SEDIMENTARY_LIMESTONE,
    "granite": CLASS_GRANITE,
}

# priority_rank tiers -- Nico's own scheme (chat, not a citation): each rank is a slot
# INSERTED into the existing volcanic > schist > greywacke > basin_fill chain. E.g.
# "rank 1 sits under volcanic and above schist" means volcanic still beats rank 1, but
# rank 1 beats schist (and everything weaker than schist, transitively). Full chain:
#   volcanic > 1 > schist > 2 > greywacke > 3 > basin_fill
# So each rank's "allowed targets" (what it may overwrite) is every base class strictly
# weaker than its own slot. Volcanic is never overwritable by any rank -- there's no
# slot above it. There's also no slot below basin_fill: a hypothetical rank 4 would be
# a no-op (basin_fill would always win), which is exactly the bug caught and fixed in
# this session's chat before this code was written (Granite South briefly at rank 4).
PRIORITY_RANK_BEATS = {
    1: {CLASS_SCHIST, CLASS_GREYWACKE, CLASS_BASIN_FILL},
    2: {CLASS_GREYWACKE, CLASS_BASIN_FILL},
    3: {CLASS_BASIN_FILL},
}


def load_authoral_zones(path: str):
    """Parse `lithology_authoral.geojson`. Deliberately NOT using
    `terrain.skeleton.build_zone_fields`: that loader's ZoneField carries Tappa-1-
    specific fields (target_elevation_m, amplitude_scale, edge_transition_km) that
    mean nothing for a rock-type layer, and its softened-edge `blend_weight` semantics
    are wrong here too -- a rock-type boundary is a hard structural fact, not something
    that should fade the way an elevation target does (same reasoning the legacy
    `classify()` above already used for schist/greywacke's own zone check). This loader
    reuses only what actually applies: matplotlib.path.Path membership on the exterior
    ring (`coordinates[0]`) -- none of the six current authored polygons have holes, so
    a hole would silently be ignored rather than erroring, not yet exercised.

    Returns a list of dicts: name, feature_type, class_code, priority_rank,
    grounded_claimed, path (matplotlib Path), bounds (minx, miny, maxx, maxy).
    """
    with open(path) as f:
        data = json.load(f)

    zones = []
    for feat in data["features"]:
        props = feat["properties"]
        ftype = props["feature_type"]
        if ftype not in AUTHORAL_FEATURE_TYPE_TO_CLASS:
            raise ValueError(
                f"unrecognized feature_type {ftype!r} in {path} (feature {props.get('name')!r}) "
                f"-- expected one of {sorted(AUTHORAL_FEATURE_TYPE_TO_CLASS)}"
            )
        rank = props["priority_rank"]
        if rank not in PRIORITY_RANK_BEATS:
            raise ValueError(
                f"priority_rank={rank} on {props.get('name')!r} has no defined tier "
                f"-- valid ranks are {sorted(PRIORITY_RANK_BEATS)} (see module docstring)"
            )
        geom = feat["geometry"]
        # Exterior ring only, per part -- holes ignored (see docstring). Handles both
        # "Polygon" (coordinates = [ring, ...holes]) and "MultiPolygon" (coordinates =
        # [[ring, ...holes], ...another part]) -- QGIS re-saved this file as
        # single-part MultiPolygon at some point in this session's editing, even
        # though every feature still has exactly one part; written generally rather
        # than assuming that stays true.
        if geom["type"] == "Polygon":
            part_rings = [geom["coordinates"][0]]
        elif geom["type"] == "MultiPolygon":
            part_rings = [part[0] for part in geom["coordinates"]]
        else:
            raise ValueError(f"unsupported geometry type {geom['type']!r} on {props.get('name')!r}")
        rings = [np.array(r, dtype=np.float64) for r in part_rings]
        all_pts = np.vstack(rings)
        minx, miny = all_pts.min(axis=0)
        maxx, maxy = all_pts.max(axis=0)
        zones.append({
            "name": props.get("name", "unnamed authoral zone"),
            "feature_type": ftype,
            "class_code": AUTHORAL_FEATURE_TYPE_TO_CLASS[ftype],
            "priority_rank": rank,
            "grounded_claimed": bool(props.get("grounded", False)),
            "paths": [Path(r) for r in rings],
            "bounds": (float(minx), float(miny), float(maxx), float(maxy)),
        })
    return zones


def apply_authoral_zones(lithology, zones, xmin, xmax, ymin, ymax, cellsize_m):
    """Composite hand-authored marble/sedimentary_limestone/granite zones on top of an
    existing DEM-native lithology array (v5's output), respecting each zone's
    `priority_rank` tier (PRIORITY_RANK_BEATS above). Zones are painted weakest-rank
    first, rank 1 last, so that if two authored zones of DIFFERENT ranks ever overlap
    each other, the higher-priority one wins deterministically -- not exercised by the
    current file (verified zero pairwise overlap among all six zones, chat QA before
    this code existed), but correct in general.

    A zone only overwrites a cell whose CURRENT class is in its own tier's allowed set
    -- it can never touch volcanic (no tier permits that), and can never touch ocean
    (CLASS_OCEAN isn't in any tier's allowed set either, so a polygon drawn partway
    into the ocean -- North Coast Limestone and Sedimentary Bay both do this on
    purpose, confirmed with Nico -- simply doesn't paint there).

    Returns (new_lithology, per_zone_stats). `per_zone_stats` carries the real numbers
    needed to judge each zone's `grounded` claim: nominal polygon area (land cells the
    polygon actually covers), the composition of what was UNDER the polygon before
    painting (the actual grounded-check evidence, not the geojson's own asserted flag),
    and how much of the nominal area actually got painted (can be smaller than nominal
    -- e.g. a rank-3 zone drawn partly over greywacke won't claim that portion).
    """
    ny, nx = lithology.shape
    y_top = ymax  # row 0 = north, matches every other raster in this pipeline
    out = lithology.copy()
    per_zone_stats = []
    cell_km2 = (cellsize_m / 1000.0) ** 2

    for zone in sorted(zones, key=lambda z: -z["priority_rank"]):  # weakest rank first, rank 1 painted last
        minx, miny, maxx, maxy = zone["bounds"]
        col0 = max(0, int((minx - xmin) / cellsize_m) - 1)
        col1 = min(nx, int((maxx - xmin) / cellsize_m) + 2)
        row0 = max(0, int((y_top - maxy) / cellsize_m) - 1)
        row1 = min(ny, int((y_top - miny) / cellsize_m) + 2)

        cols = np.arange(col0, col1)
        rows = np.arange(row0, row1)
        xs = xmin + (cols + 0.5) * cellsize_m
        ys = y_top - (rows + 0.5) * cellsize_m
        xx, yy = np.meshgrid(xs, ys)
        pts = np.column_stack([xx.ravel(), yy.ravel()])
        inside_flat = np.zeros(pts.shape[0], dtype=bool)
        for p in zone["paths"]:
            inside_flat |= p.contains_points(pts)
        inside = inside_flat.reshape(row1 - row0, col1 - col0)

        before = lithology[row0:row1, col0:col1]  # pre-authoral base, for composition reporting
        allowed = PRIORITY_RANK_BEATS[zone["priority_rank"]]
        paintable = inside & np.isin(before, list(allowed))

        vals, counts = (np.unique(before[inside], return_counts=True) if inside.any() else (np.array([]), np.array([])))
        total_inside = int(counts.sum())
        composition_pct = {
            LITHOLOGY_CLASSES.get(int(v), f"class_{v}"): float(100.0 * c / total_inside)
            for v, c in zip(vals, counts)
        } if total_inside else {}

        window = out[row0:row1, col0:col1]
        window[paintable] = zone["class_code"]
        out[row0:row1, col0:col1] = window

        per_zone_stats.append({
            "name": zone["name"],
            "feature_type": zone["feature_type"],
            "priority_rank": zone["priority_rank"],
            "grounded_claimed": zone["grounded_claimed"],
            "nominal_polygon_area_km2": total_inside * cell_km2,
            "painted_area_km2": float(paintable.sum()) * cell_km2,
            "pre_authoral_composition_pct": composition_pct,
        })

    return out, per_zone_stats


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


def place_material_pods(
    eligible_mask: np.ndarray,
    weight_field: np.ndarray | None,
    xx: np.ndarray,
    yy: np.ndarray,
    cellsize_m: float,
    n_pods: int = 10,
    min_separation_km: float = 5.0,
    radius_range_m: tuple[float, float] = (300.0, 800.0),
    candidate_pool_size: int = 500,
    seed: int = 13,
):
    """Generalizes `place_jade_pods` to any material/eligibility mask, for the four
    Vertice materials that were only ever a non-spatial per-class lookup
    (`resources.py`'s `VERTICE_MATERIALS`) until this function existed: laumontite,
    vivianite, reworked placer magnetite, native silver/copper. Jade's own
    `place_jade_pods` is left untouched (already locked/committed, no reason to risk
    changing its output) -- this is the same stochastic-pod technique, factored out
    so it isn't reimplemented four more times with copy-paste drift.

    `eligible_mask` replaces jade's `jade_eligible_mask(...)` call -- caller decides
    what "eligible ground" means for this material (e.g. the whole greywacke class,
    or basin_fill restricted to a wetness proxy -- see `run_tappa8_resource_pods.py`
    for what each of the four actually uses and why).

    `weight_field` replaces jade's `schist_grade` -- pass `None` for UNIFORM
    placement (every eligible cell equally likely) when there's no citable
    within-class gradient to weight by. Passing a fabricated field just to have
    *something* to weight by would be worse than admitting none exists -- laumontite
    is placed uniformly for exactly this reason (see the driver script).

    Same mechanics as `place_jade_pods` otherwise: percentile-style suitability zone
    already baked into `eligible_mask` by the caller, then a stochastic, weighted,
    minimum-separation pod draw -- NOT deterministic, representing "we know roughly
    where to look" rather than "here is the ore." `n_pods`/`min_separation_km`/
    `radius_range_m` are first-pass placeholders here too, same status as jade's.

    Returns (pod_mask, chosen_xy, radii_m) -- no `suitable_mask` in the return since
    the caller already has `eligible_mask` (unlike jade, which computed it inside).
    """
    ys, xs = np.where(eligible_mask)
    if len(ys) == 0:
        return np.zeros_like(eligible_mask), [], []

    rng = np.random.default_rng(seed)
    if weight_field is None:
        weights = np.ones(len(ys), dtype=np.float64)
    else:
        weights = weight_field[ys, xs].astype(np.float64)
        weights = np.clip(weights, 0.0, None)  # negative weights would break np.random.choice's p=
        if weights.sum() == 0:
            weights = np.ones(len(ys), dtype=np.float64)
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

    pod_mask = np.zeros_like(eligible_mask)
    ny, nx = eligible_mask.shape
    for idx, ((cx, cy), r_m) in enumerate(zip(chosen_xy, radii_m)):
        r_cells = int(np.ceil(r_m / cellsize_m)) + 1
        ry, rx = chosen_rc[idx]
        r0, r1 = max(0, ry - r_cells), min(ny, ry + r_cells + 1)
        c0, c1 = max(0, rx - r_cells), min(nx, rx + r_cells + 1)
        sub_xx, sub_yy = xx[r0:r1, c0:c1], yy[r0:r1, c0:c1]
        disk = (np.hypot(sub_xx - cx, sub_yy - cy) <= r_m) & eligible_mask[r0:r1, c0:c1]
        pod_mask[r0:r1, c0:c1] |= disk

    return pod_mask, chosen_xy, radii_m


# --- Overland travel friction (Tappa 6 cost-distance follow-up) --------------------
# Added for the "transport lithology-cost multiplier" item from Nico's post-Tappa-8
# follow-up list (decision doc S8f). SCOPE, confirmed with Nico before building this:
# a rock-type FRICTION WEIGHT layered onto Tappa 6's existing Tobler-hiking-function +
# Dijkstra cost graph (src/suitability/cost_distance.py), not a literal road/rail
# construction-cost model. Nico's own phrasing was about "helping/difficulting the
# construction of the transport system" -- read literally that would mean excavation
# cost, foundation stability, tunnel/bridge risk, etc. NONE of that is modelled here:
# there is no citable unit-cost data for any of it anywhere in this project, and
# building it would be a materially bigger, separate undertaking (effectively its own
# Tappa). What IS built is the closest defensible proxy available from the existing
# cost graph: a multiplier on Tobler's WALKING/TRAVEL speed, i.e. "how much does this
# rock type change how fast people and goods cross it on foot/by pack animal", which
# is what the current cost-distance model can actually represent. If Nico wants an
# actual construction-cost axis (excavatability, foundation grade, tunnelling risk)
# later, that is a distinct, uncited design task, not an extension of this one.
#
# Every multiplier below is an UNREVIEWED FIRST-PASS ESTIMATE (same status as jade's
# pod count, the silver/copper 25/75 split, etc.) -- a directional judgement call
# grounded in real-world travel-difficulty accounts for each rock type, NOT a number
# taken from a specific cited study (no such per-lithology hiking-speed dataset was
# found or is claimed to exist). Pending Nico's sign-off; NOT written to
# config/parameters.yml.
#
#   sedimentary_basin_fill  1.00  baseline/no penalty. Tobler's function is itself
#                                 calibrated on ordinary, generally-unobstructed
#                                 ground -- soft alluvial plain is close to what it
#                                 already assumes, so basin_fill is the natural anchor
#                                 (every other value is relative to this one).
#   greywacke_argillite     0.85  greywacke/schist               } NZ high-country
#   schist                  0.85  bedrock terrain               } off-track travel
#     Both get the SAME value on purpose -- no source differentiates walking speed
#     between the two specifically, and inventing a distinction would be worse than
#     admitting there isn't one. Real-world grounding: NZ tramping/mountain-safety
#     guidance consistently treats unmarked alpine/subalpine schist/greywacke terrain
#     (scree, tussock, uneven jointed rock) as materially slower than a formed track --
#     commonly a 30-50%-ish time penalty in that kind of guidance, never a single
#     canonical number. 0.85 sits at the mild end of that range deliberately, since
#     these two classes cover the majority of this world's land area (3030 km2 +
#     1773 km2 of the ~9904 km2 total, per lithology_v6) -- a more aggressive value
#     would have an outsized, poorly-justified effect on the whole network.
#   volcanic                 0.70  Real basalt/lava terrain can be dramatically slower
#                                 than ordinary ground (Hawaii/Iceland lava-field
#                                 hiking accounts describe some aa-type basalt as
#                                 barely passable, well under 1 km/h). This project
#                                 does NOT model aa vs. pahoehoe (no such distinction
#                                 exists anywhere in the lithology pipeline) -- 0.70 is
#                                 a deliberately MODERATE "jointed, blocky basalt"
#                                 value, not a worst-case one. Likely too GENTLE if any
#                                 of this world's volcanic zone is meant to be aa-style
#                                 rubble, too HARSH if it's meant to be older,
#                                 weathered/vegetated basalt -- flagged, not resolved.
#   marble                   0.60  karst dissolution terrain (sinkholes, karren,
#   sedimentary_limestone     0.60  fissured limestone pavement) is widely described in
#                                 real-world accounts (the Burren in Ireland, Mendip in
#                                 England, and this project's own karst-cave citations
#                                 already used for cave_karst) as notably hazardous/
#                                 slow off-trail going -- the same dissolution mechanism
#                                 karst_cave_candidates() already treats identically for
#                                 both classes, so they get the same travel penalty too.
#                                 This is this pipeline's most severe penalty, but both
#                                 classes are tiny (77.8 km2 + 87.7 km2) so its effect on
#                                 the overall network is necessarily local.
#   granite                  0.90  mild penalty. Real granite terrain (Sierra-Nevada-
#                                 style slab/tor country) is often noted as offering
#                                 GOOD footing on unweathered slab, with the difficulty
#                                 concentrated in jointed boulder fields rather than
#                                 the rock itself -- a light penalty, not a severe one.
#                                 Also this pipeline's smallest zone by far (13.8 km2,
#                                 a single authored polygon), so low-leverage regardless
#                                 of the exact value chosen.
LAND_TRAVEL_FRICTION = {
    CLASS_BASIN_FILL: 1.00,
    CLASS_GREYWACKE: 0.85,
    CLASS_SCHIST: 0.85,
    CLASS_VOLCANIC: 0.70,
    CLASS_MARBLE: 0.60,
    CLASS_SEDIMENTARY_LIMESTONE: 0.60,
    CLASS_GRANITE: 0.90,
}


def travel_friction_multiplier(
    lithology: np.ndarray,
    multipliers: dict = LAND_TRAVEL_FRICTION,
    default: float = 1.0,
):
    """Map each lithology class-code cell to a travel-friction multiplier
    (<=1.0 = slower than Tobler's baseline slope-only estimate would say),
    for src/suitability/cost_distance.py's `build_cost_graph(...,
    friction_multiplier=...)` hook.

    Cells whose class code isn't a key in `multipliers` -- CLASS_OCEAN (0),
    or any land_mask cell that lands on an ocean-classified lithology_v6
    pixel due to the two rasters' independent land determinations not
    perfectly agreeing at the coastline -- fall back to `default` (neutral,
    no penalty), NOT an error. Counting/logging how many land cells hit this
    fallback is the caller's job (see run_tappa8_transport_friction.py),
    since only the caller knows what "land" means for its own grid (this
    function has no land_mask of its own to check against).

    Returns a float64 array, same shape as `lithology`.
    """
    out = np.full(lithology.shape, float(default), dtype=np.float64)
    for code, mult in multipliers.items():
        out[lithology == code] = mult
    return out


# --- Excavation effort (construction-difficulty follow-up to LAND_TRAVEL_FRICTION) -
# Nico's follow-up to the transport friction multiplier (S8f): a SEPARATE relative
# index for "how much harder is it to BUILD ON/dig INTO this ground", after clarifying
# his original "construction" framing wasn't asking for a full engineering cost model
# (excavation crew-hours, blasting budgets, haulage) -- just a relative, per-class
# number, same spirit as LAND_TRAVEL_FRICTION above but a genuinely DIFFERENT physical
# property. Deliberately named "excavation effort," not "construction cost" -- this
# still only measures how hard the raw rock/sediment itself is to break/cut, NOT
# foundation bearing capacity, haulage distance, or labor/skill availability, all of
# which a real construction-cost model would also need and none of which this project
# has any data to build. See run_tappa8_excavation_effort.py / decision doc S8g for the
# full reasoning and an important caveat: this does NOT capture karst's foundation-
# STABILITY risk (sinkholes/voids under marble and sedimentary_limestone, see
# karst_cave_candidates()) -- that's a separate hazard from raw workability, and this
# index, taken alone, would make marble/limestone look like reasonably easy building
# ground, which is only true for the QUARRYING half of the question.
#
# Values are ratios relative to basin_fill = 1.0 (higher = more effort), grounded in
# real, well-documented relative hardness rather than a specific numeric source: calcite
# (limestone/marble's main mineral, Mohs ~3) is genuinely much softer than the
# quartz/feldspar that dominates greywacke/schist/granite (Mohs ~6-7), and this isn't
# just a mineralogy table exercise -- it's borne out by real historical stoneworking
# practice (ancient Egyptian/Mediterranean quarrying cut limestone with simple
# copper/bronze tools, while granite needed much harder abrasives/dolerite pounders,
# a well-attested, non-NZ-specific but genuinely general fact). UNREVIEWED first-pass
# estimates, same status as LAND_TRAVEL_FRICTION; pending Nico's sign-off; not written
# to config/parameters.yml.
#
#   sedimentary_basin_fill  1.0  baseline/anchor -- unconsolidated alluvium, diggable
#                                by hand/basic tools, no blasting.
#   sedimentary_limestone    1.3  soft carbonate rock, historically one of the easiest
#                                building stones to cut (this project's own Oamaru-stone
#                                citation, S8, is literally famous for being soft/easy
#                                to carve).
#   marble                   1.6  SAME calcite mineral as limestone, deliberately given a
#                                DIFFERENT (higher) value here, unlike LAND_TRAVEL_FRICTION
#                                where the two share a value -- metamorphic recrystallization
#                                makes marble's interlocking crystal fabric denser and
#                                measurably harder to carve than ordinary sedimentary
#                                limestone despite comparable mineral hardness, a real
#                                stoneworking distinction (this is the kind of thing a
#                                single "karst = 0.60" travel number can't capture, and
#                                doesn't need to -- it's a different physical property).
#   greywacke_argillite       2.3  greywacke/schist         } quartz-rich, well-indurated
#   schist                    2.3  bedrock                  } bedrock -- NZ engineering-
#     Same value on purpose, same reasoning as LAND_TRAVEL_FRICTION: nothing
#     differentiates them specifically. Real asymmetry NOT modeled here: schist's
#     foliation makes it genuinely easier to split ALONG cleavage planes than cut ACROSS
#     them (historically exploited -- Otago schist flagstone is a real building-stone
#     tradition) -- a single scalar can't represent that directionality, flagged rather
#     than faked.
#   volcanic                  2.6  basalt -- dense, fine-grained, tough; the standard
#                                "needs blasting" rippability-chart case.
#   granite                   3.0  hardest value here -- coarse crystalline, no foliation
#                                weakness at all (unlike schist), the canonical
#                                hardest-to-excavate common rock type in both mining/civil
#                                engineering rippability charts and historical
#                                stoneworking accounts alike.
EXCAVATION_EFFORT_MULTIPLIER = {
    CLASS_BASIN_FILL: 1.0,
    CLASS_SEDIMENTARY_LIMESTONE: 1.3,
    CLASS_MARBLE: 1.6,
    CLASS_GREYWACKE: 2.3,
    CLASS_SCHIST: 2.3,
    CLASS_VOLCANIC: 2.6,
    CLASS_GRANITE: 3.0,
}


def excavation_effort_multiplier(
    lithology: np.ndarray,
    multipliers: dict = EXCAVATION_EFFORT_MULTIPLIER,
    default: float = 1.0,
):
    """Map each lithology class-code cell to an excavation-effort multiplier
    (>=1.0 = harder/more effort than basin_fill's baseline) -- a standalone
    per-class attribute, NOT wired into cost_distance.py's graph (unlike
    `travel_friction_multiplier`, this has no "traversal" meaning -- effort
    to dig INTO a cell isn't a property of an edge between two cells).

    Same fallback convention as `travel_friction_multiplier`: an unmapped
    code (CLASS_OCEAN, or a land/lithology coastline mismatch) gets
    `default` (neutral), not an error -- meaningless for ocean cells either
    way, since nothing is ever built there.

    Returns a float64 array, same shape as `lithology`.
    """
    out = np.full(lithology.shape, float(default), dtype=np.float64)
    for code, mult in multipliers.items():
        out[lithology == code] = mult
    return out
