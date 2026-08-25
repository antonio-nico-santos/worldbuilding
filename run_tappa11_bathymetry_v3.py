"""
Tappa 11 -- Bathymetry, v3.

v1/v2 hand-coded every zone's shape in Python (rectangles, distance fields,
a few manually-placed trench polylines). Nico's request this round: draw
the shapes THEMSELVES in QGIS -- the same authored-skeleton workflow Tappa 1
already uses for the DEM (src/terrain/skeleton.py's RidgeField/ZoneField,
src/terrain/generate.py's compositing formula) -- and have this script turn
those drawings into the raster, instead of me hand-placing geometry in code.

Inputs (hand-authored in QGIS, staged from the device repo):
  data/input/bathymetry_ridges.geojson  -- 4 lines: Guardian creek (a trench,
    feature_type authored as 'ridge' after Nico's fix), North camling ridge
    1/2, South Barrier. peak_elevation_m is ADDITIVE to whatever's under it
    (matches RidgeField's own semantics -- see skeleton.py's docstring).
  data/input/bathymetry_zones.geojson -- 2 polygons: Bridge Base (plateau,
    -27m) with Underground Lake (plateau, -75m) nested fully inside it (a
    real containment, checked: 169/169 vertices inside -- this is the
    nested-zone "undersea lake" technique discussed in chat).

Design decisions made here that Nico should sanity-check against intent:
  1. build_ridge_fields()/build_zone_fields() (the actual Tappa 1 reference
     loaders in skeleton.py) are NOT used as-is: (a) they require the exact
     literal feature_type=='ridge', but 'Guardian creek' was drawn as a
     trench/valley -- a custom loader here accepts 'ridge'/'creek'/'valley'/
     'trench' as equivalent labels, since a negative peak_elevation_m
     already encodes the semantic difference; (b) they don't read
     shelf_multiplier from GeoJSON properties at all (only via an external
     {name: value} dict) -- the custom loader here reads the authored
     'shelf_muliplier' property (typo, as drawn) directly per feature.
  2. Guardian creek REPLACES v2's hand-placed trench_a/trench_b (both lived
     in the same NE-archipelago / Povo Silencioso region) rather than
     adding to them -- the whole point of authoring shapes was for Nico to
     place this by hand instead of me guessing waypoints. Povo Silencioso's
     mainland-facing depth ASYMMETRY (sheltered vs. open-sea target) has no
     authored equivalent yet, so that hardcoded gradient is kept.
  3. The ferry corridor has no authored zone/ridge yet either, so v2's
     domain-warped 4-waypoint corridor is kept unchanged.
  4. Bridge Base/Underground Lake REPLACE v2's hardcoded rectangular bridge
     zone entirely -- same role, now an authored shape. Composited using
     the exact plateau formula from Tappa 1's own generate.py: structure =
     structure*(1-w) + target_elevation_m*w, blend weight from ZoneField's
     smootherstep-on-signed-distance. Applied largest-area-first so a
     nested zone (Underground Lake) composites AFTER its container (Bridge
     Base) and correctly wins in the overlap.
  5. Multiple ridge contributions are SUMMED, not combined via np.maximum
     the way Tappa 1's mountain-building ridges are. Tappa 1's max() is
     specific to unidirectional peak-building (multiple mountains competing
     for "tallest here wins"); this domain's ridges include negative
     peaks (trenches), where max() would make a trench a no-op against any
     less-negative background. v2 already established SUM as the bathymetry
     convention for this exact reason (see its trench_a+trench_b comment) --
     kept for consistency.
  6. Background profile: SHELF_BREAK_KM moved 18->12km (Nico's request),
     AND the coast->shelf-break segment now uses an ease-out curve
     (1-(1-t)^2) instead of smootherstep -- smootherstep's zero-derivative
     start was producing unrealistic -0.7..-7.2m depths 1-2.3km offshore;
     option 1 from that chat discussion, applied here.
  7. Defense-in-depth clamp (true_ocean cells never exceed -0.5m, same as
     v1/v2) is still active. This means a ridge that manages to crest above
     sea level gets flattened to -0.5m rather than becoming new land --
     given Nico's latest edit dialed ridge heights down specifically to
     stay submerged (peak_elevation_m 35/35/75 -> 20/23/20, cutting
     surfacing risk from 44/28/69% to 8/5/0% of sampled vertices), keeping
     the clamp as a safety net rather than removing it seemed like the
     right call -- but see the clamp-hit report at the end of this run for
     exactly how much (if any) of that residual 5-8% risk actually landed
     on the clamp, since a flattened -0.5m patch can look artificial in
     QGIS if it's more than a few cells.
  8. First run of this script surfaced a real bug, not just the ridges'
     residual risk: the roughness formula had a flat 8m floor right at the
     coastline (d=0), oversized relative to how shallow the background is
     in the first ~1km (v2 likely had the same problem, just silently
     absorbed by its own early background-only clamp with no reporting).
     Result: 368,507 true_ocean cells (331.6 km2, 661 separate patches
     strung along nearly every coastline in the domain) hit the -0.5m
     clamp on the first run -- a flat "bathtub ring" artifact, not a rare
     edge case. Fixed by ramping roughness's near-shore term from 0 (at
     the coast, matching zero water depth) up to 8m over the first 1.5km
     via the same ease-out curve as point 6, instead of starting flat at
     8m. That cut the clamp area 331.6km2 -> 35.8km2 but didn't zero it --
     655 blobs remained, still scattered along the whole coastline and
     nowhere near the authored ridges/zones. Root cause #2: the coarse
     (basin-scale, 9km wavelength) noise band still had full weight right
     at the coast, where it has no physical business being -- a regional
     "high" in that slow band could align with a local fine-noise spike
     and breach sea level. Gated coarse's contribution 0->full over
     0.5-2.5km offshore. See the clamp-hit report below for the result.
  9. Sampling the actual generated raster (not the earlier pointwise
     validator, which never modeled Povo at all) along North camling ridge
     1's line found its crest at ~-287m instead of the intended ~-20..-30m
     reef depth -- a real bug, not a design surprise. Povo's box-membership
     falloff (dist_to_ne_land) measures distance to the nearest land cell
     INSIDE the NE-archipelago box with no regard for whether the query
     point itself is near the box; a point ~20km outside the box can still
     be geographically close to a land pixel just inside its edge. Fixed
     by gating on distance to the box itself (see the ridges' compositing
     section) so this only ever applies within ~3km of the box, same as
     originally intended.
"""
import sys, time, json, os
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from matplotlib.path import Path

sys.path.insert(0, "src")
from terrain.noise import Simplex2D, ridged_fbm, domain_warp  # noqa: E402

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)

XMIN, XMAX, YMIN, YMAX, CELL = -65000.0, 65000.0, -80000.0, 80000.0, 30.0
CRS_PROJ4 = "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
IN_DIR = "data/processed"
OUT_DIR = "data/processed/bathymetry"
RIDGES_PATH = "data/input/bathymetry_ridges.geojson"
ZONES_PATH = "data/input/bathymetry_zones.geojson"
os.makedirs(OUT_DIR, exist_ok=True)

log("loading DEM...")
dem = np.load(f"{IN_DIR}/dem_v3_final_30m_eroded.npy").astype(np.float32)
ny, nx = dem.shape
y_top = YMIN + ny * CELL
assert (ny, nx) == (5334, 4334)

xs = XMIN + CELL * np.arange(nx)
ys = y_top - CELL * np.arange(ny)
X, Y = np.meshgrid(xs, ys)

land = dem > 0
ocean_raw = ~land
lbl_ocean, _ = ndimage.label(ocean_raw, structure=np.ones((3, 3)))
border_labels = set(np.unique(lbl_ocean[0, :])) | set(np.unique(lbl_ocean[-1, :])) | \
                set(np.unique(lbl_ocean[:, 0])) | set(np.unique(lbl_ocean[:, -1]))
border_labels.discard(0)
true_ocean = np.isin(lbl_ocean, list(border_labels))
lake = ocean_raw & ~true_ocean
log(f"true ocean {true_ocean.mean():.4f}  lake {lake.mean():.4f}")
del ocean_raw, lbl_ocean, border_labels

lbl_land, n_land = ndimage.label(land, structure=np.ones((3, 3)))
land_sizes = ndimage.sum(np.ones_like(lbl_land), lbl_land, index=np.arange(1, n_land + 1))
main_land_label = int(np.argmax(land_sizes) + 1)
is_mainland = (lbl_land == main_land_label)
del lbl_land, land_sizes, n_land

log("distance-to-coast (EDT, native 30m)...")
dist_to_land = (ndimage.distance_transform_edt(~land) * CELL / 1000.0).astype(np.float32)
log("distance-to-mainland (EDT vs. the single largest landmass only)...")
dist_to_mainland = (ndimage.distance_transform_edt(~is_mainland) * CELL / 1000.0).astype(np.float32)
del is_mainland


def smootherstep(x, edge0, edge1):
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * t * (t * (t * 6 - 15) + 10)


def ease_out_pow2(x, edge0, edge1):
    """1-(1-t)^2 -- nonzero slope at t=0, unlike smootherstep. Only used
    for the coast->shelf-break segment (see module docstring, point 6)."""
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return 1.0 - (1.0 - t) ** 2


ROW_CHUNK = 250


def chunked_apply(func, arrays, row_chunk=ROW_CHUNK, n_outputs=1):
    ny_, nx_ = arrays[0].shape
    outs = [np.empty((ny_, nx_), dtype=np.float32) for _ in range(n_outputs)]
    for r0 in range(0, ny_, row_chunk):
        r1 = min(r0 + row_chunk, ny_)
        chunk_args = [a[r0:r1] for a in arrays]
        res = func(*chunk_args)
        if n_outputs == 1:
            outs[0][r0:r1] = res
        else:
            for i in range(n_outputs):
                outs[i][r0:r1] = res[i]
    return outs[0] if n_outputs == 1 else outs


def densify(coords, max_spacing_m):
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


# =====================================================================
# Custom loaders for the authored shapes -- see module docstring point 1
# for why skeleton.py's reference loaders aren't used as-is.
# =====================================================================
RIDGE_TYPE_ALIASES = {"ridge", "creek", "valley", "trench"}


def load_ridges(path, densify_spacing_m=250.0, default_shelf_multiplier=3.0):
    with open(path) as f:
        features = json.load(f)["features"]
    out = []
    for feat in features:
        p = feat["properties"]
        ftype = p.get("feature_type")
        if ftype not in RIDGE_TYPE_ALIASES:
            raise ValueError(f"unexpected feature_type in ridges layer: {p}")
        coords = np.array(feat["geometry"]["coordinates"], dtype=np.float64)
        dense = densify(coords, densify_spacing_m)
        name = p.get("name", "unnamed ridge")
        shelf_mult = p.get("shelf_multiplier", p.get("shelf_muliplier", default_shelf_multiplier))
        out.append({
            "name": name,
            "tree": cKDTree(dense),
            "peak_elevation_m": float(p["peak_elevation_m"]),
            "falloff_km": float(p["falloff_km"]),
            "shelf_multiplier": float(shelf_mult),
        })
    return out


def shoelace_area(ring):
    c = np.array(ring)
    x, y = c[:, 0], c[:, 1]
    return 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]) + (x[-1] * y[0] - x[0] * y[-1]))


def load_zones(path, densify_spacing_m=250.0):
    with open(path) as f:
        features = json.load(f)["features"]
    out = []
    for feat in features:
        p = feat["properties"]
        ftype = p.get("feature_type")
        if ftype not in ("plateau", "amplitude_zone"):
            raise ValueError(f"unexpected feature_type in zones layer: {p}")
        ring = np.array(feat["geometry"]["coordinates"][0], dtype=np.float64)
        dense = densify(ring, densify_spacing_m)
        name = p.get("name", "unnamed zone")
        out.append({
            "name": name,
            "feature_type": ftype,
            "target_elevation_m": p.get("target_elevation_m"),
            "amplitude_scale": float(p["amplitude_scale"]),
            "edge_transition_km": float(p["edge_transition_km"]),
            "path": Path(ring),
            "tree": cKDTree(dense),
            "area_km2": shoelace_area(ring) / 1e6,
        })
    out.sort(key=lambda z: -z["area_km2"])  # largest (outermost) first
    return out


ridges = load_ridges(RIDGES_PATH)
zones = load_zones(ZONES_PATH)
log(f"loaded {len(ridges)} authored ridges, {len(zones)} authored zones "
    f"(nesting order: {' -> '.join(z['name'] for z in zones)})")

# =====================================================================
# Two-band noise, unchanged from v2.
# =====================================================================
log("coarse ridged-fBm band (seamount/basin scale)...")
coarse_noise_gen = Simplex2D(seed=7301)
coarse = chunked_apply(
    lambda x, y: ridged_fbm(coarse_noise_gen, x, y, octaves=4, base_freq=1.0 / 9000.0,
                             lacunarity=2.0, persistence=0.55, plain_octaves=2),
    [X, Y],
)
coarse = (coarse - coarse.mean()) / coarse.std()
log(f"  coarse band stats: mean={coarse.mean():.3f} std={coarse.std():.3f}")


def fbm_value_noise(shape, seed, octaves=5, persistence=0.55):
    rng = np.random.default_rng(seed)
    total = np.zeros(shape, dtype=np.float32)
    amp, amp_sum = 1.0, 0.0
    for o in range(octaves):
        divisor = 2 ** (7 - o)
        low_shape = (max(4, shape[0] // divisor), max(4, shape[1] // divisor))
        low = rng.standard_normal(size=low_shape).astype(np.float32)
        up = ndimage.zoom(low, (shape[0] / low_shape[0], shape[1] / low_shape[1]), order=3)[: shape[0], : shape[1]]
        total += up * amp
        amp_sum += amp
        amp *= persistence
    total /= amp_sum
    return ((total - total.mean()) / total.std()).astype(np.float32)


log("fine value-noise band (texture)...")
fine = fbm_value_noise((ny, nx), seed=1104, octaves=5, persistence=0.55)

# =====================================================================
# Background: shelf -> slope -> basin. 12km break (was 18km), option-1
# ease-out curve for the coast->break segment (was smootherstep) -- see
# module docstring point 6.
# =====================================================================
log("background shelf/slope/basin (option-1 curve, near-shore roughness ramp)...")
SHELF_BREAK_KM, SHELF_DEPTH_M = 12.0, 140.0
BASIN_START_KM, BASIN_DEPTH_M = 55.0, 880.0
w_shelf = ease_out_pow2(dist_to_land, 0.0, SHELF_BREAK_KM)
shelf_profile = -SHELF_DEPTH_M * w_shelf
w_slope = smootherstep(dist_to_land, SHELF_BREAK_KM, BASIN_START_KM)
base_profile = shelf_profile * (1 - w_slope) + (-BASIN_DEPTH_M) * w_slope

rough_shelf = 8.0 * ease_out_pow2(dist_to_land, 0.0, 1.5) + 20.0 * smootherstep(dist_to_land, 0.0, SHELF_BREAK_KM)
w_toward_peak = smootherstep(dist_to_land, SHELF_BREAK_KM, 35.0)
w_past_peak = smootherstep(dist_to_land, 35.0, BASIN_START_KM)
roughness = rough_shelf * (1 - w_toward_peak) + 70.0 * w_toward_peak
roughness = roughness * (1 - w_past_peak) + 35.0 * w_past_peak
log(f"  background profile: min={base_profile[true_ocean].min():.1f} max={base_profile[true_ocean].max():.1f}")

structure = base_profile.copy()
del shelf_profile, w_shelf, w_slope, rough_shelf, w_toward_peak, w_past_peak, base_profile

# =====================================================================
# Ferry corridor -- unchanged from v2 (no authored equivalent yet, see
# module docstring point 3). Baked-in noise term dropped here; noise is
# now applied once, uniformly, at the very end via amplitude_multiplier.
# =====================================================================
log("ferry corridor (unchanged from v2, domain-warped)...")
FERRY_WAYPOINTS_KM = np.array([
    [0.82, -21.94],
    [-9.5, -30.0],
    [-19.0, -42.0],
    [-27.65, -52.46],
]) * 1000.0

ferry_dense = densify(FERRY_WAYPOINTS_KM, 250.0)
ferry_arclen = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(ferry_dense, axis=0).T))])
ferry_total_len = ferry_arclen[-1]
ferry_tree = cKDTree(ferry_dense)

Xw_f, Yw_f = chunked_apply(
    lambda x, y: domain_warp(x, y, seed=4402, wavelength_m=6000.0, amplitude_m=1800.0, octaves=3),
    [X, Y], n_outputs=2,
)
warped_xy = np.column_stack([Xw_f.ravel(), Yw_f.ravel()])
dist_to_ferry_line_m, nearest_idx = ferry_tree.query(warped_xy, k=1)
dist_to_ferry_line_km = dist_to_ferry_line_m.reshape(ny, nx) / 1000.0
t_along = (ferry_arclen[nearest_idx] / ferry_total_len).reshape(ny, nx)

halfwidth_noise_gen = Simplex2D(seed=4403)
halfwidth_raw = chunked_apply(
    lambda t: halfwidth_noise_gen.noise2(t * 6.0, np.zeros_like(t)),
    [t_along],
)
halfwidth_km = 3.75 + 0.9 * halfwidth_raw
FERRY_BLEND_KM = 3.0
w_ferry = 1.0 - smootherstep(dist_to_ferry_line_km, halfwidth_km, halfwidth_km + FERRY_BLEND_KM)

FERRY_TARGET_M = -70.0
structure = np.where(true_ocean, structure * (1 - w_ferry) + FERRY_TARGET_M * w_ferry, structure)
HALFWIDTH_KM_RANGE = (float(halfwidth_km.min()), float(halfwidth_km.max()))
log(f"  ferry core: {(w_ferry[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2, "
    f"halfwidth range {HALFWIDTH_KM_RANGE[0]:.2f}-{HALFWIDTH_KM_RANGE[1]:.2f}km")
del (Xw_f, Yw_f, warped_xy, dist_to_ferry_line_m, nearest_idx, dist_to_ferry_line_km, t_along,
     halfwidth_raw, halfwidth_km, w_ferry, ferry_dense, ferry_tree, ferry_arclen)

# =====================================================================
# Povo Silencioso -- mainland-facing asymmetry kept from v2 (no authored
# equivalent yet); trench_a/trench_b DROPPED (replaced by the authored
# "Guardian creek" ridge, applied in the general ridge pass below -- see
# module docstring point 2).
# =====================================================================
log("Povo Silencioso asymmetry (trenches now come from authored ridges)...")
NE_X0, NE_X1 = 27000.0, XMAX + 20.0
NE_Y0, NE_Y1 = 58000.0, YMAX + 20.0
ne_box = (X >= NE_X0) & (X <= NE_X1) & (Y >= NE_Y0) & (Y <= NE_Y1)
land_in_ne = land & ne_box
dist_to_ne_land = (ndimage.distance_transform_edt(~land_in_ne) * CELL / 1000.0).astype(np.float32)

# BUG FOUND while checking North camling ridge 1/2 against the actual generated
# raster (not just the pointwise validator, which never modeled Povo at all):
# dist_to_ne_land measures distance to the nearest LAND CELL INSIDE the box --
# it does NOT care whether the QUERY point is itself anywhere near the box. A
# point 20km west of NE_X0 can still be geographically close (here, ~7km) to
# some land pixel that happens to sit just inside the box's edge, and that was
# enough for w_povo to reach ~0.92 there -- overriding the plain shelf
# background almost entirely with Povo's -320m sheltered target, ~20km outside
# the archipelago this zone was meant to cover. Confirmed directly: this put
# North camling ridge 1's crest at ~-287m instead of the intended ~-20..-30m
# reef depth. Gate on distance to the BOX ITSELF (not to land within it) so
# the leak is killed beyond a short buffer past the box edge, while behavior
# fully inside the box (where this was always intended to apply) is unchanged.
dx_out = np.maximum(np.maximum(NE_X0 - X, X - NE_X1), 0.0)
dy_out = np.maximum(np.maximum(NE_Y0 - Y, Y - NE_Y1), 0.0)
dist_to_ne_box_km = np.sqrt(dx_out**2 + dy_out**2) / 1000.0
box_gate = 1.0 - smootherstep(dist_to_ne_box_km, 0.0, 3.0)
del dx_out, dy_out, dist_to_ne_box_km

POVO_FRINGE_KM = 0.3
POVO_NEAR_KM = 2.5
POVO_FAR_KM = 22.0
w_rise = smootherstep(dist_to_ne_land, POVO_FRINGE_KM, POVO_NEAR_KM)
w_fall = 1.0 - smootherstep(dist_to_ne_land, POVO_NEAR_KM, POVO_FAR_KM)
w_povo = w_rise * w_fall * box_gate
w_povo = np.where(ne_box | (dist_to_ne_land <= POVO_FAR_KM), w_povo, 0.0)
del box_gate

SHELTERED_TARGET_M = -320.0
OPEN_SEA_TARGET_M = -1400.0
open_sea_factor = smootherstep(dist_to_mainland, 2.0, 15.0)
povo_target = SHELTERED_TARGET_M + (OPEN_SEA_TARGET_M - SHELTERED_TARGET_M) * open_sea_factor
log(f"  open-sea factor stats in box: min={open_sea_factor[ne_box].min():.2f} "
    f"max={open_sea_factor[ne_box].max():.2f} mean={open_sea_factor[ne_box].mean():.2f}")

structure = np.where(true_ocean, structure * (1 - w_povo) + povo_target * w_povo, structure)
log(f"  povo core (near-shore ring, full effect): "
    f"{((w_povo > 0.9) & true_ocean).sum() * (CELL/1000)**2:.2f} km2")
del ne_box, land_in_ne, dist_to_ne_land, w_rise, w_fall, w_povo, povo_target, open_sea_factor

# =====================================================================
# Authored ridges -- summed contributions (see module docstring point 5),
# applied to `structure` (i.e. relative/additive to whatever's under
# them: general background, ferry corridor, or Povo asymmetry alike).
# =====================================================================
log("authored ridges (Guardian creek, North camling 1/2, South Barrier)...")


def ridge_contribution(X_, Y_, ridge):
    tree = ridge["tree"]
    peak_m = ridge["peak_elevation_m"]
    falloff_m = ridge["falloff_km"] * 1000.0
    shelf_m = falloff_m * ridge["shelf_multiplier"]
    k = np.log(2.0)

    def _chunk(x, y):
        xy = np.column_stack([x.ravel(), y.ravel()])
        dist_m, _ = tree.query(xy, k=1)
        dist_m = dist_m.reshape(x.shape)
        base = peak_m * np.exp(-k * (dist_m / falloff_m) ** 2)
        t = np.clip((shelf_m - dist_m) / (shelf_m - falloff_m), 0.0, 1.0)
        taper = t * t * t * (t * (t * 6 - 15) + 10)
        taper = np.where(dist_m <= falloff_m, 1.0, taper)
        return base * taper

    return chunked_apply(_chunk, [X_, Y_])


ridge_total = np.zeros((ny, nx), dtype=np.float32)
for r in ridges:
    contrib = ridge_contribution(X, Y, r)
    ridge_total += contrib
    log(f"  {r['name']}: peak={r['peak_elevation_m']:.0f}m falloff={r['falloff_km']}km "
        f"shelf_mult={r['shelf_multiplier']} -> contribution range "
        f"[{contrib.min():.1f}, {contrib.max():.1f}]m")
    del contrib

structure = np.where(true_ocean, structure + ridge_total, structure)
del ridge_total

# =====================================================================
# Authored zones -- exact plateau formula from Tappa 1's generate.py:
# structure = structure*(1-w) + target_elevation_m*w; amplitude_multiplier
# blended the same way from each zone's amplitude_scale. Largest-area-first
# order (see load_zones) so Underground Lake composites after, and wins
# inside, Bridge Base.
# =====================================================================
log("authored zones (Bridge Base, Underground Lake)...")
amplitude_multiplier = np.ones((ny, nx), dtype=np.float32)

for z in zones:
    def _blend_chunk(x, y, z=z):
        xy = np.column_stack([x.ravel(), y.ravel()])
        inside = z["path"].contains_points(xy).reshape(x.shape)
        dist_m, _ = z["tree"].query(xy, k=1)
        dist_m = dist_m.reshape(x.shape)
        signed = np.where(inside, dist_m, -dist_m)
        band_m = z["edge_transition_km"] * 1000.0
        t = np.clip((signed + band_m) / (2 * band_m), 0.0, 1.0)
        return t * t * t * (t * (t * 6 - 15) + 10)

    w = chunked_apply(lambda x, y: _blend_chunk(x, y), [X, Y])
    if z["feature_type"] == "plateau":
        structure = np.where(true_ocean, structure * (1 - w) + z["target_elevation_m"] * w, structure)
    amplitude_multiplier = amplitude_multiplier * (1 - w) + z["amplitude_scale"] * w
    log(f"  {z['name']}: area={z['area_km2']:.2f}km2 target={z['target_elevation_m']}m "
        f"amplitude_scale={z['amplitude_scale']} edge_transition={z['edge_transition_km']}km "
        f"core coverage={(w[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f}km2")
    del w

del X, Y

# =====================================================================
# Noise applied once, uniformly: fine+coarse scaled by distance-based
# roughness AND each zone's own amplitude_scale (default 1.0 outside any
# zone -- i.e. unchanged background/ferry/povo roughness).
# =====================================================================
log("applying noise (uniform recipe, roughness x amplitude_multiplier, coarse gated out of the surf zone)...")
# coarse is seamount/basin-SCALE structure -- it has no business appearing right
# at the coastline, and its long (9km) wavelength means a single regional "high"
# can align with local fine-noise spikes to breach sea level in patches strung
# along the whole coast (confirmed directly: this was the actual driver of the
# residual clamp-hit blobs after the near-shore roughness-ramp fix, none of
# which sat anywhere near an authored ridge/zone). Gated 0->full over 0.5-2.5km,
# same distance the shelf itself is still mostly flat.
coarse_gate = smootherstep(dist_to_land, 0.5, 2.5)
result = structure + fine * roughness * 0.75 * amplitude_multiplier + coarse * roughness * 0.5 * amplitude_multiplier * coarse_gate
del structure, fine, coarse, roughness, amplitude_multiplier, coarse_gate, dist_to_land, dist_to_mainland

# =====================================================================
# Defense-in-depth clamp + recombine (same as v1/v2) -- report how many
# true_ocean cells actually hit the -0.5m ceiling, since that's the
# residual surfacing risk from the authored ridges landing somewhere.
# =====================================================================
pre_clamp_max = result[true_ocean].max()
clamped_cells = ((result > -0.5) & true_ocean).sum()
log(f"pre-clamp true_ocean max: {pre_clamp_max:.2f}m ; "
    f"{clamped_cells:,} true_ocean cells ({clamped_cells * (CELL/1000)**2:.3f} km2) hit the -0.5m clamp")

result = np.where(true_ocean, np.minimum(result, -0.5), result)
final = np.where(land, dem, np.where(true_ocean, result, dem)).astype(np.float32)

assert not np.isnan(final).any() and not np.isinf(final).any()
assert np.array_equal(final[land], dem[land])
assert np.array_equal(final[lake], dem[lake])
assert final[true_ocean].max() < 0.0
log("sanity checks passed")
log(f"FINAL true-ocean depth range: {final[true_ocean].min():.1f} .. {final[true_ocean].max():.1f} m")

# =====================================================================
# Export (identical conventions to v1/v2)
# =====================================================================
log("exporting...")
np.save(f"{OUT_DIR}/bathymetry_v3_30m.npy", final)


def write_envi_raw(path_stem, array, xmin, ymin, cellsize, description, dtype="i2", nodata=-9999.0):
    _codes = {"i2": (2, "<i2"), "f4": (4, "<f4")}
    code, np_dtype = _codes[dtype]
    a2 = array[None, ...] if array.ndim == 2 else array
    nbands, h, w = a2.shape
    out = np.where(np.isnan(a2), nodata, a2)
    if dtype == "i2":
        out = np.clip(np.round(out), -32768, 32767)
    out.astype(np_dtype).tofile(path_stem + ".bin")
    y_top_local = ymin + h * cellsize
    hdr = (
        "ENVI\n"
        f"description = {{{description}}}\n"
        f"samples = {w}\nlines = {h}\nbands = {nbands}\n"
        "header offset = 0\nfile type = ENVI Standard\n"
        f"data type = {code}\ninterleave = bsq\nbyte order = 0\n"
        f"map info = {{Arbitrary, 1, 1, {xmin}, {y_top_local}, {cellsize}, {cellsize}, units=Meters}}\n"
        f"data ignore value = {nodata}\n"
    )
    with open(path_stem + ".hdr", "w") as f:
        f.write(hdr)


write_envi_raw(f"{OUT_DIR}/bathymetry_v3_30m", final, XMIN, YMIN, CELL,
                description="Tappa 11 bathymetry v3 (authored shapes: hand-drawn ridges+zones "
                            "from QGIS, 12km shelf break with ease-out nearshore curve) -- "
                            "Fictional World LCC domain",
                dtype="i2")
with open(f"{OUT_DIR}/bathymetry_v3_30m.prj", "w") as f:
    f.write(CRS_PROJ4.strip() + "\n")

meta = {
    "grid": {"shape": [ny, nx], "xmin": XMIN, "ymin": YMIN, "resolution_m": CELL, "crs_proj4": CRS_PROJ4},
    "changes_from_v2": [
        "Bridge corridor: now the authored 'Bridge Base' plateau zone (target -27m) with 'Underground "
        "Lake' plateau (-75m) nested inside it, replacing the hardcoded rectangle",
        "Povo Silencioso trenches: dropped the two hardcoded trench_a/trench_b lines, replaced by the "
        "authored 'Guardian creek' ridge (feature_type fixed from 'creek' to 'ridge' by Nico); the "
        "mainland-facing sheltered/open-sea depth asymmetry is unchanged (no authored equivalent yet)",
        "Two new authored reef ridges with no v2 predecessor: North camling ridge 1/2",
        "New authored ridge: South Barrier (near the ferry corridor's western approach)",
        "Background: shelf break 18km->12km, coast->break segment switched from smootherstep to an "
        "ease-out (1-(1-t)^2) curve -- smootherstep's flat start was giving -0.7..-7.2m depths "
        "1-2.3km offshore; new curve gives ~-22m at 1km, ~-43m at 2km",
        "Noise application unified: one recipe (fine*0.75 + coarse*0.5, scaled by distance-based "
        "roughness and each zone's amplitude_scale) applied once at the end, instead of each zone "
        "baking its own bespoke noise weights",
        "Near-shore roughness: removed a flat 8m floor at the coastline (d=0), replaced with a ramp "
        "0m->8m over the first 1.5km (ease-out curve) -- the flat floor was oversized relative to "
        "background depth in the first ~1km and produced a 331.6km2 'bathtub ring' clamp artifact "
        "along nearly every coastline on the first v3 run; cut clamp area to 35.8km2",
        "Coarse (basin-scale) noise band gated out of the immediate surf zone (0->full weight over "
        "0.5-2.5km offshore) -- it was still breaching sea level in 655 scattered coastal patches "
        "after the roughness-ramp fix alone; see module docstring point 8",
        "Fixed a real bug: Povo Silencioso's box-membership falloff was leaking ~92% weight onto "
        "North camling ridge 1 (~20km outside the NE-archipelago box) via proximity to a land pixel "
        "just inside the box edge, putting its crest at -287m instead of the intended -20..-30m reef "
        "depth. Gated on distance to the box itself, not to land within it; see module docstring point 9",
    ],
    "ridges": [
        {"name": r["name"], "peak_elevation_m": r["peak_elevation_m"], "falloff_km": r["falloff_km"],
         "shelf_multiplier": r["shelf_multiplier"]} for r in ridges
    ],
    "zones": [
        {"name": z["name"], "feature_type": z["feature_type"], "target_elevation_m": z["target_elevation_m"],
         "amplitude_scale": z["amplitude_scale"], "edge_transition_km": z["edge_transition_km"],
         "area_km2": z["area_km2"]} for z in zones
    ],
    "ferry": {"waypoints_km": FERRY_WAYPOINTS_KM.tolist(), "halfwidth_km_range": list(HALFWIDTH_KM_RANGE),
              "target_m": FERRY_TARGET_M},
    "povo_silencioso": {"sheltered_target_m": SHELTERED_TARGET_M, "open_sea_target_m": OPEN_SEA_TARGET_M},
    "shelf_break_km": SHELF_BREAK_KM, "shelf_depth_m": SHELF_DEPTH_M,
    "clamp_hit_true_ocean_cells": int(clamped_cells),
    "clamp_hit_true_ocean_km2": float(clamped_cells * (CELL / 1000.0) ** 2),
    "pre_clamp_true_ocean_max_m": float(pre_clamp_max),
    "final_true_ocean_depth_range_m": [float(final[true_ocean].min()), float(final[true_ocean].max())],
}
with open(f"{OUT_DIR}/bathymetry_v3_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
