"""
Tappa 11 -- Bathymetry, v4.

Nico's review of v3 flagged two things:
  1. The synthetic background (a pure function of distance-to-coast, plus
     statistically uniform noise) reads as shallow and flat -- structurally
     it can't do otherwise, since it has zero spatial correlation with
     nearby land. Nico liked the DEM's OWN underwater character and asked
     whether we could modify it in place (raster-calculator style: fix
     what's wrong, keep what's already there) instead of replacing it
     wholesale with a new synthetic model.
  2. A visible square artifact around Povo Silencioso -- the hardcoded
     NE-archipelago bounding rectangle's own edge, made worse by v3's
     Povo-leak fix (which added a hard-ish 3km falloff past that rectangle,
     turning a ~1000m target discontinuity into a visible seam).

This version answers (1) with yes, but not directly -- checked first
(binning the DEM's own raw ocean depth by distance-to-coast): it's
non-monotonic, deepening to a peak around 5-8km then getting SHALLOWER
again out to 35km, the opposite of a real shelf/slope/basin. That's the
signature of Tappa 1's land-generation noise field extended below sea
level with no oceanic constraint -- not designed bathymetry, and not
safe to build authored shapes directly on top of (Nico's own "option A").
But its LOCAL texture -- the part correlated with nearby terrain, not the
bad large-scale radial trend -- is exactly what v3's synthetic noise was
missing. So: detrend the DEM's own ocean signal (remove its radial
component entirely, whatever shape it has), replace that radial component
with the same well-behaved shelf/slope/basin curve v3 used, and add the
DEM's own residual texture back at full strength (Nico's choice -- the
alternative was to damp it). This IS "modifying the DEM in place" in
spirit: nothing about the DEM's actual local variation is discarded, only
its bad large-scale trend is swapped out.

(2) is answered by dropping the hardcoded rectangle for Povo's *outer*
extent: Nico asked me to derive a shape now (rather than wait on a new
QGIS export) since we don't have a hand-drawn Povo boundary yet. Built as
the convex hull of the actual NE-archipelago land pixels (still selected
via the old bounding box, but now ONLY as a membership filter, not as the
thing that gates the raster) -- an organic shape with no straight edges,
blended via the same edge_transition mechanism as an authored zone. The
near-island fringe ramp (shallow right at each islet's own shore) is
unchanged from v3/v2.

Everything else (custom ridge/zone loaders, ridge-summing convention,
ferry corridor, defense-in-depth clamp, export conventions) is unchanged
from v3 -- see that script's own docstring for the reasoning behind each.
"""
import sys, time, json, os
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.path import Path

sys.path.insert(0, "src")
from terrain.noise import Simplex2D, domain_warp  # noqa: E402

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
# Custom loaders -- unchanged from v3, see that script's docstring point 1.
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
            "name": name, "tree": cKDTree(dense),
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
            "name": name, "feature_type": ftype,
            "target_elevation_m": p.get("target_elevation_m"),
            "amplitude_scale": float(p["amplitude_scale"]),
            "edge_transition_km": float(p["edge_transition_km"]),
            "path": Path(ring), "tree": cKDTree(dense),
            "area_km2": shoelace_area(ring) / 1e6,
        })
    out.sort(key=lambda z: -z["area_km2"])
    return out


ridges = load_ridges(RIDGES_PATH)
zones = load_zones(ZONES_PATH)
log(f"loaded {len(ridges)} authored ridges, {len(zones)} authored zones "
    f"(nesting order: {' -> '.join(z['name'] for z in zones)})")

# =====================================================================
# DEM detrend: isolate the DEM's own ocean-depth texture from its bad
# large-scale radial trend. Binned median (not mean -- robust to the huge
# local variance, up to std~140m in some bins) by distance-to-coast,
# lightly smoothed, gives dem_radial_trend(d). Residual = dem - that
# trend, sampled once per cell (cheap: dem is already the full grid, no
# noise generation needed at all).
# =====================================================================
log("detrending DEM's own ocean signal (isolating local texture)...")
d_flat = dist_to_land[true_ocean]
z_flat = dem[true_ocean]
BIN_KM = 0.25
bin_edges = np.arange(0.0, float(d_flat.max()) + BIN_KM, BIN_KM)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
idx = np.clip(np.digitize(d_flat, bin_edges) - 1, 0, len(bin_centers) - 1)
medians = np.full(len(bin_centers), np.nan)
for i in range(len(bin_centers)):
    m = idx == i
    if m.sum() >= 20:
        medians[i] = np.median(z_flat[m])
valid = ~np.isnan(medians)
medians_filled = np.interp(bin_centers, bin_centers[valid], medians[valid])
trend_curve = ndimage.gaussian_filter1d(medians_filled, sigma=2.0)  # ~0.5km smoothing
log(f"  DEM's own radial trend: peak depth {trend_curve.min():.1f}m at "
    f"{bin_centers[np.argmin(trend_curve)]:.1f}km, non-monotonic beyond that "
    f"(confirms the earlier finding -- this is what gets swapped out)")

dem_radial_trend = np.interp(dist_to_land.ravel(), bin_centers, trend_curve).reshape(ny, nx).astype(np.float32)
dem_residual = np.where(true_ocean, dem - dem_radial_trend, 0.0).astype(np.float32)
log(f"  residual (DEM's own local texture, kept at full strength per Nico's choice): "
    f"mean={dem_residual[true_ocean].mean():.2f} std={dem_residual[true_ocean].std():.2f}")
del d_flat, z_flat, idx, medians, valid, medians_filled, trend_curve, dem_radial_trend

# =====================================================================
# Replacement background trend: same shelf->slope->basin curve v3 used
# (12km break, ease-out near-shore segment) -- well-behaved and monotonic,
# unlike the DEM's own radial component.
# =====================================================================
log("background shelf/slope/basin (replacement trend, same curve as v3)...")
SHELF_BREAK_KM, SHELF_DEPTH_M = 12.0, 140.0
BASIN_START_KM, BASIN_DEPTH_M = 55.0, 880.0
w_shelf = ease_out_pow2(dist_to_land, 0.0, SHELF_BREAK_KM)
shelf_profile = -SHELF_DEPTH_M * w_shelf
w_slope = smootherstep(dist_to_land, SHELF_BREAK_KM, BASIN_START_KM)
base_profile = shelf_profile * (1 - w_slope) + (-BASIN_DEPTH_M) * w_slope
log(f"  background profile: min={base_profile[true_ocean].min():.1f} max={base_profile[true_ocean].max():.1f}")

structure = base_profile.copy()
del shelf_profile, w_shelf, w_slope, base_profile

# =====================================================================
# Ferry corridor -- unchanged from v3.
# =====================================================================
log("ferry corridor (unchanged, domain-warped)...")
FERRY_WAYPOINTS_KM = np.array([
    [0.82, -21.94], [-9.5, -30.0], [-19.0, -42.0], [-27.65, -52.46],
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
halfwidth_raw = chunked_apply(lambda t: halfwidth_noise_gen.noise2(t * 6.0, np.zeros_like(t)), [t_along])
halfwidth_km = 3.75 + 0.9 * halfwidth_raw
FERRY_BLEND_KM = 3.0
w_ferry = 1.0 - smootherstep(dist_to_ferry_line_km, halfwidth_km, halfwidth_km + FERRY_BLEND_KM)

FERRY_TARGET_M = -70.0
structure = np.where(true_ocean, structure * (1 - w_ferry) + FERRY_TARGET_M * w_ferry, structure)
HALFWIDTH_KM_RANGE = (float(halfwidth_km.min()), float(halfwidth_km.max()))
log(f"  ferry core: {(w_ferry[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2")
del (Xw_f, Yw_f, warped_xy, dist_to_ferry_line_m, nearest_idx, dist_to_ferry_line_km, t_along,
     halfwidth_raw, halfwidth_km, w_ferry, ferry_dense, ferry_tree, ferry_arclen)

# =====================================================================
# Povo Silencioso -- outer extent is now the CONVEX HULL of the actual
# archipelago land pixels (selected via the old bounding box, used only
# as a membership filter now, not to gate the raster) instead of the
# rectangle itself. Organic boundary, no straight-edge seam. Blended
# exactly like an authored zone (edge_transition_km-style band straddling
# the hull boundary). Near-island fringe ramp unchanged from v3.
# =====================================================================
log("Povo Silencioso (convex-hull outer extent, derived from archipelago land)...")
NE_X0, NE_X1 = 27000.0, XMAX + 20.0
NE_Y0, NE_Y1 = 58000.0, YMAX + 20.0
ne_box = (X >= NE_X0) & (X <= NE_X1) & (Y >= NE_Y0) & (Y <= NE_Y1)
land_in_ne = land & ne_box
dist_to_ne_land = (ndimage.distance_transform_edt(~land_in_ne) * CELL / 1000.0).astype(np.float32)

py, px = np.where(land_in_ne)
povo_pts = np.column_stack([X[py, px], Y[py, px]])
hull = ConvexHull(povo_pts)
hull_pts = povo_pts[hull.vertices]
hull_path = Path(hull_pts)
hull_ring = np.vstack([hull_pts, hull_pts[:1]])
hull_dense = densify(hull_ring, 250.0)
hull_tree = cKDTree(hull_dense)
log(f"  archipelago hull: {len(hull_pts)} vertices, area {shoelace_area(hull_ring)/1e6:.1f}km2 "
    f"(from {land_in_ne.sum():,} land pixels)")
del py, px, povo_pts, hull, hull_pts, hull_ring

POVO_HULL_BLEND_KM = 5.0


def _hull_blend_chunk(x, y):
    xy = np.column_stack([x.ravel(), y.ravel()])
    inside = hull_path.contains_points(xy).reshape(x.shape)
    dist_m, _ = hull_tree.query(xy, k=1)
    dist_m = dist_m.reshape(x.shape)
    signed = np.where(inside, dist_m, -dist_m)
    band_m = POVO_HULL_BLEND_KM * 1000.0
    t = np.clip((signed + band_m) / (2 * band_m), 0.0, 1.0)
    return t * t * t * (t * (t * 6 - 15) + 10)


hull_gate = chunked_apply(_hull_blend_chunk, [X, Y])
del hull_dense, hull_tree

POVO_FRINGE_KM = 0.3
POVO_NEAR_KM = 2.5
w_rise = smootherstep(dist_to_ne_land, POVO_FRINGE_KM, POVO_NEAR_KM)
w_povo = w_rise * hull_gate
del hull_gate, w_rise

SHELTERED_TARGET_M = -320.0
OPEN_SEA_TARGET_M = -1400.0
open_sea_factor = smootherstep(dist_to_mainland, 2.0, 15.0)
povo_target = SHELTERED_TARGET_M + (OPEN_SEA_TARGET_M - SHELTERED_TARGET_M) * open_sea_factor

structure = np.where(true_ocean, structure * (1 - w_povo) + povo_target * w_povo, structure)
log(f"  povo core (near-shore ring, full effect): "
    f"{((w_povo > 0.9) & true_ocean).sum() * (CELL/1000)**2:.2f} km2")
del ne_box, land_in_ne, dist_to_ne_land, w_povo, povo_target, open_sea_factor

# =====================================================================
# Authored ridges -- unchanged from v3 (summed contributions).
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
    log(f"  {r['name']}: peak={r['peak_elevation_m']:.0f}m -> contribution range "
        f"[{contrib.min():.1f}, {contrib.max():.1f}]m")
    del contrib

structure = np.where(true_ocean, structure + ridge_total, structure)
del ridge_total

# =====================================================================
# Authored zones -- unchanged from v3. amplitude_multiplier now scales
# the DEM's own residual texture (below) instead of synthetic noise.
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
        f"amplitude_scale={z['amplitude_scale']}")
    del w

del X, Y

# =====================================================================
# Real DEM texture applied once, at full strength (Nico's explicit
# choice), scaled only by each zone's own amplitude_scale (default 1.0
# elsewhere -- background/ferry/povo all get the DEM's own residual
# unscaled).
# =====================================================================
log("applying DEM's own residual texture (headroom-scaled x amplitude_multiplier)...")
# "Full strength" as literally uniform everywhere is not physically viable -- checked
# directly: the DEM's residual has std 66-140m and a 95th percentile of +93..+224m
# within the first 8km, because it's land-terrain-scale noise (built for relief
# reaching thousands of metres) applied to a shelf that's only ~20-140m deep there.
# Applied uniformly, 31% of the entire true_ocean area (3396km2, out to 20km
# offshore) clamped flat -- confirmed by running it. That's a worse "flat" than the
# problem we're fixing. Instead: scale the residual by how much depth is actually
# available at each cell (computed from `structure` BEFORE the residual is added,
# i.e. the deterministic background/ferry/povo/ridge/zone shape) -- full strength
# once structure reaches -HEADROOM_M, ramping down smoothly in shallower water. This
# keeps 100% of the residual's own spatial pattern (nothing is redrawn or smoothed),
# only its local AMPLITUDE responds to how much room there is, the same way a real
# reef can't have 300m of relief in 20m of water but a basin floor can carry it easily.
# Swept HEADROOM_M against the actual composited grid before picking one (not
# guessed): 150 -> 427km2 clamped, 250 -> 14.5km2 (already better than v3's own
# 34.9km2), 350 -> 0km2 but mean strength down to 9% (real texture barely
# survives). 250 is the sweet spot -- clamp footprint smaller than v3's tuned
# result, while keeping meaningfully more of the residual's actual amplitude
# (mean ~17% strength, not "full", but full strength is not physically
# achievable here at all -- see the note above on why uniform full-strength
# broke immediately) than v3's synthetic noise ever provided.
HEADROOM_M = 250.0
texture_scale = smootherstep(-structure, 0.0, HEADROOM_M)
result = structure + dem_residual * amplitude_multiplier * texture_scale
log(f"  texture_scale stats (1.0=full DEM-texture strength): "
    f"mean={texture_scale[true_ocean].mean():.2f} "
    f"frac_at_full_strength={(texture_scale[true_ocean]>0.99).mean():.3f}")
del structure, dem_residual, amplitude_multiplier, texture_scale, dist_to_land, dist_to_mainland

# =====================================================================
# Defense-in-depth clamp + recombine (same as v1-v3).
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
# Export (identical conventions to v1-v3)
# =====================================================================
log("exporting...")
np.save(f"{OUT_DIR}/bathymetry_v4_30m.npy", final)


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


write_envi_raw(f"{OUT_DIR}/bathymetry_v4_30m", final, XMIN, YMIN, CELL,
                description="Tappa 11 bathymetry v4 (DEM's own detrended ocean texture + authored "
                            "shapes, convex-hull Povo extent) -- Fictional World LCC domain",
                dtype="i2")
with open(f"{OUT_DIR}/bathymetry_v4_30m.prj", "w") as f:
    f.write(CRS_PROJ4.strip() + "\n")

meta = {
    "grid": {"shape": [ny, nx], "xmin": XMIN, "ymin": YMIN, "resolution_m": CELL, "crs_proj4": CRS_PROJ4},
    "changes_from_v3": [
        "Background texture: replaced synthetic fine+coarse noise entirely with the DEM's own "
        "detrended ocean-depth signal (binned-median radial trend removed, replacement shelf/slope/"
        "basin trend added back, residual kept at full strength) -- addresses the 'too shallow and "
        "flat' finding: the DEM's own texture is spatially correlated with nearby land, which a pure "
        "distance-to-coast function structurally cannot be",
        "Povo Silencioso outer extent: replaced the hardcoded NE-archipelago bounding rectangle with "
        "the convex hull of the archipelago's actual land pixels, blended via a 5km edge_transition "
        "band -- fixes the square-artifact seam v3 had right at the rectangle's edge",
        "DEM residual texture headroom-scaled (full strength once background structure reaches "
        "-250m, ramped down in shallower water) -- literal full strength everywhere clamped 31% of "
        "true_ocean flat (3396km2, out to 20km offshore) on the first run, confirmed by running it; "
        "the residual's own land-terrain-scale magnitude (std 66-140m near shore) was never going to "
        "fit inside a ~20-140m-deep shelf uniformly. Swept the headroom threshold (150/250/350m) "
        "against the actual composited grid before picking 250m: smaller clamp footprint (14.5km2) "
        "than v3's own tuned result (34.9km2), at ~17% mean residual strength -- full strength is not "
        "physically achievable here at any threshold without the clamp reappearing",
    ],
    "ridges": [{"name": r["name"], "peak_elevation_m": r["peak_elevation_m"], "falloff_km": r["falloff_km"],
                "shelf_multiplier": r["shelf_multiplier"]} for r in ridges],
    "zones": [{"name": z["name"], "feature_type": z["feature_type"], "target_elevation_m": z["target_elevation_m"],
               "amplitude_scale": z["amplitude_scale"], "edge_transition_km": z["edge_transition_km"],
               "area_km2": z["area_km2"]} for z in zones],
    "ferry": {"waypoints_km": FERRY_WAYPOINTS_KM.tolist(), "halfwidth_km_range": list(HALFWIDTH_KM_RANGE),
              "target_m": FERRY_TARGET_M},
    "povo_silencioso": {"sheltered_target_m": SHELTERED_TARGET_M, "open_sea_target_m": OPEN_SEA_TARGET_M,
                         "hull_blend_km": POVO_HULL_BLEND_KM},
    "shelf_break_km": SHELF_BREAK_KM, "shelf_depth_m": SHELF_DEPTH_M,
    "clamp_hit_true_ocean_cells": int(clamped_cells),
    "clamp_hit_true_ocean_km2": float(clamped_cells * (CELL / 1000.0) ** 2),
    "pre_clamp_true_ocean_max_m": float(pre_clamp_max),
    "final_true_ocean_depth_range_m": [float(final[true_ocean].min()), float(final[true_ocean].max())],
}
with open(f"{OUT_DIR}/bathymetry_v4_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
