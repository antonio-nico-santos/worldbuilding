"""
Tappa 11 -- Bathymetry, v2.

v1 (run_tappa11_bathymetry.py) used pure distance-fields (distance-to-coast,
distance-to-line, distance-to-nearest-islet) for all three zones. Nico's
review of v1 in QGIS found three real problems, all traceable to that
choice:

1. Povo Silencioso read as flat/uniformly shallow -- a single radial
   distance-to-nearest-islet falloff with one target depth has no reason to
   look like a real, structured seafloor, and doesn't distinguish the open-
   sea-facing side (should be dramatic) from the mainland-facing side
   (should be calmer).
2. The ferry corridor read as a geometric "stripe" -- a constant-halfwidth
   buffer around a straight 2-point line is exactly a stripe, by
   construction.
3. The bridge corridor (and by extension anywhere depth was meant to read
   as "shallow but real") came out shallower than intended (-13..-22m
   against a -20..-40m ask).

Fix: switch to the SAME authored-shape toolkit Tappa 1 actually uses
(`src/terrain/skeleton.py`'s `RidgeField`/`ZoneField`, `src/terrain/
noise.py`'s `domain_warp` + `Simplex2D`/`ridged_fbm`) instead of ad hoc
distance fields -- hand-placed trench polylines (Gaussian decay off a line,
exactly like a Tappa 1 ridge but with a negative "peak"), domain-warped
query coordinates (turns a straight-line buffer into an organic channel,
the same fix Tappa 1 used for the coastline), and two-band noise (a coarse
band for real seamount/basin-scale structure, a fine band for texture --
same "this isn't one clean scale" lesson as Tappa 1 SS10a's spectral
analysis).
"""
import sys, time, json, os
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, "src")
from terrain.noise import Simplex2D, ridged_fbm, domain_warp  # noqa: E402

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)

XMIN, XMAX, YMIN, YMAX, CELL = -65000.0, 65000.0, -80000.0, 80000.0, 30.0
CRS_PROJ4 = "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
IN_DIR = "data/processed"
OUT_DIR = "data/processed/bathymetry"
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
del ocean_raw, lbl_ocean, border_labels  # this sandbox has a hard ~6GB memcg limit; free ASAP

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


# This sandbox's memcg has a hard ~6GB limit. `Simplex2D.noise2` (src/terrain/
# noise.py) force-casts to float64 and holds ~20 full-grid temporaries alive
# at once inside a single call; on the full 5334x4334 grid (23.1M cells,
# ~185MB per float64 array) that peaks well past 6GB and gets OOM-killed --
# observed directly, the unchunked coarse ridged_fbm call below was killed
# by the kernel OOM killer on the first attempt. `chunked_apply` calls a
# noise-evaluating function (ridged_fbm / domain_warp / a raw .noise2 call)
# on horizontal row-bands instead of the whole grid at once -- same math,
# each row still evaluated at its true (x, y) -- capping transient memory to
# one band's worth (~250 rows -> ~7MB/array instead of ~185MB/array).
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


# =====================================================================
# Two-band noise, same "one scale isn't enough" lesson as Tappa 1 SS10a.
# Coarse band: real ridged_fbm (Tappa 1's own noise family) at a long
# wavelength -- gives genuine seamount/basin-scale structure, not just
# grain. Fine band: the simpler multi-octave value-noise fBm from v1, kept
# for high-frequency seafloor texture.
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
# Background: shelf -> slope -> basin (unchanged shape from v1 -- this
# wasn't the part Nico flagged), now fed by both noise bands.
# =====================================================================
log("background shelf/slope/basin...")
SHELF_BREAK_KM, SHELF_DEPTH_M = 18.0, 140.0
BASIN_START_KM, BASIN_DEPTH_M = 55.0, 880.0
w_shelf = smootherstep(dist_to_land, 0.0, SHELF_BREAK_KM)
shelf_profile = -SHELF_DEPTH_M * w_shelf
w_slope = smootherstep(dist_to_land, SHELF_BREAK_KM, BASIN_START_KM)
base_profile = shelf_profile * (1 - w_slope) + (-BASIN_DEPTH_M) * w_slope

rough_shelf = 8.0 + 20.0 * smootherstep(dist_to_land, 0.0, SHELF_BREAK_KM)
w_toward_peak = smootherstep(dist_to_land, SHELF_BREAK_KM, 35.0)
w_past_peak = smootherstep(dist_to_land, 35.0, BASIN_START_KM)
roughness = rough_shelf * (1 - w_toward_peak) + 70.0 * w_toward_peak
roughness = roughness * (1 - w_past_peak) + 35.0 * w_past_peak

background = base_profile + fine * roughness * 0.75 + coarse * roughness * 0.5
background = np.minimum(background, -0.5)
log(f"  background: min={background[true_ocean].min():.1f} max={background[true_ocean].max():.1f}")
result = background.copy()
del base_profile, shelf_profile, w_shelf, w_slope, rough_shelf, w_toward_peak, w_past_peak, background

# =====================================================================
# ZONE: bridge corridor -- retargeted per Nico's review to -20..-40m
# (was -18m target / -13..-22m actual in v1, too shallow).
# =====================================================================
log("zone: bridge corridor...")
BRIDGE_X0, BRIDGE_X1 = -5700.0, 12200.0
BRIDGE_Y0, BRIDGE_Y1 = 40800.0, 43700.0
BRIDGE_TARGET_M = -30.0
BRIDGE_CLAMP = (-40.0, -20.0)
BRIDGE_BLEND_KM = 2.0

dx_out = np.maximum(np.maximum(BRIDGE_X0 - X, X - BRIDGE_X1), 0.0)
dy_out = np.maximum(np.maximum(BRIDGE_Y0 - Y, Y - BRIDGE_Y1), 0.0)
dist_outside_bridge_km = np.sqrt(dx_out**2 + dy_out**2) / 1000.0
w_bridge = 1.0 - smootherstep(dist_outside_bridge_km, 0.0, BRIDGE_BLEND_KM)

bridge_value = BRIDGE_TARGET_M + fine * 4.0
bridge_value = np.clip(bridge_value, *BRIDGE_CLAMP)
result = np.where(true_ocean, result * (1 - w_bridge) + bridge_value * w_bridge, result)
log(f"  bridge core: {(w_bridge[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2")
del dx_out, dy_out, dist_outside_bridge_km, w_bridge, bridge_value

# =====================================================================
# ZONE: ferry corridor -- authored as a 4-waypoint smooth polyline (not a
# straight 2-point line) with domain-warped query coordinates (so the
# corridor edges are wavy, not a geometric stripe) and a width that
# breathes along its length instead of a constant halfwidth.
# =====================================================================
log("zone: ferry corridor (authored, warped)...")
# Waypoints: mainland shore -> two interior bends -> island shore. The
# interior bends are placed off the straight line by a few km, in the
# direction of this strait's own genuinely deeper water (checked in v1:
# the straight line's deepest existing water sat slightly west of its
# midpoint), so the authored curve tracks a plausible real channel rather
# than being an arbitrary wiggle.
FERRY_WAYPOINTS_KM = np.array([
    [0.82, -21.94],
    [-9.5, -30.0],
    [-19.0, -42.0],
    [-27.65, -52.46],
]) * 1000.0

# densify the waypoint polyline (piecewise-linear through the 4 points)
def densify(coords, max_spacing_m):
    out = [coords[0]]
    for i in range(len(coords) - 1):
        p0, p1 = coords[i], coords[i + 1]
        seg_len = np.hypot(*(p1 - p0))
        n_steps = max(1, int(np.ceil(seg_len / max_spacing_m)))
        for s in range(1, n_steps + 1):
            out.append(p0 + (p1 - p0) * (s / n_steps))
    return np.array(out)

ferry_dense = densify(FERRY_WAYPOINTS_KM, 250.0)
ferry_arclen = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(ferry_dense, axis=0).T))])
ferry_total_len = ferry_arclen[-1]
ferry_tree = cKDTree(ferry_dense)

# domain-warp the query coordinates (Tappa 1's own function, new seed) --
# this is what breaks the corridor out of a perfect constant-width capsule.
log("  applying domain warp to ferry query coordinates...")
Xw_f, Yw_f = chunked_apply(
    lambda x, y: domain_warp(x, y, seed=4402, wavelength_m=6000.0, amplitude_m=1800.0, octaves=3),
    [X, Y], n_outputs=2,
)
warped_xy = np.column_stack([Xw_f.ravel(), Yw_f.ravel()])
dist_to_ferry_line_m, nearest_idx = ferry_tree.query(warped_xy, k=1)
dist_to_ferry_line_km = dist_to_ferry_line_m.reshape(ny, nx) / 1000.0
t_along = (ferry_arclen[nearest_idx] / ferry_total_len).reshape(ny, nx)

# halfwidth breathes along the route: 3-4.5km, smooth, seeded so it's
# reproducible and not correlated with the warp noise above
halfwidth_noise_gen = Simplex2D(seed=4403)
halfwidth_raw = chunked_apply(
    lambda t: halfwidth_noise_gen.noise2(t * 6.0, np.zeros_like(t)),
    [t_along],
)
halfwidth_km = 3.75 + 0.9 * halfwidth_raw
FERRY_BLEND_KM = 3.0
w_ferry = 1.0 - smootherstep(dist_to_ferry_line_km, halfwidth_km, halfwidth_km + FERRY_BLEND_KM)

FERRY_TARGET_M = -70.0
FERRY_CLAMP = (-35.0, -150.0)
ferry_value = FERRY_TARGET_M + fine * roughness * 0.4 + coarse * 20.0
ferry_value = np.clip(ferry_value, FERRY_CLAMP[1], FERRY_CLAMP[0])
result = np.where(true_ocean, result * (1 - w_ferry) + ferry_value * w_ferry, result)
HALFWIDTH_KM_RANGE = (float(halfwidth_km.min()), float(halfwidth_km.max()))
log(f"  ferry core: {(w_ferry[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2, "
    f"halfwidth range {HALFWIDTH_KM_RANGE[0]:.2f}-{HALFWIDTH_KM_RANGE[1]:.2f}km")
del (Xw_f, Yw_f, warped_xy, dist_to_ferry_line_m, nearest_idx, dist_to_ferry_line_km, t_along,
     halfwidth_raw, halfwidth_km, w_ferry, ferry_value, roughness, ferry_dense, ferry_tree, ferry_arclen)

# =====================================================================
# ZONE: Povo Silencioso -- three parts, all authored:
#  (a) a smooth open-sea/mainland-facing ASYMMETRY driven by real
#      distance-to-mainland (not distance-to-nearest-islet alone) --
#      shallower on the sheltered/mainland side, deep on the open-sea side.
#  (b) two hand-placed TRENCH lines (exact RidgeField math, negative
#      "peak") on the genuinely open-sea-facing arc -- a real, structured
#      deep channel, not just a bigger noise number.
#  (c) the two-band noise (coarse seamount/basin structure + fine texture)
#      already built above, weighted into the zone.
# =====================================================================
log("zone: Povo Silencioso (NE archipelago, authored)...")
NE_X0, NE_X1 = 27000.0, XMAX + 20.0
NE_Y0, NE_Y1 = 58000.0, YMAX + 20.0
ne_box = (X >= NE_X0) & (X <= NE_X1) & (Y >= NE_Y0) & (Y <= NE_Y1)
land_in_ne = land & ne_box
dist_to_ne_land = (ndimage.distance_transform_edt(~land_in_ne) * CELL / 1000.0).astype(np.float32)

POVO_FRINGE_KM = 0.3     # thin shallow fringe right at an islet's shore
POVO_NEAR_KM = 2.5       # full effect reached this close -- steep drop-off
POVO_FAR_KM = 22.0       # tapers back to background by here
w_rise = smootherstep(dist_to_ne_land, POVO_FRINGE_KM, POVO_NEAR_KM)
w_fall = 1.0 - smootherstep(dist_to_ne_land, POVO_NEAR_KM, POVO_FAR_KM)
w_povo = w_rise * w_fall
w_povo = np.where(ne_box | (dist_to_ne_land <= POVO_FAR_KM), w_povo, 0.0)

# (a) mainland-facing asymmetry: shallow near the mainland, deep far from
# it. Checked directly (see decision doc): dist_to_mainland inside this
# box ranges 0.03-17.4km, mean 6.3km -- a real, usable gradient, not a
# degenerate one.
SHELTERED_TARGET_M = -320.0
OPEN_SEA_TARGET_M = -1400.0
open_sea_factor = smootherstep(dist_to_mainland, 2.0, 15.0)
povo_base_target = SHELTERED_TARGET_M + (OPEN_SEA_TARGET_M - SHELTERED_TARGET_M) * open_sea_factor
log(f"  open-sea factor stats in box: min={open_sea_factor[ne_box].min():.2f} "
    f"max={open_sea_factor[ne_box].max():.2f} mean={open_sea_factor[ne_box].mean():.2f}")

# (b) two hand-placed trench lines, exact RidgeField-style Gaussian decay
# (see src/terrain/skeleton.py RidgeField.contribution) but with a
# NEGATIVE peak -- reused directly rather than reimplemented differently.
# Chunked over row-bands (via chunked_apply) for the same ~6GB-memcg reason
# as the noise calls above: column_stack-ing the full-grid (X, Y) into one
# (23.1M, 2) query array, plus the cKDTree distance/index temporaries, is
# several hundred MB of transient float64 that isn't worth holding all at
# once when it's this cheap to band instead.
def trench_contribution(X_, Y_, waypoints_m, peak_m, falloff_km, shelf_multiplier=2.5):
    dense = densify(np.array(waypoints_m, dtype=np.float64), 250.0)
    tree = cKDTree(dense)
    k = np.log(2.0)
    falloff_m = falloff_km * 1000.0
    shelf_m = falloff_m * shelf_multiplier

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

# Trench A -- "East approach": hugs the domain's open-ocean edge, checked
# directly to be the genuinely farthest-from-mainland water at every
# latitude band spanning the archipelago (see decision doc).
TRENCH_A_KM = [[63.5, 60.0], [64.8, 68.0], [63.0, 75.0], [64.5, 80.0]]
trench_a = trench_contribution(X, Y, (np.array(TRENCH_A_KM) * 1000.0).tolist(),
                                peak_m=-750.0, falloff_km=3.0, shelf_multiplier=3.0)

# Trench B -- "inter-island channel": threads the genuinely open water
# between the two largest N islands (labels 21 and 27 in the connected-
# components pass, x~43-51/y~74-80 and x~56-65/y~70-80 respectively).
TRENCH_B_KM = [[48.0, 72.0], [53.0, 75.0], [58.0, 78.0]]
trench_b = trench_contribution(X, Y, (np.array(TRENCH_B_KM) * 1000.0).tolist(),
                                peak_m=-550.0, falloff_km=2.2, shelf_multiplier=2.8)

trench_total = trench_a + trench_b  # both negative, sum deepens further where they overlap
del trench_a, trench_b

povo_target = povo_base_target + trench_total
POVO_CLAMP = (-150.0, -1900.0)
povo_value = povo_target + fine * (40.0 + 60.0 * open_sea_factor) + coarse * (60.0 + 220.0 * open_sea_factor)
povo_value = np.clip(povo_value, POVO_CLAMP[1], POVO_CLAMP[0])
result = np.where(true_ocean, result * (1 - w_povo) + povo_value * w_povo, result)
log(f"  povo core (near-shore ring, full effect): "
    f"{((w_povo > 0.9) & true_ocean).sum() * (CELL/1000)**2:.2f} km2")
del (ne_box, land_in_ne, dist_to_ne_land, w_rise, w_fall, w_povo, open_sea_factor,
     povo_base_target, trench_total, povo_target, povo_value, X, Y, coarse, fine,
     dist_to_land, dist_to_mainland)

# =====================================================================
# Recombine + defense-in-depth clamp (same as v1)
# =====================================================================
final = np.where(land, dem, np.where(true_ocean, result, dem)).astype(np.float32)
final = np.where(true_ocean, np.minimum(final, -0.5), final).astype(np.float32)

assert not np.isnan(final).any() and not np.isinf(final).any()
assert np.array_equal(final[land], dem[land])
assert np.array_equal(final[lake], dem[lake])
assert final[true_ocean].max() < 0.0
log("sanity checks passed")
log(f"FINAL true-ocean depth range: {final[true_ocean].min():.1f} .. {final[true_ocean].max():.1f} m")

# =====================================================================
# Export (identical conventions to v1)
# =====================================================================
log("exporting...")
np.save(f"{OUT_DIR}/bathymetry_v2_30m.npy", final)

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

write_envi_raw(f"{OUT_DIR}/bathymetry_v2_30m", final, XMIN, YMIN, CELL,
                description="Tappa 11 bathymetry v2 (authored: RidgeField-style trenches, "
                            "domain-warped ferry corridor, mainland-distance asymmetric Povo "
                            "Silencioso, two-band noise) -- Fictional World LCC domain",
                dtype="i2")
with open(f"{OUT_DIR}/bathymetry_v2_30m.prj", "w") as f:
    f.write(CRS_PROJ4.strip() + "\n")

meta = {
    "grid": {"shape": [ny, nx], "xmin": XMIN, "ymin": YMIN, "resolution_m": CELL, "crs_proj4": CRS_PROJ4},
    "changes_from_v1": [
        "Bridge corridor retargeted -30m (-40..-20 clamp), was -18m (-35..-3)",
        "Ferry corridor: 4-waypoint authored centerline (was 2-point straight line) + domain-warped "
        "query coordinates (Tappa 1's own domain_warp, seed 4402) + breathing halfwidth 2.85-4.65km "
        "(was constant 4km) -- replaces the v1 'stripe' with an organic channel",
        "Povo Silencioso: target depth now driven by real distance-to-mainland (asymmetric: "
        "-320m sheltered side, -1400m open-sea side, smootherstep 2-15km), plus two hand-placed "
        "RidgeField-style trench lines (-750m/-550m peaks) on the genuinely open-sea-facing arc, "
        "plus a coarse ridged-fBm noise band (seamount/basin-scale, wavelength 9km) layered under "
        "the existing fine-texture band -- was a single uniform target + one noise band",
    ],
    "bridge": {"target_m": BRIDGE_TARGET_M, "clamp_m": list(BRIDGE_CLAMP)},
    "ferry": {"waypoints_km": FERRY_WAYPOINTS_KM.tolist(), "halfwidth_km_range": list(HALFWIDTH_KM_RANGE),
              "target_m": FERRY_TARGET_M, "clamp_m": list(FERRY_CLAMP)},
    "povo_silencioso": {
        "sheltered_target_m": SHELTERED_TARGET_M, "open_sea_target_m": OPEN_SEA_TARGET_M,
        "trench_a_waypoints_km": TRENCH_A_KM, "trench_a_peak_m": -750.0,
        "trench_b_waypoints_km": TRENCH_B_KM, "trench_b_peak_m": -550.0,
        "clamp_m": list(POVO_CLAMP),
    },
    "final_true_ocean_depth_range_m": [float(final[true_ocean].min()), float(final[true_ocean].max())],
}
with open(f"{OUT_DIR}/bathymetry_v2_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
