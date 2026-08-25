"""
Tappa 11 -- Bathymetry.

Builds a bathymetry raster on the EXACT same grid/CRS/extent as the Tappa 1
DEM (data/processed/dem_v3_final_30m_eroded.npy), per the locked handoff
(tappa11_bathymetry_prompt.md, combining a Tappa 9 base spec + two Scenario-
chat lore additions for the ferry corridor and Povo Silencioso).

Key finding checked directly before designing anything: the EXISTING
sub-sea-level values in dem_v3_final_30m_eroded.npy are not a designed
bathymetry at all -- they're just the same structure()+noise() terrain
function from Tappa 1 (01_tappa1_terrain.md SS3) evaluated below zero. That
function has no shelf/slope/abyssal logic and was never validated below sea
level (Tappa 1's own doc only discusses erosion, which is fluvial-only and
runs on land). Confirmed empirically: mean ocean depth is *shallower* 30-40km
from shore than 5-10km from shore (-357m vs -449m) -- the opposite of a real
continental margin. So this script does not just carve three zones into an
existing bathymetry -- it replaces ALL ocean-cell elevations with a designed
bathymetry (shelf -> slope -> basin, distance-to-coast driven, same
"structure + noise" pattern Tappa 1 used for land), then composites the
three lore-driven zones on top via the same smootherstep-blended-zone-
override pattern Tappa 8 used for its authored lithology zones. Land cells
and lake cells (non-border-connected "ocean") are left completely untouched.
"""
import sys, time, json
import numpy as np
from scipy import ndimage

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)

# ---- domain / CRS (config/parameters.yml, same as every other Tappa) ----
XMIN, XMAX, YMIN, YMAX, CELL = -65000.0, 65000.0, -80000.0, 80000.0, 30.0
CRS_PROJ4 = "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

IN_DIR = "data/processed"
OUT_DIR = "data/processed/bathymetry"
import os
os.makedirs(OUT_DIR, exist_ok=True)

log("loading DEM...")
dem = np.load(f"{IN_DIR}/dem_v3_final_30m_eroded.npy").astype(np.float32)
ny, nx = dem.shape
y_top = YMIN + ny * CELL
assert (ny, nx) == (5334, 4334)

def xy_to_rc(x, y):
    return int(round((y_top - y) / CELL)), int(round((x - XMIN) / CELL))

# ---- masks ----
land = dem > 0
ocean_raw = ~land
log(f"land fraction (unchanged input): {land.mean():.4f}")

log("labeling ocean components (border-connected = true ocean, rest = lakes)...")
lbl_ocean, n_ocean = ndimage.label(ocean_raw, structure=np.ones((3, 3)))
border_labels = set(np.unique(lbl_ocean[0, :])) | set(np.unique(lbl_ocean[-1, :])) | \
                set(np.unique(lbl_ocean[:, 0])) | set(np.unique(lbl_ocean[:, -1]))
border_labels.discard(0)
true_ocean = np.isin(lbl_ocean, list(border_labels))
lake = ocean_raw & ~true_ocean
log(f"true ocean fraction: {true_ocean.mean():.4f}  lake fraction: {lake.mean():.4f} "
    f"(lakes are left byte-for-byte unchanged -- this stage is marine bathymetry only)")

# ---- distance to nearest land cell, in km, native 30m resolution ----
log("computing distance-to-coast (EDT, native 30m)...")
dist_to_land = ndimage.distance_transform_edt(~land) * CELL / 1000.0  # km
log(f"  done, max dist-to-coast in domain: {dist_to_land.max():.2f} km")


def smootherstep(x, edge0, edge1):
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * t * (t * (t * 6 - 15) + 10)


# ---- multi-octave value-noise fBm (no external deps, same spirit as
# Tappa 1's ridged_fbm but far simpler -- this doesn't need Tappa 1's
# ridged/warp machinery, just plausible seafloor roughness) ----
def fbm_noise(shape, seed, octaves=5, persistence=0.55):
    rng = np.random.default_rng(seed)
    total = np.zeros(shape, dtype=np.float32)
    amp = 1.0
    amp_sum = 0.0
    for o in range(octaves):
        divisor = 2 ** (7 - o)  # octave 0 -> /128 (coarse), octave4 -> /8 (finer)
        low_shape = (max(4, shape[0] // divisor), max(4, shape[1] // divisor))
        low = rng.standard_normal(size=low_shape).astype(np.float32)
        zoom_factors = (shape[0] / low_shape[0], shape[1] / low_shape[1])
        up = ndimage.zoom(low, zoom_factors, order=3)[: shape[0], : shape[1]]
        total += up * amp
        amp_sum += amp
        amp *= persistence
    total /= amp_sum
    total = (total - total.mean()) / total.std()
    return total.astype(np.float32)


log("generating seafloor roughness noise field (5 octaves)...")
noise = fbm_noise((ny, nx), seed=1104, octaves=5, persistence=0.55)
log(f"  noise stats: mean={noise.mean():.3f} std={noise.std():.3f}")

# =====================================================================
# 1. BACKGROUND bathymetry: shelf -> slope -> basin, distance-to-coast driven
# =====================================================================
log("building background shelf/slope/basin profile...")
SHELF_BREAK_KM, SHELF_DEPTH_M = 18.0, 140.0
BASIN_START_KM, BASIN_DEPTH_M = 55.0, 880.0

w_shelf = smootherstep(dist_to_land, 0.0, SHELF_BREAK_KM)
shelf_profile = -SHELF_DEPTH_M * w_shelf
w_slope = smootherstep(dist_to_land, SHELF_BREAK_KM, BASIN_START_KM)
base_profile = shelf_profile * (1 - w_slope) + (-BASIN_DEPTH_M) * w_slope

# roughness envelope: calm right at the coast, roughest mid-slope
# (real continental slopes are where canyons/scarps concentrate), a bit
# calmer again on the abyssal-basin floor.
rough_shelf = 8.0 + 20.0 * smootherstep(dist_to_land, 0.0, SHELF_BREAK_KM)
rough_slope_peak = 70.0
w_toward_slope_peak = smootherstep(dist_to_land, SHELF_BREAK_KM, 35.0)
w_past_slope_peak = smootherstep(dist_to_land, 35.0, BASIN_START_KM)
roughness = rough_shelf * (1 - w_toward_slope_peak) + rough_slope_peak * w_toward_slope_peak
roughness = roughness * (1 - w_past_slope_peak) + 35.0 * w_past_slope_peak

background = base_profile + noise * roughness
# Hard floor: a true-ocean cell must never come out at/above sea level.
# Right at the coast (dist_to_land -> 0) base_profile -> 0 and the noise
# term alone can push a few cells positive -- checked directly (first run
# of this script produced a max of +35.7 m before this clip existed).
# Same "independent hard cap as defense in depth" pattern Tappa 1's
# erosion step uses (01_tappa1_terrain.md SS4) rather than trying to
# perfectly shape the roughness envelope to avoid it.
background = np.minimum(background, -0.5)
log(f"  background depth stats over true_ocean: min={background[true_ocean].min():.1f} "
    f"max={background[true_ocean].max():.1f} mean={background[true_ocean].mean():.1f}")

result = background.copy()

# =====================================================================
# 2. ZONE: bridge corridor (B_35k <-> C_25k rail crossing, provisional)
#    Handoff bbox: x in [-4.7, 11.2] km, y in [42.0, 42.5] km
#    -> shallow shelf, calm, avoids the volcanic zone (already checked:
#       the whole volcanic lithology class is the SW island, ~50+ km away)
# =====================================================================
log("zone: bridge corridor...")
BRIDGE_X0, BRIDGE_X1 = -5700.0, 12200.0   # bbox + 1km buffer either side
BRIDGE_Y0, BRIDGE_Y1 = 40800.0, 43700.0   # bbox + buffer (captures full water body)
BRIDGE_TARGET_M = -18.0
BRIDGE_CLAMP = (-35.0, -3.0)
BRIDGE_BLEND_KM = 2.0  # smootherstep taper distance back to background, beyond the box

xs = XMIN + CELL * np.arange(nx)
ys = y_top - CELL * np.arange(ny)
X, Y = np.meshgrid(xs, ys)

# distance outside the core box (0 inside, growing outside), in km
dx_out = np.maximum(np.maximum(BRIDGE_X0 - X, X - BRIDGE_X1), 0.0)
dy_out = np.maximum(np.maximum(BRIDGE_Y0 - Y, Y - BRIDGE_Y1), 0.0)
dist_outside_box_km = np.sqrt(dx_out**2 + dy_out**2) / 1000.0
w_bridge = 1.0 - smootherstep(dist_outside_box_km, 0.0, BRIDGE_BLEND_KM)

bridge_value = BRIDGE_TARGET_M + noise * 4.0  # heavily flattened relief (firm, buildable shelf)
bridge_value = np.clip(bridge_value, *BRIDGE_CLAMP)
result = np.where(true_ocean, result * (1 - w_bridge) + bridge_value * w_bridge, result)
log(f"  bridge zone core footprint: {int((w_bridge[true_ocean] > 0.99).sum())} true-ocean cells "
    f"({(w_bridge[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2)")

# =====================================================================
# 3. ZONE: ferry corridor (mainland Circulo_E3_2k <-> SW volcanic island
#    Circulo_D_20k, already the project's established cheapest boat link,
#    Tappa 9 meta: 10.8757h / 57.7km). Straight line E3_2k->D_20k crosses
#    a single clean ~41.7km ocean stretch from (0.82,-21.94) to
#    (-27.65,-52.46) km -- used as the corridor centerline.
#    HONEST FLAG (checked directly, not assumed): a land-fraction-within-
#    12km "shelter index" along this line is 0.45-0.53 near both shores
#    but drops to 0.00-0.02 (open-ocean level) across the middle third
#    (~t=0.29-0.71, roughly x[-7,-23] y[-31,-48] km) -- this crossing is
#    NOT a tight, fully sheltered strait the way the bridge corridor is.
#    See the decision doc for the numbers and the recommendation this
#    raises for Nico/Scenario chat. Depth is still kept moderate here
#    (never dropping into the existing ~-700m trough), which is the one
#    lever bathymetry actually has -- but that is a different property
#    from wave/current exposure and this script does not pretend
#    otherwise.
# =====================================================================
log("zone: ferry corridor...")
FERRY_P0 = np.array([820.0, -21940.0])
FERRY_P1 = np.array([-27650.0, -52460.0])
FERRY_TARGET_M = -70.0
FERRY_CLAMP = (-35.0, -150.0)
FERRY_HALFWIDTH_KM = 4.0
FERRY_BLEND_KM = 3.0

seg = FERRY_P1 - FERRY_P0
seg_len2 = float(seg @ seg)
PX = X - FERRY_P0[0]
PY = Y - FERRY_P0[1]
t = np.clip((PX * seg[0] + PY * seg[1]) / seg_len2, 0.0, 1.0)
proj_x = FERRY_P0[0] + t * seg[0]
proj_y = FERRY_P0[1] + t * seg[1]
dist_to_ferry_line_km = np.sqrt((X - proj_x) ** 2 + (Y - proj_y) ** 2) / 1000.0
w_ferry = 1.0 - smootherstep(dist_to_ferry_line_km, FERRY_HALFWIDTH_KM, FERRY_HALFWIDTH_KM + FERRY_BLEND_KM)

ferry_value = FERRY_TARGET_M + noise * roughness * 0.5
ferry_value = np.clip(ferry_value, FERRY_CLAMP[1], FERRY_CLAMP[0])
result = np.where(true_ocean, result * (1 - w_ferry) + ferry_value * w_ferry, result)
log(f"  ferry zone core footprint (within {FERRY_HALFWIDTH_KM}km of centerline): "
    f"{int((w_ferry[true_ocean] > 0.99).sum())} true-ocean cells "
    f"({(w_ferry[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2)")

# =====================================================================
# 4. ZONE: Povo Silencioso (NE archipelago) -- deep, hazardous water,
#    steep drop-offs close to shore. Footprint: true-ocean cells within
#    20km of any land cell inside the NE bounding box (this archipelago
#    is a dense scatter of 300+ islets, not a handful of named islands --
#    checked directly rather than hand-picking a few), which cleanly
#    covers inter-island channels and the open-sea approach alike.
# =====================================================================
log("zone: Povo Silencioso (NE archipelago)...")
NE_X0, NE_X1 = 27000.0, XMAX + 20.0
NE_Y0, NE_Y1 = 58000.0, YMAX + 20.0
POVO_TARGET_M = -1200.0
POVO_CLAMP = (-150.0, -1900.0)
POVO_NEAR_KM = 4.0    # full target reached this close to shore -> steep drop-off
POVO_FAR_KM = 20.0    # tapers back to background by this far from any island

ne_box = (X >= NE_X0) & (X <= NE_X1) & (Y >= NE_Y0) & (Y <= NE_Y1)
land_in_ne = land & ne_box
log(f"  land cells inside NE box: {land_in_ne.sum()} ({land_in_ne.sum()*(CELL/1000)**2:.2f} km2 of islets)")
dist_to_ne_land = ndimage.distance_transform_edt(~land_in_ne) * CELL / 1000.0
w_povo = 1.0 - smootherstep(dist_to_ne_land, POVO_NEAR_KM, POVO_FAR_KM)
w_povo = np.where(ne_box | (dist_to_ne_land <= POVO_FAR_KM), w_povo, 0.0)

povo_value = POVO_TARGET_M + noise * 160.0  # amplified roughness -> real steep drop-offs/channels
povo_value = np.clip(povo_value, POVO_CLAMP[1], POVO_CLAMP[0])
result = np.where(true_ocean, result * (1 - w_povo) + povo_value * w_povo, result)
log(f"  povo zone core footprint (within {POVO_NEAR_KM}km of an islet): "
    f"{int((w_povo[true_ocean] > 0.99).sum())} true-ocean cells "
    f"({(w_povo[true_ocean] > 0.99).sum() * (CELL/1000)**2:.2f} km2)")

# =====================================================================
# 5. Recombine: land + lakes untouched, true ocean = designed bathymetry
# =====================================================================
final = np.where(land, dem, np.where(true_ocean, result, dem)).astype(np.float32)
# final defense-in-depth clamp: no true-ocean cell may be >= sea level
final = np.where(true_ocean, np.minimum(final, -0.5), final).astype(np.float32)

assert not np.isnan(final).any() and not np.isinf(final).any(), "NaN/Inf in output!"
assert final[true_ocean].max() < 0.0, "a true-ocean cell came out >= sea level!"
assert np.array_equal(final[land], dem[land]), "land cells were modified!"
assert np.array_equal(final[lake], dem[lake]), "lake cells were modified!"
log("sanity checks passed: no NaN/Inf, land+lake cells byte-identical to input DEM")

log(f"FINAL true-ocean depth range: {final[true_ocean].min():.1f} .. {final[true_ocean].max():.1f} m "
    f"(input DEM ocean range was {dem[ocean_raw].min():.1f} .. {dem[ocean_raw].max():.1f} m)")

# =====================================================================
# Export -- same conventions as dem_v3 (row0=north, xmin/ymin anchor, i2 ENVI)
# =====================================================================
log("exporting...")
np.save(f"{OUT_DIR}/bathymetry_v1_30m.npy", final)


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


write_envi_raw(f"{OUT_DIR}/bathymetry_v1_30m", final, XMIN, YMIN, CELL,
                description="Tappa 11 bathymetry v1, 30m: land+lakes unchanged from dem_v3, "
                            "true-ocean cells replaced with designed shelf/slope/basin bathymetry "
                            "+ bridge/ferry/Povo Silencioso authored zones (Fictional World LCC domain)",
                dtype="i2")
with open(f"{OUT_DIR}/bathymetry_v1_30m.prj", "w") as f:
    f.write(CRS_PROJ4.strip() + "\n")

# also export an ocean-only depth layer (land/lake = NaN) for anyone who
# wants to load just the bathymetry without re-deriving the ocean mask
depth_only = np.where(true_ocean, final, np.nan).astype(np.float32)
np.save(f"{OUT_DIR}/bathymetry_depth_only_30m.npy", depth_only)

meta = {
    "grid": {"shape": [ny, nx], "xmin": XMIN, "xmax": XMAX, "ymin": YMIN, "ymax": YMAX,
             "resolution_m": CELL, "crs_proj4": CRS_PROJ4,
             "note": "Identical grid/extent/CRS to data/processed/dem_v3_final_30m_eroded.npy "
                     "(row 0 = north, col 0 = west, same convention -- see 01_tappa1_terrain.md)."},
    "method": "Land and lake cells are byte-identical to dem_v3_final_30m_eroded.npy, unchanged. "
              "True-ocean cells (border-connected component of the <=0 mask -- same 'true ocean vs "
              "lake' distinction Tappa 9/10 established) got a newly designed bathymetry: a "
              "distance-to-coast-driven shelf(0-18km)/slope(18-55km)/basin(>55km) depth profile "
              "via smootherstep blending (same blending function Tappa 1 used for its zone "
              "overrides), textured with a 5-octave value-noise fBm field (seed=1104), then three "
              "lore-driven zones (bridge/ferry/Povo Silencioso) composited on top via the same "
              "smootherstep-weighted override pattern Tappa 8 used for its authored lithology "
              "zones. This replaces the input DEM's existing sub-sea-level values entirely for "
              "true-ocean cells -- those values were never a designed bathymetry (see this "
              "script's module docstring for the direct check that established that).",
    "background_profile": {"shelf_break_km": SHELF_BREAK_KM, "shelf_depth_m": SHELF_DEPTH_M,
                            "basin_start_km": BASIN_START_KM, "basin_depth_m": BASIN_DEPTH_M},
    "zones": {
        "bridge_corridor": {
            "source": "Handoff SS1 (from Tappa 9 base spec)",
            "bbox_km": {"x": [BRIDGE_X0/1000, BRIDGE_X1/1000], "y": [BRIDGE_Y0/1000, BRIDGE_Y1/1000]},
            "target_depth_m": BRIDGE_TARGET_M, "clamp_m": list(BRIDGE_CLAMP),
        },
        "ferry_corridor": {
            "source": "Handoff SS2 (Scenario chat addition, scenario_reference.md SS18)",
            "centerline_km": {"p0": (FERRY_P0/1000).tolist(), "p1": (FERRY_P1/1000).tolist()},
            "halfwidth_km": FERRY_HALFWIDTH_KM, "target_depth_m": FERRY_TARGET_M,
            "clamp_m": list(FERRY_CLAMP),
        },
        "povo_silencioso": {
            "source": "Handoff SS3 (Scenario chat addition, scenario_reference.md SS19)",
            "ne_box_km": {"x": [NE_X0/1000, NE_X1/1000], "y": [NE_Y0/1000, NE_Y1/1000]},
            "near_km": POVO_NEAR_KM, "far_km": POVO_FAR_KM, "target_depth_m": POVO_TARGET_M,
            "clamp_m": list(POVO_CLAMP),
        },
    },
    "true_ocean_fraction": float(true_ocean.mean()),
    "lake_fraction": float(lake.mean()),
    "final_true_ocean_depth_range_m": [float(final[true_ocean].min()), float(final[true_ocean].max())],
    "input_dem_ocean_depth_range_m": [float(dem[ocean_raw].min()), float(dem[ocean_raw].max())],
}
with open(f"{OUT_DIR}/bathymetry_v1_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
