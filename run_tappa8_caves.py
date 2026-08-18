import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np
from scipy import ndimage

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import CLASS_VOLCANIC, identify_landmasses
from geomorphology.caves import (
    resample_nearest, lava_tube_candidates, vent_proximity_weight,
    talus_pseudokarst_candidates, glacier_moulin_candidates,
    sea_cave_candidates, distance_to_ocean_km, compute_slope,
)

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- cave candidates ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]

OUTPUT_DIR = "data/processed/geomorphology"
os.makedirs(OUTPUT_DIR, exist_ok=True)

log("loading DEM + lithology + stream_mask...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask = dem > 0
lithology = np.load(f"{OUTPUT_DIR}/lithology.npy")
stream_mask = np.load("data/processed/hydrology/stream_mask.npy")
ny, nx = dem.shape

log("identifying mainland (excludes SW Island + every other islet)...")
labeled, mainland_label, sw_island_label, areas_km2 = identify_landmasses(land_mask, cellsize_m=RES_M)
mainland_mask = labeled == mainland_label
log(f"mainland={areas_km2[mainland_label]:.1f} km2, sw_island={areas_km2[sw_island_label]:.1f} km2")

log("computing slope (native 30m)...")
slope_pct = compute_slope(dem, RES_M)
mainland_land_slope = slope_pct[mainland_mask]
steep_threshold_pct = float(np.percentile(mainland_land_slope, 75.0))
log(f"steep threshold (p75 of mainland land slope) = {steep_threshold_pct:.2f}%")

log("computing distance-to-stream (native 30m)...")
dist_to_stream_km = ndimage.distance_transform_edt(~stream_mask, sampling=(RES_M, RES_M)) / 1000.0

log("computing distance-to-ocean (native 30m)...")
dist_to_ocean_km = distance_to_ocean_km(land_mask, RES_M)

log("resampling permanent_snow_mask (120m climate grid -> native 30m)...")
# permanent_snow_mask ships only as ENVI .bin/.hdr (no .npy on disk) -- read
# it directly rather than relying on a file that isn't there.
with open("data/processed/climate/permanent_snow_mask.hdr") as f:
    hdr_txt = f.read()
def hdr_val(key):
    for line in hdr_txt.splitlines():
        if line.strip().startswith(key):
            return line.split("=", 1)[1].strip()
    raise KeyError(key)
psm_samples = int(hdr_val("samples"))
psm_lines = int(hdr_val("lines"))
psm_mapinfo = hdr_val("map info").strip("{}").split(",")
psm_xmin = float(psm_mapinfo[3])
psm_ytop = float(psm_mapinfo[4])
psm_cellsize = float(psm_mapinfo[5])
psm = np.fromfile("data/processed/climate/permanent_snow_mask.bin", dtype="<i2").reshape(psm_lines, psm_samples)
log(f"permanent_snow_mask: shape=({psm_lines},{psm_samples}) xmin={psm_xmin} ytop={psm_ytop} cellsize={psm_cellsize} sum={int((psm>0).sum())}")

psm_native = resample_nearest(
    psm, coarse_xmin=psm_xmin, coarse_ymax=psm_ytop, coarse_cellsize=psm_cellsize,
    ny=ny, nx=nx, xmin=domain["xmin"], ymax=domain["ymax"], cellsize=RES_M,
)

log("computing the four cave-candidate layers...")
lava_tubes = lava_tube_candidates(lithology, land_mask)
talus = talus_pseudokarst_candidates(slope_pct, dist_to_stream_km, mainland_mask, steep_threshold_pct, stream_buffer_km=0.5)
glacier_moulin = glacier_moulin_candidates(psm_native, land_mask)
sea_caves = sea_cave_candidates(land_mask, dist_to_ocean_km, slope_pct, steep_threshold_pct, coastal_buffer_km=0.5)

log("computing lava-tube vent-proximity weight (addition beyond the literal decision text)...")
import json as _json
with open("data/input/geothermal.geojson") as f:
    geo = _json.load(f)
vents = [(feat["geometry"]["coordinates"][0], feat["geometry"]["coordinates"][1], feat["properties"]["falloff_km"]) for feat in geo["features"]]
y = domain["ymax"] - (np.arange(ny) + 0.5) * RES_M
x = domain["xmin"] + (np.arange(nx) + 0.5) * RES_M
xx, yy = np.meshgrid(x, y)
vent_weight = vent_proximity_weight(xx, yy, vents, RES_M)

cell_km2 = (RES_M / 1000.0) ** 2
areas = {
    "lava_tube_candidate_km2": float(lava_tubes.sum() * cell_km2),
    "talus_pseudokarst_candidate_km2": float(talus.sum() * cell_km2),
    "glacier_moulin_candidate_km2": float(glacier_moulin.sum() * cell_km2),
    "sea_cave_candidate_km2": float(sea_caves.sum() * cell_km2),
}
log(f"areas: {areas}")

log("exporting rasters...")
for name, arr in [
    ("cave_lava_tube", lava_tubes),
    ("cave_talus_pseudokarst", talus),
    ("cave_glacier_moulin", glacier_moulin),
    ("cave_sea_cave", sea_caves),
]:
    np.save(f"{OUTPUT_DIR}/{name}.npy", arr)
    write_envi_raw(
        f"{OUTPUT_DIR}/{name}", arr.astype(np.int16),
        xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
        description=f"Tappa 8 cave candidate eligibility mask: {name}",
        dtype="i2",
    )
    write_prj(f"{OUTPUT_DIR}/{name}.prj", CRS_PROJ4)

np.save(f"{OUTPUT_DIR}/lava_tube_vent_weight.npy", vent_weight)
np.save(f"{OUTPUT_DIR}/slope_pct_30m.npy", slope_pct.astype(np.float32))
np.save(f"{OUTPUT_DIR}/dist_to_stream_km_30m.npy", dist_to_stream_km.astype(np.float32))
np.save(f"{OUTPUT_DIR}/dist_to_ocean_km_30m.npy", dist_to_ocean_km.astype(np.float32))

meta = {
    "resolution_m": RES_M,
    "steep_threshold_pct_p75_mainland": steep_threshold_pct,
    "stream_buffer_km": 0.5,
    "coastal_buffer_km": 0.5,
    "areas_km2": areas,
    "mainland_area_km2": float(areas_km2[mainland_label]),
    "sw_island_area_km2": float(areas_km2[sw_island_label]),
    "vent_count": len(vents),
    "schist_fracture_caves": "NOT built, parked pending Nico's visual review (per 07_tappa7 S2)",
}
with open(f"{OUTPUT_DIR}/tappa8_caves_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
