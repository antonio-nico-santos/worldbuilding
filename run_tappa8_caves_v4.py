import sys, time, json
sys.path.insert(0, 'src')
import numpy as np
from scipy import ndimage

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.caves import karst_cave_candidates
from geomorphology.lithology import CLASS_MARBLE, CLASS_SEDIMENTARY_LIMESTONE

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- cave candidates v4 (NEW: karst/dissolution caves, unlocked by lithology v6's "
    "marble + sedimentary_limestone classes; the other four cave types are UNCHANGED -- none of "
    "them take marble/sedimentary_limestone as an input, so they are not re-run here) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUTPUT_DIR = "data/processed/geomorphology"

log("loading DEM + lithology v6 + stream_mask...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask = dem > 0
lithology_v6 = np.load(f"{OUTPUT_DIR}/lithology_v6.npy")
stream_mask = np.load("data/processed/hydrology/stream_mask.npy").astype(bool)

log("computing dist_to_stream_km (same method as the v1 caves script: distance_transform_edt "
    "on ~stream_mask -- not committed as a static file per S6's continuous-field policy, "
    "cheap to regenerate)...")
dist_to_stream_km = ndimage.distance_transform_edt(~stream_mask, sampling=(RES_M, RES_M)) / 1000.0

STREAM_BUFFER_KM = 0.5  # reuses talus/pseudokarst's own convention (Tappa 6 water_gentle_km), not a new number
log(f"karst eligibility = (marble OR sedimentary_limestone) AND dist_to_stream_km <= {STREAM_BUFFER_KM}km "
    "(no slope condition -- see karst_cave_candidates docstring for why)...")
karst = karst_cave_candidates(lithology_v6, land_mask, dist_to_stream_km, stream_buffer_km=STREAM_BUFFER_KM)

cell_km2 = (RES_M / 1000.0) ** 2
marble_km2 = float((lithology_v6 == CLASS_MARBLE).sum() * cell_km2)
limestone_km2 = float((lithology_v6 == CLASS_SEDIMENTARY_LIMESTONE).sum() * cell_km2)
karst_km2 = float(karst.sum() * cell_km2)
karst_in_marble_km2 = float((karst & (lithology_v6 == CLASS_MARBLE)).sum() * cell_km2)
karst_in_limestone_km2 = float((karst & (lithology_v6 == CLASS_SEDIMENTARY_LIMESTONE)).sum() * cell_km2)
log(f"soluble lithology: marble={marble_km2:.2f} km2, sedimentary_limestone={limestone_km2:.2f} km2")
log(f"karst candidates: {karst_km2:.2f} km2 total ({karst_in_marble_km2:.2f} km2 in marble = "
    f"{100*karst_in_marble_km2/marble_km2:.1f}% of marble, {karst_in_limestone_km2:.2f} km2 in "
    f"sedimentary_limestone = {100*karst_in_limestone_km2/limestone_km2:.1f}% of that class)")

log("exporting raster...")
np.save(f"{OUTPUT_DIR}/cave_karst.npy", karst)
write_envi_raw(f"{OUTPUT_DIR}/cave_karst", karst.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description=("Tappa 8 cave candidate eligibility mask: karst/dissolution (marble or "
                  f"sedimentary_limestone, within {STREAM_BUFFER_KM}km of a stream)"),
    dtype="u1")
write_prj(f"{OUTPUT_DIR}/cave_karst.prj", CRS_PROJ4)

meta = {
    "method": ("soluble lithology (marble OR sedimentary_limestone, from lithology_v6) AND "
               "dist_to_stream_km <= stream_buffer_km. No slope condition -- unlike the other four "
               "cave types, real NZ karst (Waitomo, Harwoods Hole, Punakaiki) doesn't consistently "
               "correlate with steep relief, so adding one would be an invented constraint."),
    "stream_buffer_km": STREAM_BUFFER_KM,
    "soluble_lithology_km2": {"marble": marble_km2, "sedimentary_limestone": limestone_km2},
    "karst_candidate_km2": karst_km2,
    "karst_in_marble_km2": karst_in_marble_km2,
    "karst_in_sedimentary_limestone_km2": karst_in_limestone_km2,
}
with open(f"{OUTPUT_DIR}/tappa8_caves_v4_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
