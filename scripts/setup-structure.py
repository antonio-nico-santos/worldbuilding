"""
Sets up the folder structure for the World-building GIS project.

"""
 
import os
 
FOLDERS = [
    "config",
    "src/terrain",
    "src/climate",
    "src/hydrology",
    "src/biomes",
    "src/suitability",
    "notebooks",
    "data/raw",
    "data/processed",
    "data/exports",
    "qgis",
    "docs/decisions",
]
 
GITIGNORE = """\
# Generated / regenerable data — never versioned, always reproducible from src/ + config/
data/processed/
data/raw/
 
# Python
__pycache__/
*.pyc
.venv/
venv/
 
# Jupyter
.ipynb_checkpoints/
 
# QGIS
*.qgz~
*.gpkg-shm
*.gpkg-wal
 
# OS
.DS_Store
Thumbs.db
"""
 
REQUIREMENTS = """\
numpy
rasterio
geopandas
opensimplex
pyyaml
grass-session
jupyter
"""
 
PARAMETERS_YML = """\
# Single source of truth for pipeline parameters.
# Every script in src/ should read from this file rather than hardcoding values,
# so the whole pipeline stays reproducible from one place.
 
domain:
  height_km: 210
  width_km: 100
  resolution_m: 30
 
crs:
  projection: "Lambert Conformal Conic"
  reference_point:
    lat: -44.0
    lon: 42.0
  # standard parallels for the LCC still to be chosen once the exact
  # north-south extent of the domain (in degrees) is set
 
validation_reference:
  region: "South Island, New Zealand"
  transect: null  # TBD: full coast-to-coast vs west-side-only vs Fiordland-only
  data_sources:
    - "NIWA"
    - "WorldClim / CHELSA"
 
terrain:
  seed: null            # pick any integer once you start Tappa 1
  noise_octaves: null
  noise_frequency: null
  erosion_iterations: null
  spine_path: null       # list of [x, y] control points defining the mountain skeleton (can branch)
 
climate:
  lapse_rate_c_per_1000m: -6.5
  wind_direction_deg: null   # TBD before Tappa 2
  continentality_decay: null
 
# Tappa 6 (suitability) and Tappa 7 (urban zoom) intentionally left empty —
# criteria not yet defined.
"""
 
README = """\
# World-building GIS — Portfolio Case Study 2
 
Procedurally generated fictional terrain, climate, hydrology and biomes,
built as a from-scratch GIS pipeline (not a pre-built world generator) to
demonstrate technical GIS skills for a freelance portfolio.
 
Real-world reference used for climate model validation: South Island, New
Zealand (Southern Alps).
 
See `docs/decisions/` for the per-stage planning summaries and
`config/parameters.yml` for the current pipeline parameters.
 
This repository is intentionally kept separate from the portfolio site
repository (`gis-portfolio`) — only final lightweight exports (simplified
GeoJSON, pre-rendered images) are copied into the site's `public/data/`.
"""
 
 
def main():
    for folder in FOLDERS:
        os.makedirs(folder, exist_ok=True)
        print(f"created: {folder}/")
 
    files = {
        ".gitignore": GITIGNORE,
        "requirements.txt": REQUIREMENTS,
        "config/parameters.yml": PARAMETERS_YML,
        "README.md": README,
    }
 
    for path, content in files.items():
        if os.path.exists(path):
            print(f"skipped (already exists): {path}")
            continue
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"created: {path}")
 
    print("\nDone. Next: fill in the still-open values in config/parameters.yml as you settle them.")
 
 
if __name__ == "__main__":
    main()