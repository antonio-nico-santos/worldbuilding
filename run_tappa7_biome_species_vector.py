"""
Tappa 7 -- vectorize the biome raster and join the regular/resident-bucket
species x biome percent table (locked in docs/decisions/07_tappa7_regional_
scenario.md sec.7, delivered as tappa7_fauna_biome_percent.xlsx) as attribute
columns, split by landmass (mainland vs. the SW island vs. smaller
unanalyzed landmasses -- same biome_id values are reused across landmasses,
but percent values differ, so a plain per-biome dissolve would be wrong).

Percent values below are hand-transcribed from tappa7_fauna_biome_percent.xlsx
("Species x Biome %" and "Mainland -> Island %" sheets) -- that workbook is
the source of truth; re-transcribe here if it changes.

Output -> data/processed/fauna/biome_species_percent.geojson
  One polygon per (biome, landmass) combination, every mainland/island
  species as its own _pct column, plus n_species_<category> and
  n_species_total species-richness columns (added this pass) for a quick
  biodiversity-distribution check per feature. High-threat and domesticated-
  animal buckets are NOT included here, in either the _pct or the n_species_
  columns -- those are raster suitability composites, not per-biome
  percentages, and don't fit this column model (see the other
  data/processed/fauna/*.npy outputs for those). Because of that,
  n_species_total is bucket-scoped (max 17), not a whole-world species count
  -- see the SPECIES_CATEGORY/CATEGORIES comment block below for the caveat
  in full, including the >0 presence threshold this counting uses.

CRS caveat (project-wide, not new to this file): the geometry is in the
project's LCC proj4 CRS (embedded as a legacy GeoJSON 'crs' member, matching
data/processed/suitability/circulo_candidate_sites.geojson's convention).
Modern GDAL/geopandas ignores that legacy member on read and reports
EPSG:4326 -- confirmed this happens on circulo_candidate_sites.geojson too,
so it's a pre-existing project convention limitation, not something this
script introduces. In QGIS: use Layer Properties -> assign the CRS manually
from the proj4 string in the 'crs' member if it doesn't auto-detect.
"""
import json
import os

import numpy as np
import rasterio
import geopandas as gpd
from rasterio.features import shapes
from scipy import ndimage

BASE = 'data/processed'
OUT = 'data/processed/fauna'
CRS_PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
             "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

BIOME_NAMES = {
    1: 'Permanent Snow & Ice', 2: 'Alpine Fellfield', 3: 'Alpine Tundra',
    4: 'Subalpine Wet Forest', 5: 'Subalpine Woodland', 6: 'Subalpine Dry Scrub',
    7: 'Temperate Forest', 8: 'Woodland / Shrubland', 9: 'Lowland Steppe / Grassland',
}

MAINLAND_PCT = {
    'Grassmothers':  {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0,'Subalpine Woodland':0,'Subalpine Dry Scrub':0.15,'Temperate Forest':0,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':1},
    'Blacknose':     {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0,'Subalpine Dry Scrub':0,'Temperate Forest':1,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':0},
    'Tailstand':     {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0,'Temperate Forest':0.15,'Woodland / Shrubland':1,'Lowland Steppe / Grassland':0.15},
    'Flashfrog':     {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0.15,'Subalpine Wet Forest':1,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0,'Temperate Forest':0.15,'Woodland / Shrubland':0,'Lowland Steppe / Grassland':0},
    'Scattermouse':  {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0,'Temperate Forest':0.15,'Woodland / Shrubland':1,'Lowland Steppe / Grassland':0.15},
    'Snaketail':     {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0,'Subalpine Dry Scrub':0.15,'Temperate Forest':1,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':1},
    'Furypack':      {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0,'Subalpine Dry Scrub':0.15,'Temperate Forest':1,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':1},
    'Cryburrow':     {'Permanent Snow & Ice':0.15,'Alpine Fellfield':1,'Alpine Tundra':1,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0.15,'Temperate Forest':0,'Woodland / Shrubland':0,'Lowland Steppe / Grassland':0},
    'Deergoat':      {'Permanent Snow & Ice':0.15,'Alpine Fellfield':1,'Alpine Tundra':1,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0.15,'Temperate Forest':0,'Woodland / Shrubland':0,'Lowland Steppe / Grassland':0},
    'Treefox':       {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0.15,'Temperate Forest':1,'Woodland / Shrubland':1,'Lowland Steppe / Grassland':1},
    'Mudlizard':     {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0.15,'Subalpine Wet Forest':1,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0,'Temperate Forest':1,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':0},
    'Quillhog':      {'Permanent Snow & Ice':0,'Alpine Fellfield':0.8,'Alpine Tundra':0.8,'Subalpine Wet Forest':0.8,'Subalpine Woodland':0.8,'Subalpine Dry Scrub':0.8,'Temperate Forest':0.8,'Woodland / Shrubland':0.8,'Lowland Steppe / Grassland':0.8},
    'Farsmell':      {'Permanent Snow & Ice':0,'Alpine Fellfield':0.8,'Alpine Tundra':0.8,'Subalpine Wet Forest':0.8,'Subalpine Woodland':0.8,'Subalpine Dry Scrub':0.8,'Temperate Forest':0.8,'Woodland / Shrubland':0.8,'Lowland Steppe / Grassland':0.8},
    'Meatcleaner':   {'Permanent Snow & Ice':0,'Alpine Fellfield':0.7,'Alpine Tundra':0.7,'Subalpine Wet Forest':0.7,'Subalpine Woodland':0.7,'Subalpine Dry Scrub':0.7,'Temperate Forest':0.7,'Woodland / Shrubland':0.7,'Lowland Steppe / Grassland':0.7},
    'Trinketbird':   {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0,'Subalpine Woodland':0,'Subalpine Dry Scrub':0,'Temperate Forest':0.15,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':1},
    'Rustowl (mainland ecotype)': {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0.15,'Subalpine Woodland':0.15,'Subalpine Dry Scrub':0.15,'Temperate Forest':1,'Woodland / Shrubland':1,'Lowland Steppe / Grassland':1},
    'Clicksnake':    {'Permanent Snow & Ice':0,'Alpine Fellfield':0,'Alpine Tundra':0,'Subalpine Wet Forest':0,'Subalpine Woodland':0,'Subalpine Dry Scrub':0,'Temperate Forest':0,'Woodland / Shrubland':0.15,'Lowland Steppe / Grassland':1},
}

ISLAND_PCT = {
    'Grassmothers':  {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Blacknose':     {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Flashfrog':     {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Scattermouse':  {'Woodland / Shrubland': 0.50, 'Lowland Steppe / Grassland': 0.08},
    'Snaketail':     {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Furypack':      {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Treefox':       {'Woodland / Shrubland': 0.25, 'Lowland Steppe / Grassland': 0.25},
    'Quillhog':      {'Woodland / Shrubland': 0.20, 'Lowland Steppe / Grassland': 0.20},
    'Cryburrow':     {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Deergoat':      {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Farsmell':      {'Woodland / Shrubland': 0.72, 'Lowland Steppe / Grassland': 0.72},
    'Meatcleaner':   {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Mudlizard':     {'Woodland / Shrubland': 0.80, 'Lowland Steppe / Grassland': 0.65},
    'Trinketbird':   {'Woodland / Shrubland': 0,    'Lowland Steppe / Grassland': 0},
    'Clicksnake':    {'Woodland / Shrubland': 0.15, 'Lowland Steppe / Grassland': 1.0},  # parity w/ mainland
    'Rustowl (island ecotype)': {'Woodland / Shrubland': 1.0, 'Lowland Steppe / Grassland': 1.0},
}

SHARED_SPECIES = ['Grassmothers', 'Blacknose', 'Flashfrog', 'Scattermouse', 'Snaketail', 'Furypack',
                   'Cryburrow', 'Deergoat', 'Treefox', 'Mudlizard', 'Quillhog', 'Farsmell', 'Meatcleaner',
                   'Trinketbird', 'Clicksnake']

# Category = the niche each species was designed under in docs/decisions/07_tappa7_regional_
# scenario.md sec.7 ("Status" list) / scenario_reference.md, NOT a formal taxonomic grouping --
# e.g. "bird" here means the bird-niche design pass, "mesopredator" the mesopredator-niche pass.
# Only the regular/resident bucket's 17 species are covered (see module docstring) -- Reaper and
# Twinshadows (high-threat bucket, tundra/alpine-apex-predator and mesopredator niches
# respectively) are DELIBERATELY excluded from these counts, same as they're excluded from the
# _pct columns above. A "total_species_count" column here is NOT total world biodiversity --
# it's biodiversity within this one bucket's data model. Flagging this explicitly rather than
# letting the column name imply more than it delivers.
SPECIES_CATEGORY = {
    'Grassmothers': 'grassland_large_grazer',
    'Blacknose': 'mid_size_browser', 'Tailstand': 'mid_size_browser',
    'Flashfrog': 'reptile_amphibian', 'Clicksnake': 'reptile_amphibian', 'Mudlizard': 'reptile_amphibian',
    'Scattermouse': 'small_mammal', 'Quillhog': 'small_mammal', 'Snaketail': 'small_mammal',
    'Treefox': 'mesopredator', 'Furypack': 'mesopredator',
    'Trinketbird': 'bird', 'Rustowl': 'bird',
    'Cryburrow': 'alpine_prey_base', 'Deergoat': 'alpine_prey_base',
    'Farsmell': 'decomposer_scavenger', 'Meatcleaner': 'decomposer_scavenger',
}
CATEGORIES = ['grassland_large_grazer', 'mid_size_browser', 'reptile_amphibian', 'small_mammal',
              'mesopredator', 'bird', 'alpine_prey_base', 'decomposer_scavenger']

# "Present" = pct > 0, matching this table's own vocabulary (0% is a deliberate, real "confirmed
# absent" value in the archetype method, not a placeholder) -- see the Specialist/Named-multi-
# biome/generalist archetypes in 07_tappa7_regional_scenario.md sec.7. This is a real modeling
# choice, not the only defensible one: it counts marginal 15%-neighbor presence the same as a
# 100%-core biome, so it reads as "species detectable here at all", not "core range here". A
# stricter threshold (e.g. >=50%) would give a materially different, more conservative map --
# worth knowing which question you're actually asking before reading the counts.


def main():
    biome_id = np.load(f'{BASE}/biomes/biome_id.npy').astype(np.int32)
    land = biome_id != 0
    res_x, res_y = 120.0, 120.0
    xmin, ymax = -65000.0, 80000.0
    transform = rasterio.transform.from_origin(xmin, ymax, res_x, res_y)

    labeled, n = ndimage.label(land, structure=np.ones((3, 3), dtype=int))
    sizes = ndimage.sum(land, labeled, index=range(1, n + 1))
    order = np.argsort(sizes)[::-1]
    mainland_label, island_label = order[0] + 1, order[1] + 1
    landmass = np.zeros_like(biome_id, dtype=np.int32)
    landmass[labeled == mainland_label] = 1
    landmass[labeled == island_label] = 2
    key = np.where(land, biome_id * 10 + landmass, 0).astype(np.int32)

    geoms = [{'properties': {'key': int(v)}, 'geometry': geom}
             for geom, v in shapes(key, mask=(key != 0), transform=transform)]
    gdf = gpd.GeoDataFrame.from_features(geoms, crs=CRS_PROJ4)
    dissolved = gdf.dissolve(by='key').reset_index()
    dissolved['biome_id'] = dissolved['key'] // 10
    dissolved['landmass_code'] = dissolved['key'] % 10
    dissolved['biome_name'] = dissolved['biome_id'].map(BIOME_NAMES)
    dissolved['landmass'] = dissolved['landmass_code'].map(
        {0: 'other (unanalyzed minor landmass)', 1: 'mainland', 2: 'SW island'})
    dissolved['area_km2'] = (dissolved.geometry.area / 1e6).round(2)

    for sp in SHARED_SPECIES:
        col = f'{sp}_pct'
        vals = []
        for _, row in dissolved.iterrows():
            if row['landmass_code'] == 1:
                v = MAINLAND_PCT.get(sp, {}).get(row['biome_name'])
            elif row['landmass_code'] == 2:
                v = ISLAND_PCT.get(sp, {}).get(row['biome_name'])
            else:
                v = None
            vals.append(v)
        dissolved[col] = vals

    # Tailstand: mainland only -- island presence is authorial (human introduction),
    # not this join. Left null on island rows rather than 0 (0 would wrongly read as
    # "confirmed absent via this model").
    vals, notes = [], []
    for _, row in dissolved.iterrows():
        if row['landmass_code'] == 1:
            vals.append(MAINLAND_PCT.get('Tailstand', {}).get(row['biome_name']))
            notes.append(None)
        elif row['landmass_code'] == 2:
            vals.append(None)
            notes.append('Not this join -- island presence is authorial (human introduction by ferry), see 07_tappa7 doc')
        else:
            vals.append(None)
            notes.append(None)
    dissolved['Tailstand_pct'] = vals
    dissolved['Tailstand_note'] = notes

    # Rustowl: mainland and island are separate ecotypes/populations, kept as two columns.
    vals_m, vals_i = [], []
    for _, row in dissolved.iterrows():
        vals_m.append(MAINLAND_PCT.get('Rustowl (mainland ecotype)', {}).get(row['biome_name']) if row['landmass_code'] == 1 else None)
        vals_i.append(ISLAND_PCT.get('Rustowl (island ecotype)', {}).get(row['biome_name']) if row['landmass_code'] == 2 else None)
    dissolved['Rustowl_mainland_ecotype_pct'] = vals_m
    dissolved['Rustowl_island_ecotype_pct'] = vals_i

    dissolved.loc[dissolved['landmass_code'] == 0, 'note'] = (
        'Smaller landmass distinct from the already-locked SW island (794.6 km2) -- surfaced by this '
        'vectorization, never analyzed for mainland/island percentages before. Species columns left '
        'null, not assumed.'
    )

    # Species-richness columns, per category and total. "Present" = pct > 0 (see CATEGORIES note
    # above). For landmass_code == 0 ("other", unanalyzed minor landmass) every species column is
    # null by design -- counting those as 0 would silently read as "confirmed zero biodiversity"
    # when the truth is "not modeled." Left null there too, same policy as the _pct columns.
    def species_col(sp, landmass_code):
        if sp == 'Rustowl':
            if landmass_code == 1:
                return 'Rustowl_mainland_ecotype_pct'
            elif landmass_code == 2:
                return 'Rustowl_island_ecotype_pct'
            return None
        return f'{sp}_pct'

    for cat in CATEGORIES:
        cat_species = [sp for sp, c in SPECIES_CATEGORY.items() if c == cat]
        counts = []
        for _, row in dissolved.iterrows():
            if row['landmass_code'] == 0:
                counts.append(None)
                continue
            n = 0
            for sp in cat_species:
                col = species_col(sp, row['landmass_code'])
                if col is None:
                    continue
                v = row.get(col)
                if v is not None and not (isinstance(v, float) and np.isnan(v)) and v > 0:
                    n += 1
            counts.append(n)
        dissolved[f'n_species_{cat}'] = counts

    dissolved['n_species_total'] = dissolved[[f'n_species_{c}' for c in CATEGORIES]].sum(
        axis=1, skipna=False).astype('Int64')
    # skipna=False so a row with any null category (the landmass_code==0 rows, where every
    # category is null) totals to null too, rather than silently summing to 0.

    dissolved = dissolved.drop(columns=['key'])

    os.makedirs(OUT, exist_ok=True)
    out_path = f'{OUT}/biome_species_percent.geojson'
    dissolved.to_file(out_path, driver='GeoJSON')

    # Re-inject the legacy 'crs' member -- GDAL's GeoJSON writer drops it per RFC7946,
    # but the rest of this project's GeoJSON files carry it (see module docstring caveat).
    d = json.load(open(out_path))
    d['crs'] = {"type": "proj4", "properties": {"proj4": CRS_PROJ4}}
    d['name'] = 'tappa7_biome_species_percent'
    json.dump(d, open(out_path, 'w'), indent=2)

    print(f'Saved {out_path}, {len(dissolved)} features')
    print(dissolved[['biome_name', 'landmass', 'area_km2']].sort_values(['landmass', 'biome_name']).to_string())


if __name__ == '__main__':
    main()
