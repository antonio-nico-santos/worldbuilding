import sys
sys.path.insert(0, 'src')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from terrain.skeleton import load_geojson, build_zone_fields
from geomorphology.lithology import _grid_xy
from params import load_params

OUT = "data/processed/geomorphology"
dom = load_params("config/parameters.yml")["domain"]

lith_v4 = np.load(f"{OUT}/lithology_v4.npy")
lith_v5 = np.load(f"{OUT}/lithology_v5.npy")
schist_grade_v5 = np.load(f"{OUT}/schist_grade_v5.npy")
jade_suit_v5 = np.load(f"{OUT}/jade_suitable_v5.npy")
jade_pods_v5 = np.load(f"{OUT}/jade_pods_v5.npy")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
ny, nx = dem.shape

colors = ["#0b3d63", "#d9c48a", "#8a8a5c", "#7a4a8a", "#b23a2f"]
cmap = ListedColormap(colors)

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
for ax, lith, title in zip(axes, [lith_v4, lith_v5],
                           ["v4: DEM ridge network + falloff + noise", "v5: DEM-native (elevation + relief), no ridge geometry"]):
    ax.imshow(lith, cmap=cmap, vmin=0, vmax=4, interpolation="nearest")
    ax.set_title(title, fontsize=13)
    ax.axis("off")
plt.tight_layout()
plt.savefig("/tmp/wb_lithology_v4_vs_v5.png", dpi=110)
plt.close()

# overlay authored zone outlines to visually confirm the SE-plains fix
zones = build_zone_fields(load_geojson("data/input/terrain_zones.geojson"))
fig, ax = plt.subplots(figsize=(12, 14))
ax.imshow(lith_v5, cmap=cmap, vmin=0, vmax=4, interpolation="nearest")
for zone in zones:
    if zone.feature_type != "plateau":
        continue
    pts = zone.path.vertices
    cols = (pts[:, 0] - dom["xmin"]) / dom["resolution_m"]
    rows = (dom["ymax"] - pts[:, 1]) / dom["resolution_m"]
    ax.plot(cols, rows, color="cyan", linewidth=1.5)
    cx, cy = cols.mean(), rows.mean()
    ax.text(cx, cy, zone.name, color="white", fontsize=9, ha="center",
            bbox=dict(facecolor="black", alpha=0.5, pad=1))
ax.set_title("v5 with authored plateau/plains zone outlines overlaid (QA only, not a classification input)", fontsize=12)
ax.axis("off")
plt.tight_layout()
plt.savefig("/tmp/wb_lithology_v5_zones_overlay.png", dpi=110)
plt.close()

fig, ax = plt.subplots(figsize=(11, 11))
hillshade_bg = np.where(dem > 0, dem, np.nan)
ax.imshow(hillshade_bg, cmap="gray", alpha=0.35)
grade_masked = np.ma.masked_where(lith_v5 != 3, schist_grade_v5)
ax.imshow(grade_masked, cmap="Purples", vmin=0, vmax=1, alpha=0.7)
suit_overlay = np.ma.masked_where(~jade_suit_v5, jade_suit_v5)
ax.imshow(suit_overlay, cmap=ListedColormap(["#00e0e0"]), alpha=0.5)
pods_overlay = np.ma.masked_where(~jade_pods_v5, jade_pods_v5)
ax.imshow(pods_overlay, cmap=ListedColormap(["#00ff00"]), alpha=0.95)
ax.set_title("Jade v5: schist grade (purple) / suitability zone p80+ (cyan) / 10 discrete pods (green)", fontsize=12)
ax.axis("off")
plt.tight_layout()
plt.savefig("/tmp/wb_jade_v5_overview.png", dpi=110)
plt.close()

# island-focused zoom panel -- Nico's request to extend elevation+relief logic to the island.
# island occupies roughly rows 3773:5268, cols 183:1538 in this raster (identify_landmasses bbox).
r0, r1, c0, c1 = 3700, 5334, 100, 1620
fig, axes = plt.subplots(1, 2, figsize=(16, 11))
axes[0].imshow(lith_v5[r0:r1, c0:c1], cmap=cmap, vmin=0, vmax=4, interpolation="nearest")
axes[0].set_title("Island, v5 lithology (volcanic core / basin_fill margin, island-own relief p25 threshold)", fontsize=11)
axes[0].axis("off")
relief_crop = relief_2km_for_island = np.load(f"{OUT}/relief_2km.npy")[r0:r1, c0:c1]
dem_crop = dem[r0:r1, c0:c1]
relief_masked = np.ma.masked_where(dem_crop <= 0, relief_crop)
im = axes[1].imshow(relief_masked, cmap="viridis")
axes[1].set_title("Island local relief (2km window) -- basin_fill = below island's own p25", fontsize=11)
axes[1].axis("off")
plt.colorbar(im, ax=axes[1], fraction=0.04, label="relief (m)")

zones_i = build_zone_fields(load_geojson("data/input/terrain_zones.geojson"))
for zone in zones_i:
    if zone.name != "Island plateau":
        continue
    pts = zone.path.vertices
    cols_z = (pts[:, 0] - dom["xmin"]) / dom["resolution_m"] - c0
    rows_z = (dom["ymax"] - pts[:, 1]) / dom["resolution_m"] - r0
    for ax in axes:
        ax.plot(cols_z, rows_z, color="cyan", linewidth=1.5)
plt.tight_layout()
plt.savefig("/tmp/wb_lithology_v5_island_zoom.png", dpi=110)
plt.close()

print("class areas v4:", {c: int((lith_v4==c).sum()) for c in range(5)})
print("class areas v5:", {c: int((lith_v5==c).sum()) for c in range(5)})
print("done")
