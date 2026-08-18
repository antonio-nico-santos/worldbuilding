import sys
sys.path.insert(0, 'src')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

OUT = "data/processed/geomorphology"
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land = dem > 0

lava_v1 = np.load(f"{OUT}/cave_lava_tube.npy")
lava_v2 = np.load(f"{OUT}/cave_lava_tube_v2.npy")
glacier_v1 = np.load(f"{OUT}/cave_glacier_moulin.npy")
glacier_v2 = np.load(f"{OUT}/cave_glacier_moulin_v2.npy")
talus = np.load(f"{OUT}/cave_talus_pseudokarst.npy")
sea = np.load(f"{OUT}/cave_sea_cave.npy")

bg = np.where(land, dem, np.nan)

fig, axes = plt.subplots(2, 2, figsize=(16, 16))

def panel(ax, mask_v1, mask_v2, title):
    ax.imshow(bg, cmap="gray", alpha=0.35)
    if mask_v1 is not None:
        ov1 = np.ma.masked_where(~mask_v1, mask_v1)
        ax.imshow(ov1, cmap=ListedColormap(["#ff5555"]), alpha=0.35)
    ov2 = np.ma.masked_where(~mask_v2, mask_v2)
    ax.imshow(ov2, cmap=ListedColormap(["#00ff00"]), alpha=0.85)
    ax.set_title(title, fontsize=12)
    ax.axis("off")

panel(axes[0,0], lava_v1, lava_v2, f"Lava tube: v1 red={lava_v1.sum()*0.0009:.1f}km2 -> v2 green={lava_v2.sum()*0.0009:.1f}km2")
panel(axes[0,1], glacier_v1, glacier_v2, f"Glacier/moulin: v1 red={glacier_v1.sum()*0.0009:.1f}km2 -> v2 green={glacier_v2.sum()*0.0009:.1f}km2")
panel(axes[1,0], None, talus, f"Talus/pseudokarst (unchanged): {talus.sum()*0.0009:.1f}km2")
panel(axes[1,1], None, sea, f"Sea caves (unchanged): {sea.sum()*0.0009:.1f}km2")

plt.tight_layout()
plt.savefig("/tmp/wb_caves_v2_overview.png", dpi=110)
plt.close()

# --- lava tube v3 (re-masked against lithology_v5's island volcanic/basin_fill split) vs v2 ---
lava_v3 = np.load(f"{OUT}/cave_lava_tube_v3.npy")
lith_v5 = np.load(f"{OUT}/lithology_v5.npy")

r0, r1, c0, c1 = 3700, 5334, 100, 1620  # island bbox, see identify_landmasses row/col range check
bg_island = np.where(land, dem, np.nan)[r0:r1, c0:c1]
colors = ["#0b3d63", "#d9c48a", "#8a8a5c", "#7a4a8a", "#b23a2f"]
lith_cmap = ListedColormap(colors)

fig, axes = plt.subplots(1, 3, figsize=(22, 10))
axes[0].imshow(lith_v5[r0:r1, c0:c1], cmap=lith_cmap, vmin=0, vmax=4, interpolation="nearest")
axes[0].set_title("Island lithology v5 (volcanic core / basin_fill margin)", fontsize=11)
axes[0].axis("off")

panel(axes[1], lava_v2[r0:r1, c0:c1], lava_v3[r0:r1, c0:c1],
      f"Lava tube candidates: v2/lithology_v2 red={lava_v2.sum()*0.0009:.1f}km2 -> "
      f"v3/lithology_v5 green={lava_v3.sum()*0.0009:.1f}km2")

diff_lost = lava_v2[r0:r1, c0:c1] & ~lava_v3[r0:r1, c0:c1]
diff_gained = lava_v3[r0:r1, c0:c1] & ~lava_v2[r0:r1, c0:c1]
axes[2].imshow(bg_island, cmap="gray", alpha=0.35)
axes[2].imshow(np.ma.masked_where(~diff_lost, diff_lost), cmap=ListedColormap(["#ff2222"]), alpha=0.8)
axes[2].imshow(np.ma.masked_where(~diff_gained, diff_gained), cmap=ListedColormap(["#22ff22"]), alpha=0.8)
axes[2].set_title("Delta only: red=lost (now outside volcanic zone), green=gained (recalibrated slope threshold)", fontsize=10)
axes[2].axis("off")

plt.tight_layout()
plt.savefig("/tmp/wb_caves_v3_lava_tube_island.png", dpi=110)
plt.close()

print("done")
