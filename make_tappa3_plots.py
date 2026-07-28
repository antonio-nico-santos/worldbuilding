#!/usr/bin/env python3
"""Validation/overview plots for Tappa 3 — see docs/decisions/03_tappa3_snow.md."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CDIR = Path("data/processed/climate")
OUT = Path("docs/decisions/assets/03_tappa3_snow")
OUT.mkdir(parents=True, exist_ok=True)

elev = np.load(CDIR / "surface_elevation_m.npy")
land = np.load(CDIR / "land_mask.npy")
n_snow = np.load(CDIR / "months_with_snow.npy")
season = np.load(CDIR / "seasonality_index_c.npy")
balance = np.load(CDIR / "mass_balance_mm.npy")
perm = np.load(CDIR / "permanent_snow_mask.npy")
naive_perm = (np.load(CDIR / "temperature_monthly_c.npy").max(axis=0) < 0.0) & land
meta = json.loads((CDIR / "tappa3_snow_meta.json").read_text())

sea_mask = ~land
hill = np.ma.masked_where(sea_mask, elev)

# --- Fig 1: overview maps ------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 6))

ax = axes[0]
im = ax.imshow(np.ma.masked_where(sea_mask, n_snow), cmap="Blues", vmin=0, vmax=12)
ax.imshow(np.ma.masked_where(land, elev), cmap="ocean", alpha=0.5)
plt.colorbar(im, ax=ax, fraction=0.046, label="months with modelled snow")
ax.set_title("Months with snow (precip-aware)")
ax.axis("off")

ax = axes[1]
im = ax.imshow(np.ma.masked_where(sea_mask, season), cmap="magma", vmin=7.5, vmax=11)
ax.imshow(np.ma.masked_where(land, elev), cmap="ocean", alpha=0.5)
plt.colorbar(im, ax=ax, fraction=0.046, label="seasonality index (C)")
ax.set_title("Seasonality (warmest - coldest month)")
ax.axis("off")

ax = axes[2]
comp = np.zeros(elev.shape, dtype=np.uint8)
comp[land] = 1
comp[naive_perm] = 2
comp[perm & ~naive_perm] = 3
# distinct hues by design: sea (cool grey) vs. land (near-white) vs. the two
# permanent-snow categories in warm colors that pop against both backgrounds
cmap = matplotlib.colors.ListedColormap(["#9ecae1", "#f7fbff", "#6baed6", "#e6550d"])
im = ax.imshow(comp, cmap=cmap, vmin=0, vmax=3)
ax.set_title("Permanent snow: naive (T-only) vs mass-balance")
handles = [
    matplotlib.patches.Patch(color="#9ecae1", label="sea"),
    matplotlib.patches.Patch(color="#f7fbff", label="land, not permanent"),
    matplotlib.patches.Patch(color="#6baed6", label="naive AND mass-balance"),
    matplotlib.patches.Patch(color="#e6550d", label="mass-balance ONLY (new)"),
]
ax.legend(handles=handles, loc="lower left", fontsize=8, framealpha=0.9)
ax.axis("off")

plt.tight_layout()
fig.savefig(OUT / "03_overview.png", dpi=150)
plt.close(fig)

# --- Fig 2: balance vs elevation, windward vs leeward ---------------------
fig, ax = plt.subplots(figsize=(8, 6))
for side, color in [("wet", "tab:blue"), ("dry", "tab:red")]:
    table = meta["elevation_balance_bands"][side]
    z = [row[0] for row in table]
    b = [row[1] for row in table]
    ax.plot(b, z, color=color, label=f"{'windward (wet tercile)' if side=='wet' else 'leeward (dry tercile)'}")
ax.axvline(0, color="k", lw=0.8)
ela_w = meta["summary"]["ela_windward_m"]
ela_l = meta["summary"]["ela_leeward_m"]
if ela_w:
    ax.axhline(ela_w, color="tab:blue", ls="--", lw=1, label=f"ELA windward = {ela_w:.0f} m")
if ela_l:
    ax.axhline(ela_l, color="tab:red", ls="--", lw=1, label=f"ELA leeward = {ela_l:.0f} m")
ax.set_xlabel("annual mass balance, mm w.e. (binned mean)")
ax.set_ylabel("elevation (m)")
ax.set_title("Mass balance vs elevation: windward/leeward ELA")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
plt.tight_layout()
fig.savefig(OUT / "03_ela_bands.png", dpi=150)
plt.close(fig)

# --- Fig 3: sensitivity comparison ----------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
labels = ["mass-balance\n(this stage)", "naive\n(T-only, Tappa 2)"]
locked = [meta["summary"]["permanent_snow_area_km2"], meta["summary"]["naive_permanent_snow_area_km2"]]
sens = [meta["sensitivity_lapse_seasonal_amplitude_0"]["permanent_snow_area_km2"],
        meta["sensitivity_lapse_seasonal_amplitude_0"]["naive_permanent_snow_area_km2"]]
x = np.arange(2)
w = 0.35
ax.bar(x - w/2, locked, w, label="lapse_seasonal_amplitude = 0.3 (locked)")
ax.bar(x + w/2, sens, w, label="lapse_seasonal_amplitude = 0.0 (sensitivity)")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("permanent snow area (km2)")
ax.set_title("Sensitivity to lapse_seasonal_amplitude_c_per_km")
ax.legend(fontsize=8)
for i, (lv, sv) in enumerate(zip(locked, sens)):
    ax.annotate(f"{lv:.0f}", (i - w/2, lv), ha="center", va="bottom", fontsize=8)
    ax.annotate(f"{sv:.0f}", (i + w/2, sv), ha="center", va="bottom", fontsize=8)
plt.tight_layout()
fig.savefig(OUT / "03_sensitivity.png", dpi=150)
plt.close(fig)

print("wrote", list(OUT.glob("*.png")))
