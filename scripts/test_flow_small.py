import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import numpy as np
from hydrology.flow import priority_flood_d8, accumulate_flow

# --- 1. synthetic correctness tests -----------------------------------
print("=== synthetic tests ===")

# a) single pit in the middle of a bowl, sea on the border
dem = np.array([
    [0, 0, 0, 0, 0],
    [0, 5, 5, 5, 0],
    [0, 5, 1, 5, 0],   # pit at center (1), surrounded by rim of 5
    [0, 5, 5, 5, 0],
    [0, 0, 0, 0, 0],
], dtype=np.float64)
seed = dem <= 0
filled, recv, order = priority_flood_d8(dem, seed, epsilon=0.01)
print("pit test filled:\n", filled)
center = 2 * 5 + 2
assert filled[2, 2] > 5.0, "pit should be filled above its rim"
# check strictly decreasing along path from center to a seed
i = center
path = [i]
while recv[i] != i:
    path.append(recv[i])
    i = recv[i]
elevs = [filled.ravel()[p] for p in path]
assert all(elevs[k] > elevs[k + 1] for k in range(len(elevs) - 1)), f"not strictly decreasing: {elevs}"
print("  OK: pit filled, strictly-decreasing path to sea:", elevs)

# b) flat plateau (no epsilon tie-break issue: must not get stuck / must
#    still get a single consistent direction for every plateau cell)
dem2 = np.full((7, 7), 10.0)
dem2[0, :] = 0.0  # north edge = sea
seed2 = np.zeros_like(dem2, dtype=bool)
seed2[0, :] = True
filled2, recv2, order2 = priority_flood_d8(dem2, seed2, epsilon=0.01)
# every cell must reach the sea in a finite number of hops with strictly
# decreasing elevation (guards against a receiver cycle on the flat)
for r in range(7):
    for c in range(7):
        i = r * 7 + c
        steps = 0
        while recv2[i] != i and steps < 100:
            i = recv2[i]
            steps += 1
        assert recv2[i] == i, f"cell ({r},{c}) never reached a sink (cycle?)"
print("  OK: flat plateau -- every cell reaches sea, no cycles")

# c) weighted accumulation sanity: total accumulated weight at the sea
#    outlets must equal total weight of the whole domain (mass conservation)
w = np.ones_like(dem2)
accum = accumulate_flow(recv2, order2, w)
outlet_total = sum(accum[i] for i in range(49) if recv2[i] == i)
print("  weighted accum: total domain weight =", w.sum(), " sum at outlets =", outlet_total)
assert abs(outlet_total - w.sum()) < 1e-6, "mass not conserved at outlets"
print("  OK: mass conservation holds")

# --- 2. real crop: timing + qualitative check -------------------------
print("\n=== real DEM crop ===")
dem_full = np.load("/tmp/wb/data/processed/dem_v3_final_30m_eroded.npy")
print("full DEM shape", dem_full.shape)

for size in (300, 600, 1000):
    crop = dem_full[1500:1500 + size, 1500:1500 + size].astype(np.float64)
    seed = np.zeros_like(crop, dtype=bool)
    seed |= crop <= 0.0
    seed[0, :] = seed[-1, :] = seed[:, 0] = seed[:, -1] = True
    t0 = time.time()
    filled, recv, order = priority_flood_d8(crop, seed, epsilon=1e-4)
    dt = time.time() - t0
    n = crop.size
    print(f"  {size}x{size} ({n:>9,} cells): {dt:6.2f}s  ({dt/n*1e6:.3f} us/cell)")

    # sanity: land fraction, no NaN, filled >= original everywhere
    assert np.all(filled >= crop - 1e-9), "filled dem went below original somewhere"
    assert not np.any(np.isnan(filled))
    print(f"    land fraction in crop: {(crop>0).mean():.3f}, "
          f"max fill raise: {(filled-crop).max():.2f} m")
