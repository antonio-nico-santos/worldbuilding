"""
Tappa 8 -- one-time local decompression step.

The device-commit path this data was pushed through (Claude's remote-device
bridge, not a project tool) rejects files over 20 MB and calls over 100 MB
total. lithology_v5, jade_pods_v5, jade_suitable_v5, basin_fill_grounded_v5,
and all four cave_* rasters are ~22 MB each as raw .npy/.bin (5334x4334,
1 byte/cell) -- just over that per-file limit even at the smallest lossless
dtype (uint8, see src/terrain/raster_io.py's "u1" ENVI option, added this
stage). They compress to under 1 MB apiece (categorical/boolean data, huge
contiguous regions -- gzip -9 ratios of 0.1-2.8%), so they were shipped as
plain gzip (stdlib, nothing project-specific) instead.

Run this ONCE after pulling this commit:
    python scripts/decompress_tappa8_data.py

It finds every *.gz under data/processed/geomorphology/, decompresses it
back to the original filename (lossless -- gzip, not a lossy re-encode),
and deletes the .gz afterward. Safe to re-run (skips a target that already
exists and matches the .gz's uncompressed content).
"""

import gzip
import shutil
import sys
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "processed" / "geomorphology"


def main():
    gz_files = sorted(DATA_DIR.glob("*.gz"))
    if not gz_files:
        print(f"No .gz files found under {DATA_DIR} -- nothing to do.")
        return

    for gz_path in gz_files:
        target = gz_path.with_suffix("")  # strip .gz
        print(f"decompressing {gz_path.name} -> {target.name} ...", end=" ", flush=True)
        with gzip.open(gz_path, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
        gz_path.unlink()
        print(f"done ({target.stat().st_size / 1e6:.2f} MB)")

    print(f"\n{len(gz_files)} file(s) decompressed and .gz originals removed.")


if __name__ == "__main__":
    sys.exit(main())
