"""
Tappa 8 -- basin fill internal zonation (Arable / Wetland-backswamp /
Estuarine-coastal), implementing as an actual raster the LOGIC already
authored in scenario_reference.md S22.4. That section explicitly left
"exact area/threshold percentiles for all three sub-uses... for Tappa 8 to
calibrate against the actual rasters" -- this module is that calibration,
not a new design. The logic itself (three sub-uses, which fields drive
them, coastal claims priority over the relief split) is NOT this module's
invention; only the specific threshold values are.

Two fields, per S22.4, both already-existing Tappa 8 machinery -- no new
spatial primitive, matching S22.4's own framing ("reusing spatial fields
Tappa 8's own pipeline already computes"):
- relief_2km (terrain_relief.compute_local_relief, the SAME field
  lithology v5 already uses for schist/greywacke/basin_fill) -- a
  drainage/local-relief proxy. Higher relief = better-drained = Arable;
  lower relief = poorly-drained = Wetland/backswamp.
- dist_to_ocean_km (caves.distance_to_ocean_km, the SAME field already
  used for sea caves and placer magnetite's coastal weighting) -- coastal
  proximity. Below the coastal threshold, regardless of relief, =
  Estuarine/coastal margin.

Classification order matters and mirrors S22.4's own prose exactly: the
coastal check runs FIRST and claims cells unconditionally, because S22.4
defines wetland as "poorly-drained, low local relief, INLAND (not
coastal-adjacent)" -- that definition presupposes the coastal band has
already been carved out. The relief split (Arable vs. Wetland) only
applies to whatever basin_fill remains after the coastal claim.

Derived from a DEM/distance field, not hand-authored as geometry -- the
same "DEM-native beats line/polygon-distance" reasoning S0-S2 of the
decision doc already established for lithology v5 itself (a hand-drawn
zone can't track real, locally-irregular drainage/coastal geometry the
way a computed field does), applied here to the sub-zonation question
S22.19 explicitly left open ("whether Tappa 8 actually authors it as real
geometry").
"""
import numpy as np

BASIN_SUBZONE_CLASSES = {
    0: "not_basin_fill",
    1: "arable",
    2: "wetland_backswamp",
    3: "estuarine_coastal",
}
CLASS_ARABLE = 1
CLASS_WETLAND_BACKSWAMP = 2
CLASS_ESTUARINE_COASTAL = 3


def classify_basin_fill_zones(
    basin_fill_mask: np.ndarray,
    relief_2km: np.ndarray,
    dist_to_ocean_km: np.ndarray,
    coastal_threshold_km: float = 1.5,
    wetland_relief_percentile: float = 25.0,
):
    """Returns `(subzone, stats)`. `subzone` is a uint8 array of
    BASIN_SUBZONE_CLASSES codes, 0 (not_basin_fill) everywhere
    `basin_fill_mask` is False.

    `coastal_threshold_km`: an ABSOLUTE distance, not a percentile --
    unlike relief (no natural absolute scale), "how far does
    estuarine/tidal influence plausibly reach" is a physically meaningful
    km quantity. Closer in spirit to karst/talus's `stream_buffer_km` and
    sea_cave's `coastal_buffer_km` (both absolute, both 0.5 km) than to a
    percentile threshold -- scaled up from those 0.5 km point-feature
    eligibility tests because this defines a broad LANDSCAPE ZONE, not a
    specific cave-formation trigger. UNREVIEWED first-pass value.

    `wetland_relief_percentile`: computed over the NON-COASTAL basin_fill
    population's OWN relief_2km distribution -- not mainland's. basin_fill
    was already excluded from schist/greywacke's higher-relief thresholds
    by `classify_from_terrain`, so re-percentiling within basin_fill
    itself (rather than reusing a mainland-wide percentile) is the
    locally-consistent choice, same spirit as every other percentile
    threshold in this project being calibrated against the population it
    actually applies to. 25.0 (bottom quartile of relief = wetland) is
    sized so Arable stays basin_fill's "clear majority", per S22.4's own
    prose -- a sanity-checked choice, not an independently derived one.

    `stats` records both the percentile used and the raw relief value (m)
    it resolved to against this world's actual data, for the decision-doc
    record -- same "record the calibration, not just the knob" convention
    used throughout this project.
    """
    subzone = np.zeros(basin_fill_mask.shape, dtype=np.uint8)

    coastal = basin_fill_mask & (dist_to_ocean_km <= coastal_threshold_km)
    subzone[coastal] = CLASS_ESTUARINE_COASTAL

    remaining = basin_fill_mask & ~coastal
    relief_remaining = relief_2km[remaining]
    wetland_relief_threshold_m = float(np.percentile(relief_remaining, wetland_relief_percentile))
    wetland = remaining & (relief_2km <= wetland_relief_threshold_m)
    subzone[wetland] = CLASS_WETLAND_BACKSWAMP

    arable = remaining & ~wetland
    subzone[arable] = CLASS_ARABLE

    stats = {
        "coastal_threshold_km": coastal_threshold_km,
        "wetland_relief_percentile": wetland_relief_percentile,
        "wetland_relief_threshold_m": wetland_relief_threshold_m,
        "n_basin_fill_cells": int(basin_fill_mask.sum()),
        "n_coastal_cells": int(coastal.sum()),
        "n_wetland_cells": int(wetland.sum()),
        "n_arable_cells": int(arable.sum()),
    }
    return subzone, stats
