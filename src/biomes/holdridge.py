"""
Holdridge Life Zone classification (Holdridge 1967) -- a general-purpose,
world-agnostic implementation. Nothing in this module is specific to this
project's domain; it takes any (12, ny, nx) monthly temperature stack (deg C)
and any (ny, nx) annual precipitation field (mm) and returns the standard
Holdridge fields. The project-specific fine-tuning (which zones this
particular world actually populates, what to call them, the permanent-snow
override) lives in `world_biomes.py`, deliberately kept separate so this
module could be dropped into an unrelated project unchanged.

Why Holdridge over the more commonly tutorial-cited Whittaker diagram: this
project's stated two-axis brief (00_pre_project_planning.md) is temperature +
precipitation. Holdridge honours that -- its third input, potential
evapotranspiration (PET), is *derived* from temperature by a fixed constant,
not an independent measurement -- while being fully formulaic: every boundary
below is a citable number, not a hand-drawn polygon vertex digitized off a
scanned textbook figure (which is what adopting the Whittaker diagram would
require, and which this project already declined to do once for a comparable
reason -- see docs/decisions/04_tappa4_hydrology.md S2 on the rejected
area x slope^2 channel-initiation threshold).

Formulas and constants, and what's actually verified
-----------------------------------------------------
* **Biotemperature** (`biotemperature`): the annual mean of monthly mean
  temperatures, with every month's value clipped to [0, 30] deg C before
  averaging -- Holdridge's own reasoning is that most plants are dormant
  both below freezing and above ~30 deg C, so temperature outside that
  band should not contribute to (or be double-penalized against) the
  biological "growing" signal the index is meant to capture. Confirmed
  directly against two independent secondary sources this session: the
  Wikipedia summary of Holdridge (1967) and Aguirre et al.'s 2025
  *Scientific Data* paper on global climate-classification uncertainty
  ("those months with a mean temperature above 30.0 C and below 0.0 C are
  considered as 30.0 C and 0.0 C, respectively").
* **Potential evapotranspiration**: `PET_mm = 58.93 * biotemperature_c`.
  The 58.93 constant is Holdridge's own empirical calibration (an annual
  PET, in mm, per degree of biotemperature) -- confirmed verbatim in the
  same 2025 *Scientific Data* paper's methods section. **Known weak
  spot, not silently ignored**: this constant was fit mostly against
  tropical/subtropical stations; several secondary sources note the
  resulting PET ratio becomes numerically fragile at low biotemperature
  (below roughly 3 deg C) precisely because PET itself collapses towards
  zero there. This world's land-mean biotemperature sits low (see
  `world_biomes.py`), so this is a live caveat here, not a hypothetical
  one -- restated in the decision doc.
* **Humidity province** (`humidity_province`): named by the ratio
  `PET_mm / annual_precip_mm`, in 8 bands spaced by powers of 2 --
  Superarid (>32), Perarid (16-32), Arid (8-16), Semiarid (4-8), Subhumid
  (2-4), Humid (1-2), Perhumid (0.5-1), Superhumid (<0.5). The log2
  spacing and the 0.125-32 total axis range are corroborated by multiple
  secondary sources this session; Holdridge's original 1967 monograph
  (the ultimate primary source) was not itself accessed.
* **Belt** (`belt`): named by biotemperature alone, in 7 bands with
  boundaries at 1.5 / 3 / 6 / 12 / 18 / 24 deg C (Polar, Subpolar, Boreal,
  Cool Temperate, Warm Temperate, Subtropical, Tropical). The 1.5-3-6-12-24
  sequence (a clean doubling each step) is directly confirmed by an
  EPA/ORNL life-zones report; the same report's own summary table
  collapses Warm Temperate and Subtropical into one row, so the interior
  18 deg C split is reproduced here from secondary literature consensus,
  not independently re-derived.
* **Latitudinal vs. altitudinal naming**: Holdridge's original diagram
  uses different NAMES for the same biotemperature bands depending on
  whether the coldness comes from latitude (Polar/Subpolar/Boreal/...) or
  elevation on a warmer-based mountain (Nival/Alpine/Subalpine/Montane/
  Premontane/...). The numeric boundaries are identical either way -- only
  the label changes, and which label set is "correct" additionally depends
  on the biotemperature of the base region the mountain rises from (a
  cold-temperate-based massif and a tropical-based one use different
  numbers of intermediate belt names for the same physical temperature
  drop). Reproducing that context-dependent renaming correctly requires
  the primary diagram and is out of scope for a reusable module; `belt()`
  below always returns the latitudinal name set, and `world_biomes.py`
  does its own renaming appropriate to this world's specific base
  biotemperature, rather than this module guessing.
* **No official life-zone names** (e.g. "Boreal Wet Forest"): Holdridge's
  diagram assigns roughly 30-38 traditionally-named zones to specific
  belt x humidity-province cells, read directly off the printed diagram.
  Without pixel-accurate access to that diagram this session, inventing
  those compound names risks getting them wrong for cells only weakly
  attested in secondary sources. This module instead returns the belt and
  humidity-province codes as two independently-verified axes; compound
  official-style names are NOT fabricated here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "PET_CONSTANT",
    "BELT_BOUNDARIES_C",
    "BELT_NAMES",
    "HUMIDITY_BOUNDARIES",
    "HUMIDITY_NAMES",
    "HoldridgeResult",
    "biotemperature",
    "potential_evapotranspiration_mm",
    "pet_ratio",
    "belt_index",
    "humidity_index",
    "classify",
]

# Holdridge (1967); PET constant and biotemperature clipping confirmed
# against Aguirre et al. (2025, Scientific Data) -- see module docstring.
PET_CONSTANT = 58.93  # mm PET per year, per degree C of biotemperature

# Ascending biotemperature boundaries (deg C) between the 7 latitudinal
# belts -- Polar | Subpolar | Boreal | Cool Temperate | Warm Temperate |
# Subtropical | Tropical. 1.5/3/6/12/24 confirmed directly (EPA/ORNL life
# zones report); the interior 18 boundary is secondary-literature
# consensus, not independently re-derived this session (see docstring).
BELT_BOUNDARIES_C = (1.5, 3.0, 6.0, 12.0, 18.0, 24.0)
BELT_NAMES = (
    "Polar", "Subpolar", "Boreal", "Cool Temperate",
    "Warm Temperate", "Subtropical", "Tropical",
)

# Ascending PET-ratio boundaries between the 8 humidity provinces --
# Superhumid | Perhumid | Humid | Subhumid | Semiarid | Arid | Perarid |
# Superarid. Powers of 2, per module docstring.
HUMIDITY_BOUNDARIES = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
HUMIDITY_NAMES = (
    "Superhumid", "Perhumid", "Humid", "Subhumid",
    "Semiarid", "Arid", "Perarid", "Superarid",
)


@dataclass
class HoldridgeResult:
    biotemperature_c: np.ndarray
    pet_mm: np.ndarray
    pet_ratio: np.ndarray
    belt_idx: np.ndarray        # 0..6 into BELT_NAMES
    humidity_idx: np.ndarray    # 0..7 into HUMIDITY_NAMES


def biotemperature(temp_c_monthly: np.ndarray) -> np.ndarray:
    """Annual biotemperature from a (12, ny, nx) (or (12,) or (12, n))
    monthly-mean-temperature stack: each month clipped to [0, 30] deg C,
    then averaged. See module docstring for the citation."""
    clipped = np.clip(temp_c_monthly, 0.0, 30.0)
    return clipped.mean(axis=0).astype(np.float32)


def potential_evapotranspiration_mm(biotemp_c: np.ndarray) -> np.ndarray:
    """PET (mm/yr) = 58.93 * biotemperature. See module docstring."""
    return (PET_CONSTANT * biotemp_c).astype(np.float32)


def pet_ratio(biotemp_c: np.ndarray, precip_mm_annual: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """PET / annual precipitation. `eps` guards cells with ~zero precipitation
    (not expected on land here, but this module makes no assumption about
    the caller's domain)."""
    pet = potential_evapotranspiration_mm(biotemp_c)
    return (pet / np.maximum(precip_mm_annual, eps)).astype(np.float32)


def belt_index(biotemp_c: np.ndarray) -> np.ndarray:
    """Index (0..6) into BELT_NAMES via BELT_BOUNDARIES_C (searchsorted:
    boundary values fall into the WARMER of the two adjacent belts, i.e.
    right-open bins [b_i, b_{i+1}) -- consistent with `humidity_index`)."""
    return np.searchsorted(BELT_BOUNDARIES_C, biotemp_c, side="right").astype(np.int8)


def humidity_index(ratio: np.ndarray) -> np.ndarray:
    """Index (0..7) into HUMIDITY_NAMES via HUMIDITY_BOUNDARIES."""
    return np.searchsorted(HUMIDITY_BOUNDARIES, ratio, side="right").astype(np.int8)


def classify(temp_c_monthly: np.ndarray, precip_mm_annual: np.ndarray) -> HoldridgeResult:
    """Full Holdridge classification from a monthly temperature stack and
    an annual precipitation field of matching spatial shape."""
    tb = biotemperature(temp_c_monthly)
    pet = potential_evapotranspiration_mm(tb)
    ratio = (pet / np.maximum(precip_mm_annual, 1e-6)).astype(np.float32)
    return HoldridgeResult(
        biotemperature_c=tb,
        pet_mm=pet,
        pet_ratio=ratio,
        belt_idx=belt_index(tb),
        humidity_idx=humidity_index(ratio),
    )
