"""
Monthly near-surface air temperature (Tappa 2).

Model
-----
Temperature is built as a *sea-level* field with its own seasonal cycle, and
only then reduced to the actual land surface by a lapse rate:

    T(x, y, m) = T_sl(x, y, m) - gamma(m) * z(x, y)

    T_sl(x, y, m) = T_ref
                    + dT_dlat * (lat(y) - lat_ref)
                    + A(x, y) * cos(2*pi * (m - m_peak(x, y)) / 12)

Applying the lapse rate to the whole column (rather than adding a seasonal
term on top of an already-reduced annual mean) is the physically correct
order and is what lets the lapse rate itself vary by month.

Term hierarchy, measured on this domain
---------------------------------------
* elevation        23.4 C sea level -> 3596 m peak at 6.5 C/km  (18.0 C at 5.0)
* latitude          ~1.0 C across the full 1.437 deg N-S span
* continentality    ~1-1.5 C of extra annual *range* at the most inland point

Elevation outruns everything else by more than an order of magnitude. The
latitude and continentality terms are kept because they are cheap and
methodologically expected, not because they drive the map — see the decision
doc, which says so explicitly rather than implying a three-way balance.

Continentality
--------------
Max distance-to-coast on this island is 29.9 km (mean 8.3 km). That is about
a third of South Island NZ's width, so this world is strongly maritime by
construction and continentality can only ever be a second-order term here.
It is modelled the standard way — as a function of distance to the sea — but
it acts on the *annual amplitude* and on the *phase lag*, not on the annual
mean. Continentality does not make a place colder on average; it makes its
summers hotter, its winters colder, and its seasonal peak arrive earlier.

Lapse rate
----------
`config/parameters.yml` originally carried -6.5 C/km, which is the ICAO
standard-atmosphere / global-average environmental lapse rate. It is the
wrong value for this world's stated validation reference: Norton (1985)
derived 5.0 C/km as the mean annual lapse rate from 301 New Zealand climate
stations, and it has been confirmed since at 5.0 (Franz Josef, Anderson et
al. 2006) and 5.4 with strong winter inversions (Ben Ohau, Doughty 2013).
Maritime mid-latitude ranges are shallower than the global average because
the air crossing them is close to moist-adiabatic. The default here is
therefore 5.0 C/km, with a seasonal modulation reproducing the observed
pattern (steeper in summer, shallower in winter when inversions set in).
At the 3596 m peak the two choices differ by 5.4 C, which is the difference
between a permanent snowline and no permanent snowline — not a detail.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

__all__ = ["TemperatureParams", "distance_to_coast_km", "temperature_month", "temperature_year"]


@dataclass
class TemperatureParams:
    # --- sea-level annual mean ---
    t_ref_sea_level_c: float = 11.5
    """Annual mean sea-level temperature at `lat_ref_deg`. NZ coastal stations
    near this latitude sit around 11-12 C."""

    lat_ref_deg: float = -44.0
    """Latitude the reference temperature is quoted at (domain centre)."""

    dt_dlat_c_per_deg: float = 0.70
    """Meridional gradient, C per degree of latitude *towards the equator*.
    Positive: moving north (towards -43.28) gets warmer."""

    # --- lapse rate ---
    lapse_rate_c_per_km: float = 5.0
    """Annual mean. See module docstring for why this is 5.0 and not 6.5."""

    lapse_seasonal_amplitude_c_per_km: float = 0.3
    """Half-range of the seasonal lapse-rate cycle.

    Deliberately conservative, and worth understanding before changing. A
    single column lapse rate cannot represent both of the things the
    observations show: winter *valley* inversions (a boundary-layer effect
    near the ground) and free-atmosphere summit temperatures (which are about
    as seasonal as the surface). Any positive value here trades the second
    for the first, because it damps the annual range aloft by
    2 * amplitude * z: at 0.7 C/km the 3.6 km summit would end up with only
    3.0 C of annual range against the coast's 8.0, which is too flat to
    defend. At 0.3 the summit keeps ~5.8 C.

    Set to 0.0 for a strictly constant lapse rate. Tappa 3's snow metrics are
    the place this choice actually bites — see the sensitivity note in the
    decision doc."""

    lapse_peak_month: float = 2.0
    """Month of the steepest (most negative) lapse rate. February = SH late
    summer, when surface heating is strongest and inversions are absent."""

    # --- seasonal cycle ---
    amplitude_coastal_c: float = 4.0
    """Half the annual temperature range right at the shore."""

    amplitude_inland_c: float = 7.0
    """Asymptotic half-range far inland. Never fully reached here: at the most
    inland cell (29.9 km) the decay function only gets to ~0.70."""

    continentality_scale_km: float = 30.0
    """e-folding distance of the continentality decay, 1 - exp(-d/L).

    Kept at a physically defensible maritime-to-continental transition scale
    rather than shrunk to make the term look important on a small island.
    The consequence, measured: at low elevation the annual range goes from
    8.05 C within 2 km of the shore to 10.74 C at 16-30 km inland. Real, but
    a fifth of what elevation does over the same distance.

    Note this effect is nearly invisible if you bin the surface field by
    distance-to-coast alone, because on this island "inland" and "high" are
    the same cells, and the seasonal lapse term damps the annual range with
    height by almost exactly as much as continentality raises it. The two
    have to be separated by holding elevation fixed."""

    peak_month_coastal: float = 2.2
    """Warmest month at the coast. Maritime thermal inertia puts the peak two
    months after the solstice."""

    peak_month_shift_inland: float = 0.9
    """How much earlier (in months) the peak arrives far inland."""

    # --- validation bookkeeping ---
    notes: dict = field(default_factory=dict)


def distance_to_coast_km(land: np.ndarray, res_m: float) -> np.ndarray:
    """Euclidean distance from every land cell to the nearest sea cell, in km.

    Replaces GRASS `r.grow.distance` from the original Tappa 2 sketch — no
    GRASS bindings are available in this environment, and
    `scipy.ndimage.distance_transform_edt` computes the identical quantity.
    Sea cells get 0.
    """
    return (ndimage.distance_transform_edt(land, sampling=res_m) / 1000.0).astype(np.float32)


def _continentality(dist_km: np.ndarray, scale_km: float) -> np.ndarray:
    return 1.0 - np.exp(-dist_km / scale_km)


def lapse_rate_c_per_km(month: float, p: TemperatureParams) -> float:
    """Month-dependent lapse rate (positive number, C per km of ascent)."""
    return p.lapse_rate_c_per_km + p.lapse_seasonal_amplitude_c_per_km * np.cos(
        2.0 * np.pi * (month - p.lapse_peak_month) / 12.0
    )


def sea_level_temperature(
    month: float, lat_deg: np.ndarray, dist_km: np.ndarray, p: TemperatureParams
) -> np.ndarray:
    c = _continentality(dist_km, p.continentality_scale_km)
    annual_mean = p.t_ref_sea_level_c + p.dt_dlat_c_per_deg * (lat_deg - p.lat_ref_deg)
    amplitude = p.amplitude_coastal_c + (p.amplitude_inland_c - p.amplitude_coastal_c) * c
    peak = p.peak_month_coastal - p.peak_month_shift_inland * c
    return annual_mean + amplitude * np.cos(2.0 * np.pi * (month - peak) / 12.0)


def temperature_month(
    month: float,
    elevation_m: np.ndarray,
    lat_deg: np.ndarray,
    dist_km: np.ndarray,
    p: TemperatureParams | None = None,
) -> np.ndarray:
    """Mean near-surface temperature (C) for one month.

    `elevation_m` is the land surface; below sea level it is clipped to 0 so
    that ocean cells report a sea-surface temperature rather than a spurious
    warm anomaly from a negative lapse-rate term.
    """
    p = p or TemperatureParams()
    z_km = np.maximum(elevation_m, 0.0) / 1000.0
    return (sea_level_temperature(month, lat_deg, dist_km, p) - lapse_rate_c_per_km(month, p) * z_km).astype(
        np.float32
    )


def temperature_year(
    elevation_m: np.ndarray,
    lat_deg: np.ndarray,
    dist_km: np.ndarray,
    p: TemperatureParams | None = None,
) -> np.ndarray:
    """(12, ny, nx) stack of monthly mean temperatures, month index 0 = January."""
    return np.stack(
        [temperature_month(m, elevation_m, lat_deg, dist_km, p) for m in range(1, 13)]
    )
