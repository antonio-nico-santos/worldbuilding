"""
Tappa 6 -- solar exposure / insolation proxy. The one criterion with no
GRASS `r.sun` (or any GIS solar-radiation tool) available in this sandbox, so
it is reimplemented from published, citable formulas rather than a from-
scratch radiative-transfer model -- same posture as noise.py/erosion.py/
flow.py in the prior stages.

Method, in order:

1. Solar geometry (Duffie & Beckman, "Solar Engineering of Thermal
   Processes"): declination, per-cell solar elevation/azimuth from
   latitude + hour angle. Verified against hand-calculated equinox/solstice
   values at this world's -44 deg latitude before use (see
   `docs/decisions/` Tappa 6 notes) -- e.g. equinox noon gives elevation =
   90-|lat| exactly, June (day 172) gives a LOW noon sun (southern winter),
   December (day 355) a HIGH one (southern summer): the formula is
   hemisphere-agnostic given a signed latitude, no special-casing needed.

2. Clear-sky horizontal radiation via FAO-56 (Allen et al. 1998, "Crop
   evapotranspiration", Irrigation and Drainage Paper 56): Ra (extra-
   terrestrial daily radiation) from the standard closed-form daily
   integral, then Rso = (0.75 + 2e-5*z) * Ra -- elevation-corrected clear
   sky radiation, thinner atmosphere at altitude passing more direct beam.
   This is the ONLY citable/analytic piece; it supplies the atmospheric
   attenuation magnitude but says nothing about terrain geometry.

3. Terrain geometry is handled separately and numerically, because an
   analytic all-day formula can't incorporate per-hour horizon shading.
   For each representative day, step through hour angles and compute, at
   each step, the incidence angle on the actual tilted/oriented surface
   (Duffie & Beckman's cos(theta) formula) and whether the sun is blocked
   by surrounding terrain (see `horizon_angles`). The ratio of the
   (shaded, tilted) numerical sum to the (unshaded, horizontal) numerical
   sum over the same hour steps is then applied to Rso -- this keeps the
   atmospheric-attenuation number analytic/citable while letting the
   numerical integration carry only the geometric redistribution.

4. Diffuse sky radiation: a fixed 15% clear-sky diffuse fraction (a
   documented placeholder, not independently sourced -- same honesty
   posture as Tappa 3's sigma_day_c) scaled by an approximate sky-view
   factor: (1+cos(slope))/2 (Liu & Jordan 1960 isotropic-sky model) further
   discounted by the mean horizon angle across all 16 directions, as a
   crude stand-in for surrounding terrain blocking part of the sky dome
   (not a rigorous sky-view integral).

HONEST LIMITATIONS (stated up front, not discovered later):
  - No cloud cover / real atmospheric turbidity -- this is a CLEAR-SKY
    proxy, useful for ranking sites against each other, not a bankable
    energy-yield estimate.
  - Horizon shading uses nearest-integer-pixel ray marching in 16
    directions out to 10 km, at the elevation grid's own 120 m resolution
    -- it will miss shading from terrain features narrower than the
    120 m cell (irrelevant at settlement-siting scale) and slightly
    mis-times shadow onset/offset versus true azimuth (bounded by
    360/16 = 22.5 deg direction spacing).
  - No equation-of-time or per-cell solar-noon longitude correction --
    defensible given the whole domain is ~1.2 deg of longitude wide.
  - The 0.5 hr integration step and 15% diffuse fraction are both
    resolution/parameter choices, not derived constants.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "fao56_declination_rad",
    "eccentricity_correction",
    "extraterrestrial_radiation_MJ_m2_day",
    "solar_elevation_azimuth_deg",
    "horizon_angles_deg",
    "monthly_insolation_MJ_m2_day",
]

_GSC = 0.0820  # MJ m-2 min-1, FAO-56's solar constant
_DIFFUSE_FRACTION = 0.15
_MID_MONTH_DOY = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
_DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def fao56_declination_rad(day_of_year):
    return 0.409 * np.sin(2 * np.pi / 365.0 * day_of_year - 1.39)


def eccentricity_correction(day_of_year):
    return 1 + 0.033 * np.cos(2 * np.pi / 365.0 * day_of_year)


def extraterrestrial_radiation_MJ_m2_day(lat_rad, dec_rad, day_of_year):
    """FAO-56 eq. 21 -- daily extraterrestrial radiation on a horizontal
    surface, plus the sunset hour angle `ws` (both needed downstream)."""
    ws = np.arccos(np.clip(-np.tan(lat_rad) * np.tan(dec_rad), -1, 1))
    dr = eccentricity_correction(day_of_year)
    ra = (
        (24 * 60 / np.pi)
        * _GSC
        * dr
        * (ws * np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(dec_rad) * np.sin(ws))
    )
    return ra, ws


def solar_elevation_azimuth_deg(lat_deg, dec_deg, hour_angle_deg):
    """Duffie & Beckman solar position. Azimuth is measured clockwise from
    North (0-360), NOT the northern-hemisphere "from South" convention --
    this keeps the incidence-angle formula hemisphere-agnostic. Verified
    at lat=-44: equinox noon -> elev=46 az=0 (due N); day 172 (S winter)
    noon -> elev~22.5 (low); day 355 (S summer) noon -> elev~69.5 (high);
    morning/afternoon azimuths mirror exactly around North."""
    lat, dec, H = np.radians(lat_deg), np.radians(dec_deg), np.radians(hour_angle_deg)
    sin_elev = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(H)
    elev = np.arcsin(np.clip(sin_elev, -1, 1))
    cos_az = (np.sin(dec) - np.sin(lat) * np.sin(elev)) / (np.cos(lat) * np.cos(elev) + 1e-12)
    az = np.arccos(np.clip(cos_az, -1, 1))
    az = np.where(H > 0, 2 * np.pi - az, az)
    return np.degrees(elev), np.degrees(az)


def horizon_angles_deg(elevation_m, cellsize_m, n_dirs=16, max_dist_m=10000):
    """Per-cell horizon angle (degrees) in `n_dirs` azimuth bins (0=N,
    clockwise), by marching outward at nearest-integer-pixel steps to
    `max_dist_m` and keeping the steepest elevation angle seen. Verified
    against a synthetic flat-plain-plus-wall case (known geometry) before
    use on the real DEM -- see Tappa 6 planning notes.

    Returns (horizon_deg[n_dirs, ny, nx], azimuths_deg[n_dirs]).
    """
    ny, nx = elevation_m.shape
    max_steps = int(max_dist_m / cellsize_m)
    azimuths = np.linspace(0, 360, n_dirs, endpoint=False)
    horizon = np.zeros((n_dirs, ny, nx), dtype=np.float32)
    yy, xx = np.mgrid[0:ny, 0:nx]
    for di, az in enumerate(azimuths):
        rad = np.radians(az)
        dy_unit, dx_unit = -np.cos(rad), np.sin(rad)  # row0=North convention
        running_max_tan = np.full((ny, nx), -10.0, dtype=np.float32)
        for step in range(1, max_steps + 1):
            dist_m = step * cellsize_m
            oy = np.round(yy + dy_unit * step).astype(np.int64)
            ox = np.round(xx + dx_unit * step).astype(np.int64)
            valid = (oy >= 0) & (oy < ny) & (ox >= 0) & (ox < nx)
            oy_c, ox_c = np.clip(oy, 0, ny - 1), np.clip(ox, 0, nx - 1)
            tan_ang = (elevation_m[oy_c, ox_c] - elevation_m) / dist_m
            running_max_tan = np.maximum(running_max_tan, np.where(valid, tan_ang, -10.0))
        horizon[di] = np.degrees(np.arctan(running_max_tan))
    return horizon, azimuths


def _horizon_at_azimuth(horizon, horizon_az_deg, sun_az_deg):
    """Circular linear interpolation of the per-direction horizon array at
    each cell's own sun azimuth (sun_az_deg has the grid's shape)."""
    n = len(horizon_az_deg)
    step = 360.0 / n
    pos = (sun_az_deg % 360) / step
    i0 = np.floor(pos).astype(int) % n
    i1 = (i0 + 1) % n
    frac = pos - np.floor(pos)
    h0 = np.take_along_axis(horizon, i0[None, :, :], axis=0)[0]
    h1 = np.take_along_axis(horizon, i1[None, :, :], axis=0)[0]
    return h0 * (1 - frac) + h1 * frac


def monthly_insolation_MJ_m2_day(
    lat_deg,
    elevation_m,
    slope_deg,
    aspect_deg,
    horizon_deg,
    horizon_az_deg,
    hour_step=0.5,
    diffuse_fraction=_DIFFUSE_FRACTION,
):
    """Returns (12, ny, nx) float32 -- representative-mid-month daily total
    insolation (direct + diffuse, MJ/m2/day) under this clear-sky proxy,
    per the module's documented method."""
    lat_rad = np.radians(lat_deg)
    slope_rad = np.radians(slope_deg)
    ny, nx = elevation_m.shape

    svf_slope = (1 + np.cos(slope_rad)) / 2.0
    mean_horizon_deg = horizon_deg.mean(axis=0)
    svf_terrain_extra = np.clip(1 - mean_horizon_deg / 90.0, 0, 1)
    svf = svf_slope * svf_terrain_extra

    monthly = np.zeros((12, ny, nx), dtype=np.float32)
    for mi, day_of_year in enumerate(_MID_MONTH_DOY):
        dec = fao56_declination_rad(day_of_year)
        ra, ws = extraterrestrial_radiation_MJ_m2_day(lat_rad, dec, day_of_year)
        rso = (0.75 + 2e-5 * np.maximum(elevation_m, 0)) * ra

        max_ws_hours = (np.degrees(ws) / 15.0).max()
        hour_offsets = np.arange(-max_ws_hours, max_ws_hours + 1e-9, hour_step)

        num_tilted = np.zeros((ny, nx))
        den_horizontal = np.zeros((ny, nx))
        for dh in hour_offsets:
            hour_angle_deg = 15.0 * dh
            active = np.abs(hour_angle_deg) <= np.degrees(ws)
            elev_deg, sun_az_deg = solar_elevation_azimuth_deg(np.degrees(lat_rad), np.degrees(dec), hour_angle_deg)
            elev_rad = np.radians(elev_deg)
            zen_rad = np.pi / 2 - elev_rad
            cos_zen = np.maximum(np.cos(zen_rad), 0.0)

            daz = np.radians(sun_az_deg - aspect_deg)
            cos_inc = np.maximum(
                np.cos(slope_rad) * np.cos(zen_rad) + np.sin(slope_rad) * np.sin(zen_rad) * np.cos(daz), 0.0
            )
            horizon_here = _horizon_at_azimuth(horizon_deg, horizon_az_deg, sun_az_deg)
            unshaded = elev_deg > horizon_here
            contrib = active & unshaded & (elev_deg > 0)

            num_tilted += np.where(contrib, cos_inc, 0.0)
            den_horizontal += np.where(active & (elev_deg > 0), cos_zen, 0.0)

        ratio = np.where(den_horizontal > 1e-6, num_tilted / np.maximum(den_horizontal, 1e-6), 0.0)
        monthly[mi] = (rso * (1 - diffuse_fraction) * ratio + rso * diffuse_fraction * svf).astype(np.float32)
    return monthly


def annual_insolation_MJ_m2(monthly_daily_MJ_m2_day):
    return sum(monthly_daily_MJ_m2_day[i] * _DAYS_IN_MONTH[i] for i in range(12))


def normalize_suitability(arr, land_mask, pctl=(2, 98)):
    p_lo, p_hi = np.percentile(arr[land_mask], pctl)
    suit = np.clip((arr - p_lo) / (p_hi - p_lo), 0.0, 1.0)
    return np.where(land_mask, suit, np.nan).astype(np.float32), (float(p_lo), float(p_hi))
