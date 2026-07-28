"""
Orographic precipitation via Smith & Barstad (2004) linear theory (Tappa 2).

Reference
---------
Smith, R.B. & Barstad, I. (2004), "A linear theory of orographic
precipitation", *Journal of the Atmospheric Sciences* 61(12), 1377-1391.

Why this and not an upslope heuristic
-------------------------------------
The usual worldbuilding recipe is `P ~ dot(wind, grad(h))` plus an ad-hoc
leeward decay. It produces a picture, but it has three defects that matter
for a piece meant to be defended technically:

1. Precipitation lands exactly on the slope that produced the uplift. Real
   condensate is advected downwind while it forms and while it falls, so the
   maximum sits *upwind of* — or beyond — the crest depending on wind speed,
   never pinned to the steepest slope.
2. The rain shadow has to be bolted on as a separate rule, with its own
   free parameters, rather than emerging from the same physics.
3. There is no way to justify any particular constant.

Smith & Barstad solves the steady-state advection of vertically integrated
cloud water and hydrometeors over terrain, linearised about a uniform wind,
in Fourier space:

    P_hat(k,l) = Cw * i * sigma * h_hat(k,l)
                 / [ (1 - i*m*Hw) * (1 + i*sigma*tau_c) * (1 + i*sigma*tau_f) ]

    sigma = u*k + v*l                      (intrinsic frequency)
    m     = sqrt( (Nm^2 - sigma^2)/sigma^2 * (k^2 + l^2) )   (vertical wavenumber)

Windward enhancement, the downwind displacement of the maximum, and leeward
drying all fall out of that one expression. Every constant has a name, a
unit, and a measurable real-world range, which is the whole point.

Coupling to the temperature model
---------------------------------
`Cw` (uplift sensitivity) and `Hw` (moisture scale height) are not tuned
independently here — both are derived from the month's sea-level temperature
through Clausius-Clapeyron and the moist adiabat (`thermo_from_temperature`).
Warm months hold more moisture and rain harder for the same wind; cold months
less. So the seasonal precipitation cycle is a consequence of the seasonal
temperature cycle plus the seasonal wind cycle, not a third hand-tuned knob.

Known limits (stated here so the decision doc does not have to claim more
than the model delivers)
------------------------
* Linear: valid while the flow goes *over* the barrier rather than around or
  blocked. For a 3.4 km ridge at U ~ 12 m/s and Nm = 0.005 the non-dimensional
  mountain height Nm*H/U is ~1.4, i.e. at the edge of where blocking starts to
  matter. The model will therefore overstate crest precipitation somewhat and
  understate upwind-flank precipitation.
* No moisture budget: an air parcel is never depleted by an upwind range, so a
  second range downwind is not starved the way it should be. Over a 130 km
  island this is minor; it would not be over a continent.
* Steady state, one wind and one stability per month, so this is a monthly
  *mean* field, not a synoptic sequence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:  # scipy's pocketfft is ~2x numpy's here, but numpy is a fine fallback
    from scipy import fft as _fft
except ImportError:  # pragma: no cover
    from numpy import fft as _fft

__all__ = [
    "PrecipParams",
    "orographic_precip",
    "precip_month",
    "wind_speed",
    "saturation_vapour_density",
    "moist_adiabatic_lapse_rate",
    "thermo_from_temperature",
    "wind_vector",
]

_RV = 461.5      # J kg-1 K-1, gas constant for water vapour
_RD = 287.05     # J kg-1 K-1, gas constant for dry air
_LV = 2.5e6      # J kg-1, latent heat of vaporisation
_CP = 1005.0     # J kg-1 K-1
_G = 9.80665     # m s-2


# --------------------------------------------------------------------------
# thermodynamics
# --------------------------------------------------------------------------

def saturation_vapour_pressure(t_c: float | np.ndarray) -> np.ndarray:
    """Saturation vapour pressure (Pa) over liquid water, Bolton (1980)."""
    t_c = np.asarray(t_c, dtype=float)
    return 611.2 * np.exp(17.67 * t_c / (t_c + 243.5))


def saturation_vapour_density(t_c: float | np.ndarray) -> np.ndarray:
    """Saturation water-vapour density (kg m-3) — the `rho_Sref` of S&B."""
    t_k = np.asarray(t_c, dtype=float) + 273.15
    return saturation_vapour_pressure(t_c) / (_RV * t_k)


def moist_adiabatic_lapse_rate(t_c: float | np.ndarray, p_pa: float = 100_000.0) -> np.ndarray:
    """Saturated adiabatic lapse rate (K m-1) at temperature `t_c`, pressure `p_pa`."""
    t_k = np.asarray(t_c, dtype=float) + 273.15
    es = saturation_vapour_pressure(t_c)
    r = 0.622 * es / (p_pa - es)              # saturation mixing ratio
    num = 1.0 + _LV * r / (_RD * t_k)
    den = 1.0 + _LV**2 * r / (_RV * _CP * t_k**2)
    return _G / _CP * num / den


def thermo_from_temperature(t_sea_level_c: float, env_lapse_c_per_km: float) -> tuple[float, float]:
    """Return ``(Cw, Hw)`` for a given sea-level temperature.

    * ``Cw = rho_Sref * Gamma_m / gamma``  [kg m-3] — how much condensate a
      unit of uplift produces. S&B eq. (12).
    * ``Hw = Rv * T^2 / (L * gamma)``      [m] — the depth over which moisture
      is distributed, from Clausius-Clapeyron with a constant lapse rate.

    Note both scale *inversely* with the environmental lapse rate: choosing
    5.0 C/km rather than the standard-atmosphere 6.5 raises Hw by ~30%, which
    is why the lapse-rate decision in `temperature.py` is not confined to the
    temperature model.
    """
    gamma = env_lapse_c_per_km / 1000.0                       # K m-1
    gamma_m = float(moist_adiabatic_lapse_rate(t_sea_level_c))
    cw = float(saturation_vapour_density(t_sea_level_c)) * gamma_m / gamma
    t_k = t_sea_level_c + 273.15
    hw = _RV * t_k**2 / (_LV * gamma)
    return cw, hw


def wind_vector(from_deg: float, speed: float) -> tuple[float, float]:
    """Convert a meteorological "wind from" bearing to a (u, v) map vector.

    `from_deg` is clockwise from north, meteorological convention: 250 = a
    west-south-westerly, i.e. air arriving from 250 deg and travelling
    towards 70 deg (east-north-east). Returns eastward and northward
    components in the projected CRS.
    """
    toward = np.radians((from_deg + 180.0) % 360.0)
    return float(speed * np.sin(toward)), float(speed * np.cos(toward))


# --------------------------------------------------------------------------
# the model
# --------------------------------------------------------------------------

@dataclass
class PrecipParams:
    wind_from_deg: float = 250.0
    """Prevailing wind bearing. Locked at 250 (WSW) in Tappa 2: the domain
    sits in the 44 deg S westerly belt, and within the western quadrant 250
    maximises orographic effectiveness against this particular skeleton
    (0.77 overall, no ridge below 0.67) — see the decision doc."""

    wind_speed_ms: float = 12.0
    """Annual mean cross-barrier wind speed."""

    wind_speed_seasonal_amplitude_ms: float = 2.0
    """Half-range of the seasonal wind cycle. The SH mid-latitude jet shifts
    equatorward in winter, so 44 S — on the equatorward flank of the belt —
    sees its strongest westerlies then.

    This single number sets the whole seasonal precipitation cycle, because
    the two seasonal drivers pull against each other: winter has more wind
    but colder, drier air. Measured domain-mean winter:summer ratios are
    0.76 at amplitude 0 (moisture wins outright), 1.02 at 1.5, 1.12 at 2.0
    and 1.36 at 3.0. 2.0 is chosen to land near the reference region's
    behaviour — West Coast rainfall is close to uniform through the year,
    Canterbury has a mild winter/autumn maximum — without either the flat
    cycle or a pronounced monsoon-like one."""

    wind_peak_month: float = 8.0
    """Month of strongest westerlies (August, SH late winter)."""

    nm_moist_stability: float = 0.005
    """Moist Brunt-Vaisala frequency, s-1. S&B's tested range is 0 (neutral,
    convective) to 0.01 (stably stratified). 0.005 is their mid-range value
    for saturated mid-latitude flow."""

    tau_c_s: float = 1000.0
    """Cloud-water -> hydrometeor conversion time, s."""

    tau_f_s: float = 1000.0
    """Hydrometeor fallout time, s. Together tau_c + tau_f set how far
    downwind of the uplift the rain actually lands: U*(tau_c+tau_f) = 24 km
    at 12 m/s, comparable to the width of the Spine itself."""

    orographic_duty_cycle: float = 0.07
    """Fraction of the month during which the modelled saturated upslope flow
    is actually happening. THIS IS THE SINGLE CALIBRATED PARAMETER of the
    precipitation model, and it exists because Smith & Barstad is an *event*
    model, not a climatology.

    Run at face value the model returns the instantaneous rate of a fully
    saturated 12 m/s airstream lifted over the Spine: ~43 mm/hr on the
    windward flank, which is a perfectly reasonable heavy-rain rate and an
    absurd monthly mean (~170,000 mm/yr). Real climatology is that rate
    times the fraction of hours it is occurring. Because the model is linear
    in Cw, applying the duty cycle is identical to using an effective
    Cw = Cw_thermo * duty, which is how it is implemented.

    0.07 is calibrated (not assumed) against South Island NZ: it puts the
    crest maximum at ~12,600 mm/yr against the Southern Alps' ~11,000-13,000,
    and the windward coastal strip at ~5,000-7,000 against Milford Sound's
    6,545. See the decision doc's validation table."""

    background_precip_mm_per_month: float = 60.0
    """Non-orographic (synoptic / frontal) precipitation — the rain a flat
    island at this latitude would get anyway. 60 mm/month = 720 mm/yr, in
    line with flat maritime NZ sites (Christchurch 612, Invercargill 1035).
    Linear theory makes no attempt at this component; it is additive and
    stated separately rather than folded into a fitted constant."""

    lee_floor_fraction: float = 0.45
    """Deepest achievable lee suppression, as a fraction of the background.

    Left to itself the linear solution drives the lee to -230,000 mm/yr
    equivalent, and the usual `max(P, 0)` truncation then sets 27% of the
    island's land area to *exactly zero* rainfall — which no maritime
    mid-latitude place on Earth records. The unbounded drying is an artefact
    of the model having no moisture budget: it can remove more water than the
    airstream ever carried. Capping suppression at 45% of background puts the
    driest cells at ~324 mm/yr, against Alexandra (Central Otago, the deepest
    real rain shadow in the reference region) at ~300 mm/yr. The cap is
    applied as a softplus rather than a hard `max` so it does not stamp a
    flat-valued plateau across the whole lee."""

    pad_km: float = 40.0
    """Zero-taper margin added around the domain before the FFT. Without it
    the transform's implicit periodicity would let the eastern lee wrap
    round and rain on the western windward coast."""


def _taper_pad(h: np.ndarray, pad: int) -> np.ndarray:
    """Pad with zeros and blend the terrain down to 0 with a raised cosine.

    A hard zero edge is itself a step function and would inject broadband
    spectral energy — the same failure mode as the single-frequency domain
    warp in Tappa 1 (decision doc 01, S3d). The blend runs over the padded
    margin, so no real terrain is modified.
    """
    ny, nx = h.shape
    out = np.zeros((ny + 2 * pad, nx + 2 * pad), dtype=float)
    out[pad : pad + ny, pad : pad + nx] = h
    ramp = 0.5 * (1.0 - np.cos(np.pi * (np.arange(pad) + 0.5) / pad))
    wy = np.concatenate([ramp, np.ones(ny), ramp[::-1]])
    wx = np.concatenate([ramp, np.ones(nx), ramp[::-1]])
    return out * wy[:, None] * wx[None, :]


def orographic_precip(
    h: np.ndarray,
    res_m: float,
    *,
    u: float,
    v: float,
    cw: float,
    hw: float,
    nm: float,
    tau_c: float,
    tau_f: float,
    seconds: float,
    pad_km: float = 40.0,
    return_anomaly: bool = False,
) -> np.ndarray:
    """Precipitation (mm over `seconds`) for one steady wind/moisture state.

    `h` is the uplift-driving surface: ``max(dem, 0)`` area-averaged to the
    working grid, north-up (row 0 = north). Bathymetry must already be
    clipped — air flows over the sea surface, not the sea floor.
    """
    pad = int(round(pad_km * 1000.0 / res_m))
    # work south-up so that increasing row index = increasing northing, which
    # makes the l-wavenumber sign convention match the map's v component.
    hp = _taper_pad(h, pad)[::-1, :]
    ny, nx = hp.shape

    k = 2.0 * np.pi * _fft.fftfreq(nx, d=res_m)[None, :]   # eastward
    l = 2.0 * np.pi * _fft.fftfreq(ny, d=res_m)[:, None]   # northward
    sigma = u * k + v * l

    with np.errstate(divide="ignore", invalid="ignore"):
        m_sq = (nm**2 - sigma**2) / sigma**2 * (k**2 + l**2)
        m = np.sqrt(m_sq.astype(complex))
        # Radiation condition: for vertically propagating waves (sigma^2 < Nm^2)
        # the vertical wavenumber must carry energy upward, which fixes
        # sign(m) = sign(sigma). Getting this branch wrong flips the whole
        # field's upwind/downwind asymmetry — the rain lands on the lee side.
        prop = (m_sq > 0) & np.isfinite(m_sq)
        m = np.where(prop & (sigma < 0), -m, m)

        h_hat = _fft.fft2(hp)
        denom = (1.0 - 1j * m * hw) * (1.0 + 1j * sigma * tau_c) * (1.0 + 1j * sigma * tau_f)
        p_hat = cw * 1j * sigma * h_hat / denom

    p_hat[~np.isfinite(p_hat)] = 0.0        # sigma == 0: flow parallel to the
                                            # feature does no work, no precip
    p = np.real(_fft.ifft2(p_hat))[::-1, :]  # back to north-up
    p = p[pad : pad + h.shape[0], pad : pad + h.shape[1]]

    return (p * seconds).astype(np.float32)   # kg m-2 s-1 == mm s-1


def _softplus_floor(x: np.ndarray, floor: float, knee: float) -> np.ndarray:
    """Smooth asymptote to `floor` from above; identity well above it."""
    d = (x - floor) / knee
    return floor + knee * np.where(d > 30.0, d, np.log1p(np.exp(np.minimum(d, 30.0))))


def wind_speed(month: float, p: PrecipParams) -> float:
    return p.wind_speed_ms + p.wind_speed_seasonal_amplitude_ms * np.cos(
        2.0 * np.pi * (month - p.wind_peak_month) / 12.0
    )


def precip_month(
    uplift_h: np.ndarray,
    res_m: float,
    month: float,
    t_sea_level_c: float,
    env_lapse_c_per_km: float,
    p: PrecipParams | None = None,
    days: float = 365.25 / 12.0,
    reference_flux: float | None = None,
) -> tuple[np.ndarray, dict]:
    """Total monthly precipitation (mm) plus the state used to produce it.

    `t_sea_level_c` is that month's domain-mean sea-level temperature, taken
    straight from the temperature model — this is the coupling that makes the
    seasonal precipitation cycle a *consequence* of the seasonal temperature
    and wind cycles rather than an independent set of knobs.

    `reference_flux` is the annual-mean moisture flux `U * rho_s(T)`, used to
    scale the non-orographic background with the season the same way the
    orographic term scales. Pass the value returned in the annual mean run,
    or leave it None to keep the background constant through the year.
    """
    p = p or PrecipParams()
    u_mag = wind_speed(month, p)
    u, v = wind_vector(p.wind_from_deg, u_mag)
    cw, hw = thermo_from_temperature(t_sea_level_c, env_lapse_c_per_km)

    anomaly = orographic_precip(
        uplift_h,
        res_m,
        u=u,
        v=v,
        cw=cw * p.orographic_duty_cycle,
        hw=hw,
        nm=p.nm_moist_stability,
        tau_c=p.tau_c_s,
        tau_f=p.tau_f_s,
        seconds=days * 86400.0,
        pad_km=p.pad_km,
    )

    flux = u_mag * float(saturation_vapour_density(t_sea_level_c))
    bg = p.background_precip_mm_per_month * (days / (365.25 / 12.0))
    if reference_flux:
        bg *= flux / reference_flux
    floor = bg * p.lee_floor_fraction

    total = _softplus_floor(bg + anomaly, floor, max(floor * 0.5, 1e-6))
    state = {
        "month": month,
        "wind_speed_ms": u_mag,
        "t_sea_level_c": t_sea_level_c,
        "cw_kg_m3": cw,
        "cw_effective_kg_m3": cw * p.orographic_duty_cycle,
        "hw_m": hw,
        "background_mm": bg,
        "moisture_flux": flux,
    }
    return total.astype(np.float32), state
