"""
Snow, seasonality, and permanent-snow (Tappa 3 — derived climate metrics).

This module answers the three questions `00_pre_project_planning.md` set out
for this stage: months-with-snow (count), a seasonality/escursion index, and
a permanent snow line. All three are built directly on Tappa 2's monthly
temperature and precipitation stacks — no new climate physics, only
derived post-processing of what Tappa 2 already produced.

Why this isn't just a temperature threshold
--------------------------------------------
`02_tappa2_climate.md` S5/S8 flagged the naive definition used there
("permanent snow = all 12 monthly means below 0 C") as a known
overestimate of elevation (equivalently, an underestimate of area): real
Southern Alps glaciers persist well below the elevation where the warmest
month's mean temperature crosses 0 C, precisely because snowfall
accumulation there is large enough to outlast the melt season. A pure
temperature threshold cannot represent that; it has no notion of how much
snow actually falls.

The fix here has two parts, both needed together:

1. **Rain/snow partition** (`snow_fraction`): monthly precipitation is
   split into snow and rain by monthly-mean temperature, not counted as
   "all-or-nothing" snow below 0 C. This distinguishes a cold-and-dry month
   from a cold-and-wet one, which the original count could not.
2. **A minimal mass-balance model** (`annual_mass_balance`): monthly snow
   accumulation is weighed against monthly melt (a positive-degree-day
   model), and "permanent snow" is redefined as *cells where the annual
   balance does not go negative* — i.e. where accumulation keeps up with
   melt over a full year — rather than a temperature threshold at all.

Both pieces share one modelling device, used nowhere else in this project:
treating each month's daily temperature as **normally distributed around
the monthly mean**, with an assumed day-to-day standard deviation
`sigma_day_c`. This is the honest way to ask two questions that a monthly
mean cannot answer on its own — "what fraction of the month's days (and,
by assumption, precipitation mass) fell below the rain/snow threshold?"
and "how many degree-days above 0 C did this month actually accumulate?"
— without pretending the model has daily data it does not have. Both
reduce to closed-form expressions of a standard Normal, derived below.

Citations and honest limitations for every constant are in
`docs/decisions/03_tappa3_snow.md`; the summary:

* Rain/snow 50%-threshold `t50_snow_rain_c = 1.0` C, dry-bulb — Jennings et
  al. (2018), *Nature Communications* 9:1148, ~17,000 NH stations, 95% of
  stations between -0.4 and 2.4 C. This is the *between-station* spread of
  where the local threshold sits; it does not by itself set the
  within-month day-to-day sharpness of the transition (see next point).
* Day-to-day standard deviation `sigma_day_c = 3.0` C — **not** independently
  verified for this world's real-world validation reference (no NIWA figure
  for West Coast South Island daily variability was found); a documented
  estimate bridging a global generic daily-temperature std. dev. (~2 C) and
  a common generic ice-sheet-model default (~5 C). Open to recalibration —
  see the decision doc.
* Degree-day factor `degree_day_factor_mm_per_c_day = 4.5` — Hock (2003),
  *J. Hydrol.* 282, Table 1: observed snow DDFs cluster 3.2-7.6 mm w.e. per
  C per day (full reported range 2.5-16.9); 4.5 sits in the clustered range.

What this model does NOT do, by design, kept simple deliberately: no ice
dynamics/flow, no multi-year snowpack carry-over beyond "does this month's
accumulation exceed this month's melt", no refreezing of meltwater, no
rain-on-snow melt enhancement, no avalanche/wind redistribution of snow.
"Permanent snow" here is a single-year steady-state proxy for where an
equilibrium line altitude (ELA) would sit under a repeating climatology,
not a real glacier mass-balance model — see the decision doc for how this
is validated against real Southern Alps ELA figures despite that gap.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

__all__ = [
    "SnowParams",
    "snow_fraction",
    "expected_positive_degree_days",
    "monthly_snow_rain",
    "months_with_snow",
    "seasonality_index",
    "annual_mass_balance",
]


@dataclass
class SnowParams:
    t50_snow_rain_c: float = 1.0
    """Dry-bulb monthly-mean temperature at which precipitation is modelled
    as 50% snow / 50% rain. Jennings et al. (2018) — see module docstring."""

    sigma_day_c: float = 3.0
    """Assumed day-to-day (synoptic) standard deviation of daily temperature
    around the monthly mean. Drives both the sharpness of the rain/snow
    transition and the positive-degree-day estimate. See module docstring
    for why this is a documented estimate, not a verified figure."""

    degree_day_factor_mm_per_c_day: float = 4.5
    """Snow melt, mm water-equivalent per positive-degree-day. Hock (2003)
    — see module docstring."""

    min_snow_month_mm: float = 5.0
    """Minimum modelled monthly snowfall (mm w.e.) for a month to count
    towards `months_with_snow`. Filters out months where the snow fraction
    is technically nonzero (the Normal-CDF tail never truly reaches zero)
    but negligible — a smootherstep-style "effectively zero" cutoff, same
    spirit as Tappa 1's shelf taper (01_tappa1_terrain.md S3a)."""


def snow_fraction(t_mean_c: np.ndarray, p: SnowParams | None = None) -> np.ndarray:
    """Fraction of a month's precipitation mass falling as snow.

    Modelled as the fraction of the month's days with daily temperature
    below the rain/snow threshold, under the assumption that daily
    temperature is Normal(t_mean_c, sigma_day_c) and that precipitation
    intensity is uncorrelated with which specific days are warm or cold
    within the month (stated simplification, not measured).

    fraction = P(T_daily < t50) = Phi((t50 - t_mean_c) / sigma_day_c)

    Monotonic in t_mean_c by construction: 1.0 in the deep cold, 0.0 in the
    warm limit, and it is `sigma_day_c` — not a separately fitted curve
    width — that sets how gradual the transition is. A colder monthly mean
    always gives a fraction at least as high as a warmer one.
    """
    p = p or SnowParams()
    return stats.norm.cdf((p.t50_snow_rain_c - t_mean_c) / p.sigma_day_c)


def expected_positive_degree_days(t_mean_c: np.ndarray, days_in_month: float, p: SnowParams | None = None) -> np.ndarray:
    """Expected sum of daily max(T, 0) over a month, i.e. monthly PDD (C.day).

    Closed form for E[max(X, 0)] when X ~ Normal(mu, sigma^2):

        E[max(X, 0)] = mu * Phi(mu / sigma) + sigma * phi(mu / sigma)

    This is the standard rectified-Gaussian expectation (the same
    expression that prices a financial call option under a Normal, or gives
    the mean of a half-truncated Normal). It matters here because a month
    whose *mean* temperature sits just below 0 C still has some warm days
    contributing real melt — a plain "PDD = max(mean, 0) * days" would
    silently zero out melt in every such month, understating ablation right
    around the elevation the permanent-snow boundary actually sits at.
    """
    p = p or SnowParams()
    mu, sigma = t_mean_c, p.sigma_day_c
    z = mu / sigma
    pdd_per_day = mu * stats.norm.cdf(z) + sigma * stats.norm.pdf(z)
    return pdd_per_day * days_in_month


def monthly_snow_rain(precip_mm: np.ndarray, temp_c: np.ndarray, p: SnowParams | None = None):
    """Split a (12, ny, nx) monthly precipitation stack into (snow_mm, rain_mm)."""
    p = p or SnowParams()
    frac = snow_fraction(temp_c, p)
    snow_mm = precip_mm * frac
    rain_mm = precip_mm * (1.0 - frac)
    return snow_mm.astype(np.float32), rain_mm.astype(np.float32)


def months_with_snow(snow_mm: np.ndarray, p: SnowParams | None = None) -> np.ndarray:
    """Count of months (0-12) whose modelled snowfall clears `min_snow_month_mm`."""
    p = p or SnowParams()
    return (snow_mm >= p.min_snow_month_mm).sum(axis=0).astype(np.int16)


def seasonality_index(temp_c: np.ndarray) -> np.ndarray:
    """Escursion/seasonality index: warmest-month mean minus coldest-month mean, per cell.

    No methodological choice here beyond what Tappa 2 already computed —
    this is a direct read of the monthly temperature stack's own range.
    """
    return (temp_c.max(axis=0) - temp_c.min(axis=0)).astype(np.float32)


def annual_mass_balance(
    precip_mm: np.ndarray,
    temp_c: np.ndarray,
    days_in_month: np.ndarray,
    p: SnowParams | None = None,
):
    """Annual snow accumulation, melt, and net balance per cell.

    Returns (accum_mm, melt_mm, balance_mm), each (ny, nx), all annual sums.
    `balance_mm >= 0` is this project's permanent-snow / perennial-firn
    proxy: accumulation keeps pace with melt over a full repeating annual
    cycle, so snow set down this year is not fully removed before the next
    accumulation season — the discrete, one-year analogue of an
    equilibrium-line altitude. It is not a substitute for a real
    multi-year glacier mass-balance model (see module docstring for what
    is deliberately left out).
    """
    p = p or SnowParams()
    snow_mm, _ = monthly_snow_rain(precip_mm, temp_c, p)
    accum_mm = snow_mm.sum(axis=0)

    n_months = temp_c.shape[0]
    melt_mm = np.zeros_like(accum_mm)
    for i in range(n_months):
        pdd = expected_positive_degree_days(temp_c[i], days_in_month[i], p)
        melt_mm += p.degree_day_factor_mm_per_c_day * pdd

    balance_mm = accum_mm - melt_mm
    return accum_mm.astype(np.float32), melt_mm.astype(np.float32), balance_mm.astype(np.float32)
