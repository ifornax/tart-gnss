"""Primary visibility-domain direction estimator.

Two stages, both fitting the complex visibilities directly (no wrapped
phase differences are ever subtracted):

1. **Coarse search.** The profiled objective ``chi2(l, m)`` (complex amplitude
   eliminated analytically, see :mod:`visibility_model`) is evaluated on a
   grid of directions uniformly spaced on the sphere inside a search cap of
   angular radius ``search_radius_deg``. The grid step is a fraction
   (``1/grid_oversample``) of the array resolution ``lambda / B_max``. The cap
   centre is the TLE-predicted direction when available; it only bounds the
   search and never enters the objective. ``blind=True`` searches the whole
   visible hemisphere instead.
2. **Local refinement.** All four parameters ``(l, m, Re C, Im C)`` are refined
   with :func:`scipy.optimize.least_squares` (``trf``, bounded ``-1 <= l, m <= 1``)
   on the stacked real residual vector ``[Re r, Im r]``, ``r = sqrt(w)(V - C e^{-i phi})``,
   with a ``linear``, ``huber`` or ``soft_l1`` loss. For the robust losses the
   loss scale is ``f_scale = huber_k * sigma_robust`` where
   ``sigma_robust = 1.4826 * MAD`` of the residual components of the coarse
   solution (otherwise the loss would act on arbitrary visibility units).

The solution is flagged ``BOUNDARY_SOLUTION`` (and never replaced by the cap
centre) if it lies on the cap edge, if the coarse optimum was on the edge, or
if it is a sidelobe of a brighter source just outside the cap: when the
profiled ``chi2`` in a guard band ``GUARD_BAND_RESOLUTIONS x lambda/B_max``
wide outside the cap is lower than anywhere inside, that source is subtracted
and the cap searched again; if the in-cap optimum then moves by more than half
a resolution element it was a sidelobe. (A genuine in-cap source with a
brighter, unrelated satellite nearby stays put and is not flagged.) Sources
further outside than the guard band cannot be detected this way, a known
limitation of any bounded search.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np

from tart_gnss.stage3.errors import DirectionFitError, InsufficientBaselinesError
from tart_gnss.stage3.geometry import (
    SPEED_OF_LIGHT_M_S,
    angular_separation_deg,
    azel_to_unit,
    lm_to_azel,
    offset_directions,
)
from tart_gnss.stage3.visibility_model import (
    VisibilitySamples,
    predict_visibilities,
    profiled_chi2,
    residual_jacobian,
    residuals,
)

logger = logging.getLogger(__name__)

ROBUST_LOSSES = ("linear", "huber", "soft_l1")
#: Normal-consistency factor for the median absolute deviation.
MAD_TO_SIGMA = 1.4826
#: Width of the guard band checked outside the search cap, in resolution elements.
GUARD_BAND_RESOLUTIONS = 2.0
#: Largest grid evaluated at once (directions); bounds memory use.
GRID_CHUNK = 4096


@dataclass
class SearchRegion:
    """A spherical cap of angular radius ``radius_deg`` about ENU direction ``centre``."""

    centre: np.ndarray
    radius_deg: float
    centre_source: str

    @classmethod
    def around_azel(cls, az_deg: float, el_deg: float, radius_deg: float,
                    source: str = "TLE prediction") -> "SearchRegion":
        """Cap centred on an azimuth/elevation."""
        return cls(azel_to_unit(az_deg, el_deg), float(radius_deg), source)

    @classmethod
    def hemisphere(cls) -> "SearchRegion":
        """The whole sky above the horizon (blind search)."""
        return cls(np.array([0.0, 0.0, 1.0]), 90.0, "blind (zenith, full hemisphere)")

    def separation_deg(self, s) -> np.ndarray:
        """Angular distance of direction(s) ``s`` from the cap centre."""
        return angular_separation_deg(s, self.centre)

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        from tart_gnss.stage3.geometry import unit_to_azel

        az, el = unit_to_azel(self.centre)
        return {"centre_az_deg": float(az), "centre_el_deg": float(el), "radius_deg": self.radius_deg,
                "centre_source": self.centre_source}


@dataclass
class EstimatorOptions:
    """Tunable estimator settings (all recorded in the Stage-3 manifest)."""

    robust_loss: str = "huber"
    search_radius_deg: float = 10.0
    blind: bool = False
    grid_oversample: float = 5.0
    huber_k: float = 1.345
    max_nfev: int = 200
    min_elevation_deg: float = 0.0

    def __post_init__(self) -> None:
        if self.robust_loss not in ROBUST_LOSSES:
            raise ValueError(f"robust_loss must be one of {ROBUST_LOSSES}")
        if not 0 < self.search_radius_deg <= 90:
            raise ValueError("search radius must be in (0, 90] degrees")

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        return dict(self.__dict__)


@dataclass
class DirectionCovariance:
    """Local covariance of a fitted direction."""

    sigma_l: float
    sigma_m: float
    cov_lm: float
    sigma_az_deg: float
    sigma_el_deg: float
    cov_az_el: float
    zenith_az_uncertain: bool = False
    method: str = ""


@dataclass
class DirectionEstimate:
    """Result of one visibility-domain fit (one station, one time bin)."""

    station_id: str
    l: float
    m: float
    az_deg: float
    el_deg: float
    complex_amplitude: complex
    chi2: float
    reduced_chi2: float
    rms_complex_residual: float
    median_abs_residual: float
    max_abs_residual: float
    sigma_vis: float
    snr: float
    residual_tail_ratio: float
    n_samples: int
    n_baselines_used: int
    success: bool
    optimizer_status: str
    boundary: bool
    invalid_geometry: bool
    f_scale: float
    robust_loss: str
    grid_step_deg: float
    n_grid: int
    separation_from_centre_deg: float
    norad_id: Optional[int] = None
    time_utc: Optional[datetime] = None
    covariance: Optional[DirectionCovariance] = None
    quality_flags: List[str] = field(default_factory=list)

    @property
    def unit_vector(self) -> np.ndarray:
        """ENU unit vector of the estimate."""
        n = np.sqrt(max(0.0, 1.0 - self.l ** 2 - self.m ** 2))
        return np.array([self.l, self.m, n])


def resolution_rad(samples: VisibilitySamples) -> float:
    """Array resolution ``lambda_min / B_max`` (radians)."""
    b = samples.max_baseline_m()
    if b <= 0:
        raise InsufficientBaselinesError("no non-zero baselines")
    return SPEED_OF_LIGHT_M_S / float(samples.frequency_hz.max()) / b


def search_grid(region: SearchRegion, step_rad: float, min_elevation_deg: float = 0.0) -> np.ndarray:
    """Directions on a square tangent-plane grid (exact great-circle offsets) inside the cap.

    Directions below ``min_elevation_deg`` (or the horizon) are dropped.
    """
    radius = np.radians(region.radius_deg)
    n = int(np.ceil(radius / step_rad))
    x = np.arange(-n, n + 1) * step_rad
    xx, yy = np.meshgrid(x, x)
    keep = np.hypot(xx, yy) <= radius + 1e-12
    s = offset_directions(region.centre, xx[keep], yy[keep])
    floor = np.sin(np.radians(max(min_elevation_deg, 0.0))) + 1e-9
    return s[s[:, 2] > floor]


def robust_sigma(values: np.ndarray) -> float:
    """``1.4826 * median(|x - median(x)|)``."""
    v = np.asarray(values, dtype=float)
    return float(MAD_TO_SIGMA * np.median(np.abs(v - np.median(v))))


def coarse_search(samples: VisibilitySamples, grid: np.ndarray):
    """Evaluate the profiled ``chi2`` on ``grid``; return (best index, chi2 array, C_hat array)."""
    chi2 = np.empty(len(grid))
    amp = np.empty(len(grid), dtype=complex)
    for i in range(0, len(grid), GRID_CHUNK):
        c2, a = profiled_chi2(samples, grid[i:i + GRID_CHUNK])
        chi2[i:i + GRID_CHUNK], amp[i:i + GRID_CHUNK] = c2, a
    return int(np.argmin(chi2)), chi2, amp


def _stack(r: np.ndarray) -> np.ndarray:
    return np.concatenate([r.real, r.imag])


def _clip_lm(l: float, m: float):
    r = np.hypot(l, m)
    if r >= 1.0:
        scale = (1.0 - 1e-12) / r
        return l * scale, m * scale
    return l, m


def refine(
    samples: VisibilitySamples,
    l0: float,
    m0: float,
    c0: complex,
    robust_loss: str,
    f_scale: float,
    max_nfev: int = 200,
):
    """Refine ``(l, m, Re C, Im C)`` with :func:`scipy.optimize.least_squares`.

    Returns the scipy ``OptimizeResult``. ``l^2 + m^2 > 1`` is prevented by
    radial clipping inside the residual function (the caller flags solutions
    on the unit circle).
    """
    from scipy.optimize import least_squares

    def fun(p):
        l, m = _clip_lm(p[0], p[1])
        return _stack(residuals(samples, l, m, complex(p[2], p[3])))

    def jac(p):
        l, m = _clip_lm(p[0], p[1])
        return residual_jacobian(samples, l, m, complex(p[2], p[3]))

    big = np.inf
    lo, hi = [-1.0, -1.0, -big, -big], [1.0, 1.0, big, big]
    p0 = np.array([np.clip(l0, -1 + 1e-9, 1 - 1e-9), np.clip(m0, -1 + 1e-9, 1 - 1e-9), c0.real, c0.imag])
    kwargs = {} if robust_loss == "linear" else {"f_scale": f_scale}
    return least_squares(fun, p0, jac=jac, bounds=(lo, hi), method="trf", loss=robust_loss,
                         x_scale="jac", max_nfev=max_nfev, **kwargs)


def fit_statistics(samples: VisibilitySamples, l: float, m: float, c: complex) -> Dict[str, float]:
    """Goodness-of-fit statistics of a fitted point-source model.

    * ``chi2 = sum w |V - M|^2`` and ``reduced_chi2 = chi2 / (2N - 4)``. These
      have a probabilistic chi-square meaning only if the weights are inverse
      noise variances (for TART ``WEIGHT = 1``, so they are diagnostics in
      squared visibility units).
    * ``rms_complex_residual = sqrt(sum w |r|^2 / sum w)``; median/max of ``|V - M|``.
    * ``sigma_vis``: robust per-component noise ``1.4826 * MAD`` of the stacked
      real and imaginary parts of ``sqrt(w)(V - M)``.
    * ``snr = |C_hat| sqrt(sum w / 2) / sigma_vis``: the amplitude estimate
      divided by its standard error (the variance of ``C_hat`` is
      ``2 sigma^2 / sum w`` for per-component noise ``sigma / sqrt(w)``), i.e.
      the matched-filter detection significance of the point source.
    * ``residual_tail_ratio = sqrt(chi2 / (2N - 4)) / sigma_vis``: ~1 for
      Gaussian residuals; larger values indicate outliers or unmodelled
      structure (e.g. other satellites).
    """
    r = samples.vis - predict_visibilities(samples.baseline_enu_m, samples.frequency_hz, l, m, c)
    w = samples.weight
    chi2 = float(np.sum(w * np.abs(r) ** 2))
    dof = max(2 * samples.n - 4, 1)
    sigma = robust_sigma(_stack(np.sqrt(w) * r))
    sum_w = float(w.sum())
    snr = abs(c) * np.sqrt(sum_w / 2.0) / sigma if sigma > 0 else float("inf")
    return {
        "chi2": chi2,
        "dof": dof,
        "reduced_chi2": chi2 / dof,
        "rms_complex_residual": float(np.sqrt(chi2 / sum_w)),
        "median_abs_residual": float(np.median(np.abs(r))),
        "max_abs_residual": float(np.abs(r).max()),
        "sigma_vis": sigma,
        "snr": float(snr),
        "residual_tail_ratio": float(np.sqrt(chi2 / dof) / sigma) if sigma > 0 else float("inf"),
    }


def _sidelobe_of_outside_source(samples, region, grid, best, chi2_grid, step, res_rad, min_el) -> bool:
    """True if the in-cap optimum is a sidelobe of a brighter source just outside the cap.

    The guard band (``GUARD_BAND_RESOLUTIONS`` resolution elements outside the
    cap) is searched. If it holds a lower ``chi2`` than the cap, that point
    source is subtracted and the cap is searched again: a genuine in-cap source
    stays put, a sidelobe moves by more than half a resolution element.
    """
    guard_deg = min(region.radius_deg + GUARD_BAND_RESOLUTIONS * np.degrees(res_rad), 90.0)
    wide = search_grid(SearchRegion(region.centre, guard_deg, region.centre_source), step, min_el)
    band = wide[region.separation_deg(wide) > region.radius_deg]
    if not len(band):
        return False
    b_best, b_chi2, b_amp = coarse_search(samples, band)
    if not b_chi2[b_best] < chi2_grid[best]:
        return False
    s_out = band[b_best]
    cleaned = VisibilitySamples(
        samples.baseline_enu_m, samples.frequency_hz,
        samples.vis - predict_visibilities(samples.baseline_enu_m, samples.frequency_hz,
                                           s_out[0], s_out[1], b_amp[b_best]),
        samples.weight, samples.antenna1, samples.antenna2, samples.time_mjd_seconds,
    )
    best2, _, _ = coarse_search(cleaned, grid)
    moved = float(angular_separation_deg(grid[best2], grid[best]))
    logger.debug("brighter source %.2f deg outside the cap; in-cap optimum moves %.3f deg after subtracting it",
                 float(region.separation_deg(s_out)) - region.radius_deg, moved)
    return moved > 0.5 * np.degrees(res_rad)


def estimate_direction(
    samples: VisibilitySamples,
    initial_az_deg: Optional[float] = None,
    initial_el_deg: Optional[float] = None,
    search_radius_deg: float = 10.0,
    robust_loss: str = "huber",
    options: Optional[EstimatorOptions] = None,
    station_id: str = "",
    n_baselines_used: Optional[int] = None,
) -> DirectionEstimate:
    """Estimate source direction from complex interferometric visibilities.

    ``initial_az_deg``/``initial_el_deg`` (e.g. a TLE prediction) only centre
    the coarse search cap. Without them, or with ``options.blind``, the whole
    hemisphere is searched. Uncertainties are added separately
    (:func:`tart_gnss.stage3.uncertainty.estimate_direction_uncertainty`).

    Raises
    ------
    InsufficientBaselinesError
        If there are too few samples to constrain 4 real parameters.
    DirectionFitError
        If the search region contains no direction above the horizon.
    """
    opts = options or EstimatorOptions(robust_loss=robust_loss, search_radius_deg=search_radius_deg)
    if samples.n < 3:
        raise InsufficientBaselinesError(f"{samples.n} valid visibilities; at least 3 are needed")
    if opts.blind or initial_az_deg is None or initial_el_deg is None:
        region = SearchRegion.hemisphere()
    else:
        region = SearchRegion.around_azel(initial_az_deg, initial_el_deg, opts.search_radius_deg)
    res_rad = resolution_rad(samples)
    step = res_rad / opts.grid_oversample
    grid = search_grid(region, step, opts.min_elevation_deg)
    if len(grid) == 0:
        raise DirectionFitError("search region contains no direction above the horizon")
    best, chi2_grid, amp_grid = coarse_search(samples, grid)
    s0, c0 = grid[best], amp_grid[best]
    coarse_edge = bool(region.radius_deg < 90.0 and
                       region.separation_deg(s0) > region.radius_deg - 1.5 * np.degrees(step))
    sidelobe = False
    if region.radius_deg < 90.0:
        sidelobe = _sidelobe_of_outside_source(samples, region, grid, best, chi2_grid, step, res_rad,
                                               opts.min_elevation_deg)

    lin = residuals(samples, s0[0], s0[1], c0)
    f_scale = opts.huber_k * robust_sigma(_stack(lin))
    if not f_scale > 0:
        f_scale = 1.0
    try:
        res = refine(samples, s0[0], s0[1], c0, opts.robust_loss, f_scale, opts.max_nfev)
        l, m = _clip_lm(res.x[0], res.x[1])
        c = complex(res.x[2], res.x[3])
        success = bool(res.success) and np.all(np.isfinite(res.x))
        status = f"{res.status}: {res.message}"
    except (ValueError, np.linalg.LinAlgError) as exc:
        logger.debug("refinement failed: %s", exc)
        l, m, c, success, status = s0[0], s0[1], c0, False, f"exception: {exc}"
    r2 = l * l + m * m
    n = np.sqrt(max(0.0, 1.0 - r2))
    s_hat = np.array([l, m, n])
    sep = float(region.separation_deg(s_hat))
    on_horizon = r2 >= 1.0 - 1e-9
    boundary = bool(coarse_edge or sidelobe or on_horizon or
                    (region.radius_deg < 90.0 and sep > region.radius_deg - 0.5 * np.degrees(step)))
    az, el = lm_to_azel(l, m)
    stats = fit_statistics(samples, l, m, c)
    logger.debug("%s fit: l=%.5f m=%.5f |C|=%.4g snr=%.1f nfev=%s status=%s", station_id, l, m, abs(c),
                 stats["snr"], getattr(locals().get("res"), "nfev", None), status)
    return DirectionEstimate(
        station_id=station_id, l=float(l), m=float(m), az_deg=float(az), el_deg=float(el),
        complex_amplitude=c, chi2=stats["chi2"], reduced_chi2=stats["reduced_chi2"],
        rms_complex_residual=stats["rms_complex_residual"], median_abs_residual=stats["median_abs_residual"],
        max_abs_residual=stats["max_abs_residual"], sigma_vis=stats["sigma_vis"], snr=stats["snr"],
        residual_tail_ratio=stats["residual_tail_ratio"], n_samples=samples.n,
        n_baselines_used=int(n_baselines_used if n_baselines_used is not None
                             else np.unique(samples.baseline_ids()).size),
        success=bool(success), optimizer_status=status, boundary=boundary,
        invalid_geometry=bool(on_horizon or el < opts.min_elevation_deg), f_scale=float(f_scale),
        robust_loss=opts.robust_loss, grid_step_deg=float(np.degrees(step)), n_grid=int(len(grid)),
        separation_from_centre_deg=sep,
    )
