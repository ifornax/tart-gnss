"""Direction uncertainties: local (Jacobian) covariance and baseline bootstrap.

Local covariance
----------------
At the solution ``p = (l, m, Re C, Im C)`` the analytic Jacobian ``J`` of the
stacked weighted residuals ``[Re r, Im r]`` gives

    Sigma_p = sigma_hat^2 (J^T J)^-1,

where the weights are already inside ``J`` (``J^T W J`` in the usual notation)
and ``sigma_hat^2`` is the per-component residual variance:
``chi2 / (2N - 4)`` for the linear loss, ``(1.4826 MAD)^2`` of the residual
components for robust losses (insensitive to the outliers the robust loss
down-weights). The ``(l, m)`` block of ``Sigma_p`` is the direction covariance
(the amplitude is marginalised, not fixed).

It is propagated to azimuth/elevation with the analytic Jacobian
:func:`tart_gnss.stage3.geometry.azel_jacobian`. That linearisation fails near
the zenith, where azimuth is undefined: when the zenith distance
``r = sqrt(l^2 + m^2)`` is below ``3 x`` the largest ``(l, m)`` standard
deviation, the estimate is flagged ``ZENITH_AZ_UNCERTAIN`` and ``sigma_az`` is
reported as NaN (Stage 4 should use the ``(l, m)`` covariance, which stays valid).

Bootstrap
---------
Baselines (antenna pairs, with all their samples in the bin) are resampled
with replacement and the fit is refined again from the point estimate.
Azimuths are unwrapped about the point estimate before percentiles are taken.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

from tart_gnss.stage3.errors import CoordinateTransformError
from tart_gnss.stage3.estimator import (
    DirectionCovariance,
    DirectionEstimate,
    _clip_lm,
    refine,
    robust_sigma,
)
from tart_gnss.stage3.geometry import azel_jacobian, lm_to_azel, wrap_azimuth_near
from tart_gnss.stage3.visibility_model import VisibilitySamples, residual_jacobian, residuals

ZENITH_SIGMA_FACTOR = 3.0


def estimate_direction_uncertainty(result: DirectionEstimate, samples: VisibilitySamples) -> DirectionCovariance:
    """Estimate the local covariance of a fitted direction (see module docstring).

    Returns a covariance whose entries are NaN if ``J^T J`` is singular.
    """
    c = result.complex_amplitude
    jac = residual_jacobian(samples, result.l, result.m, c)
    r = residuals(samples, result.l, result.m, c)
    stacked = np.concatenate([r.real, r.imag])
    dof = max(stacked.size - 4, 1)
    if result.robust_loss == "linear":
        sigma2 = float(stacked @ stacked) / dof
        method = "sigma^2 (J^T J)^-1, sigma^2 = chi2/(2N-4)"
    else:
        sigma2 = robust_sigma(stacked) ** 2
        method = "sigma^2 (J^T J)^-1, sigma = 1.4826 MAD of residual components"
    try:
        jtj = jac.T @ jac
        if np.linalg.cond(jtj) > 1e14:
            raise np.linalg.LinAlgError("ill-conditioned")
        cov = sigma2 * np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        nan = float("nan")
        return DirectionCovariance(nan, nan, nan, nan, nan, nan, False, method + " (singular)")
    cov_lm = cov[:2, :2]
    return propagate_to_azel(result.l, result.m, cov_lm, method)


def propagate_to_azel(l: float, m: float, cov_lm: np.ndarray, method: str = "") -> DirectionCovariance:
    """Propagate an ``(l, m)`` covariance to azimuth/elevation (degrees)."""
    sig_l, sig_m = float(np.sqrt(cov_lm[0, 0])), float(np.sqrt(cov_lm[1, 1]))
    r = float(np.hypot(l, m))
    max_sigma = float(np.sqrt(np.max(np.linalg.eigvalsh(cov_lm))))
    zenith = r < ZENITH_SIGMA_FACTOR * max_sigma
    try:
        g = azel_jacobian(l, m)
        cov_ae = np.degrees(1.0) ** 2 * (g @ cov_lm @ g.T)
        s_az, s_el, c_ae = float(np.sqrt(cov_ae[0, 0])), float(np.sqrt(cov_ae[1, 1])), float(cov_ae[0, 1])
    except CoordinateTransformError:
        s_az = s_el = c_ae = float("nan")
        zenith = zenith or r < 1e-12
    if zenith:
        s_az, c_ae = float("nan"), float("nan")
    return DirectionCovariance(sig_l, sig_m, float(cov_lm[0, 1]), s_az, s_el, c_ae, bool(zenith), method)


@dataclass
class BootstrapResult:
    """Empirical percentiles from a baseline bootstrap."""

    n_requested: int
    n_success: int
    az_p16: float
    az_p50: float
    az_p84: float
    el_p16: float
    el_p50: float
    el_p84: float

    def to_dict(self) -> Dict[str, float]:
        """Plain dictionary (for table rows)."""
        return dict(self.__dict__)


def bootstrap_direction(
    result: DirectionEstimate,
    samples: VisibilitySamples,
    n_boot: int,
    rng: np.random.Generator,
    max_nfev: int = 100,
) -> Optional[BootstrapResult]:
    """Baseline bootstrap of a fitted direction (see module docstring).

    Returns ``None`` if ``n_boot <= 0``.
    """
    if n_boot <= 0:
        return None
    ids = samples.baseline_ids()
    uniq, inverse = np.unique(ids, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(uniq.size)]
    az, el = [], []
    for _ in range(n_boot):
        pick = rng.integers(0, uniq.size, uniq.size)
        index = np.concatenate([members[i] for i in pick])
        try:
            res = refine(samples.subset(index), result.l, result.m, result.complex_amplitude,
                         result.robust_loss, result.f_scale, max_nfev)
        except (ValueError, np.linalg.LinAlgError):
            continue
        if not res.success:
            continue
        l, m = _clip_lm(res.x[0], res.x[1])
        a, e = lm_to_azel(l, m)
        az.append(float(a))
        el.append(float(e))
    if not az:
        nan = float("nan")
        return BootstrapResult(n_boot, 0, nan, nan, nan, nan, nan, nan)
    az_u = wrap_azimuth_near(np.array(az), result.az_deg)
    p_az = np.percentile(az_u, [16, 50, 84])
    p_el = np.percentile(el, [16, 50, 84])
    return BootstrapResult(n_boot, len(az), *map(float, p_az), *map(float, p_el))
