"""Descriptive time-series statistics of the residuals.

Residual epochs are *not* independent: the Stage-5 smoothed trajectory is
correlated from epoch to epoch by construction. Summaries are therefore
descriptive (mean, median, RMS, MAD scale, percentiles, extremes); no
p-value or chi-square is formed over epochs.

Trend
    ``dr_j(t) = a_j + b_j (t - t0)`` by weighted least squares with the
    Stage-5 variance of the component (projected into RTN) as weights. The
    slope is a mismatch *rate*, not an orbit correction. Its formal sigma
    assumes independent epochs and is reported only for reference.
    ``STRONG_RESIDUAL_TREND`` is descriptive: the fitted change over the
    arc, ``|b| * span``, exceeds 3 times the median Stage-5 sigma of that
    component.

Autocorrelation
    Lag-1 autocorrelation of each mean-removed component, and the
    corresponding effective sample size ``n (1 - rho) / (1 + rho)``.
    ``RESIDUAL_AUTOCORRELATED`` uses the white-noise bound ``3 / sqrt(n)``;
    a smoothed trajectory nearly always exceeds it, which is expected and
    does not by itself indicate a problem.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

MAD_TO_SIGMA = 1.4826


def describe(x) -> Dict[str, float]:
    """Mean, median, RMS, MAD-based scale, 16/84th percentiles, min, max and count of finite values."""
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    if not v.size:
        nan = float("nan")
        return {"n": 0, "mean": nan, "median": nan, "rms": nan, "mad_sigma": nan, "p16": nan, "p84": nan,
                "min": nan, "max": nan}
    med = float(np.median(v))
    return {
        "n": int(v.size), "mean": float(np.mean(v)), "median": med, "rms": float(np.sqrt(np.mean(v * v))),
        "mad_sigma": float(MAD_TO_SIGMA * np.median(np.abs(v - med))),
        "p16": float(np.percentile(v, 16)), "p84": float(np.percentile(v, 84)),
        "min": float(np.min(v)), "max": float(np.max(v)),
    }


def trend(t_s, y, sigma: Optional[np.ndarray] = None) -> Dict[str, float]:
    """Weighted linear fit ``y = a + b (t - t0)``; ``t0`` is the first finite epoch."""
    t = np.asarray(t_s, dtype=float)
    y = np.asarray(y, dtype=float)
    s = np.ones_like(y) if sigma is None else np.asarray(sigma, dtype=float)
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    nan = float("nan")
    if ok.sum() < 2 or np.ptp(t[ok]) == 0:
        return {"intercept": nan, "slope": nan, "slope_sigma_formal": nan, "span_s": 0.0, "change_over_span": nan,
                "weighted": sigma is not None}
    t0 = t[ok][0]
    a_mat = np.vstack([np.ones(ok.sum()), t[ok] - t0]).T
    w = 1.0 / s[ok] ** 2
    n_mat = a_mat.T @ (a_mat * w[:, None])
    coef = np.linalg.solve(n_mat, a_mat.T @ (w * y[ok]))
    cov = np.linalg.inv(n_mat)
    span = float(np.ptp(t[ok]))
    return {"intercept": float(coef[0]), "slope": float(coef[1]), "slope_sigma_formal": float(np.sqrt(cov[1, 1])),
            "span_s": span, "change_over_span": float(coef[1] * span), "weighted": sigma is not None}


def lag1_autocorrelation(x) -> float:
    """Lag-1 autocorrelation of the mean-removed finite sequence."""
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return float("nan")
    v = v - v.mean()
    den = float(np.sum(v * v))
    return float(np.sum(v[:-1] * v[1:]) / den) if den > 0 else float("nan")


def effective_sample_size(n: int, rho: float) -> float:
    """``n (1 - rho) / (1 + rho)`` for an AR(1)-like sequence (``n`` if rho is not positive)."""
    if not np.isfinite(rho) or rho <= 0:
        return float(n)
    return float(n * (1 - rho) / (1 + rho))
