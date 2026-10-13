"""Innovation statistics: TLE-independent model-consistency diagnostics.

For a correct model and correct covariances the innovations ``nu_k`` are
zero-mean, white, with covariance ``S_k``; the NIS ``nu^T S^-1 nu`` is
``chi^2_3`` (mean 3, median 2.37). Statistics are taken over the measurements
the quality policy offered to the filter (``attempted``), gated ones included,
after a burn-in of :data:`BURN_IN` measurements: the filter starts from an
almost uninformative prior, so the first two innovations (which fix the six
state parameters) have ``S`` dominated by the prior and NIS near zero.

Whiteness uses the whitened innovations ``e_k = L_k^-1 nu_k``
(``S_k = L_k L_k^T``): the lag-``j`` autocorrelation
``rho_j = sum e_k . e_k+j / sum |e_k|^2`` over consecutive attempted
measurements has standard deviation ``~ 1 / sqrt(3 n_pairs)`` when white.

Run-level flags
    ``HIGH_NIS``          either the mean NIS of the innovations inside the
                          gate exceeds the ``nis_probability`` quantile of
                          ``chi^2_{3n} / n``, or more innovations exceed the
                          gate than the ``nis_probability`` quantile of
                          ``Binomial(n, 1 - gate probability)`` allows
                          (covariances too small or model wrong). Isolated
                          gated outliers therefore do not raise it; the
                          reported ``mean_nis`` includes them.
    ``DYNAMICS_MISMATCH`` lag-1 autocorrelation above ``3 / sqrt(3 n_pairs)``
                          (correlated innovations: unmodelled dynamics)
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from tart_gnss.stage5.filter import FilterResult

BURN_IN = 2
MIN_FOR_FLAGS = 5
DEFAULT_NIS_PROBABILITY = 0.9973


def whitened(innovation: np.ndarray, innovation_cov: np.ndarray) -> np.ndarray:
    """``L^-1 nu`` per row (NaN where unavailable)."""
    out = np.full_like(innovation, np.nan)
    for k, (nu, s) in enumerate(zip(innovation, innovation_cov)):
        if np.all(np.isfinite(nu)) and np.all(np.isfinite(s)):
            try:
                out[k] = np.linalg.solve(np.linalg.cholesky(s), nu)
            except np.linalg.LinAlgError:
                pass
    return out


def statistics_mask(fr: FilterResult) -> np.ndarray:
    """Measurements entering the statistics: attempted, finite NIS, after the burn-in."""
    m = fr.attempted & np.isfinite(fr.nis)
    idx = np.flatnonzero(m)
    m[idx[:BURN_IN]] = False
    return m


def autocorrelation(white: np.ndarray, mask: np.ndarray, max_lag: int = 3) -> Dict[int, float]:
    """Lag-``j`` autocorrelation of the whitened innovations over the masked sequence."""
    e = white[mask]
    den = float(np.sum(e * e))
    out = {}
    for lag in range(1, max_lag + 1):
        if len(e) > lag and den > 0:
            out[lag] = float(np.sum(e[:-lag] * e[lag:]) / den)
        else:
            out[lag] = float("nan")
    return out


def summarise(fr: FilterResult, nis_probability: float = DEFAULT_NIS_PROBABILITY) -> Dict[str, object]:
    """Innovation summary and run-level flags (see the module docstring)."""
    from scipy.stats import chi2

    mask = statistics_mask(fr)
    nis = fr.nis[mask]
    n = int(mask.sum())
    white = whitened(fr.innovation, fr.innovation_cov)
    acf = autocorrelation(white, mask)
    n_pairs = max(n - 1, 0)
    acf_limit = 3.0 / np.sqrt(3 * n_pairs) if n_pairs else float("nan")
    nis_limit = float(chi2.ppf(nis_probability, 3 * n) / n) if n else float("nan")  # reference only
    inn = fr.innovation[mask]
    flags: List[str] = []
    from scipy.stats import binom

    mean_nis = float(np.mean(nis)) if n else float("nan")
    inside = nis[nis <= fr.gate_threshold]
    mean_inside = float(np.mean(inside)) if inside.size else float("nan")
    inside_limit = float(chi2.ppf(nis_probability, 3 * inside.size) / inside.size) if inside.size else float("nan")
    n_above = int(np.sum(nis > fr.gate_threshold))
    above_limit = int(binom.ppf(nis_probability, n, 1.0 - fr.config.gate_probability)) if n else 0
    if n >= MIN_FOR_FLAGS and (not inside.size or mean_inside > inside_limit or n_above > above_limit):
        flags.append("HIGH_NIS")
    if n >= MIN_FOR_FLAGS and acf[1] > acf_limit:
        flags.append("DYNAMICS_MISMATCH")
    att = int(fr.attempted.sum())
    return {
        "n_attempted": att,
        "n_accepted": int(fr.accepted.sum()),
        "n_outliers": int(fr.outlier.sum()),
        "accepted_fraction": float(fr.accepted.sum() / att) if att else float("nan"),
        "burn_in": BURN_IN,
        "n_statistics": n,
        "mean_nis": mean_nis,
        "median_nis": float(np.median(nis)) if n else float("nan"),
        "expected_mean_nis": 3.0,
        "expected_median_nis": float(chi2.ppf(0.5, 3)),
        "mean_nis_limit": nis_limit,
        "mean_nis_inside_gate": mean_inside,
        "mean_nis_inside_gate_limit": inside_limit,
        "n_above_gate": n_above,
        "n_above_gate_limit": above_limit,
        "fraction_above_gate": float(np.mean(nis > fr.gate_threshold)) if n else float("nan"),
        "gate_threshold": fr.gate_threshold,
        "innovation_rms_m": float(np.sqrt(np.mean(np.sum(inn ** 2, axis=1)))) if n else float("nan"),
        "innovation_rms_xyz_m": [float(v) for v in np.sqrt(np.mean(inn ** 2, axis=0))] if n else [],
        "autocorrelation": {str(k): v for k, v in acf.items()},
        "autocorrelation_limit_lag1": acf_limit,
        "flags": flags,
    }


def dynamics_diagnostics(summary: Dict[str, object]) -> Dict[str, object]:
    """The TLE-independent subset used to compare dynamics models."""
    return {k: summary[k] for k in ("innovation_rms_m", "mean_nis", "median_nis", "accepted_fraction")} | {
        "lag1_autocorrelation": summary["autocorrelation"].get("1")}
