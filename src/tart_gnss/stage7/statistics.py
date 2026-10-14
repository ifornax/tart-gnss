"""Rank correlation with time-trend control and an autocorrelation-aware effective sample size.

Diagnostics in one pass share smooth time trends (geometry, residuals and
metadata all drift during a pass), and two monotone series have a rank
correlation near 1 with no causal link. Every correlation is therefore
reported twice:

* ``coefficient``: the plain rank (or Pearson/Kendall) coefficient;
* ``partial_coefficient``: the rank correlation after removing a linear
  trend in time from the ranks of both variables. Decisions use this one.

Its p-value uses an effective sample size
``n_eff = n (1 - r_x r_y) / (1 + r_x r_y)`` (``r`` the lag-1 autocorrelations
of the two detrended sequences, Bretherton et al. 1999), so smooth,
autocorrelated series do not produce tiny p-values. P-values are descriptive.
"""

from __future__ import annotations

from typing import Dict

import numpy as np

METHODS = ("spearman", "pearson", "kendall")


def _ranks(x):
    from scipy.stats import rankdata

    return rankdata(x)


def _detrend(v, t):
    a = np.vstack([np.ones_like(t), t - t.mean()]).T
    coef, *_ = np.linalg.lstsq(a, v, rcond=None)
    return v - a @ coef


def _lag1(v):
    v = v - v.mean()
    den = float(v @ v)
    return float(v[:-1] @ v[1:] / den) if den > 0 and len(v) > 2 else 0.0


def correlation(x, y, t=None, method: str = "spearman") -> Dict[str, float]:
    """Coefficient, time-detrended partial coefficient, p-values (naive and n_eff), n and n_eff."""
    from scipy import stats

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.arange(len(x), dtype=float) if t is None else np.asarray(t, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(t)
    x, y, t = x[ok], y[ok], t[ok]
    n = int(ok.sum())
    nan = float("nan")
    out = {"n_samples": n, "coefficient": nan, "p_value_naive": nan, "partial_coefficient": nan,
           "p_value": nan, "n_effective": nan, "x_unique": int(len(np.unique(x))) if n else 0,
           "x_minority": int(np.min(np.unique(x, return_counts=True)[1])) if n else 0}
    if n < 4 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return out
    if method == "spearman":
        r, p = stats.spearmanr(x, y)
        xr, yr = _ranks(x), _ranks(y)
    elif method == "kendall":
        r, p = stats.kendalltau(x, y)
        xr, yr = _ranks(x), _ranks(y)
    elif method == "pearson":
        r, p = stats.pearsonr(x, y)
        xr, yr = x, y
    else:
        raise ValueError(f"unknown correlation method {method!r}")
    order = np.argsort(t)
    xd, yd = _detrend(xr[order], t[order]), _detrend(yr[order], t[order])
    if np.std(xd) == 0 or np.std(yd) == 0:
        rp = 0.0
    else:
        rp = float(np.corrcoef(xd, yd)[0, 1])
    rho = _lag1(xd) * _lag1(yd)
    neff = float(np.clip(n * (1 - rho) / (1 + rho), 3.0, n)) if rho > 0 else float(n)
    tstat = rp * np.sqrt(max(neff - 2.0, 1.0) / max(1 - rp * rp, 1e-12))
    pe = float(2 * stats.t.sf(abs(tstat), max(neff - 2.0, 1.0)))
    out.update({"coefficient": float(r), "p_value_naive": float(p), "partial_coefficient": rp, "p_value": pe,
                "n_effective": neff})
    return out


def holm(p_values) -> np.ndarray:
    """Holm-Bonferroni adjusted p-values (NaN kept)."""
    p = np.asarray(p_values, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if not len(ok):
        return out
    order = ok[np.argsort(p[ok])]
    m = len(order)
    running = 0.0
    for i, j in enumerate(order):
        running = max(running, min(1.0, (m - i) * p[j]))
        out[j] = running
    return out


def runs(mask, times_s):
    """Consecutive runs of True: list of ``(start_index, end_index, length, duration_s)``."""
    mask = np.asarray(mask, bool)
    out, start = [], None
    for i, v in enumerate(list(mask) + [False]):
        if v and start is None:
            start = i
        elif not v and start is not None:
            out.append((start, i - 1, i - start, float(times_s[i - 1] - times_s[start])))
            start = None
    return out
