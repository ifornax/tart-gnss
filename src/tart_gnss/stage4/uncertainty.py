"""Position uncertainty: summaries of the linear covariance, bootstrap and Monte Carlo checks.

The linear covariance (:mod:`tart_gnss.stage4.triangulation`) assumes the
angular errors are small compared with the crossing angles between the sight
lines. The parametric bootstrap here makes no such assumption: each station
direction is redrawn from its Stage-3 covariance (Gaussian in the tangent
plane) and the triangulation is repeated. :func:`monte_carlo_check` compares
the empirical scatter of many such realisations with the predicted covariance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import numpy as np

from tart_gnss.stage4.errors import Stage4Error
from tart_gnss.stage4.geometry import tangent_covariance
from tart_gnss.stage4.triangulation import TriangulationResult, triangulate_lines, triangulate_weighted


def triangulator(method: str) -> Callable[..., TriangulationResult]:
    """The triangulation function for ``method`` ('weighted' or 'unweighted')."""
    if method == "weighted":
        return triangulate_weighted
    if method == "unweighted":
        return triangulate_lines
    raise ValueError(f"unknown method {method!r}")


def covariance_components(cov: np.ndarray) -> Dict[str, float]:
    """``cov_xx`` ... ``cov_zz`` and ``sigma_x/y/z_m`` from a 3x3 covariance."""
    c = np.asarray(cov, dtype=float)
    names = ("x", "y", "z")
    out = {f"cov_{names[i]}{names[j]}": float(c[i, j]) for i in range(3) for j in range(i, 3)}
    for i, n in enumerate(names):
        out[f"sigma_{n}_m"] = float(np.sqrt(c[i, i])) if np.isfinite(c[i, i]) and c[i, i] >= 0 else float("nan")
    return out


def perturb_directions(sight_vectors: np.ndarray, angular_covariances: np.ndarray,
                       rng: np.random.Generator) -> np.ndarray:
    """Redraw each unit vector from its tangent-plane Gaussian covariance."""
    out = np.empty_like(np.asarray(sight_vectors, dtype=float))
    for i, s in enumerate(sight_vectors):
        e, c = tangent_covariance(s, angular_covariances[i])
        v = s + e @ rng.multivariate_normal(np.zeros(2), c)
        out[i] = v / np.linalg.norm(v)
    return out


@dataclass
class BootstrapResult:
    """Percentiles of bootstrap positions (m) and stability diagnostics."""

    n_requested: int
    n_success: int
    percentiles: Dict[str, float] = field(default_factory=dict)
    std_m: np.ndarray = field(default_factory=lambda: np.full(3, np.nan))
    failure_fraction: float = float("nan")
    max_sigma_ratio: float = float("nan")
    max_median_shift_sigma: float = float("nan")


def bootstrap_position(
    station_positions_ecef_m: np.ndarray,
    sight_vectors_ecef: np.ndarray,
    angular_covariances: np.ndarray,
    n: int,
    rng: np.random.Generator,
    method: str = "weighted",
    linear_sigma_m: Optional[np.ndarray] = None,
    point_estimate_m: Optional[np.ndarray] = None,
) -> BootstrapResult:
    """Parametric bootstrap of the triangulated position.

    Returns the 16th/50th/84th percentiles of x, y, z, the fraction of
    resamples that failed (exception, non-finite, or a non-positive range),
    the largest ratio between the bootstrap half-width ``(p84 - p16) / 2``
    and the linear ``sigma`` (either way round), and the largest shift of the
    bootstrap median from ``point_estimate_m`` in units of the linear sigma
    (a sign of estimator bias, i.e. of the non-linear regime).
    """
    fn = triangulator(method)
    pts = []
    for _ in range(n):
        s = perturb_directions(sight_vectors_ecef, angular_covariances, rng)
        try:
            res = fn(station_positions_ecef_m, s, angular_covariances)
        except Stage4Error:
            continue
        if res.success and np.all(res.station_ranges_m > 0):
            pts.append(res.position_ecef_m)
    out = BootstrapResult(n_requested=n, n_success=len(pts))
    out.failure_fraction = 1.0 - len(pts) / n if n else float("nan")
    if len(pts) >= 2:
        p = np.percentile(np.array(pts), [16, 50, 84], axis=0)
        for i, ax in enumerate("xyz"):
            for j, q in enumerate((16, 50, 84)):
                out.percentiles[f"{ax}_p{q}"] = float(p[j, i])
        out.std_m = np.std(np.array(pts), axis=0, ddof=1)
        if linear_sigma_m is not None and np.all(np.isfinite(linear_sigma_m)) and np.all(linear_sigma_m > 0):
            half = 0.5 * (p[2] - p[0])
            ratio = half / linear_sigma_m
            out.max_sigma_ratio = float(np.max(np.maximum(ratio, 1.0 / np.maximum(ratio, 1e-300))))
            if point_estimate_m is not None:
                out.max_median_shift_sigma = float(np.max(np.abs(p[1] - point_estimate_m) / linear_sigma_m))
    return out


@dataclass
class MonteCarloCheck:
    """Empirical scatter vs predicted covariance over noisy realisations."""

    n: int
    predicted_sigma_m: np.ndarray
    empirical_sigma_m: np.ndarray
    mean_error_m: np.ndarray
    mean_mahalanobis2: float
    los_axis_predicted_sigma_m: float
    los_axis_empirical_sigma_m: float

    @property
    def sigma_ratio(self) -> np.ndarray:
        """Empirical / predicted sigma for x, y, z."""
        return self.empirical_sigma_m / self.predicted_sigma_m


def monte_carlo_check(
    truth_ecef_m: np.ndarray,
    station_positions_ecef_m: np.ndarray,
    angular_covariances: np.ndarray,
    n: int,
    rng: np.random.Generator,
    method: str = "weighted",
) -> MonteCarloCheck:
    """Compare empirical position scatter with the predicted covariance.

    Exact sight vectors are computed from ``truth_ecef_m``, perturbed ``n``
    times by their angular covariances and triangulated. The predicted
    covariance is the mean of the per-realisation covariances. Also returns
    the mean squared Mahalanobis distance (expected 3) and the scatter along
    the mean line of sight (the least constrained axis).
    """
    r = np.asarray(station_positions_ecef_m, dtype=float)
    d = np.asarray(truth_ecef_m, dtype=float)[None] - r
    s0 = d / np.linalg.norm(d, axis=1)[:, None]
    fn = triangulator(method)
    errs, covs, m2 = [], [], []
    for _ in range(n):
        res = fn(r, perturb_directions(s0, angular_covariances, rng), angular_covariances)
        e = res.position_ecef_m - truth_ecef_m
        errs.append(e)
        covs.append(res.covariance_ecef_m2)
        m2.append(float(e @ np.linalg.solve(res.covariance_ecef_m2, e)))
    errs, cov = np.array(errs), np.mean(covs, axis=0)
    los = s0.mean(axis=0)
    los /= np.linalg.norm(los)
    return MonteCarloCheck(
        n=n, predicted_sigma_m=np.sqrt(np.diag(cov)), empirical_sigma_m=np.std(errs, axis=0, ddof=1),
        mean_error_m=errs.mean(axis=0), mean_mahalanobis2=float(np.mean(m2)),
        los_axis_predicted_sigma_m=float(np.sqrt(los @ cov @ los)),
        los_axis_empirical_sigma_m=float(np.std(errs @ los, ddof=1)),
    )
