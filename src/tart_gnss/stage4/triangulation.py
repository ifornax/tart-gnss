"""Least-squares intersection of station sight lines (ITRS/ECEF, metres).

Station ``i`` at ``r_i`` sees the satellite along the unit vector ``s_i``, so
``x = r_i + rho_i s_i`` with range ``rho_i > 0``. With the projector
``P_i = I - s_i s_i^T`` the perpendicular offset of a point ``x`` from line
``i`` is ``d_i = P_i (x - r_i)``.

Unweighted (baseline) estimator
    minimise ``sum_i |P_i (x - r_i)|^2``, i.e. solve ``A x = b`` with
    ``A = sum_i P_i`` and ``b = sum_i P_i r_i``. Closed form: no initial
    guess and no TLE.

Weighted (production) estimator
    In the tangent basis ``E_i = [e1_i e2_i]`` of ``s_i`` the transverse
    offset is ``q_i = E_i^T (x - r_i)``. A direction error ``delta s_i``
    displaces the line by ``-rho_i delta s_i`` at the satellite, so
    ``Cov(q_i) = rho_i^2 C_i`` with ``C_i = E_i^T Sigma_theta,i E_i`` (2x2,
    rad^2). Minimise ``chi2 = sum_i q_i^T W_i q_i`` with
    ``W_i = (rho_i^2 C_i)^-1``: solve
    ``(sum_i E_i W_i E_i^T) x = sum_i E_i W_i E_i^T r_i``, update the ranges
    ``rho_i = s_i^T (x - r_i)`` and iterate (starting from the unweighted
    solution) until the position changes by less than ``tol_m``.

Position covariance (local linearisation in the direction errors)
    ``Sigma_x = N^-1 (sum_i E_i W_i (rho_i^2 C_i) W_i E_i^T) N^-1`` with
    ``N = sum_i E_i W_i E_i^T``. For the weighted estimator this reduces to
    ``N^-1``; for the unweighted one (``W_i = I``) it is the sandwich form.
    Valid while the angular errors are small compared with the crossing
    angles between the lines (see ``docs/stage4_triangulation.md``).

Point estimates are computed with :func:`numpy.linalg.solve`; when the normal
matrix is numerically rank deficient or ill-conditioned, an SVD least-squares
solution (:func:`numpy.linalg.lstsq`) is used instead and recorded in
``TriangulationResult.solver``. ``A`` is never explicitly inverted for the
point estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from tart_gnss.stage4.errors import CovarianceError, InsufficientStationsError, TriangulationError
from tart_gnss.stage4.geometry import pairwise_angles_deg, tangent_basis, tangent_covariance

#: Above this condition number the point estimate falls back to SVD least squares.
NUMERICAL_FALLBACK_CONDITION = 1e12
#: Relative eigenvalue below which a position covariance is treated as singular.
COVARIANCE_RCOND = 1e-14


@dataclass
class TriangulationResult:
    """Outcome of one sight-line intersection.

    ``chi2``/``reduced_chi2``/``normalised_residuals`` need directional
    covariances (``None``/NaN otherwise); ``dof = 2 N - 3``.
    """

    position_ecef_m: np.ndarray
    covariance_ecef_m2: np.ndarray
    station_ranges_m: np.ndarray
    line_residuals_m: np.ndarray
    normalised_residuals: np.ndarray
    condition_number: float
    matrix_rank: int
    min_crossing_angle_deg: float
    max_crossing_angle_deg: float
    chi2: Optional[float]
    reduced_chi2: Optional[float]
    dof: int
    method: str
    solver: str
    n_iterations: int
    converged: bool
    covariance_ok: bool
    success: bool
    time_utc: Optional[object] = None
    quality_flags: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def sigma_ecef_m(self) -> np.ndarray:
        """1-sigma uncertainties of x, y, z (NaN when the covariance is unavailable)."""
        return np.sqrt(np.clip(np.diag(self.covariance_ecef_m2), 0.0, None)) if self.covariance_ok \
            else np.full(3, np.nan)


def _validate(station_positions_ecef_m, sight_vectors_ecef) -> Tuple[np.ndarray, np.ndarray]:
    r = np.asarray(station_positions_ecef_m, dtype=float)
    s = np.asarray(sight_vectors_ecef, dtype=float)
    if r.ndim != 2 or r.shape[1] != 3 or s.shape != r.shape:
        raise TriangulationError("station positions and sight vectors must both have shape (N, 3)")
    if len(r) < 2:
        raise InsufficientStationsError(f"{len(r)} sight line(s): at least 2 are needed to intersect")
    if not (np.all(np.isfinite(r)) and np.all(np.isfinite(s))):
        raise TriangulationError("non-finite station position or sight vector")
    norms = np.linalg.norm(s, axis=1)
    if np.any(norms < 1e-12):
        raise TriangulationError("zero-length sight vector")
    if np.any(np.abs(norms - 1.0) > 1e-6):
        raise TriangulationError("sight vectors must be unit vectors")
    return r, s / norms[:, None]


def projectors(sight_vectors: np.ndarray) -> np.ndarray:
    """``P_i = I - s_i s_i^T`` for each row of ``sight_vectors`` (N x 3 x 3)."""
    s = np.asarray(sight_vectors, dtype=float)
    return np.eye(3)[None] - s[:, :, None] * s[:, None, :]


def geometry_metrics(sight_vectors: np.ndarray) -> Tuple[float, int, float, float]:
    """Condition number and rank of ``A = sum P_i`` and min/max pairwise crossing angle (deg)."""
    a = projectors(sight_vectors).sum(axis=0)
    sv = np.linalg.svd(a, compute_uv=False)
    cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf")
    rank = int(np.sum(sv > sv[0] * 1e-12))
    ang = pairwise_angles_deg(sight_vectors)
    return cond, rank, float(ang.min()), float(ang.max())


def _solve(n: np.ndarray, rhs: np.ndarray) -> Tuple[np.ndarray, str]:
    """Solve ``n x = rhs`` with ``solve``; SVD least squares if ill-conditioned."""
    sv = np.linalg.svd(n, compute_uv=False)
    if sv[-1] > 0 and sv[0] / sv[-1] < NUMERICAL_FALLBACK_CONDITION:
        try:
            return np.linalg.solve(n, rhs), "solve"
        except np.linalg.LinAlgError:
            pass
    x, *_ = np.linalg.lstsq(n, rhs, rcond=None)
    return x, "lstsq"


def _covariance(n: np.ndarray, meat: Optional[np.ndarray]) -> Tuple[np.ndarray, bool]:
    """``N^-1 meat N^-1`` (or ``N^-1``); returns ``(cov, ok)`` with ``ok`` false if singular/non-PD."""
    sv = np.linalg.svd(n, compute_uv=False)
    if not (np.all(np.isfinite(n)) and sv[-1] > sv[0] * COVARIANCE_RCOND):
        return np.full((3, 3), np.nan), False
    ninv = np.linalg.inv(n)
    cov = ninv if meat is None else ninv @ meat @ ninv
    cov = 0.5 * (cov + cov.T)
    ok = bool(np.all(np.isfinite(cov)) and np.min(np.linalg.eigvalsh(cov)) > 0.0)
    return cov, ok


def _tangent_terms(s: np.ndarray, angular_covariances: Optional[np.ndarray]):
    """Tangent bases ``E_i`` (N x 3 x 2) and angular covariances ``C_i`` (N x 2 x 2, or None)."""
    es, cs = [], []
    for i in range(len(s)):
        if angular_covariances is None:
            es.append(np.column_stack(tangent_basis(s[i])))
        else:
            e, c = tangent_covariance(s[i], angular_covariances[i])
            if not np.all(np.isfinite(c)) or np.min(np.linalg.eigvalsh(c)) <= 0.0:
                raise CovarianceError(f"sight line {i}: angular covariance is not positive definite")
            es.append(e)
            cs.append(c)
    return np.array(es), (np.array(cs) if cs else None)


def _finish(r, s, x, n_norm, meat, es, cs, method, solver, n_iter, converged, extra_notes=()) -> TriangulationResult:
    rho = np.einsum("ij,ij->i", s, x[None] - r)
    q = np.einsum("nij,ni->nj", es, x[None] - r)
    resid = np.linalg.norm(q, axis=1)
    cond, rank, amin, amax = geometry_metrics(s)
    dof = 2 * len(r) - 3
    chi2 = red = None
    norm_res = np.full(len(r), np.nan)
    if cs is not None:
        rho_w = np.maximum(np.abs(rho), 1.0)
        norm_res = np.array([np.sqrt(q[i] @ np.linalg.solve(rho_w[i] ** 2 * cs[i], q[i])) for i in range(len(r))])
        chi2 = float(np.sum(norm_res ** 2))
        red = chi2 / dof if dof > 0 else None
    cov, cov_ok = _covariance(n_norm, meat) if (meat is not None or method == "weighted") else \
        (np.full((3, 3), np.nan), False)
    success = bool(np.all(np.isfinite(x)) and rank == 3 and converged)
    return TriangulationResult(
        position_ecef_m=x, covariance_ecef_m2=cov, station_ranges_m=rho, line_residuals_m=resid,
        normalised_residuals=norm_res, condition_number=cond, matrix_rank=rank,
        min_crossing_angle_deg=amin, max_crossing_angle_deg=amax, chi2=chi2, reduced_chi2=red, dof=dof,
        method=method, solver=solver, n_iterations=n_iter, converged=converged, covariance_ok=cov_ok,
        success=success, notes=list(extra_notes),
    )


def triangulate_lines(
    station_positions_ecef_m: np.ndarray,
    sight_vectors_ecef: np.ndarray,
    angular_covariances: Optional[np.ndarray] = None,
) -> TriangulationResult:
    """Estimate the 3-D closest-intersection point of multiple sight lines (unweighted).

    Parameters
    ----------
    station_positions_ecef_m, sight_vectors_ecef
        ``(N, 3)`` station positions (m) and unit sight vectors, ECEF.
    angular_covariances
        Optional ``(N, 3, 3)`` ECEF tangent covariances of the sight vectors
        (rad^2). Not used for the point estimate; when given, the position
        covariance (sandwich form) and ``chi2`` are also returned.

    Raises
    ------
    InsufficientStationsError
        Fewer than two lines.
    TriangulationError
        Malformed or non-finite inputs.
    CovarianceError
        A supplied angular covariance is not positive definite.
    """
    r, s = _validate(station_positions_ecef_m, sight_vectors_ecef)
    p = projectors(s)
    a, b = p.sum(axis=0), np.einsum("nij,nj->i", p, r)
    x, solver = _solve(a, b)
    es, cs = _tangent_terms(s, angular_covariances)
    meat = None
    if cs is not None:
        rho = np.maximum(np.abs(np.einsum("ij,ij->i", s, x[None] - r)), 1.0)
        meat = np.einsum("nia,nab,njb->ij", es, rho[:, None, None] ** 2 * cs, es)
    return _finish(r, s, x, a, meat, es, cs, "unweighted", solver, 1, True)


def triangulate_weighted(
    station_positions_ecef_m: np.ndarray,
    sight_vectors_ecef: np.ndarray,
    angular_covariances: np.ndarray,
    max_iterations: int = 50,
    tol_m: float = 1e-3,
) -> TriangulationResult:
    """Estimate the 3-D position using directional uncertainty weights.

    Iteratively reweighted: start from :func:`triangulate_lines`, weight each
    line by ``(rho_i^2 C_i)^-1`` with the current ranges, re-solve, and repeat
    until the position moves by less than ``tol_m``. Ranges are floored at
    1 m (in absolute value) when forming weights; a negative range is reported
    in the result, never hidden.

    Parameters
    ----------
    angular_covariances
        ``(N, 3, 3)`` ECEF tangent covariances of the sight vectors (rad^2),
        e.g. from :func:`tart_gnss.stage4.geometry.lm_covariance_to_ecef`.

    Raises
    ------
    CovarianceError
        A directional covariance is singular or not positive definite.
    """
    r, s = _validate(station_positions_ecef_m, sight_vectors_ecef)
    if angular_covariances is None:
        raise CovarianceError("weighted triangulation needs angular covariances")
    cov3 = np.asarray(angular_covariances, dtype=float)
    if cov3.shape != (len(r), 3, 3):
        raise CovarianceError("angular covariances must have shape (N, 3, 3)")
    es, cs = _tangent_terms(s, cov3)
    proj = np.einsum("nia,nja->nij", es, es)
    x, solver = _solve(proj.sum(axis=0), np.einsum("nij,nj->i", proj, r))
    solvers = {solver}
    converged = False
    n_iter = 0
    n_mat = None
    notes = []
    for n_iter in range(1, max_iterations + 1):
        rho = np.maximum(np.abs(np.einsum("ij,ij->i", s, x[None] - r)), 1.0)
        w = np.linalg.inv(rho[:, None, None] ** 2 * cs)
        ew = np.einsum("nia,nab,njb->nij", es, w, es)
        n_mat = ew.sum(axis=0)
        # Scale to O(1) for the solve; the covariance uses the unscaled matrix.
        scale = 1.0 / np.trace(n_mat)
        x_new, solver = _solve(n_mat * scale, np.einsum("nij,nj->i", ew, r) * scale)
        solvers.add(solver)
        step = float(np.linalg.norm(x_new - x))
        x = x_new
        if not np.all(np.isfinite(x)):
            break
        if step < tol_m:
            converged = True
            break
    if not converged:
        notes.append(f"weighted iteration did not converge in {max_iterations} iterations")
    rho = np.maximum(np.abs(np.einsum("ij,ij->i", s, x[None] - r)), 1.0)
    w = np.linalg.inv(rho[:, None, None] ** 2 * cs)
    n_mat = np.einsum("nia,nab,njb->ij", es, w, es)
    return _finish(r, s, x, n_mat, None, es, cs, "weighted", "lstsq" if "lstsq" in solvers else "solve",
                   n_iter, converged, notes)
