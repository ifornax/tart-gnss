"""``observed - TLE`` residuals: Cartesian (GCRS), RTN, angular, range and Mahalanobis.

Permanent sign convention: **observed minus TLE**, for every residual.

Measurement-conditional distances use only the Stage-5 covariance and are
conditional on the TLE state being treated as fixed::

    D2_pos   = dr^T C_r^-1 dr          (position_mahalanobis_obs_cov)
    D2_state = dx^T P^-1 dx            (state_mahalanobis_obs_cov)

Reference-aware distances add the *empirical reference dispersion* of a TLE
ensemble (never a formal orbit-determination covariance)::

    C_total = C_obs + C_TLE,empirical  (position/state_mahalanobis_total_cov)

Every distance is per epoch. Stage-5 epochs are strongly correlated, so the
distances are never summed into a chi-square over epochs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from tart_gnss.stage6.rtn import covariance_to_rtn, rtn_bases

MAX_CONDITION = 1e14


@dataclass
class StateTrajectory:
    """States ``(N, 6)`` at ``times_s`` with optional covariances ``(N, 6, 6)`` (GCRS, SI units)."""

    times_s: np.ndarray
    states: np.ndarray
    covariances: Optional[np.ndarray] = None


@dataclass
class ResidualTable:
    """Per-epoch residual arrays (all ``observed - reference``)."""

    times_s: np.ndarray
    dr: np.ndarray
    dv: np.ndarray
    bases: np.ndarray
    dr_rtn: np.ndarray
    dv_rtn: np.ndarray
    cov_r_rtn: np.ndarray
    cov_v_rtn: np.ndarray
    angular_deg: np.ndarray
    range_m: np.ndarray
    along_track_time_s: np.ndarray
    d2_pos_obs: np.ndarray
    d2_state_obs: np.ndarray
    d2_pos_total: np.ndarray
    d2_state_total: np.ndarray
    covariance_ok: np.ndarray
    notes: List[str] = field(default_factory=list)

    @property
    def position_norm_m(self) -> np.ndarray:
        """``|dr|``."""
        return np.linalg.norm(self.dr, axis=1)

    @property
    def velocity_norm_m_s(self) -> np.ndarray:
        """``|dv|``."""
        return np.linalg.norm(self.dv, axis=1)


def mahalanobis2(delta: np.ndarray, cov: np.ndarray) -> float:
    """``delta^T C^-1 delta`` by Cholesky; NaN if ``C`` is not usable (non-finite, not PD, ill-conditioned)."""
    if not (np.all(np.isfinite(delta)) and np.all(np.isfinite(cov))):
        return float("nan")
    c = 0.5 * (cov + cov.T)
    try:
        w = np.linalg.eigvalsh(c)
        if w[0] <= 0 or w[-1] / w[0] > MAX_CONDITION:
            return float("nan")
        lo = np.linalg.cholesky(c)
    except np.linalg.LinAlgError:
        return float("nan")
    z = np.linalg.solve(lo, delta)
    return float(z @ z)


def covariance_usable(cov: np.ndarray) -> bool:
    """Finite, symmetric within 1e-6 relative, positive definite, condition <= MAX_CONDITION."""
    if not np.all(np.isfinite(cov)):
        return False
    scale = float(np.max(np.abs(cov))) or 1.0
    if float(np.max(np.abs(cov - cov.T))) > 1e-6 * scale:
        return False
    w = np.linalg.eigvalsh(0.5 * (cov + cov.T))
    return bool(w[0] > 0 and w[-1] / w[0] <= MAX_CONDITION)


def angular_separation_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Geocentric angle between position vectors (atan2 form, accurate at small angles)."""
    cr = np.linalg.norm(np.cross(a, b), axis=1)
    dt = np.sum(a * b, axis=1)
    return np.degrees(np.arctan2(cr, dt))


def compute_residuals(observed: StateTrajectory, reference: StateTrajectory,
                      reference_dispersion: Optional[np.ndarray] = None) -> ResidualTable:
    """Calculate observed-minus-reference Cartesian and RTN residuals.

    ``reference`` supplies the RTN axes (its own state). ``reference_dispersion``
    (``(N, 6, 6)``, empirical) enables the total-covariance distances.
    """
    xo, xr = np.asarray(observed.states, float), np.asarray(reference.states, float)
    if xo.shape != xr.shape:
        raise ValueError("observed and reference trajectories must have the same epochs")
    n = len(xo)
    p = observed.covariances if observed.covariances is not None else np.full((n, 6, 6), np.nan)
    dx = xo - xr
    dr, dv = dx[:, :3], dx[:, 3:]
    q = rtn_bases(xr[:, :3], xr[:, 3:])
    dr_rtn = np.einsum("nij,nj->ni", q, dr)
    dv_rtn = np.einsum("nij,nj->ni", q, dv)
    cov_r = np.array([covariance_to_rtn(qk, pk[:3, :3]) for qk, pk in zip(q, p)])
    cov_v = np.array([covariance_to_rtn(qk, pk[3:, 3:]) for qk, pk in zip(q, p)])
    ok = np.array([covariance_usable(pk[:3, :3]) for pk in p])
    d2p = np.array([mahalanobis2(dr[k], p[k][:3, :3]) if ok[k] else np.nan for k in range(n)])
    d2s = np.array([mahalanobis2(dx[k], p[k]) if ok[k] else np.nan for k in range(n)])
    d2pt = np.full(n, np.nan)
    d2st = np.full(n, np.nan)
    if reference_dispersion is not None:
        for k in range(n):
            if ok[k]:
                tot = p[k] + reference_dispersion[k]
                d2pt[k] = mahalanobis2(dr[k], tot[:3, :3])
                d2st[k] = mahalanobis2(dx[k], tot)
    speed = np.linalg.norm(xr[:, 3:], axis=1)
    notes = []
    if np.any(ok & ~np.isfinite(d2s)):
        notes.append(f"{int(np.sum(ok & ~np.isfinite(d2s)))} epoch(s): 6x6 covariance too ill-conditioned for "
                     "the state distance (position distance still given)")
    return ResidualTable(
        times_s=np.asarray(observed.times_s, float), dr=dr, dv=dv, bases=q, dr_rtn=dr_rtn, dv_rtn=dv_rtn,
        cov_r_rtn=cov_r, cov_v_rtn=cov_v, angular_deg=angular_separation_deg(xo[:, :3], xr[:, :3]),
        range_m=np.linalg.norm(xo[:, :3], axis=1) - np.linalg.norm(xr[:, :3], axis=1),
        along_track_time_s=dr_rtn[:, 1] / speed, d2_pos_obs=d2p, d2_state_obs=d2s, d2_pos_total=d2pt,
        d2_state_total=d2st, covariance_ok=ok, notes=notes,
    )


def rows(res: ResidualTable) -> List[Dict[str, object]]:
    """Per-epoch residual fields (column names of the residual table)."""
    out = []
    for k in range(len(res.times_s)):
        cr, cv = res.cov_r_rtn[k], res.cov_v_rtn[k]

        def sd(c, i):
            return float(np.sqrt(c[i, i])) if np.isfinite(c[i, i]) and c[i, i] >= 0 else float("nan")

        out.append({
            "dr_x_m": res.dr[k, 0], "dr_y_m": res.dr[k, 1], "dr_z_m": res.dr[k, 2],
            "position_residual_3d_m": float(np.linalg.norm(res.dr[k])),
            "dv_x_m_s": res.dv[k, 0], "dv_y_m_s": res.dv[k, 1], "dv_z_m_s": res.dv[k, 2],
            "velocity_residual_3d_m_s": float(np.linalg.norm(res.dv[k])),
            "dr_radial_m": res.dr_rtn[k, 0], "dr_along_track_m": res.dr_rtn[k, 1],
            "dr_cross_track_m": res.dr_rtn[k, 2],
            "dv_radial_m_s": res.dv_rtn[k, 0], "dv_along_track_m_s": res.dv_rtn[k, 1],
            "dv_cross_track_m_s": res.dv_rtn[k, 2],
            "sigma_radial_obs_m": sd(cr, 0), "sigma_along_track_obs_m": sd(cr, 1),
            "sigma_cross_track_obs_m": sd(cr, 2),
            "cov_radial_along_obs_m2": cr[0, 1], "cov_radial_cross_obs_m2": cr[0, 2],
            "cov_along_cross_obs_m2": cr[1, 2],
            "sigma_dv_radial_obs_m_s": sd(cv, 0), "sigma_dv_along_track_obs_m_s": sd(cv, 1),
            "sigma_dv_cross_track_obs_m_s": sd(cv, 2),
            "geocentric_angular_residual_deg": res.angular_deg[k], "geocentric_range_residual_m": res.range_m[k],
            "equivalent_along_track_time_offset_s": res.along_track_time_s[k],
            "position_mahalanobis_obs_cov": res.d2_pos_obs[k], "state_mahalanobis_obs_cov": res.d2_state_obs[k],
            "position_mahalanobis_total_cov": res.d2_pos_total[k],
            "state_mahalanobis_total_cov": res.d2_state_total[k],
        })
    return out
