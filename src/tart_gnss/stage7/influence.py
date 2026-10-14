"""Leave-one-station-out and pairwise re-triangulation, and station direction residuals.

All solutions reuse the Stage-4 code: sight lines from the Stage-3 rows with
:func:`tart_gnss.stage4.loader.sight_line` and the frozen Stage-4 station
coordinates, intersected with :func:`tart_gnss.stage4.triangulation.triangulate_weighted`
(no triangulation mathematics is duplicated). The full solution is first
recomputed and compared with the Stage-4 table; the largest difference is
recorded as a reproduction check.

Leave-one-out (LOO), per epoch with >= 3 stations and per station ``i``
    ``shift = x_full - x_-i`` (ECEF; RTN via the Stage-6 TLE basis), and
    ``loo_mahalanobis = shift^T (C_-i - C_full)^+ shift``. For nested
    least-squares estimates with consistent covariances the covariance of
    the difference is exactly ``C_-i - C_full``; it has rank <= 2 (one
    station carries two angular constraints), so its pseudo-inverse is used
    and ``loo_dof`` records the rank (chi^2_dof under the null). If the
    matrix is not positive semidefinite only the physical shift is reported. Removing a station
    from a 3-station network leaves two lines: ``degraded_geometry`` is set
    when fewer than 3 stations remain or the condition number exceeds
    ``max_condition``. ``loo_d2_tle`` is the LOO solution's distance from the
    TLE position with ``C_-i`` (used by the reference-based station test).

Pairwise, per epoch and station pair
    the two-line weighted solution (``pair_position``), the closest-approach
    separation of the two lines ``pair_line_separation_m``, their angle, and
    the separation normalised by its predicted sigma
    ``pair_separation_z = sep / sqrt(n^T (rho_1^2 C_1 + rho_2^2 C_2) n)``
    (``n`` the common normal, ``C`` the tangent angular covariances).

Station direction residuals
    the Stage-5 smoothed state, moved to the station's own sample time
    (``r + v dt``), rotated to ECEF and seen from the station, compared with
    the measured direction: angle, ``(l, m)`` direction-cosine differences in
    the station's own East/North frame (well defined at the zenith, unlike
    azimuth), Az/El differences, and the components along / across the
    predicted apparent motion. All residuals are measured - predicted.
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Optional, Tuple

import numpy as np

from tart_gnss.stage4.errors import CovarianceError, TriangulationError
from tart_gnss.stage4.geometry import ecef_to_azel, ecef_to_enu_vector
from tart_gnss.stage4.loader import sight_line
from tart_gnss.stage4.triangulation import triangulate_lines, triangulate_weighted
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices
from tart_gnss.stage7.alignment import d2

DEFAULT_MAX_CONDITION = 1e8


def lines_for_epoch(station_rows, stations) -> List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    """``(station_id, r, s, cov3)`` for the contributing stations with a linked Stage-3 row."""
    out = []
    for _, row in station_rows.iterrows():
        sid = row["station_id"]
        if sid not in stations or not row.get("link_stage3", False):
            continue
        try:
            s, cov3, _ = sight_line(row, stations[sid])
        except (TriangulationError, CovarianceError):
            continue
        out.append((sid, stations[sid].ecef_m, s, cov3))
    return out


def solve(lines):
    """Weighted Stage-4 intersection (unweighted if the covariances are unusable)."""
    r = np.array([x[1] for x in lines])
    s = np.array([x[2] for x in lines])
    c = np.array([x[3] for x in lines])
    try:
        return triangulate_weighted(r, s, c)
    except (CovarianceError, TriangulationError):
        return triangulate_lines(r, s)


def _rtn(shift_ecef, m_e2g, q):
    if not (np.all(np.isfinite(q)) and np.all(np.isfinite(shift_ecef))):
        return np.full(3, np.nan)
    return q @ (m_e2g @ shift_ecef)


def loo_d2(shift, c_loo, c_full, rel_tol: float = 1e-6) -> Tuple[float, int]:
    """``(D^2, dof)`` of a leave-one-out shift with the covariance of the difference ``C_-i - C_full``.

    Removing one station removes at most two angular constraints, so the
    difference covariance has rank <= 2: the distance uses its pseudo-inverse
    on the eigenvectors with eigenvalue above ``rel_tol`` x the largest, and
    ``dof`` is that rank (chi^2_dof under the null). NaN if the matrix is not
    finite or has a significantly negative eigenvalue.
    """
    if not (np.all(np.isfinite(c_loo)) and np.all(np.isfinite(c_full)) and np.all(np.isfinite(shift))):
        return float("nan"), 0
    c = c_loo - c_full
    w, v = np.linalg.eigh(0.5 * (c + c.T))
    if w[-1] <= 0 or w[0] < -1e-3 * w[-1]:
        return float("nan"), 0
    keep = w > rel_tol * w[-1]
    z = v[:, keep].T @ shift
    return float(np.sum(z * z / w[keep])), int(keep.sum())


def pair_diagnostic(a, b) -> Dict[str, object]:
    """Closest approach of two sight lines and its normalised separation."""
    (sa, ra, ua, ca), (sb, rb, ub, cb) = a, b
    w = ra - rb
    aa, bb, cc = ua @ ua, ua @ ub, ub @ ub
    dd, ee = ua @ w, ub @ w
    den = aa * cc - bb * bb
    angle = float(np.degrees(np.arccos(np.clip(abs(bb), -1, 1))))
    if den <= 1e-18:
        return {"pair_line_separation_m": np.nan, "pair_geometry_angle_deg": angle, "pair_separation_z": np.nan,
                "pair_position": np.full(3, np.nan)}
    ta = (bb * ee - cc * dd) / den
    tb = (aa * ee - bb * dd) / den
    pa, pb = ra + ta * ua, rb + tb * ub
    n = np.cross(ua, ub)
    n /= np.linalg.norm(n)
    sep = float(abs((pa - pb) @ n))
    var = ta ** 2 * (n @ ca @ n) + tb ** 2 * (n @ cb @ n)
    z = sep / np.sqrt(var) if np.isfinite(var) and var > 0 else np.nan
    try:
        pos = solve([a, b]).position_ecef_m
    except (TriangulationError, CovarianceError, np.linalg.LinAlgError):
        pos = 0.5 * (pa + pb)
    return {"pair_line_separation_m": sep, "pair_geometry_angle_deg": angle, "pair_separation_z": float(z),
            "pair_position": pos}


def compute_leave_one_out(epochs, stations_tbl, stations, leave_one_out: bool = True, pairwise: bool = True,
                          max_condition: float = DEFAULT_MAX_CONDITION):
    """Recompute localisation while excluding each station in turn (and per station pair).

    Returns ``(influence_rows, pair_rows, reproduction_max_diff_m)``.
    """
    influence, pairs = [], []
    max_diff = 0.0
    if stations_tbl is None or len(stations_tbl) == 0:
        return influence, pairs, float("nan")
    for ei, ep in epochs.iterrows():
        rows = stations_tbl[stations_tbl["epoch_index"] == ei]
        lines = lines_for_epoch(rows, stations)
        if len(lines) < 2:
            continue
        full = solve(lines)
        x4 = ep["stage4_position_ecef"]
        if np.all(np.isfinite(x4)):
            max_diff = max(max_diff, float(np.linalg.norm(full.position_ecef_m - x4)))
        q = ep["rtn_basis_gcrs"]
        m = ep["ecef_to_gcrs"]
        if leave_one_out and len(lines) >= 3:
            for i, (sid, *_r) in enumerate(lines):
                sub = [x for j, x in enumerate(lines) if j != i]
                loo = solve(sub)
                shift = full.position_ecef_m - loo.position_ecef_m
                degraded = len(sub) < 3 or not (loo.condition_number <= max_condition)
                influence.append({
                    "epoch_index": ei, "time_utc": ep["time_utc"], "station_id": sid, "n_stations": len(lines),
                    "n_stations_remaining": len(sub),
                    "loo_x_ecef_m": loo.position_ecef_m[0], "loo_y_ecef_m": loo.position_ecef_m[1],
                    "loo_z_ecef_m": loo.position_ecef_m[2], "loo_shift_3d_m": float(np.linalg.norm(shift)),
                    **dict(zip(("loo_shift_radial_m", "loo_shift_along_track_m", "loo_shift_cross_track_m"),
                               _rtn(shift, m, q))),
                    **dict(zip(("loo_mahalanobis", "loo_dof"),
                               loo_d2(shift, loo.covariance_ecef_m2, full.covariance_ecef_m2))),
                    "loo_condition_number": loo.condition_number,
                    "loo_min_crossing_angle_deg": loo.min_crossing_angle_deg, "degraded_geometry": bool(degraded),
                    "loo_d2_tle": d2(loo.position_ecef_m - ep["tle_ecef"], loo.covariance_ecef_m2),
                    "full_d2_tle": d2(full.position_ecef_m - ep["tle_ecef"], full.covariance_ecef_m2),
                })
        if pairwise:
            for a, b in combinations(lines, 2):
                pd_ = pair_diagnostic(a, b)
                pairs.append({"epoch_index": ei, "time_utc": ep["time_utc"], "pair_station_1": a[0],
                              "pair_station_2": b[0], **pd_, "n_stations": len(lines)})
    return influence, pairs, max_diff


def compute_station_direction_residuals(epochs, stations_tbl, stations) -> List[Dict[str, object]]:
    """Compare each station direction with the network-derived Stage-5 trajectory."""
    out = []
    if stations_tbl is None or len(stations_tbl) == 0:
        return out
    ts = stations_tbl["station_time_s"].to_numpy(float)
    mats = itrs_to_gcrs_matrices(ts) if len(ts) else np.zeros((0, 3, 3))
    mats_next = itrs_to_gcrs_matrices(ts + 1.0) if len(ts) else np.zeros((0, 3, 3))
    for k, (_, row) in enumerate(stations_tbl.iterrows()):
        rec = {"direction_residual_angle_deg": np.nan, "direction_residual_az_deg": np.nan,
               "direction_residual_el_deg": np.nan, "direction_residual_l_deg": np.nan,
               "direction_residual_m_deg": np.nan, "direction_residual_along_motion_deg": np.nan,
               "direction_residual_across_motion_deg": np.nan, "predicted_el_deg": np.nan}
        ei = int(row["epoch_index"])
        sid = row["station_id"]
        if ei < 0 or sid not in stations or not row.get("link_stage3", False):
            out.append(rec)
            continue
        state = epochs.iloc[ei]["stage5_state_gcrs"]
        if not np.all(np.isfinite(state)):
            out.append(rec)
            continue
        st = stations[sid]
        try:
            s_meas, _, _ = sight_line(row, st)
        except (TriangulationError, CovarianceError):
            out.append(rec)
            continue
        dt = float(row["station_time_s"] - epochs.iloc[ei]["time_s"])

        def unit_at(extra, mat):
            r = state[:3] + state[3:] * (dt + extra)
            d = mat.T @ r - st.ecef_m
            return d / np.linalg.norm(d)

        u = unit_at(0.0, mats[k])
        u1 = unit_at(1.0, mats_next[k])
        enu_p = ecef_to_enu_vector(u, st.latitude_deg, st.longitude_deg)
        enu_p1 = ecef_to_enu_vector(u1, st.latitude_deg, st.longitude_deg)
        lm_meas = np.array([float(row["l"]), float(row["m"])])
        if not np.all(np.isfinite(lm_meas)):
            enu_m = ecef_to_enu_vector(s_meas, st.latitude_deg, st.longitude_deg)
            lm_meas = enu_m[:2]
        res = lm_meas - enu_p[:2]  # measured - predicted direction cosines (l = East, m = North)
        motion = enu_p1[:2] - enu_p[:2]
        mn = np.linalg.norm(motion)
        az_p, el_p, _ = ecef_to_azel(st.ecef_m + u, st.ecef_m, st.latitude_deg, st.longitude_deg)
        daz = (float(row.get("az_obs_deg", np.nan)) - az_p + 180.0) % 360.0 - 180.0
        perp = np.array([-motion[1], motion[0]]) / mn if mn > 0 else np.full(2, np.nan)
        rec.update({
            "direction_residual_angle_deg": float(np.degrees(np.arccos(np.clip(s_meas @ u, -1, 1)))),
            "direction_residual_az_deg": daz,
            "direction_residual_el_deg": float(row.get("el_obs_deg", np.nan)) - el_p,
            "direction_residual_l_deg": float(np.degrees(res[0])), "direction_residual_m_deg": float(
                np.degrees(res[1])),
            "direction_residual_along_motion_deg": float(np.degrees(res @ motion / mn)) if mn > 0 else np.nan,
            "direction_residual_across_motion_deg": float(np.degrees(res @ perp)) if mn > 0 else np.nan,
            "predicted_el_deg": el_p,
        })
        out.append(rec)
    return out


def station_bias(stations_tbl) -> Dict[str, Dict[str, float]]:
    """Robust per-station direction bias (median, MAD scale) and descriptive trends vs elevation and time."""
    from tart_gnss.stage6.statistics import MAD_TO_SIGMA, trend

    out = {}
    if stations_tbl is None or len(stations_tbl) == 0:
        return out
    for sid, g in stations_tbl.groupby("station_id"):
        rec = {"n": int(g["direction_residual_angle_deg"].notna().sum())}
        for c in ("direction_residual_az_deg", "direction_residual_el_deg", "direction_residual_angle_deg",
                  "direction_residual_along_motion_deg"):
            v = g[c].to_numpy(float)
            v = v[np.isfinite(v)]
            med = float(np.median(v)) if v.size else float("nan")
            rec[f"median_{c}"] = med
            rec[f"mad_{c}"] = float(MAD_TO_SIGMA * np.median(np.abs(v - med))) if v.size else float("nan")
        for c, x in (("time", "station_time_s"), ("elevation", "predicted_el_deg")):
            tr = trend(g[x].to_numpy(float), g["direction_residual_angle_deg"].to_numpy(float))
            rec[f"angle_slope_vs_{c}"] = tr["slope"]
        out[str(sid)] = rec
    return out


def tangent_vector_deg(row) -> Optional[np.ndarray]:
    """``(l, m)`` direction residual vector (deg) of a station row (None if unavailable)."""
    v = np.array([row.get("direction_residual_l_deg", np.nan), row.get("direction_residual_m_deg", np.nan)], float)
    return v if np.all(np.isfinite(v)) else None
