"""Time alignment of the stage outputs (never by row number).

Epoch table (one row per Stage-6 row)
    Stage-6 residual row -> Stage-5 trajectory row (same ``time_utc``) ->
    Stage-5 innovation row (same time) -> Stage-4 position row (same time;
    Stage-5 epochs are Stage-4 epochs) -> TLE state (``tle_trajectory``).
    Matches use the nearest time within ``tolerance_s`` (default 0.01 s); a
    missing match is recorded in ``link_stage5`` / ``link_stage4``.

Station table (one row per Stage-4 station contribution)
    keyed by (reference ``time_utc``, ``station_id``) and linked to its
    Stage-3 direction row by (``station_id``, ``station_time_utc``) and to the
    Stage-2 native file covering the station time (calibration).

Derived per epoch: the Stage-4 position minus the TLE position (ECEF and
RTN) and its squared Mahalanobis distance with the Stage-4 covariance,
``stage4_d2_tle``. Stage-4 positions are epoch-by-epoch estimates, so unlike
the smoothed Stage-5/6 quantities these are close to independent between
epochs; they are what persistence is measured on.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from tart_gnss.stage4.timing import time_seconds
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices
from tart_gnss.stage7.loader import DiagnosticInputs, calibration_at, stage2_calibration_table

DEFAULT_TOLERANCE_S = 0.01
NUMERIC_EPOCH_COLUMNS = (
    "position_residual_3d_m", "dr_radial_m", "dr_along_track_m", "dr_cross_track_m", "position_mahalanobis_obs_cov",
    "position_mahalanobis_total_cov", "tle_age_days", "stage5_nis", "stage4_condition_number",
    "stage4_min_crossing_angle_deg", "stage4_max_crossing_angle_deg", "stage4_n_stations", "stage4_reduced_chi2",
)


def nearest(sorted_t: np.ndarray, t: float, tol: float) -> Optional[int]:
    """Index of the nearest element of ``sorted_t`` to ``t`` within ``tol`` (None if none)."""
    if len(sorted_t) == 0 or not np.isfinite(t):
        return None
    i = int(np.searchsorted(sorted_t, t))
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(sorted_t) and abs(sorted_t[j] - t) <= tol and (
                best is None or abs(sorted_t[j] - t) < abs(sorted_t[best] - t)):
            best = j
    return best


def _index(df, col="time_utc"):
    if df is None or len(df) == 0:
        return np.zeros(0), np.zeros(0, int)
    t = time_seconds(df[col])
    order = np.argsort(t, kind="stable")
    return t[order], order


def _cov3(row, prefix="cov_"):
    xx, xy, xz = (float(row.get(f"{prefix}xx", np.nan)), float(row.get(f"{prefix}xy", np.nan)),
                  float(row.get(f"{prefix}xz", np.nan)))
    yy, yz, zz = float(row.get(f"{prefix}yy", np.nan)), float(row.get(f"{prefix}yz", np.nan)), float(
        row.get(f"{prefix}zz", np.nan))
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def d2(delta: np.ndarray, cov: np.ndarray) -> float:
    """Squared Mahalanobis distance (NaN if the covariance is unusable)."""
    from tart_gnss.stage6.residuals import mahalanobis2

    return mahalanobis2(np.asarray(delta, float), np.asarray(cov, float))


def build_epoch_table(inp: DiagnosticInputs, tolerance_s: float = DEFAULT_TOLERANCE_S):
    """One row per Stage-6 epoch with linked Stage-5, Stage-4 and TLE quantities."""
    import pandas as pd

    r6 = inp.residuals.reset_index(drop=True)
    t6 = time_seconds(r6["time_utc"])
    t5, o5 = _index(inp.trajectory)
    ti, oi = _index(inp.innovations)
    t4, o4 = _index(inp.positions)
    tt, ot = _index(inp.tle_trajectory)
    rows = []
    for k in range(len(r6)):
        a = r6.iloc[k]
        t = float(t6[k])
        rec: Dict[str, object] = {
            "time_utc": a["time_utc"], "time_s": t, "norad_id": inp.norad_id,
            "stage6_valid": bool(a.get("valid", False)),
            "position_residual_3d_m": a.get("position_residual_3d_m"),
            "dr_radial_m": a.get("dr_radial_m"), "dr_along_track_m": a.get("dr_along_track_m"),
            "dr_cross_track_m": a.get("dr_cross_track_m"),
            "position_mahalanobis_obs_cov": a.get("position_mahalanobis_obs_cov"),
            "position_mahalanobis_total_cov": a.get("position_mahalanobis_total_cov"),
            "tle_age_days": a.get("tle_age_days"), "precision_validation": bool(a.get("precision_validation", False)),
            "stage6_flags": list(a.get("quality_flags") or []),
        }
        j5 = nearest(t5, t, tolerance_s)
        rec["link_stage5"] = j5 is not None
        if j5 is not None:
            b = inp.trajectory.iloc[o5[j5]]
            rec.update({"stage5_nis": b.get("nis"), "stage5_flags": list(b.get("quality_flags") or []),
                        "stage5_valid": bool(b.get("valid", False))})
            rec["stage5_state_gcrs"] = np.array([b.get(c) for c in (
                "x_gcrs_smoothed_m", "y_gcrs_smoothed_m", "z_gcrs_smoothed_m", "vx_gcrs_smoothed_m_s",
                "vy_gcrs_smoothed_m_s", "vz_gcrs_smoothed_m_s")], dtype=float)
        else:
            rec.update({"stage5_nis": np.nan, "stage5_flags": [], "stage5_valid": False,
                        "stage5_state_gcrs": np.full(6, np.nan)})
        ji = nearest(ti, t, tolerance_s)
        if ji is not None and "nis_first_pass" in inp.innovations.columns:
            rec["stage5_nis"] = inp.innovations.iloc[oi[ji]].get("nis_first_pass", rec["stage5_nis"])
        j4 = nearest(t4, t, tolerance_s)
        rec["link_stage4"] = j4 is not None
        if j4 is not None:
            c = inp.positions.iloc[o4[j4]]
            rec.update({"stage4_condition_number": c.get("condition_number"),
                        "stage4_min_crossing_angle_deg": c.get("min_crossing_angle_deg"),
                        "stage4_max_crossing_angle_deg": c.get("max_crossing_angle_deg"),
                        "stage4_n_stations": c.get("n_stations"), "stage4_valid": bool(c.get("valid", False)),
                        "stage4_flags": list(c.get("quality_flags") or []),
                        "stage4_reduced_chi2": c.get("reduced_chi2"),
                        "stage4_position_ecef": np.array([c.get("x_ecef_m"), c.get("y_ecef_m"), c.get("z_ecef_m")],
                                                         dtype=float),
                        "stage4_cov_ecef": _cov3(c)})
        else:
            rec.update({"stage4_condition_number": np.nan, "stage4_min_crossing_angle_deg": np.nan,
                        "stage4_max_crossing_angle_deg": np.nan, "stage4_n_stations": 0, "stage4_valid": False,
                        "stage4_flags": [], "stage4_reduced_chi2": np.nan, "stage4_position_ecef": np.full(3, np.nan),
                        "stage4_cov_ecef": np.full((3, 3), np.nan)})
        jt = nearest(tt, t, tolerance_s)
        if jt is not None:
            tr = inp.tle_trajectory.iloc[ot[jt]]
            rec["tle_ecef"] = np.array([tr["x_ecef_m"], tr["y_ecef_m"], tr["z_ecef_m"]], dtype=float)
            rec["rtn_basis_gcrs"] = np.array([[tr[f"{v}_hat_{c}"] for c in "xyz"] for v in ("r", "t", "n")],
                                             dtype=float) if "r_hat_x" in tr else np.full((3, 3), np.nan)
        else:
            rec["tle_ecef"] = np.full(3, np.nan)
            rec["rtn_basis_gcrs"] = np.full((3, 3), np.nan)
        rows.append(rec)
    df = pd.DataFrame(rows)
    for col in NUMERIC_EPOCH_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype(float)
    if len(df):
        mats = itrs_to_gcrs_matrices(df["time_s"].to_numpy(float))
        df["ecef_to_gcrs"] = list(mats)
        d4, d4r, s4 = [], [], []
        for _, r in df.iterrows():
            d = r["stage4_position_ecef"] - r["tle_ecef"]
            d4.append(d2(d, r["stage4_cov_ecef"]))
            q = r["rtn_basis_gcrs"]
            d4r.append(q @ (r["ecef_to_gcrs"] @ d) if np.all(np.isfinite(q)) else np.full(3, np.nan))
            c = r["stage4_cov_ecef"]
            s4.append(float(np.sqrt(np.trace(c))) if np.all(np.isfinite(c)) else np.nan)
        df["stage4_d2_tle"] = d4
        df["stage4_residual_3d_m"] = [float(np.linalg.norm(r["stage4_position_ecef"] - r["tle_ecef"]))
                                      for _, r in df.iterrows()]
        df["stage4_dr_rtn_m"] = d4r
        df["stage4_sigma_3d_m"] = s4
    return df


def build_station_table(inp: DiagnosticInputs, epochs, tolerance_s: float = DEFAULT_TOLERANCE_S):
    """One row per Stage-4 station contribution linked to Stage 3 and Stage 2."""
    import pandas as pd

    con = inp.contributions
    if con is None or len(con) == 0:
        return pd.DataFrame()
    cal = stage2_calibration_table(inp.stage2)
    te = epochs["time_s"].to_numpy(float) if len(epochs) else np.zeros(0)
    order_e = np.argsort(te)
    te_sorted = te[order_e]
    d3 = inp.directions
    idx3: Dict[str, tuple] = {}
    if d3 is not None and len(d3):
        for sid, g in d3.groupby("station_id"):
            t3 = time_seconds(g["time_utc"])
            o = np.argsort(t3)
            idx3[str(sid)] = (t3[o], g.index.to_numpy()[o])
    ref_col = "reference_time_utc" if "reference_time_utc" in con.columns else "time_utc"
    t_ref = time_seconds(con[ref_col])
    t_st = time_seconds(con["station_time_utc"])
    rows = []
    for k in range(len(con)):
        c = con.iloc[k]
        sid = str(c["station_id"])
        je = nearest(te_sorted, float(t_ref[k]), tolerance_s)
        rec: Dict[str, object] = {
            "time_utc": c[ref_col], "time_s": float(t_ref[k]), "station_id": sid, "norad_id": inp.norad_id,
            "station_time_s": float(t_st[k]), "station_time_offset_s": c.get("station_time_offset_s"),
            "line_residual_m": c.get("line_residual_m"), "normalised_line_residual": c.get("normalised_residual"),
            "sigma_theta_deg": c.get("sigma_theta_deg"), "range_m": c.get("range_m"),
            "stage3_flags": list(c.get("stage3_quality_flags") or []),
            "stage4_flags": list(c.get("stage4_station_flags") or []),
            "epoch_index": int(order_e[je]) if je is not None else -1,
            "time_mismatch": "TIME_MISMATCH" in list(c.get("stage4_station_flags") or []),
        }
        rec["calibration_warning"] = bool(inp.stations.get(sid) and inp.stations[sid].calibration_warnings)
        j3 = None
        if sid in idx3:
            j3 = nearest(idx3[sid][0], float(t_st[k]), tolerance_s)
        rec["link_stage3"] = j3 is not None
        if j3 is not None:
            d = d3.loc[idx3[sid][1][j3]]
            for col in ("l", "m", "sigma_l", "sigma_m", "cov_lm", "az_obs_deg", "el_obs_deg", "sigma_az_deg",
                        "sigma_el_deg", "cov_az_el", "snr", "reduced_chi2", "flagged_fraction",
                        "vis_image_separation_deg", "n_baselines_used"):
                rec[col] = d.get(col, np.nan)
        cr = calibration_at(cal, sid, float(t_st[k])) if cal else None
        rec["link_stage2"] = cr is not None
        rec["zero_gain_count"] = cr["zero_gain_count"] if cr else np.nan
        rec["gain_solution_index"] = cr["solution_index"] if cr else np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    for col in ("station_time_offset_s", "line_residual_m", "normalised_line_residual", "sigma_theta_deg", "range_m",
                "l", "m", "sigma_l", "sigma_m", "cov_lm", "az_obs_deg", "el_obs_deg", "sigma_az_deg", "sigma_el_deg",
                "cov_az_el", "snr", "reduced_chi2", "flagged_fraction", "vis_image_separation_deg", "n_baselines_used",
                "zero_gain_count", "gain_solution_index"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype(float)
    return df
