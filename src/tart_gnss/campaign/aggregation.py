"""Build the run, epoch and station population tables from loaded Stage-7 runs.

Everything is read from completed products: the Stage-7 diagnosis and
station-influence tables and summary, and (through the recorded provenance)
the Stage-6 residuals and summary, Stage-5 trajectory and manifest, Stage-4
positions and station contributions, and the Stage-3 manifest. Rows are
joined by time (nearest within 0.01 s), never by row number.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np

from tart_gnss.campaign.loader import CampaignRun
from tart_gnss.campaign.schema import EPOCH_COLUMNS, RUN_COLUMNS, STATION_COLUMNS

TOL_S = 0.01


@dataclass
class CampaignTables:
    """Run, epoch and station tables (pandas DataFrames)."""

    runs: object
    epochs: object
    stations: object


def _secs(series) -> np.ndarray:
    import pandas as pd

    t = pd.to_datetime(series, utc=True)
    return ((t - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds()).to_numpy(float)


def _lookup(df, t_query: np.ndarray, cols: List[str]) -> Dict[str, np.ndarray]:
    """Values of ``cols`` at the rows of ``df`` nearest each query time (NaN beyond tolerance)."""
    out = {c: np.full(len(t_query), np.nan) for c in cols}
    if df is None or not len(df) or not len(t_query):
        return out
    t = _secs(df["time_utc"])
    order = np.argsort(t)
    ts = t[order]
    idx = np.clip(np.searchsorted(ts, t_query), 1, len(ts) - 1) if len(ts) > 1 else np.zeros(len(t_query), int)
    best = idx.copy()
    if len(ts) > 1:
        left = np.abs(ts[idx - 1] - t_query) <= np.abs(ts[idx] - t_query)
        best = np.where(left, idx - 1, idx)
    ok = np.abs(ts[best] - t_query) <= TOL_S
    rows = order[best]
    for c in cols:
        if c in df.columns:
            v = df[c].to_numpy(float)[rows]
            out[c] = np.where(ok, v, np.nan)
    return out


def _median(v) -> float:
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else float("nan")


def _rms(v) -> float:
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(np.sqrt(np.mean(a * a))) if a.size else float("nan")


def _cov_eigen(positions) -> np.ndarray:
    """Per-row sqrt eigenvalues (descending) of the Stage-4 3x3 ECEF covariance."""
    if positions is None or not len(positions):
        return np.zeros((0, 3))
    out = []
    for _, r in positions.iterrows():
        c = np.array([[r.get("cov_xx"), r.get("cov_xy"), r.get("cov_xz")], [r.get("cov_xy"), r.get("cov_yy"),
                      r.get("cov_yz")], [r.get("cov_xz"), r.get("cov_yz"), r.get("cov_zz")]], dtype=float)
        if np.all(np.isfinite(c)):
            out.append(np.sqrt(np.clip(np.linalg.eigvalsh(c), 0, None))[::-1])
        else:
            out.append(np.full(3, np.nan))
    return np.array(out)


def run_id_for(run: CampaignRun, start_utc: Optional[str]) -> str:
    """``<YYYYMMDD>_NORAD<id>_<HHMMSS>`` from the observation start (``unknown`` parts if absent)."""
    stamp = "unknown"
    if start_utc:
        try:
            stamp = datetime.fromisoformat(str(start_utc).replace("Z", "+00:00")).strftime("%Y%m%d_%H%M%S")
        except ValueError:
            pass
    d, _, t = stamp.partition("_")
    return f"{d}_NORAD{run.norad_id}_{t or 'unknown'}"


def build_run_tables(run: CampaignRun, run_id: str, reference_frame: str = "rtn"):
    """``(run_row, epoch_rows, station_rows)`` for one loaded run."""
    diag = run.diagnosis
    dsel = diag[diag["diagnosed"].fillna(False).astype(bool)] if len(diag) else diag
    t_ep = _secs(dsel["time_utc"]) if len(dsel) else np.zeros(0)
    m6, m5 = run.stage_manifests.get(6) or {}, run.stage_manifests.get(5) or {}
    s6 = run.stage6_summary or {}
    res = _lookup(run.residuals, t_ep, ["velocity_residual_3d_m_s", "dr_x_m", "dr_y_m", "dr_z_m",
                                        "sigma_radial_obs_m", "sigma_along_track_obs_m", "sigma_cross_track_obs_m",
                                        "geocentric_angular_residual_deg"])
    tr = _lookup(run.trajectory, t_ep, ["sigma_x_smoothed_m", "sigma_y_smoothed_m", "sigma_z_smoothed_m",
                                        "sigma_x_filtered_m", "sigma_y_filtered_m", "sigma_z_filtered_m",
                                        "sigma_vx_smoothed_m_s", "sigma_vy_smoothed_m_s", "sigma_vz_smoothed_m_s",
                                        "x_ecef_smoothed_m", "y_ecef_smoothed_m", "z_ecef_smoothed_m"])
    pos = _lookup(run.positions, t_ep, ["sigma_x_m", "sigma_y_m", "sigma_z_m", "x_ecef_m", "y_ecef_m", "z_ecef_m",
                                        "max_crossing_angle_deg"])
    s5 = np.sqrt(tr["sigma_x_smoothed_m"] ** 2 + tr["sigma_y_smoothed_m"] ** 2 + tr["sigma_z_smoothed_m"] ** 2)
    s5f = np.sqrt(tr["sigma_x_filtered_m"] ** 2 + tr["sigma_y_filtered_m"] ** 2 + tr["sigma_z_filtered_m"] ** 2)
    s5v = np.sqrt(tr["sigma_vx_smoothed_m_s"] ** 2 + tr["sigma_vy_smoothed_m_s"] ** 2 + tr["sigma_vz_smoothed_m_s"] ** 2)
    s4 = np.sqrt(pos["sigma_x_m"] ** 2 + pos["sigma_y_m"] ** 2 + pos["sigma_z_m"] ** 2)
    d45 = np.sqrt((pos["x_ecef_m"] - tr["x_ecef_smoothed_m"]) ** 2 + (pos["y_ecef_m"] - tr["y_ecef_smoothed_m"]) ** 2
                  + (pos["z_ecef_m"] - tr["z_ecef_smoothed_m"]) ** 2)
    pos_rows = None
    if run.positions is not None and len(run.positions):
        tp = _secs(run.positions["time_utc"])
        keep = np.array([np.any(np.abs(t_ep - x) <= TOL_S) for x in tp]) if len(t_ep) else np.zeros(len(tp), bool)
        pos_rows = run.positions[keep]
    eig = _cov_eigen(pos_rows)
    con = run.contributions
    el_ep = np.full(len(t_ep), np.nan)
    if con is not None and len(con):
        ref_col = "reference_time_utc" if "reference_time_utc" in con.columns else "time_utc"
        tc = _secs(con[ref_col])
        elc = con["el_obs_deg"].to_numpy(float)
        for i, t in enumerate(t_ep):
            sel = np.abs(tc - t) <= TOL_S
            if sel.any():
                el_ep[i] = float(np.nanmean(elc[sel]))
    sat = run.satellite
    stations = sorted(((run.stage_manifests.get(4) or {}).get("stations") or {}).keys()
                      or (run.summary.get("station_summaries") or {}).keys())
    infl = run.influence
    if len(infl) and len(t_ep):
        ti = _secs(infl["time_utc"])
        st_infl = infl[np.array([np.any(np.abs(t_ep - x) <= TOL_S) for x in ti])]
    else:
        st_infl = infl.iloc[:0] if len(infl) else infl
    st_used = sorted(str(x) for x in st_infl["station_id"].dropna().unique()) if len(st_infl) else []
    station_set = "+".join(st_used)
    times = dsel["time_utc"] if len(dsel) else []
    start = str(min(times)) if len(times) else None
    end = str(max(times)) if len(times) else None
    precise = bool(run.summary.get("precision_validation", False))
    flags = list(run.metadata_flags) + list(run.summary.get("quality_flags") or [])
    if run.broken:
        flags.append("BROKEN_PROVENANCE")
    inn = m5.get("innovation_statistics") or {}
    dr = {c: dsel[c].to_numpy(float) if c in dsel else np.full(len(dsel), np.nan)
          for c in ("dr_radial_m", "dr_along_track_m", "dr_cross_track_m", "position_residual_3d_m",
                    "position_mahalanobis_obs_cov", "position_mahalanobis_total_cov")}
    if reference_frame == "gcrs":
        comp = {"radial": res["dr_x_m"], "along_track": res["dr_y_m"], "cross_track": res["dr_z_m"]}
    else:
        comp = {"radial": dr["dr_radial_m"], "along_track": dr["dr_along_track_m"],
                "cross_track": dr["dr_cross_track_m"]}
    loo_med = _median(infl["loo_shift_3d_m"]) if len(infl) and "loo_shift_3d_m" in infl else float("nan")
    run_row = {
        "run_id": run_id, "manifest_path": str(run.manifest_path), "manifest_sha256": run.sha256,
        "git_commit": (run.manifest.get("software") or {}).get("git_commit"),
        "pipeline_version": (run.manifest.get("software") or {}).get("tart_gnss"),
        "norad_id": run.norad_id, "satellite_name": sat.get("satellite_name"), "constellation": run.constellation,
        "observation_date": start[:10] if start else None, "start_utc": start, "end_utc": end,
        "duration_s": float(np.ptp(t_ep)) if len(t_ep) else 0.0, "stations": st_used or stations,
        "n_stations": len(st_used or stations), "station_set": station_set,
        "n_epochs": int(len(diag)), "n_valid_epochs": int(len(dsel)),
        "usable_fraction": float(len(dsel) / len(diag)) if len(diag) else 0.0,
        "cohort": "precision_validation" if precise else "descriptive_only", "precision_validation": precise,
        "tle_source": (m6.get("tle_reference") or {}).get("source_kind"),
        "tle_epoch_utc": s6.get("tle_epoch_utc") or (m6.get("tle_reference") or {}).get("epoch_utc"),
        "median_tle_age_days": _median(dsel["tle_age_days"]) if len(dsel) else float("nan"),
        "max_abs_tle_age_days": float(np.nanmax(np.abs(dsel["tle_age_days"].to_numpy(float)))) if len(dsel) else
        float("nan"),
        "median_elevation_deg": _median(el_ep),
        "median_snr": _median(st_infl["snr"]) if len(st_infl) and "snr" in st_infl else float("nan"),
        "median_flagged_fraction": _median(st_infl["flagged_fraction"]) if len(st_infl) and "flagged_fraction" in
        st_infl else float("nan"),
        "median_sigma_az_deg": _median(con["sigma_az_deg"]) if con is not None and len(con) else float("nan"),
        "median_sigma_el_deg": _median(con["sigma_el_deg"]) if con is not None and len(con) else float("nan"),
        "median_sigma_theta_deg": _median(con["sigma_theta_deg"]) if con is not None and len(con) else float("nan"),
        "median_vis_image_separation_deg": _median(st_infl["vis_image_separation_deg"]) if len(st_infl) and
        "vis_image_separation_deg" in st_infl else float("nan"),
        "median_direction_residual_deg": _median(st_infl["direction_residual_angle_deg"]) if len(st_infl) else
        float("nan"),
        "median_condition_number": _median(dsel["stage4_condition_number"]) if len(dsel) else float("nan"),
        "median_min_crossing_angle_deg": _median(dsel["stage4_min_crossing_angle_deg"]) if len(dsel) else float("nan"),
        "median_stage4_sigma_3d_m": _median(s4),
        "median_sigma_major_m": _median(eig[:, 0]) if len(eig) else float("nan"),
        "median_sigma_intermediate_m": _median(eig[:, 1]) if len(eig) else float("nan"),
        "median_sigma_minor_m": _median(eig[:, 2]) if len(eig) else float("nan"),
        "median_stage5_sigma_3d_m": _median(s5), "median_stage5_filtered_sigma_3d_m": _median(s5f),
        "median_stage5_velocity_sigma_m_s": _median(s5v), "stage5_rmse_vs_stage4_m": _rms(d45),
        "stage5_median_nis": inn.get("median_nis"), "stage5_fraction_above_gate": inn.get("fraction_above_gate"),
        "stage5_lag1_autocorrelation": (inn.get("autocorrelation") or {}).get("1"),
        "stage5_n_resets": (m5.get("filter") or {}).get("n_resets"),
        "stage5_dynamics_mismatch": "DYNAMICS_MISMATCH" in list(((m5.get("summary") or {}).get("run_flags")) or []),
        "median_position_residual_3d_m": _median(dr["position_residual_3d_m"]),
        "rms_position_residual_3d_m": _rms(dr["position_residual_3d_m"]),
        "median_radial_m": _median(comp["radial"]), "median_along_track_m": _median(comp["along_track"]),
        "median_cross_track_m": _median(comp["cross_track"]), "rms_radial_m": _rms(comp["radial"]),
        "rms_along_track_m": _rms(comp["along_track"]), "rms_cross_track_m": _rms(comp["cross_track"]),
        "median_sigma_radial_obs_m": _median(res["sigma_radial_obs_m"]),
        "median_sigma_along_track_obs_m": _median(res["sigma_along_track_obs_m"]),
        "median_sigma_cross_track_obs_m": _median(res["sigma_cross_track_obs_m"]),
        "median_velocity_residual_m_s": _median(res["velocity_residual_3d_m_s"]),
        "median_angular_residual_deg": _median(res["geocentric_angular_residual_deg"]),
        "median_position_mahalanobis_obs_cov": _median(dr["position_mahalanobis_obs_cov"]),
        "median_position_mahalanobis_total_cov": _median(dr["position_mahalanobis_total_cov"]),
        "median_leave_one_out_shift_m": loo_med,
        "overall_evidence_level": run.summary.get("overall_evidence_level"),
        "overall_diagnostic_tags": list(run.summary.get("overall_diagnostic_tags") or []),
        "diagnostic_confidence": run.summary.get("overall_confidence"),
        "most_influential_station": run.summary.get("most_influential_station"),
        "quality_flags": sorted(set(flags)),
    }
    epoch_rows = []
    for i, (_, r) in enumerate(dsel.iterrows()):
        epoch_rows.append({
            "run_id": run_id, "time_utc": str(r["time_utc"]), "norad_id": run.norad_id,
            "constellation": run.constellation, "position_residual_3d_m": r.get("position_residual_3d_m"),
            "velocity_residual_3d_m_s": res["velocity_residual_3d_m_s"][i], "dr_radial_m": r.get("dr_radial_m"),
            "dr_along_track_m": r.get("dr_along_track_m"), "dr_cross_track_m": r.get("dr_cross_track_m"),
            "dr_x_m": res["dr_x_m"][i], "dr_y_m": res["dr_y_m"][i], "dr_z_m": res["dr_z_m"][i],
            "position_mahalanobis_obs_cov": r.get("position_mahalanobis_obs_cov"),
            "position_mahalanobis_total_cov": r.get("position_mahalanobis_total_cov"),
            "stage4_sigma_3d_m": s4[i], "stage5_sigma_3d_m": s5[i],
            "stage4_condition_number": r.get("stage4_condition_number"),
            "stage4_min_crossing_angle_deg": r.get("stage4_min_crossing_angle_deg"),
            "stage4_max_crossing_angle_deg": pos["max_crossing_angle_deg"][i],
            "stage4_n_stations": r.get("stage4_n_stations"), "stage5_nis": r.get("stage5_nis"),
            "mean_elevation_deg": el_ep[i], "evidence_level": r.get("evidence_level"),
            "diagnostic_tags": list(r.get("diagnostic_tags") or []), "precision_validation": precise,
            "tle_age_days": r.get("tle_age_days"),
        })
    station_rows = []
    sums = run.summary.get("station_summaries") or {}
    for sid, ss in sums.items():
        g = st_infl[st_infl["station_id"] == sid] if len(st_infl) else st_infl
        cg = con[con["station_id"] == sid] if con is not None and len(con) else None
        if cg is not None and len(cg) and len(t_ep):
            ref_col = "reference_time_utc" if "reference_time_utc" in cg.columns else "time_utc"
            tcg = _secs(cg[ref_col])
            cg = cg[np.array([np.any(np.abs(t_ep - x) <= TOL_S) for x in tcg])]
        station_rows.append({
            "run_id": run_id, "station_id": sid, "norad_id": run.norad_id, "constellation": run.constellation,
            "n_epochs": ss.get("n_epochs"), "usable_fraction": ss.get("usable_fraction"),
            "median_direction_residual_deg": ss.get("median_direction_residual_deg"),
            "median_line_residual_m": ss.get("median_line_residual_m"),
            "median_leave_one_out_shift_m": ss.get("median_leave_one_out_shift_m"),
            "median_leave_one_out_mahalanobis": ss.get("median_leave_one_out_mahalanobis"),
            "median_snr": ss.get("median_snr"), "median_flagged_fraction": ss.get("median_flagged_fraction"),
            "zero_gain_fraction": ss.get("zero_gain_fraction"),
            "calibration_change_count": ss.get("calibration_change_count"),
            "median_time_offset_s": ss.get("timing_offset_median_s"), "mad_time_offset_s": ss.get("timing_offset_mad_s"),
            "median_elevation_deg": _median(cg["el_obs_deg"]) if cg is not None and len(cg) else float("nan"),
            "median_sigma_theta_deg": _median(cg["sigma_theta_deg"]) if cg is not None and len(cg) else float("nan"),
            "median_sigma_az_deg": _median(cg["sigma_az_deg"]) if cg is not None and len(cg) else float("nan"),
            "median_sigma_el_deg": _median(cg["sigma_el_deg"]) if cg is not None and len(cg) else float("nan"),
            "median_vis_image_separation_deg": ss.get("median_image_visibility_offset_deg"),
            "station_tags": list(ss.get("evidence_tags") or []), "health_flags": list(ss.get("health_flags") or []),
            "precision_validation": precise, "overall_evidence_level": run.summary.get("overall_evidence_level"),
        })
    return run_row, epoch_rows, station_rows


def build_campaign_tables(runs: List[CampaignRun], run_ids: List[str], reference_frame: str = "rtn") -> CampaignTables:
    """Build run, epoch and station population tables."""
    import pandas as pd

    r_rows, e_rows, s_rows = [], [], []
    for run, rid in zip(runs, run_ids):
        a, b, c = build_run_tables(run, rid, reference_frame)
        r_rows.append(a)
        e_rows += b
        s_rows += c

    def frame(rows, cols):
        df = pd.DataFrame(rows)
        for col in cols:
            if col not in df.columns:
                df[col] = None
        return df[cols]

    return CampaignTables(frame(r_rows, RUN_COLUMNS), frame(e_rows, EPOCH_COLUMNS), frame(s_rows, STATION_COLUMNS))
