"""Cross-run station stability (one summary row per station; no opaque score).

Every statistic is over the station's runs (one per-run value per run).
Stability is the run-to-run scatter of the station's per-run median
direction residual: ``std_direction_residual_deg`` (the main measure; a robust
MAD would hide a station that is biased only in some runs, which is exactly
the instability of interest) with the MAD and 84th percentile alongside. Fractions of
runs in which the station carried each Stage-7 station tag are listed
separately. No combined ranking is produced: the table is a transparent
comparison.
"""

from __future__ import annotations

import numpy as np

from tart_gnss.campaign.bootstrap import bootstrap_run_level, median_of
from tart_gnss.campaign.statistics import MAD_TO_SIGMA

TAGS = ("STATION_SPECIFIC", "CALIBRATION_CORRELATED", "TIMING_CORRELATED", "LOCALISATION_CORRELATED")


def _med(v):
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else float("nan")


def _mad(v):
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(MAD_TO_SIGMA * np.median(np.abs(a - np.median(a)))) if a.size else float("nan")


def _std(v):
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(np.std(a, ddof=1)) if a.size > 1 else float("nan")


def _p84(v):
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(np.percentile(a, 84)) if a.size else float("nan")


def summarize_stations(station_table, n_bootstrap: int = 0, seed: int = 20260625, min_runs: int = 3):
    """Compute cross-run station stability statistics."""
    import pandas as pd

    rows = []
    if station_table is None or not len(station_table):
        return pd.DataFrame()
    for sid, g in station_table.groupby("station_id"):
        n_runs = int(g["run_id"].nunique())
        rec = {
            "station_id": sid, "n_runs": n_runs, "under_sampled": n_runs < min_runs,
            "n_epochs": int(pd.to_numeric(g["n_epochs"], errors="coerce").fillna(0).sum()),
            "usable_fraction": _med(g["usable_fraction"]),
            "median_direction_residual_deg": _med(g["median_direction_residual_deg"]),
            "mad_direction_residual_deg": _mad(g["median_direction_residual_deg"]),
            "std_direction_residual_deg": _std(g["median_direction_residual_deg"]),
            "p84_direction_residual_deg": _p84(g["median_direction_residual_deg"]),
            "p84_leave_one_out_mahalanobis": _p84(g["median_leave_one_out_mahalanobis"]),
            "median_line_residual_m": _med(g["median_line_residual_m"]),
            "median_leave_one_out_shift_m": _med(g["median_leave_one_out_shift_m"]),
            "median_leave_one_out_mahalanobis": _med(g["median_leave_one_out_mahalanobis"]),
            "median_snr": _med(g["median_snr"]), "median_flagged_fraction": _med(g["median_flagged_fraction"]),
            "median_sigma_theta_deg": _med(g["median_sigma_theta_deg"]),
            "zero_gain_incidence": float(np.mean(pd.to_numeric(g["zero_gain_fraction"], errors="coerce").fillna(0) > 0)),
            "gain_change_incidence": float(np.mean(pd.to_numeric(g["calibration_change_count"], errors="coerce")
                                                   .fillna(0) > 0)),
            "timing_offset_median_s": _med(g["median_time_offset_s"]),
            "timing_offset_mad_s": _mad(g["median_time_offset_s"]),
        }
        for t in TAGS:
            rec[f"fraction_runs_{t.lower()}"] = float(np.mean([t in (x or []) for x in g["station_tags"]]))
        rec["fraction_runs_any_station_tag"] = float(np.mean([bool(set(x or []) & set(TAGS)) for x in g["station_tags"]]))
        if n_bootstrap and n_runs >= 2:
            per_run = g.groupby("run_id", as_index=False)["median_direction_residual_deg"].median()
            ci = bootstrap_run_level(per_run, median_of("median_direction_residual_deg"), n_bootstrap, seed)
            rec["median_direction_residual_ci95"] = [ci.ci_low, ci.ci_high]
        else:
            rec["median_direction_residual_ci95"] = [float("nan"), float("nan")]
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("station_id").reset_index(drop=True)


def stability_highlights(summary, min_runs: int = 3):
    """Most stable stations (lowest run-to-run MAD) and most frequently flagged ones (adequately sampled only)."""
    if summary is None or not len(summary):
        return {"most_stable": [], "most_frequently_flagged": [], "note": "no stations"}
    ok = summary[summary["n_runs"] >= min_runs]
    note = "" if len(ok) else f"no station has >= {min_runs} runs; highlights not computed"
    stable = ok.sort_values("std_direction_residual_deg")["station_id"].head(3).tolist() if len(ok) else []
    flagged = ok[ok["fraction_runs_any_station_tag"] > 0].sort_values(
        "fraction_runs_any_station_tag", ascending=False)["station_id"].head(3).tolist() if len(ok) else []
    return {"most_stable": stable, "most_frequently_flagged": flagged, "note": note,
            "criterion_stable": "lowest run-to-run standard deviation of the per-run median direction residual",
            "criterion_flagged": "highest fraction of runs with a Stage-7 station tag"}
