"""Grouped comparisons: constellation, satellite, station set, geometry, elevation, SNR and TLE age.

All grouped statistics aggregate per run first (one value per run) and
compare groups only when each has at least ``min_runs_per_group`` runs
(smaller groups are listed as ``under_sampled``). Residual statistics are
taken from the precision cohort; uncertainty and SNR statistics from all
included runs.
"""

from __future__ import annotations

from typing import Dict, Sequence

import numpy as np

from tart_gnss.campaign.statistics import binned, compare_groups, dependence, describe

CROSSING_EDGES = (0.0, 2.0, 4.0, 6.0, 10.0, 20.0, 90.0)
TLE_AGE_EDGES = (0.0, 1.0, 3.0, 7.0, 14.0, np.inf)
ELEVATION_EDGES = (20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 90.0)


def _evidence_distribution(levels) -> Dict[str, int]:
    v = [int(x) for x in levels if x is not None and x == x]
    return {str(k): v.count(k) for k in range(6)}


def group_summary(runs, by: str, science, min_runs_per_group: int = 3):
    """One row per value of ``by``: n, precision fraction, medians of performance and evidence distribution."""
    import pandas as pd

    rows = []
    if not len(runs):
        return pd.DataFrame()
    for key, g in runs.groupby(by, dropna=False):
        gs = science[science[by] == key] if len(science) else science
        rec = {by: key, "n_runs": int(len(g)), "n_precision_runs": int(g["precision_validation"].sum()),
               "precision_validation_fraction": float(g["precision_validation"].mean()),
               "under_sampled": len(g) < min_runs_per_group}
        for c in ("median_snr", "median_sigma_theta_deg", "median_stage4_sigma_3d_m", "median_stage5_sigma_3d_m",
                  "median_min_crossing_angle_deg", "median_condition_number", "usable_fraction"):
            d = describe(g[c].astype(float))
            rec.update({c: d["median"], f"{c}_mad": d["mad"], f"{c}_p16": d["p16"], f"{c}_p84": d["p84"]})
        for c in ("median_position_residual_3d_m", "median_radial_m", "median_along_track_m", "median_cross_track_m"):
            d = describe(gs[c].astype(float)) if len(gs) else describe([])
            rec.update({f"tle_relative_{c}": d["median"], f"tle_relative_{c}_mad": d["mad"]})
        rec["evidence_levels"] = _evidence_distribution(g["overall_evidence_level"])
        rows.append(rec)
    return pd.DataFrame(rows)


def compare(runs, by: str, column: str, min_runs_per_group: int = 3) -> Dict[str, object]:
    """Mann-Whitney / Kruskal-Wallis of per-run ``column`` across ``by`` groups (with effect sizes)."""
    groups = {str(k): g[column].astype(float).tolist() for k, g in runs.groupby(by)} if len(runs) else {}
    return {"by": by, "metric": column, **compare_groups(groups, min_runs_per_group)}


def satellite_summary(runs, science, min_runs_per_group: int = 3):
    """Per NORAD ID: runs, time span, median TLE-relative residual and run-to-run scatter, uncertainty, levels."""
    import pandas as pd

    rows = []
    for nid, g in runs.groupby("norad_id"):
        gs = science[science["norad_id"] == nid] if len(science) else science
        res = describe(gs["median_position_residual_3d_m"].astype(float)) if len(gs) else describe([])
        rows.append({
            "norad_id": nid, "satellite_name": g["satellite_name"].iloc[0], "constellation": g["constellation"].iloc[0],
            "n_runs": int(len(g)), "n_precision_runs": int(len(gs)), "first_start_utc": g["start_utc"].min(),
            "last_start_utc": g["start_utc"].max(), "median_tle_relative_residual_3d_m": res["median"],
            "run_to_run_scatter_mad_m": res["mad"],
            "run_to_run_scatter_std_m": float(np.std(gs["median_position_residual_3d_m"].astype(float), ddof=1))
            if len(gs) > 1 else float("nan"),
            "median_stage4_sigma_3d_m": describe(g["median_stage4_sigma_3d_m"].astype(float))["median"],
            "median_stage5_sigma_3d_m": describe(g["median_stage5_sigma_3d_m"].astype(float))["median"],
            "evidence_levels": _evidence_distribution(g["overall_evidence_level"]),
            "under_sampled": len(g) < min_runs_per_group,
        })
    return pd.DataFrame(rows)


def geometry_summary(runs, science):
    """Per-run medians binned by minimum crossing angle, plus per-run dependence statistics."""
    tab = binned(runs, "median_min_crossing_angle_deg", ["median_stage4_sigma_3d_m", "median_sigma_major_m",
                                                         "median_leave_one_out_shift_m"], CROSSING_EDGES, " deg")
    if len(science):
        tab_s = binned(science, "median_min_crossing_angle_deg", ["median_position_residual_3d_m"], CROSSING_EDGES,
                       " deg")
        tab["tle_relative_residual_3d_m_median"] = tab_s["median_position_residual_3d_m_median"]
        tab["n_precision_runs"] = tab_s["n_runs"]
    dep = {
        "log10_stage4_sigma_vs_crossing_angle": dependence(runs["median_min_crossing_angle_deg"],
                                                           np.log10(runs["median_stage4_sigma_3d_m"].astype(float))),
        "log10_stage4_sigma_vs_log10_condition_number": dependence(
            np.log10(runs["median_condition_number"].astype(float)),
            np.log10(runs["median_stage4_sigma_3d_m"].astype(float))),
        "loo_shift_vs_crossing_angle": dependence(runs["median_min_crossing_angle_deg"],
                                                  runs["median_leave_one_out_shift_m"]),
        "residual_vs_crossing_angle": dependence(science["median_min_crossing_angle_deg"],
                                                 science["median_position_residual_3d_m"]) if len(science) else None,
    }
    # Station count changes the geometry far more than the crossing angle; compare like with like.
    dep["log10_stage4_sigma_vs_crossing_angle_by_n_stations"] = {
        str(int(k)): dependence(g["median_min_crossing_angle_deg"], np.log10(g["median_stage4_sigma_3d_m"].astype(float)))
        for k, g in runs.groupby("n_stations")}
    dep["median_stage4_sigma_3d_m_by_n_stations"] = {
        str(int(k)): describe(g["median_stage4_sigma_3d_m"].astype(float))["median"] for k, g in runs.groupby("n_stations")}
    return tab, dep


def tle_age_summary(runs):
    """Per-run residuals binned by |TLE age| (all included runs, cohort marked), with signed-age dependence."""
    r = runs.copy()
    r["abs_tle_age_days"] = r["median_tle_age_days"].astype(float).abs()
    tab = binned(r, "abs_tle_age_days", ["median_position_residual_3d_m", "median_along_track_m", "median_radial_m"],
                 TLE_AGE_EDGES, " d")
    tab["n_precision_runs"] = [int(((r["abs_tle_age_days"] >= lo) & (r["abs_tle_age_days"] < hi)
                                    & r["precision_validation"]).sum()) for lo, hi in zip(TLE_AGE_EDGES[:-1],
                                                                                         TLE_AGE_EDGES[1:])]
    dep = {"residual_vs_abs_tle_age": dependence(r["abs_tle_age_days"], r["median_position_residual_3d_m"]),
           "residual_vs_signed_tle_age": dependence(r["median_tle_age_days"], r["median_position_residual_3d_m"]),
           "note": "descriptive; includes stale (descriptive-cohort) runs by construction; no causation inferred"}
    return tab, dep


def elevation_summary(stations):
    """Station-level (run x station) direction uncertainty and residual binned by median elevation."""
    return binned(stations, "median_elevation_deg", ["median_sigma_theta_deg", "median_direction_residual_deg",
                                                     "median_line_residual_m"], ELEVATION_EDGES, " deg")


def snr_summary(stations, n_bins: int = 4):
    """Station-level metrics in SNR quantile bins (linear: SNR spans less than a decade in practice)."""
    s = stations["median_snr"].astype(float)
    s = s[np.isfinite(s)]
    if len(s) < n_bins:
        edges: Sequence[float] = (0.0, np.inf)
    else:
        q = np.unique(np.quantile(s, np.linspace(0, 1, n_bins + 1)))
        edges = list(q[:-1]) + [np.inf]
    return binned(stations, "median_snr", ["median_sigma_theta_deg", "median_direction_residual_deg",
                                           "median_leave_one_out_mahalanobis"], edges)
