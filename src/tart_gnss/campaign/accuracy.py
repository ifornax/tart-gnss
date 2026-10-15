"""Localisation, trajectory and TLE-relative performance (with careful terminology).

Terms used throughout (never "absolute accuracy"):

* **measurement uncertainty**: formal 1-sigma from the fits (Stage-3 direction
  sigma, Stage-4 position covariance, Stage-5 smoothed covariance);
* **TLE-relative residual**: Stage-6 observed - TLE, meaningful as a
  performance measure only in the ``precision_validation`` cohort;
* **network repeatability**: pass-to-pass scatter (see
  :mod:`tart_gnss.campaign.reproducibility`);
* **trajectory residual**: Stage-4 positions about the Stage-5 trajectory
  (``stage5_rmse_vs_stage4_m``), not a TLE comparison.

Every population statistic is computed over runs (one per-run median per
run) with a run-level bootstrap interval.
"""

from __future__ import annotations

from typing import Dict

from tart_gnss.campaign.bootstrap import bootstrap_run_level, median_of
from tart_gnss.campaign.statistics import describe

TERMS = {
    "measurement_uncertainty": "formal 1-sigma of the Stage-3/4/5 fits",
    "tle_relative_residual": "Stage-6 observed - TLE (precision-validation cohort only by default)",
    "network_repeatability": "pass-to-pass scatter of per-pass summaries",
    "trajectory_residual": "Stage-4 positions about the Stage-5 trajectory (RMS)",
}

DIRECTION = ["median_sigma_az_deg", "median_sigma_el_deg", "median_sigma_theta_deg",
             "median_vis_image_separation_deg", "median_direction_residual_deg"]
LOCALISATION = ["median_stage4_sigma_3d_m", "median_sigma_major_m", "median_sigma_intermediate_m",
                "median_sigma_minor_m"]
TRAJECTORY = ["median_stage5_filtered_sigma_3d_m", "median_stage5_sigma_3d_m", "median_stage5_velocity_sigma_m_s",
              "stage5_rmse_vs_stage4_m", "stage5_median_nis", "stage5_fraction_above_gate",
              "stage5_lag1_autocorrelation", "stage5_n_resets"]
RESIDUAL = ["median_position_residual_3d_m", "median_radial_m", "median_along_track_m", "median_cross_track_m",
            "median_velocity_residual_m_s", "median_angular_residual_deg"]


def performance(runs, n_bootstrap: int, seed: int) -> Dict[str, Dict[str, object]]:
    """Per-run-median population statistics with run-level bootstrap CIs of the median."""
    out: Dict[str, Dict[str, object]] = {}
    for group, cols in (("direction", DIRECTION), ("localisation_3d", LOCALISATION), ("trajectory", TRAJECTORY),
                        ("tle_relative_residual", RESIDUAL)):
        g = {}
        for c in cols:
            if c not in runs:
                continue
            d = describe(runs[c].astype(float)) if len(runs) else describe([])
            ci = bootstrap_run_level(runs, median_of(c), n_bootstrap, seed) if len(runs) else None
            g[c] = {**d, "median_ci95": [ci.ci_low, ci.ci_high] if ci else None}
        out[group] = g
    if len(runs):
        out["trajectory"]["dynamics_mismatch_fraction"] = float(runs["stage5_dynamics_mismatch"].astype(bool).mean())
    return out


def epoch_residual_distribution(epochs) -> Dict[str, Dict[str, float]]:
    """Descriptive (pooled-epoch) distributions; for display only, never for inference."""
    return {c: describe(epochs[c].astype(float)) for c in ("position_residual_3d_m", "dr_radial_m", "dr_along_track_m",
                                                          "dr_cross_track_m", "velocity_residual_3d_m_s")
            if c in epochs}
