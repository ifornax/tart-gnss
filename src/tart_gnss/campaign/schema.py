"""Canonical campaign schema: supported versions, constellation names, table columns."""

from __future__ import annotations

from typing import Optional

SUPPORTED_STAGE7_SCHEMA = (1,)
CONSTELLATIONS = ("GPS", "Galileo", "BeiDou", "GLONASS", "Unknown")
_ALIASES = {"gps": "GPS", "navstar": "GPS", "galileo": "Galileo", "gal": "Galileo", "beidou": "BeiDou",
            "bds": "BeiDou", "compass": "BeiDou", "glonass": "GLONASS", "glo": "GLONASS"}

RUN_COLUMNS = [
    "run_id", "manifest_path", "manifest_sha256", "git_commit", "pipeline_version", "norad_id", "satellite_name",
    "constellation", "observation_date", "start_utc", "end_utc", "duration_s", "stations", "n_stations",
    "station_set", "n_epochs", "n_valid_epochs", "usable_fraction", "cohort", "precision_validation",
    "tle_source", "tle_epoch_utc", "median_tle_age_days", "max_abs_tle_age_days",
    "median_elevation_deg", "median_snr", "median_flagged_fraction", "median_sigma_az_deg", "median_sigma_el_deg",
    "median_sigma_theta_deg", "median_vis_image_separation_deg", "median_direction_residual_deg",
    "median_condition_number", "median_min_crossing_angle_deg", "median_stage4_sigma_3d_m",
    "median_sigma_major_m", "median_sigma_intermediate_m", "median_sigma_minor_m",
    "median_stage5_sigma_3d_m", "median_stage5_filtered_sigma_3d_m", "median_stage5_velocity_sigma_m_s",
    "stage5_rmse_vs_stage4_m", "stage5_median_nis", "stage5_fraction_above_gate", "stage5_lag1_autocorrelation",
    "stage5_n_resets", "stage5_dynamics_mismatch",
    "median_position_residual_3d_m", "rms_position_residual_3d_m", "median_radial_m", "median_along_track_m",
    "median_cross_track_m", "rms_radial_m", "rms_along_track_m", "rms_cross_track_m",
    "median_sigma_radial_obs_m", "median_sigma_along_track_obs_m", "median_sigma_cross_track_obs_m",
    "median_velocity_residual_m_s", "median_angular_residual_deg",
    "median_position_mahalanobis_obs_cov", "median_position_mahalanobis_total_cov",
    "median_leave_one_out_shift_m", "overall_evidence_level", "overall_diagnostic_tags", "diagnostic_confidence",
    "most_influential_station", "quality_flags",
]
EPOCH_COLUMNS = [
    "run_id", "time_utc", "norad_id", "constellation", "position_residual_3d_m", "velocity_residual_3d_m_s",
    "dr_radial_m", "dr_along_track_m", "dr_cross_track_m", "dr_x_m", "dr_y_m", "dr_z_m",
    "position_mahalanobis_obs_cov", "position_mahalanobis_total_cov", "stage4_sigma_3d_m", "stage5_sigma_3d_m",
    "stage4_condition_number", "stage4_min_crossing_angle_deg", "stage4_max_crossing_angle_deg", "stage4_n_stations",
    "stage5_nis", "mean_elevation_deg", "evidence_level", "diagnostic_tags", "precision_validation", "tle_age_days",
]
STATION_COLUMNS = [
    "run_id", "station_id", "norad_id", "constellation", "n_epochs", "usable_fraction",
    "median_direction_residual_deg", "median_line_residual_m", "median_leave_one_out_shift_m",
    "median_leave_one_out_mahalanobis", "median_snr", "median_flagged_fraction", "zero_gain_fraction",
    "calibration_change_count", "median_time_offset_s", "mad_time_offset_s", "median_elevation_deg",
    "median_sigma_theta_deg", "median_sigma_az_deg", "median_sigma_el_deg", "median_vis_image_separation_deg",
    "station_tags", "health_flags", "precision_validation", "overall_evidence_level",
]
LIST_COLUMNS = ("overall_diagnostic_tags", "diagnostic_tags", "station_tags", "health_flags", "quality_flags",
                "stations")


def canonical_constellation(value: Optional[str]) -> str:
    """Canonical constellation name from structured metadata (``Unknown`` if absent or unrecognised)."""
    if not value:
        return "Unknown"
    v = str(value).strip()
    for c in CONSTELLATIONS:
        if v.lower() == c.lower():
            return c
    return _ALIASES.get(v.lower(), "Unknown")
