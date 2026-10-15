"""Canonical schema."""

from tart_gnss.campaign.schema import EPOCH_COLUMNS, RUN_COLUMNS, STATION_COLUMNS, canonical_constellation


def test_constellations():
    assert canonical_constellation("gps") == "GPS" and canonical_constellation("Galileo") == "Galileo"
    assert canonical_constellation("BDS") == "BeiDou" and canonical_constellation(None) == "Unknown"
    assert canonical_constellation("GPS BIIF-9 (PRN 26)") == "Unknown"  # names are never parsed


def test_required_columns():
    for c in ("run_id", "norad_id", "satellite_name", "constellation", "start_utc", "end_utc", "duration_s",
              "n_epochs", "n_valid_epochs", "stations", "n_stations", "precision_validation", "tle_epoch_utc",
              "median_tle_age_days", "max_abs_tle_age_days", "overall_evidence_level", "overall_diagnostic_tags",
              "diagnostic_confidence", "median_sigma_major_m", "median_sigma_intermediate_m", "median_sigma_minor_m"):
        assert c in RUN_COLUMNS
    for c in ("run_id", "time_utc", "velocity_residual_3d_m_s", "stage4_sigma_3d_m", "stage5_sigma_3d_m",
              "stage5_nis", "evidence_level", "tle_age_days"):
        assert c in EPOCH_COLUMNS
    for c in ("run_id", "station_id", "median_leave_one_out_mahalanobis", "zero_gain_fraction", "station_tags"):
        assert c in STATION_COLUMNS
