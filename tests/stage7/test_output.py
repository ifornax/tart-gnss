"""Stage-7 tables."""

import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage7.output import (  # noqa: E402
    CORRELATION_COLUMNS,
    DIAGNOSIS_COLUMNS,
    INFLUENCE_COLUMNS,
    read_table,
    write_tables,
)


def test_round_trip_and_multi_satellite_keys(tmp_path):
    diag = [{"time_utc": "2026-06-25T09:00:00Z", "norad_id": 40534, "diagnostic_tags": ["NORMAL"], "evidence_level": 0}]
    infl = [{"time_utc": "2026-06-25T09:00:00Z", "norad_id": 40534, "station_id": "a", "stage3_flags": ["X", "Y"]}]
    paths = write_tables(diag, infl, [], pd.DataFrame(), tmp_path)
    d = read_table(paths["diagnosis_csv"])
    assert list(d.columns) == DIAGNOSIS_COLUMNS and d["diagnostic_tags"].iloc[0] == ["NORMAL"]
    i = read_table(paths["station_influence_csv"])
    assert list(i.columns) == INFLUENCE_COLUMNS and i["stage3_flags"].iloc[0] == ["X", "Y"]
    assert list(pd.read_csv(paths["correlations_csv"]).columns) == CORRELATION_COLUMNS
    for cols in (DIAGNOSIS_COLUMNS, INFLUENCE_COLUMNS):
        assert "norad_id" in cols and "time_utc" in cols
    assert "station_id" in INFLUENCE_COLUMNS


def test_spec_columns():
    for c in ("max_leave_one_out_shift_m", "most_influential_station", "evidence_level", "diagnostic_confidence",
              "stage5_nis", "stage4_min_crossing_angle_deg", "tle_age_days", "precision_validation"):
        assert c in DIAGNOSIS_COLUMNS
    for c in ("loo_shift_radial_m", "loo_mahalanobis", "degraded_geometry", "direction_residual_az_deg",
              "calibration_warning", "timing_warning"):
        assert c in INFLUENCE_COLUMNS
