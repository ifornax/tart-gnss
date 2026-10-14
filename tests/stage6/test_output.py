"""Stage-6 tables."""

import pytest

pytest.importorskip("pandas")

from tart_gnss.stage6.output import RESIDUAL_COLUMNS, TLE_TRAJECTORY_COLUMNS, read_table, write_tables  # noqa: E402


def test_round_trip(tmp_path):
    row = {"time_utc": "2026-10-06T18:00:00.000000Z", "dr_radial_m": 1.5, "valid": True,
           "quality_flags": ["TLE_STALE", "NOT_PRECISION_VALIDATION"], "stage5_quality_flags": ["GOOD"]}
    trow = {"time_utc": row["time_utc"], "x_gcrs_m": 1.0, "r_hat_x": 1.0}
    paths = write_tables([row], [trow], tmp_path)
    for key in ("residuals_csv", "residuals_parquet"):
        if paths[key]:
            df = read_table(paths[key])
            assert list(df.columns) == RESIDUAL_COLUMNS
            assert df["quality_flags"].iloc[0] == ["TLE_STALE", "NOT_PRECISION_VALIDATION"]
            assert df["dr_radial_m"].iloc[0] == 1.5
    tk = paths["tle_trajectory_parquet"] or paths["tle_trajectory_csv"]
    assert list(read_table(tk).columns) == TLE_TRAJECTORY_COLUMNS


def test_spec_columns_present():
    for c in ("position_mahalanobis_obs_cov", "state_mahalanobis_total_cov", "equivalent_along_track_time_offset_s",
              "geocentric_angular_residual_deg", "precision_validation", "n_tle_ensemble", "tle_age_days"):
        assert c in RESIDUAL_COLUMNS
