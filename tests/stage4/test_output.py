"""Stage-4 tables and manifest serialisation."""

import json
import math

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage4.output import (  # noqa: E402
    CONTRIBUTION_COLUMNS,
    POSITION_COLUMNS,
    read_table,
    write_manifest,
    write_tables,
)


def test_required_columns_present():
    spec_pos = ["time_utc", "time_mjd", "norad_id", "x_ecef_m", "y_ecef_m", "z_ecef_m", "cov_xx", "cov_xy",
                "cov_xz", "cov_yy", "cov_yz", "cov_zz", "sigma_x_m", "sigma_y_m", "sigma_z_m", "n_stations",
                "rms_line_residual_m", "max_line_residual_m", "condition_number", "min_crossing_angle_deg",
                "max_crossing_angle_deg", "chi2", "reduced_chi2", "fit_success", "quality_flags", "latitude_deg",
                "longitude_deg", "altitude_m", "x_p16", "x_p50", "x_p84", "z_p84"]
    spec_con = ["time_utc", "station_id", "station_time_utc", "station_time_offset_s", "station_x_ecef_m",
                "station_y_ecef_m", "station_z_ecef_m", "az_obs_deg", "el_obs_deg", "sigma_az_deg", "sigma_el_deg",
                "cov_az_el", "range_m", "line_residual_m", "stage3_quality_flags", "stage4_station_flags"]
    assert set(spec_pos) <= set(POSITION_COLUMNS)
    assert set(spec_con) <= set(CONTRIBUTION_COLUMNS)


def test_round_trip(tmp_path):
    pos = [{"time_utc": "2026-06-25T09:14:00.000000Z", "x_ecef_m": 1.0, "quality_flags": ["HIGH_CHI2", "GOOD"],
            "valid": True},
           {"time_utc": "2026-06-25T09:14:01.000000Z", "x_ecef_m": 2.0, "quality_flags": [], "valid": False}]
    con = [{"time_utc": "2026-06-25T09:14:00.000000Z", "station_id": "a",
            "station_time_utc": "2026-06-25T09:14:00.300000Z", "stage3_quality_flags": ["LOW_SNR"],
            "stage4_station_flags": []}]
    paths = write_tables(pos, con, tmp_path)
    for key in ("positions_csv", "station_contributions_csv"):
        assert paths[key]
    for p in (paths["positions_csv"], paths["positions_parquet"]):
        if p is None:
            continue
        df = read_table(p)
        assert list(df.columns) == POSITION_COLUMNS
        assert df.loc[0, "quality_flags"] == ["HIGH_CHI2", "GOOD"] and df.loc[1, "quality_flags"] == []
    c = read_table(paths["station_contributions_csv"])
    assert c.loc[0, "stage3_quality_flags"] == ["LOW_SNR"]
    assert (c.loc[0, "station_time_utc"] - c.loc[0, "time_utc"]).total_seconds() == pytest.approx(0.3)


def test_manifest_nan_and_inf_become_null(tmp_path):
    p = write_manifest({"a": float("nan"), "b": [float("inf"), 1.0], "c": np.float64(2.5), "d": np.arange(2)},
                       tmp_path / "m.json")
    doc = json.loads(open(p).read())
    assert doc == {"a": None, "b": [None, 1.0], "c": 2.5, "d": [0, 1]}
    assert not any(isinstance(v, float) and math.isnan(v) for v in doc.values())
