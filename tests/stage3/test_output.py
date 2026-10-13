"""Table and manifest writing; round trip through CSV and Parquet."""

import json
import math

import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage3 import output  # noqa: E402


def _row(**kw):
    r = {c: math.nan for c in output.COLUMNS}
    r.update(time_utc="2026-06-25T09:14:00.000000Z", time_mjd=61216.38, station_id="alpha", tart_id="zz-alpha",
             norad_id=40534, primary_method="visibility", az_obs_deg=210.0, el_obs_deg=44.0, fit_success=True,
             optimizer_status="1: ok", usable=True, quality_flags=["LOW_SNR", "CALIBRATION_WARNING"],
             stage2_three_station_valid=False)
    r.update(kw)
    return r


def test_schema_contains_required_fields():
    required = {"time_utc", "time_mjd", "station_id", "norad_id", "l", "m", "sigma_l", "sigma_m", "cov_lm",
                "az_obs_deg", "el_obs_deg", "sigma_az_deg", "sigma_el_deg", "cov_az_el", "complex_amp_real",
                "complex_amp_imag", "apparent_amplitude", "chi2", "reduced_chi2", "rms_complex_residual",
                "median_abs_residual", "snr", "n_baselines_total", "n_baselines_used", "flagged_fraction",
                "az_image_deg", "el_image_deg", "image_snr", "vis_image_separation_deg", "fit_success",
                "optimizer_status", "quality_flags", "az_tle_deg", "el_tle_deg"}
    assert required <= set(output.COLUMNS)


def test_csv_and_parquet_round_trip(tmp_path):
    paths = output.write_tables([_row(), _row(quality_flags=[], station_id="beta")], tmp_path)
    df = output.read_directions(paths["csv"])
    assert df["quality_flags"].tolist() == [["LOW_SNR", "CALIBRATION_WARNING"], []]
    assert df["stage2_three_station_valid"].tolist() == [False, False]
    if paths["parquet"]:
        dp = output.read_directions(paths["parquet"])
        assert dp["quality_flags"].tolist() == [["LOW_SNR", "CALIBRATION_WARNING"], []]
        assert dp["az_obs_deg"].tolist() == [210.0, 210.0]


def test_manifest_nan_becomes_null(tmp_path):
    p = output.write_manifest({"a": math.nan, "b": [1.0, math.nan]}, tmp_path / "m.json")
    assert json.loads(open(p).read()) == {"a": None, "b": [1.0, None]}
