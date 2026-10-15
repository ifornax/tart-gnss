"""Synthetic campaign validation (specification section 47)."""

import numpy as np
import pandas as pd
import pytest


def _tab(d, name):
    return pd.read_csv(d / "tables" / f"{name}.csv")


def test_campaign_runs_and_cohorts(campaign):
    code, d, s, m = campaign
    assert code == 0
    assert s["n_discovered_runs"] == s["n_included_runs"] == 17 and s["n_excluded_runs"] == 0
    assert s["n_precision_validation_runs"] == 16 and s["n_descriptive_only_runs"] == 1  # the 20-day-old TLE
    assert s["n_satellites"] == 3 and s["constellations"] == ["BeiDou", "GPS", "Galileo"]
    assert s["n_unique_stations"] == 4 and s["sample_size_tier"] == "basic_inference"
    assert s["inference_unit"].startswith("run")


def test_clean_runs(campaign):  # 47.1
    _, d, _, _ = campaign
    r = pd.read_csv(d / "campaign_runs.csv")
    clean = r[r["run_id"].isin(["20261007_NORAD90002_021211"])]
    assert (clean["overall_evidence_level"] == 0).all()
    gps0 = r[r["run_id"] == "20261007_NORAD40534_021211"].iloc[0]
    assert gps0["median_position_residual_3d_m"] < 5 * gps0["median_stage5_sigma_3d_m"]


def test_unstable_station(campaign):  # 47.2
    _, d, s, _ = campaign
    st = _tab(d, "station_summary").set_index("station_id")
    assert st.loc["madagascar", "fraction_runs_station_specific"] > 0.1
    assert st.drop("madagascar")["fraction_runs_station_specific"].max() == 0
    three = st.loc[["madagascar", "mauritius", "south_africa"]]
    assert three["p84_direction_residual_deg"].idxmax() == "madagascar"
    assert three["std_direction_residual_deg"].idxmax() == "madagascar"
    assert s["most_frequently_flagged_stations"] == ["madagascar"]


def test_geometry_dependence(campaign):  # 47.3
    _, _, s, _ = campaign
    g = s["geometry_dependence"]
    two = g["log10_stage4_sigma_vs_crossing_angle_by_n_stations"]["2"]
    assert two["n_runs"] == 3 and two["spearman_r"] == pytest.approx(-1.0) and two["theil_sen_slope"] < 0
    by_n = g["median_stage4_sigma_3d_m_by_n_stations"]
    assert by_n["2"] > by_n["3"] > by_n["4"]


def test_tle_age_dependence(campaign):  # 47.4
    _, d, s, _ = campaign
    a = s["tle_age_dependence"]["residual_vs_abs_tle_age"]
    assert a["spearman_r"] > 0.4 and a["theil_sen_slope"] > 0
    # injected drift: 0.4 s/day x ~3.87 km/s = 1.55 km/day; the estimate is the right order
    assert 500 < a["theil_sen_slope"] < 3000
    t = _tab(d, "tle_age_summary")
    assert t.loc[t["bin"] == ">14 d", "n_precision_runs"].iloc[0] == 0


def test_constellation_difference(campaign):  # 47.5
    _, d, s, _ = campaign
    snr = s["constellation_comparisons"]["median_snr"]
    assert snr["test"] == "Mann-Whitney U" and snr["p_value"] < 0.01
    assert snr["effect_size"] == pytest.approx(-1.0)  # Galileo SNR below GPS in every pair
    assert "BeiDou" in snr["under_sampled"]
    c = _tab(d, "constellation_summary").set_index("constellation")
    assert c.loc["Galileo", "median_snr"] < c.loc["GPS", "median_snr"] and bool(c.loc["BeiDou", "under_sampled"])


def test_repeatability(campaign):  # 47.6
    _, _, s, _ = campaign
    gal = s["repeatability"]["90001"]
    assert gal["n_passes"] == 3
    expected = 3874.0  # |v| x 1 s: passes offset by -1, 0, +1 s
    assert gal["along_track"]["between_pass_std_m"] == pytest.approx(expected, rel=0.15)
    assert gal["along_track"]["ratio_between_to_within"] > 3


def test_bootstrap_intervals_present(campaign):
    _, _, s, _ = campaign
    for key in ("median_stage4_3d_uncertainty_m", "median_tle_relative_residual_3d_m", "median_tle_relative_radial_m"):
        lo, hi = s[key]["ci95"]
        assert lo <= s[key]["median"] <= hi and np.isfinite(lo)


def test_terminology_and_no_physical_rates(campaign):
    _, d, s, _ = campaign
    text = (d / "summary.json").read_text().lower()
    assert "absolute accuracy" not in text and "tle_relative" in text
    for bad in ("manoeuvre", "clock_failure", "hardware_failure", "rfi_event"):
        assert bad not in text
