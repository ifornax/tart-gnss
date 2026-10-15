"""Grouped comparisons and dependence tables."""

import numpy as np
import pandas as pd

from tart_gnss.campaign.constellation_metrics import compare, geometry_summary, tle_age_summary


def test_compare_respects_group_size():
    runs = pd.DataFrame({"constellation": ["GPS"] * 4 + ["Galileo"] * 4 + ["BeiDou"], "median_snr":
                         [20, 21, 19, 20, 12, 13, 11, 12, 16]})
    c = compare(runs, "constellation", "median_snr", 3)
    assert c["test"] == "Mann-Whitney U" and c["under_sampled"] == ["BeiDou"] and c["effect_size"] == -1.0
    c = compare(runs, "constellation", "median_snr", 5)
    assert c["test"].startswith("unavailable")


def test_tle_age_bins_and_dependence():
    age = np.array([0.5, 2, 5, 8, 10, 20.0])
    runs = pd.DataFrame({"run_id": list("abcdef"), "median_tle_age_days": age, "precision_validation": age < 14,
                         "median_position_residual_3d_m": 1000 * age + 100, "median_along_track_m": 1000 * age,
                         "median_radial_m": 0.0})
    tab, dep = tle_age_summary(runs)
    assert list(tab["n_runs"]) == [1, 1, 1, 2, 1] and dep["residual_vs_abs_tle_age"]["spearman_r"] == 1.0
    assert tab["n_precision_runs"].iloc[-1] == 0


def test_geometry_by_station_count():
    runs = pd.DataFrame({"run_id": list("abcdef"), "n_stations": [2, 2, 2, 3, 3, 3],
                         "median_min_crossing_angle_deg": [2, 5, 8, 2, 2.1, 2.2],
                         "median_stage4_sigma_3d_m": [50e3, 30e3, 20e3, 15e3, 15e3, 16e3],
                         "median_sigma_major_m": 1.0, "median_leave_one_out_shift_m": 1.0,
                         "median_condition_number": [1e3, 500, 300, 200, 210, 220],
                         "median_position_residual_3d_m": 1.0})
    _, dep = geometry_summary(runs, runs)
    assert dep["log10_stage4_sigma_vs_crossing_angle_by_n_stations"]["2"]["spearman_r"] == -1.0
