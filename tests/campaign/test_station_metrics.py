"""Station stability."""

import pandas as pd

from tart_gnss.campaign.station_metrics import stability_highlights, summarize_stations


def _tab():
    rows = []
    for k in range(6):
        for sid, v in (("a", 0.01), ("b", 0.01 + (0.05 if k < 2 else 0)), ("c", 0.01 + 0.001 * k)):
            rows.append({"run_id": f"r{k}", "station_id": sid, "n_epochs": 10, "usable_fraction": 1.0,
                         "median_direction_residual_deg": v, "median_line_residual_m": 1.0,
                         "median_leave_one_out_shift_m": 1.0, "median_leave_one_out_mahalanobis": 1.0,
                         "median_snr": 20, "median_flagged_fraction": 0.1, "median_sigma_theta_deg": 0.005,
                         "zero_gain_fraction": 0.0, "calibration_change_count": 0, "median_time_offset_s": 0.1,
                         "station_tags": ["STATION_SPECIFIC"] if sid == "b" and k < 2 else []})
    return pd.DataFrame(rows)


def test_intermittent_instability_is_visible():
    s = summarize_stations(_tab(), n_bootstrap=50, seed=1).set_index("station_id")
    assert s["std_direction_residual_deg"].idxmax() == "b"
    assert s.loc["b", "mad_direction_residual_deg"] < s.loc["b", "std_direction_residual_deg"]  # MAD hides it
    assert s.loc["b", "fraction_runs_station_specific"] == 1 / 3
    h = stability_highlights(s.reset_index())
    assert h["most_stable"][0] == "a" and h["most_frequently_flagged"] == ["b"]


def test_under_sampled_highlights():
    h = stability_highlights(summarize_stations(_tab()[_tab()["run_id"].isin(["r0", "r1"])]))
    assert h["most_stable"] == [] and "no station" in h["note"]
