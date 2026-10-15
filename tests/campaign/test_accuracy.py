"""Performance summaries and terminology."""

import pandas as pd

from tart_gnss.campaign.accuracy import TERMS, performance


def test_performance_groups(campaign):
    _, d, s, _ = campaign
    p = s["performance_analysis_cohort"]
    assert set(p) == {"direction", "localisation_3d", "trajectory", "tle_relative_residual"}
    loc = p["localisation_3d"]
    assert loc["median_sigma_major_m"]["median"] > loc["median_sigma_minor_m"]["median"] * 3  # anisotropic
    assert 0 <= p["trajectory"]["dynamics_mismatch_fraction"] <= 1
    assert "absolute" not in " ".join(TERMS.values())


def test_empty():
    out = performance(pd.DataFrame(columns=["run_id", "median_stage4_sigma_3d_m", "stage5_dynamics_mismatch"]), 10, 1)
    assert out["localisation_3d"]["median_stage4_sigma_3d_m"]["n"] == 0
