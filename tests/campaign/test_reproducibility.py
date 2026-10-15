"""Pass-to-pass repeatability."""

import numpy as np
import pandas as pd

from tart_gnss.campaign.reproducibility import repeatability


def test_known_scatter():
    runs = pd.DataFrame({"norad_id": [1, 1, 1, 2], "run_id": list("abcd"), "cohort": ["p"] * 4,
                         "median_radial_m": [0.0, 0, 0, 0], "median_along_track_m": [-100.0, 0, 100, 0],
                         "median_cross_track_m": [5.0, 5, 5, 0], "median_sigma_radial_obs_m": 10.0,
                         "median_sigma_along_track_obs_m": 20.0, "median_sigma_cross_track_obs_m": 10.0,
                         "overall_diagnostic_tags": [["A"], ["A"], ["A", "B"], ["A"]],
                         "overall_evidence_level": [1, 1, 1, 0]})
    st = pd.DataFrame({"run_id": list("abc"), "station_id": "x", "median_direction_residual_deg": [0.1, 0.1, 0.1],
                       "median_leave_one_out_mahalanobis": [1.0, 2.0, 3.0]})
    r = repeatability(runs, st)
    assert list(r) == ["1"]
    assert r["1"]["along_track"]["between_pass_std_m"] == 100.0
    assert r["1"]["along_track"]["ratio_between_to_within"] == 5.0
    assert abs(r["1"]["direction_bias_between_pass_std_deg"]["x"]) < 1e-12
    assert np.isclose(r["1"]["evidence_tag_jaccard_mean"], (1 + 0.5 + 0.5) / 3)
