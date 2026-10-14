"""Synthetic diagnosis scenarios through the full Stage 3-7 chain (specification sections 35-37)."""

import numpy as np
import pytest


def test_clean_network_is_level_0(runs):  # 35.1
    r = runs("clean")
    assert r.code == 0 and r.level == 0 and r.tags == {"NORMAL"}
    assert r.summary["n_anomalous_epochs"] == 0 and r.summary["overall_confidence"] == "HIGH"


def test_biased_station_is_identified(runs):  # 35.2
    r = runs("bias")
    assert "STATION_SPECIFIC" in r.tags and r.level == 2
    assert "LOCALISATION_CORRELATED" not in r.tags  # the picked station explains its own pair inconsistency
    assert r.summary["station_test"]["station_specific"] == ["madagascar"]
    assert r.summary["most_influential_station"] == "madagascar"
    assert r.summary["station_test"]["picks"] == {"madagascar": 36}
    # With three stations a leave-one-out shift is large whichever station is removed (the biased line stays in
    # the healthy-station removals), so the influence ranking alone is weak evidence; the station test is not.
    loo = r.influence.groupby("station_id")["loo_mahalanobis"].median()
    assert loo.idxmax() == "madagascar"


def test_leave_one_out_error_against_truth(runs):  # 36
    # Excluding the biased station reduces the error against the truth; excluding a healthy one does not.
    import pandas as pd

    from tart_gnss.stage7.synthetic import truth_ecef

    from .conftest import scenarios

    sc = scenarios()["bias"]
    r = runs("bias")

    def secs(col):
        return (col - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().to_numpy()

    pos = pd.read_parquet(r.chain["stage4"].parent / "positions.parquet")
    full_err = np.median([np.linalg.norm(np.array(p) - truth_ecef(sc, t)) for p, t in
                          zip(pos[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy(), secs(pd.to_datetime(
                              pos["time_utc"], utc=True)))])
    errs = {}
    for sid, g in r.influence.groupby("station_id"):
        errs[sid] = np.median([np.linalg.norm(np.array(p) - truth_ecef(sc, t)) for p, t in
                               zip(g[["loo_x_ecef_m", "loo_y_ecef_m", "loo_z_ecef_m"]].to_numpy(), secs(g["time_utc"]))])
    assert errs["madagascar"] < 0.3 * full_err
    assert errs["mauritius"] > 0.7 * full_err and errs["south_africa"] > 0.7 * full_err


def test_calibration_correlated_bias(runs):  # 35.3
    r = runs("calibration")
    assert "CALIBRATION_CORRELATED" in r.tags and r.level >= 3
    best = r.summary["strongest_calibration_correlation"]
    assert best["station_id"] == "madagascar" and best["x_variable"] == "flagged_fraction"


def test_timing_correlated(runs):  # 35.4
    import numpy as np
    from scipy.stats import spearmanr

    r = runs("timing")
    assert "TIMING_CORRELATED" in r.tags
    t = r.correlations[r.correlations["interpretation_tag"] == "TIMING_CORRELATED"]
    assert len(t) and set(t["x_variable"]) <= {"station_time_offset_s", "abs_station_time_offset_s",
                                              "max_abs_station_time_offset_s", "mean_station_time_offset_s"}
    # Direction residuals are taken at each station's own (correctly recorded) time, so the offset shows in the
    # unmatched Stage-4 solution: the line residual of the late station grows with its offset.
    mu = r.influence[r.influence["station_id"] == "mauritius"]
    assert spearmanr(mu["station_time_offset_s"], mu["line_residual_m"])[0] > 0.6
    assert mu["direction_residual_angle_deg"].median() < 0.03  # a few sigma (the Stage-5 track itself is pulled)


def test_timing_signatures_in_rtn(runs):
    # A network-common time offset is purely along-track; one station's offset maps through the geometry.
    net = runs("network").diagnosis
    assert (net["dr_along_track_m"].abs() > 5 * net["dr_radial_m"].abs()).mean() > 0.9
    assert (net["dr_along_track_m"].abs() > 5 * net["dr_cross_track_m"].abs()).mean() > 0.9


def test_poor_geometry_is_not_network_consistent(runs):  # 35.5
    r = runs("geometry")
    assert "GEOMETRY_CORRELATED" in r.tags
    assert "NETWORK_CONSISTENT" not in r.tags and r.level < 4


def test_dynamics_mismatch(runs):  # 35.6
    r = runs("dynamics")
    assert "TRAJECTORY_MODEL_CORRELATED" in r.tags and r.level < 4
    assert set(r.summary["stage5_dynamics"]["run_flags"]) & {"DYNAMICS_MISMATCH", "HIGH_NIS"}


def test_stale_reference_is_limited(runs):  # 35.7
    r = runs("stale")
    assert "REFERENCE_LIMITED" in r.tags and r.level <= 3
    assert r.summary["precision_validation"] is False
    assert r.summary["overall_confidence"] != "HIGH"


def test_network_consistent_offset(runs):  # 35.8
    r = runs("network")
    assert "NETWORK_CONSISTENT" in r.tags and "STATION_SPECIFIC" not in r.tags
    assert r.influence["loo_mahalanobis"].median() < 5


def test_unexplained_requires_persistence_and_clean_diagnostics(runs, tmp_path):  # 35.9
    r = runs("network")
    assert r.level == 5 and "UNEXPLAINED" in r.tags and r.summary["unexplained_intervals"]
    # The same data with a persistence requirement longer than the arc stop at Level 4.
    from tart_gnss.stage7 import cli

    out = tmp_path / "s7"
    assert cli.main(["--stage6-manifest", str(r.chain["stage6"]), "--output-dir", str(out), "--no-plots",
                     "--min-persistent-epochs", "100"]) == 0
    import json

    s = json.loads(next(out.glob("*/*/summary.json")).read_text())
    assert s["overall_evidence_level"] == 4 and "UNEXPLAINED" not in s["overall_diagnostic_tags"]


def test_unknown_reference_uncertainty_blocks_level_5(runs):
    # A single TLE: the anomaly statistic treats it as exact, so the TLE cannot be ruled out as the explanation.
    r = runs("noensemble")
    assert "REFERENCE_UNCERTAINTY_UNKNOWN" in r.summary["quality_flags"]
    assert r.level == 4 and "UNEXPLAINED" not in r.tags and r.summary["overall_confidence"] != "HIGH"
    assert runs("network").summary["anomaly_statistic"] == "position_mahalanobis_total_cov"


def test_missing_calibration_metadata_is_not_evidence(runs):  # 33
    r = runs("nostage2")
    assert "CALIBRATION_METADATA_MISSING" in r.summary["quality_flags"]
    assert r.level == 4 and "UNEXPLAINED" not in r.tags
    assert r.summary["overall_confidence"] != "HIGH"
    assert "CALIBRATION_CORRELATED" not in r.tags


def test_random_metadata_is_not_promoted(runs):  # 37
    r = runs("random")
    assert r.level == 0
    assert not (r.tags & {"CALIBRATION_CORRELATED", "LOCALISATION_CORRELATED"})
    assert (r.correlations["n_samples"] > 0).any()


@pytest.mark.parametrize("name", ["clean", "bias", "network"])
def test_full_solution_reproduces_stage4(runs, name):
    assert runs(name).summary["full_solution_reproduction_max_m"] < 1e-3
