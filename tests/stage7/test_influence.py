"""Leave-one-out and pairwise re-triangulation, direction residuals."""

import numpy as np
import pytest

pytest.importorskip("pandas")

from tart_gnss.stage7.influence import loo_d2, pair_diagnostic  # noqa: E402


def test_loo_d2_rank_and_value():
    c_full = np.diag([1.0, 1.0, 1.0])
    c_loo = np.diag([5.0, 1.0, 1.0])  # removing the station added variance in x only
    d, k = loo_d2(np.array([4.0, 0.0, 0.0]), c_loo, c_full)
    assert k == 1 and d == pytest.approx(4.0)
    assert np.isnan(loo_d2(np.ones(3), np.diag([0.5, 1, 1.0]), c_full)[0])


def test_pair_separation_and_z():
    r1, r2 = np.zeros(3), np.array([1e6, 0, 0])
    s1 = np.array([0.0, 0.0, 1.0])
    target = np.array([5e5, 0, 2e7])
    s2 = (target - r2) / np.linalg.norm(target - r2)
    s1 = (target + np.array([0, 3e3, 0]) - r1)
    s1 /= np.linalg.norm(s1)
    c = np.eye(3) * (1e-4) ** 2
    out = pair_diagnostic(("a", r1, s1, c), ("b", r2, s2, c))
    assert out["pair_line_separation_m"] == pytest.approx(3e3, rel=0.02)
    assert out["pair_separation_z"] == pytest.approx(3e3 / np.sqrt(2 * (2e7 * 1e-4) ** 2), rel=0.05)


def test_biased_station_shift_and_degraded_geometry(runs):
    r = runs("bias")
    inf = r.influence
    assert (inf["n_stations_remaining"] == 2).all() and inf["degraded_geometry"].all()
    assert (inf["loo_dof"] >= 1).all()
    med = inf.groupby("station_id")["direction_residual_angle_deg"].median()
    assert med.idxmax() == "madagascar"
    pairs = r.pairs
    bad = pairs[pairs["pair_separation_z"] > 4]
    assert len(bad) and bad.apply(lambda x: "madagascar" in (x.pair_station_1, x.pair_station_2), axis=1).all()


def test_clean_direction_residuals_are_small(runs):
    inf = runs("clean").influence
    assert inf["direction_residual_angle_deg"].median() < 0.01  # sigma 0.005 deg
    assert inf["loo_mahalanobis"].median() < 3
