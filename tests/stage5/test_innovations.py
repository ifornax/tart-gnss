"""Innovation statistics and run-level consistency flags."""

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("pandas")

from tart_gnss.stage5 import synthetic as syn  # noqa: E402
from tart_gnss.stage5.innovations import autocorrelation, whitened  # noqa: E402

from .helpers import run_track, times, truth  # noqa: E402


def test_whitening_gives_unit_covariance():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(3, 3))
    s = a @ a.T + np.eye(3)
    nu = rng.multivariate_normal(np.zeros(3), s, size=20000)
    w = whitened(nu, np.repeat(s[None], len(nu), axis=0))
    assert np.allclose(np.cov(w.T), np.eye(3), atol=0.05)
    mask = np.ones(len(w), bool)
    assert abs(autocorrelation(w, mask)[1]) < 0.03


def test_consistent_model_gives_chi2_3_statistics():
    rng = np.random.default_rng(1)
    t = times(60)
    res = run_track(syn.stage4_rows(t, truth(t), 2e3, 500.0, rng))
    s = res.summary
    assert 2.0 < s["mean_nis"] < 4.2
    assert s["flags"] == []
    assert s["fraction_above_gate"] < 0.05


def test_underestimated_covariance_raises_high_nis():
    rng = np.random.default_rng(2)
    t = times(60)
    rows = syn.stage4_rows(t, truth(t), 2e3, 500.0, rng)
    for r in rows:  # the Stage-4 covariance claims 3x better precision than the noise has
        for c in ("cov_xx", "cov_xy", "cov_xz", "cov_yy", "cov_yz", "cov_zz"):
            r[c] /= 9.0
    res = run_track(rows)
    assert "HIGH_NIS" in res.summary["flags"]
    assert res.summary["mean_nis"] > 9.0


def test_wrong_dynamics_raise_dynamics_mismatch():
    # Constant velocity with the two-body process noise cannot follow an orbit: innovations correlate.
    rng = np.random.default_rng(3)
    t = times(60, step=20.0)
    res = run_track(syn.stage4_rows(t, truth(t), 200.0, 200.0, rng), dynamics="constant-velocity", sigma_a=1e-4)
    assert {"HIGH_NIS", "DYNAMICS_MISMATCH"} <= set(res.summary["flags"])


def test_isolated_outlier_does_not_raise_high_nis():
    rng = np.random.default_rng(4)
    t = times(60)
    res = run_track(syn.stage4_rows(t, truth(t), 2e3, 500.0, rng, outliers={30: np.array([0.0, 5e4, 0.0])}))
    s = res.summary
    assert s["n_above_gate"] == 1 and s["mean_nis"] > 10.0  # the reported mean includes the outlier
    assert "HIGH_NIS" not in s["flags"]
