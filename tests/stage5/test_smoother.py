"""RTS smoother: reduces uncertainty, ends at the filter, equals the batch solution in the linear case."""

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("pandas")

from tart_gnss.stage5.filter import FilterConfig, TrajectoryMeasurement, ekf_filter  # noqa: E402
from tart_gnss.stage5.smoother import rts_smooth  # noqa: E402


def _cv_measurements(rng, n=15):
    t = np.sort(rng.uniform(0, 300, n))
    t[0] = 0.0
    x0 = np.array([2e7, 1e7, -1e7, 1000.0, 2000.0, -500.0])
    out = []
    for tk in t:
        a = rng.normal(size=(3, 3))
        c = a @ a.T * 1e5 + np.eye(3) * 1e5
        y = x0[:3] + x0[3:] * tk + rng.multivariate_normal(np.zeros(3), c)
        out.append(TrajectoryMeasurement(time_s=tk, position_gcrs_m=y, covariance_gcrs_m2=c))
    return out


def test_smoothed_equals_batch_generalised_least_squares():
    # Linear dynamics, no process noise: the smoother at t0 must equal the batch GLS fit of all data
    # (up to the 1e-4 information carried by the inflated prior).
    rng = np.random.default_rng(0)
    ms = _cv_measurements(rng)
    fr = ekf_filter(ms, FilterConfig(dynamics="constant-velocity", process_accel_sigma_m_s2=0.0,
                                     prior_inflation=1e8))
    sm = rts_smooth(fr)
    n_mat, b = np.zeros((6, 6)), np.zeros(6)
    for m in ms:
        h = np.hstack([np.eye(3), m.time_s * np.eye(3)])
        w = np.linalg.inv(m.covariance_gcrs_m2)
        n_mat += h.T @ w @ h
        b += h.T @ w @ m.position_gcrs_m
    x_batch = np.linalg.solve(n_mat, b)
    assert np.allclose(sm.x_smooth[0], x_batch, rtol=0, atol=np.sqrt(np.diag(np.linalg.inv(n_mat))) * 1e-3)
    assert np.allclose(sm.p_smooth[0], np.linalg.inv(n_mat), rtol=1e-4)


def test_smoother_reduces_uncertainty_and_ends_at_filter():
    rng = np.random.default_rng(1)
    ms = _cv_measurements(rng)
    fr = ekf_filter(ms, FilterConfig(dynamics="constant-velocity", process_accel_sigma_m_s2=0.01))
    sm = rts_smooth(fr)
    assert np.allclose(sm.x_smooth[-1], fr.x_filt[-1])
    assert np.allclose(sm.p_smooth[-1], fr.p_filt[-1])
    for k in range(len(ms)):
        assert np.trace(sm.p_smooth[k]) <= np.trace(fr.p_filt[k]) * (1 + 1e-9)
        assert np.linalg.eigvalsh(sm.p_smooth[k]).min() > 0
