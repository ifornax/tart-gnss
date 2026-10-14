"""observed - TLE residuals: sign, RTN projection, Mahalanobis, angles, time-shift physics."""

import numpy as np
import pytest

pytest.importorskip("skyfield")

from tart_gnss.stage5.dynamics import GM_EARTH_M3_S2  # noqa: E402
from tart_gnss.stage6.residuals import (  # noqa: E402
    StateTrajectory,
    angular_separation_deg,
    compute_residuals,
    mahalanobis2,
)

from .helpers import inject_rtn, times, truth  # noqa: E402

T = times(10)
X = truth(T)


def _res(obs, cov=None, disp=None):
    n = len(T)
    c = np.repeat(np.diag([1e6] * 3 + [1.0] * 3)[None], n, axis=0) if cov is None else cov
    return compute_residuals(StateTrajectory(T, obs, c), StateTrajectory(T, X), disp)


def test_sign_is_observed_minus_tle():  # 34.2
    off = np.array([1500.0, -700.0, 250.0])
    obs = X.copy()
    obs[:, :3] += off
    r = _res(obs)
    assert np.allclose(r.dr, off, atol=1e-6)
    assert np.allclose(r.position_norm_m, np.linalg.norm(off))


def test_rtn_injection_recovered():  # 34.3
    r = _res(inject_rtn(X, (1000.0, 2000.0, -500.0)))
    assert np.allclose(r.dr_rtn, [1000.0, 2000.0, -500.0], atol=1e-6)
    assert np.allclose(r.range_m, 1000.0, atol=1.0)  # second-order terms only
    assert np.allclose(r.along_track_time_s, 2000.0 / np.linalg.norm(X[:, 3:], axis=1))


def test_velocity_injection_recovered():  # 34.4
    r = _res(inject_rtn(X, dv_rtn=(0.3, -1.2, 0.8)))
    assert np.allclose(r.dv_rtn, [0.3, -1.2, 0.8], atol=1e-9)
    assert np.allclose(r.dr, 0.0, atol=1e-6)


def test_time_shift_signature():
    # Observed = SGP4 at t + tau: along-track |v| tau, radial velocity -GM/r^2 tau (gravity), time offset tau.
    tau = 2.0
    r = _res(truth(T + tau))
    speed = np.linalg.norm(X[:, 3:], axis=1)
    radius = np.linalg.norm(X[:, :3], axis=1)
    assert np.allclose(r.dr_rtn[:, 1], speed * tau, rtol=2e-3)
    assert np.allclose(r.along_track_time_s, tau, rtol=2e-3)
    assert np.allclose(r.dv_rtn[:, 0], -GM_EARTH_M3_S2 / radius ** 2 * tau, rtol=0.05)
    assert np.all(np.abs(r.dr_rtn[:, 2]) < 1.0)


def test_mahalanobis_values():
    obs = X.copy()
    obs[:, 0] += 3000.0
    cov = np.repeat(np.diag([1e6, 4e6, 9e6, 1, 1, 1.0])[None], len(T), axis=0)
    r = _res(obs, cov)
    assert np.allclose(r.d2_pos_obs, 9.0)
    assert np.allclose(r.d2_state_obs, 9.0)
    assert np.all(np.isnan(r.d2_pos_total))
    disp = np.repeat(np.diag([3e6, 0, 0, 0, 0, 0.0])[None], len(T), axis=0)
    r2 = _res(obs, cov, disp)
    assert np.allclose(r2.d2_pos_total, 9e6 / 4e6)


def test_unusable_covariance_gives_nan_distance():
    bad = np.repeat(np.diag([1e6, 1e6, -1.0, 1, 1, 1])[None], len(T), axis=0)
    r = _res(inject_rtn(X, (10.0, 0, 0)), bad)
    assert not r.covariance_ok.any() and np.all(np.isnan(r.d2_pos_obs))
    assert mahalanobis2(np.ones(3), np.diag([1.0, 1.0, 1e-20])) != mahalanobis2(np.ones(3), np.eye(3))


def test_angular_separation_small_angles():
    a = np.array([[2.6e7, 0, 0]])
    b = np.array([[2.6e7, 26.0, 0]])  # 1e-6 rad
    assert angular_separation_deg(a, b)[0] == pytest.approx(np.degrees(1e-6), rel=1e-6)
