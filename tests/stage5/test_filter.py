"""Forward EKF: update algebra, gating, irregular cadence, covariance checks."""

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("pandas")

from tart_gnss.stage5 import synthetic as syn  # noqa: E402
from tart_gnss.stage5.errors import FilterDivergenceError, InitialisationError  # noqa: E402
from tart_gnss.stage5.filter import (  # noqa: E402
    H,
    FilterConfig,
    TrajectoryMeasurement,
    check_covariance,
    ekf_filter,
    gate_threshold,
    kalman_update,
)

from .helpers import run_track, times, truth  # noqa: E402


def test_gate_threshold_value():
    assert gate_threshold(0.9973) == pytest.approx(14.156, abs=1e-3)


def test_joseph_update_equals_textbook_form():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(6, 6))
    p = a @ a.T + np.eye(6)
    r = np.diag([2.0, 3.0, 4.0])
    x = rng.normal(size=6)
    y = rng.normal(size=3)
    xn, pn, nu, s, nis = kalman_update(x, p, y, r)
    s_ref = H @ p @ H.T + r
    k = p @ H.T @ np.linalg.inv(s_ref)
    assert np.allclose(xn, x + k @ (y - H @ x))
    assert np.allclose(pn, (np.eye(6) - k @ H) @ p, atol=1e-10)
    assert nis == pytest.approx(nu @ np.linalg.inv(s_ref) @ nu)


def test_check_covariance():
    p = np.diag([1.0, 2, 3, 4, 5, 6])
    p[0, 1] = 1e-12
    out = check_covariance(p, "x")
    assert np.allclose(out, out.T)
    with pytest.raises(FilterDivergenceError):
        check_covariance(np.diag([1.0, -1, 1, 1, 1, 1]), "x")
    with pytest.raises(FilterDivergenceError):
        check_covariance(np.full((6, 6), np.nan), "x")


def test_noise_free_constant_velocity_is_recovered_exactly():
    t = np.arange(10) * 7.0
    x0 = np.array([2e7, 1e7, -1e7, 1000.0, 2000.0, -500.0])
    ms = [TrajectoryMeasurement(time_s=tk, position_gcrs_m=x0[:3] + x0[3:] * tk, covariance_gcrs_m2=np.eye(3) * 1e6)
          for tk in t]
    fr = ekf_filter(ms, FilterConfig(dynamics="constant-velocity", process_accel_sigma_m_s2=0.0))
    assert np.allclose(fr.x_filt[-1][:3], x0[:3] + x0[3:] * t[-1], atol=1e-3)
    assert np.allclose(fr.x_filt[-1][3:], x0[3:], atol=1e-5)
    assert fr.accepted.all()


def test_first_measurement_must_be_usable():
    m = TrajectoryMeasurement(time_s=0.0, position_gcrs_m=np.zeros(3), covariance_gcrs_m2=np.eye(3), used=False)
    with pytest.raises(InitialisationError):
        ekf_filter([m], FilterConfig())


def test_irregular_cadence_uses_actual_times():
    # Jittered epochs (+-4 s on a 10 s grid): the filter must use the true dt. Feeding it the nominal
    # grid instead (as a uniform-cadence implementation would) is visibly inconsistent.
    rng = np.random.default_rng(3)
    t = times(40, jitter=4.0, rng=rng)
    x = truth(t)
    rows = syn.stage4_rows(t, x, 300.0, 300.0, rng)
    good = run_track(rows)
    assert good.summary["mean_nis"] < 5.0
    nominal = times(40)
    for r, tn in zip(rows, nominal):
        from datetime import datetime, timezone

        r["time_utc"] = datetime.fromtimestamp(tn, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    bad = run_track(rows)
    assert bad.summary["mean_nis"] > 10 * good.summary["mean_nis"]
    err = good.smoothed.x_smooth - x
    assert np.sqrt(np.mean(np.sum(err[:, :3] ** 2, axis=1))) < 300.0


def test_robust_inflate_assimilates_outlier_with_inflated_covariance():
    rng = np.random.default_rng(4)
    t = times(30)
    rows = syn.stage4_rows(t, truth(t), 1e3, 1e3, rng, outliers={20: np.array([3e4, 0, 0])})
    none = run_track(rows)
    infl = run_track(rows, robust_update="inflate")
    j = 20
    assert none.filtered.outlier[j] and not none.filtered.accepted[j]
    assert infl.filtered.outlier[j] and infl.filtered.accepted[j]
    assert infl.filtered.inflation[j] > 1.0


def test_divergence_recovery_reinitialises_instead_of_locking_out():
    # After 8 epochs the measured track jumps by 300 km (e.g. a different source): with hard gating and
    # near-deterministic dynamics every later measurement would be gated for ever. The filter must
    # restart from the data after three consecutive gated measurements, and the smoother must not
    # smooth across the restart.
    rng = np.random.default_rng(5)
    t = times(30)
    x = truth(t)
    jump = {k: np.array([3e5, 0.0, 0.0]) for k in range(8, 30)}
    res = run_track(syn.stage4_rows(t, x, 1e3, 1e3, rng, outliers=jump))
    fr, sm = res.filtered, res.smoothed
    assert list(np.flatnonzero(fr.reset)) == [8]
    assert list(np.flatnonzero(fr.diverged)) == [8, 9, 10]
    assert fr.accepted[8:].all() and len(fr.initials) == 2
    assert all("FILTER_DIVERGENCE" in res.trajectory[k]["quality_flags"] for k in (8, 9, 10))
    assert res.trajectory[8]["valid"]
    m = np.array([res.measurements.matrices[k] for k in range(30)])
    ecef = np.einsum("nji,nj->ni", m, sm.x_smooth[:, :3])
    true_ecef = np.einsum("nji,nj->ni", m, x[:, :3])
    assert np.all(np.linalg.norm(ecef[:8] - true_ecef[:8], axis=1) < 3e4)
    assert np.all(np.linalg.norm(ecef[8:] - true_ecef[8:] - [3e5, 0, 0], axis=1) < 3e4)
    assert res.summary["n_resets"] == 1
