"""TLE-independent initial state: hierarchy, GLS covariance, prior inflation."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage5.errors import InitialisationError  # noqa: E402
from tart_gnss.stage5.initialisation import initial_state  # noqa: E402

R0 = np.array([1e7, 2e7, -1e7])
V0 = np.array([1000.0, -2000.0, 3000.0])
C = np.diag([1e6, 4e6, 9e6])


def _pts(ts):
    return [R0 + V0 * t for t in ts], [C] * len(ts)


def test_linear_fit_recovers_velocity_exactly():
    ts = [0.0, 10.0, 20.0, 30.0]
    y, c = _pts(ts)
    init = initial_state(ts, y, c, window_s=30.0, dynamics="constant-velocity")
    assert init.method == "linear_fit" and init.n_points == 4
    assert np.allclose(init.state, np.concatenate([R0, V0]))
    # GLS covariance of the slope for equal weights: sigma^2 / sum (t - tbar)^2 per axis.
    assert init.fit_covariance[3, 3] == pytest.approx(1e6 / 500.0)
    assert np.allclose(init.prior_covariance, 1e4 * init.fit_covariance)
    assert not init.weak


def test_window_excludes_later_points():
    ts = [0.0, 10.0, 20.0, 100.0]
    y, c = _pts(ts)
    init = initial_state(ts, y, c, window_s=25.0, dynamics="constant-velocity")
    assert init.n_points == 3


def test_finite_difference_fallback_and_weak_flag():
    ts = [0.0, 40.0, 80.0]
    y, c = _pts(ts)
    init = initial_state(ts, y, c, window_s=30.0, dynamics="constant-velocity")
    assert init.method == "finite_difference" and init.weak
    assert np.allclose(init.state[3:], V0)
    assert init.fit_covariance[3, 3] == pytest.approx(2e6 / 40.0 ** 2)


def test_zero_velocity_last_resort():
    y, c = _pts([0.0])
    init = initial_state([0.0], y, c)
    assert init.method == "zero_velocity" and init.weak
    assert np.allclose(init.state[3:], 0.0)
    assert np.sqrt(init.fit_covariance[3, 3]) == pytest.approx(1e4)


def test_large_velocity_sigma_is_weak():
    ts = [0.0, 1.0, 2.0]
    y, c = _pts(ts)
    init = initial_state(ts, y, [C * 1e8] * 3, window_s=30.0, dynamics="constant-velocity")
    assert init.method == "linear_fit" and init.weak


def test_gravity_correction_for_two_body():
    from tart_gnss.stage5.synthetic import example_initial_state, kepler_propagate

    x0 = example_initial_state()
    ts = [0.0, 10.0, 20.0, 30.0]
    y = [kepler_propagate(x0, t)[:3] for t in ts]
    init = initial_state(ts, y, [C] * 4, dynamics="two-body")
    assert np.linalg.norm(init.state[3:] - x0[3:]) < 0.05  # gravity known, not fitted
    init_cv = initial_state(ts, y, [C] * 4, dynamics="constant-velocity")
    assert np.linalg.norm(init_cv.state[3:] - x0[3:]) > 1.0


def test_no_measurement_raises():
    with pytest.raises(InitialisationError):
        initial_state([], [], [])
