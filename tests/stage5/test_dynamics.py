"""Dynamics: constants, two-body integration vs analytic Kepler, STM, process noise, J2."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage5 import dynamics as d  # noqa: E402
from tart_gnss.stage5.errors import DynamicsError  # noqa: E402
from tart_gnss.stage5.synthetic import example_initial_state, kepler_propagate  # noqa: E402

X0 = example_initial_state()


def test_constants_match_library_values():
    from sgp4.earth_gravity import wgs84
    from skyfield.constants import ERAD

    assert d.GM_EARTH_M3_S2 == pytest.approx(wgs84.mu * 1e9, rel=2e-7)
    assert d.R_EARTH_M == ERAD
    assert d.J2_EARTH == pytest.approx(wgs84.j2, rel=1e-5)


@pytest.mark.parametrize("dt", [10.0, -45.0, 3600.0])
def test_two_body_matches_analytic_kepler(dt):
    x, _ = d.propagate_two_body(X0, dt)
    k = kepler_propagate(X0, dt)
    assert np.allclose(x[:3], k[:3], atol=1e-3)
    assert np.allclose(x[3:], k[3:], atol=1e-6)


def test_conserves_energy_and_angular_momentum():
    x, _ = d.propagate_two_body(X0, 7200.0)

    def energy(s):
        return 0.5 * s[3:] @ s[3:] - d.GM_EARTH_M3_S2 / np.linalg.norm(s[:3])

    assert energy(x) == pytest.approx(energy(X0), rel=1e-11)
    assert np.allclose(np.cross(x[:3], x[3:]), np.cross(X0[:3], X0[3:]), rtol=1e-11)


def test_state_transition_matrix_matches_finite_differences():
    dt = 300.0
    _, phi = d.propagate_two_body(X0, dt)
    fd = np.zeros((6, 6))
    for j in range(6):
        h = np.zeros(6)
        h[j] = 1.0 if j < 3 else 1e-3
        fd[:, j] = (kepler_propagate(X0 + h, dt) - kepler_propagate(X0 - h, dt)) / (2 * h[j])
    assert np.max(np.abs(fd - phi)) / np.max(np.abs(phi)) < 1e-7


def test_zero_dt_is_identity():
    for dyn in d.DYNAMICS:
        x, f = d.propagate(X0, 0.0, dyn)
        assert np.array_equal(x, X0)
        assert np.array_equal(f, np.eye(6))


def test_constant_velocity_is_exact():
    x, f = d.propagate(X0, 12.5, "constant-velocity")
    assert np.allclose(x[:3], X0[:3] + 12.5 * X0[3:])
    assert np.allclose(x[3:], X0[3:])
    assert np.allclose(f[:3, 3:], 12.5 * np.eye(3))


def test_process_noise_formula_and_psd():
    q = d.process_noise(10.0, 2.0)
    assert q[0, 0] == pytest.approx(4.0 * 1e4 / 4)
    assert q[0, 3] == pytest.approx(4.0 * 1e3 / 2)
    assert q[3, 3] == pytest.approx(4.0 * 1e2)
    assert q[0, 1] == 0.0
    assert np.linalg.eigvalsh(q).min() >= -1e-9 * q.max()
    # Longer gaps give larger noise (gaps must widen the prediction).
    assert d.process_noise(100.0, 2.0)[0, 0] > q[0, 0]


def test_j2_magnitude_and_gradient():
    a = d.j2_acceleration(X0[:3])
    assert 1e-5 < np.linalg.norm(a) < 1e-4  # ~5e-5 m/s^2 at GNSS radius
    g = d.gravity_gradient(X0[:3], "j2")
    fd = np.zeros((3, 3))
    for j in range(3):
        h = np.zeros(3)
        h[j] = 10.0
        fd[:, j] = (d.acceleration(X0[:3] + h, "j2") - d.acceleration(X0[:3] - h, "j2")) / 20.0
    assert np.allclose(g, fd, rtol=1e-6, atol=1e-15)
    x2, _ = d.propagate(X0, 600.0, "j2")
    xk = kepler_propagate(X0, 600.0)
    assert 0.1 < np.linalg.norm(x2[:3] - xk[:3]) < 100.0  # J2 displaces by ~0.5 a t^2 ~ 10 m


def test_errors():
    with pytest.raises(DynamicsError):
        d.propagate(np.array([1.0, 0, 0, 0, 0, 0]), 10.0, "two-body")
    with pytest.raises(DynamicsError):
        d.propagate(X0, 10.0, "kepler")
