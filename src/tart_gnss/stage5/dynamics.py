"""Dynamical models, state transition matrices and process noise (GCRS, SI units).

State ``x = [x, y, z, vx, vy, vz]`` (m, m/s) in the GCRS (inertial) frame.

Models
    ``constant-velocity``  ``r(t+dt) = r + v dt``, ``v(t+dt) = v``; exact,
                           ``F = [[I, dt I], [0, I]]``. Diagnostic baseline.
    ``two-body``           ``r'' = -GM r / |r|^3`` (default, production).
    ``j2``                 two-body plus the Earth's J2 zonal term about the
                           true rotation pole (optional).

Orbital models are integrated with ``scipy.integrate.solve_ivp`` (DOP853,
``rtol = 1e-10``, ``atol = 1e-6``) together with the variational equations
``Phi' = A(t) Phi``, ``A = [[0, I], [G(r), 0]]``, so the state-transition
matrix ``Phi = d x(t+dt) / d x(t)`` used by the EKF and the RTS smoother is the
exact Jacobian of the integrated flow. ``G = d a / d r`` is analytic for two
body, ``G = GM / r^3 (3 r_hat r_hat^T - I)``; the J2 part is a central
difference of the analytic J2 acceleration.

Constants (IERS Conventions 2010, Table 1.1; GM also equals the WGS84 value of
NIMA TR8350.2): ``GM = 3.986004418e14 m^3 s^-2`` (geocentric, including the
atmosphere), ``J2 = 1.0826359e-3``, equatorial radius ``a_E = 6378136.6 m``
(identical to ``skyfield.constants.ERAD``). The tests cross-check GM against
``sgp4.earth_gravity.wgs84.mu`` (as a constant only; no TLE is involved).

Process noise
    Unmodelled accelerations are white acceleration noise of standard
    deviation ``sigma_a`` per axis (``--process-accel-sigma``, m/s^2), held
    constant over each interval ``dt`` (the discrete white-noise-acceleration
    form given in the Stage-5 specification)::

        Q(dt) = sigma_a^2 [[dt^4/4 I, dt^3/2 I], [dt^3/2 I, dt^2 I]]

    (as specified for constant velocity). The same ``Q`` is used for the
    orbital models: over the short TART intervals (seconds to minutes, much
    shorter than the ~12 h orbital period) the gravity gradient changes
    ``Q`` only at second order in ``dt``. ``Q`` grows with ``dt``, so long
    gaps automatically widen the predicted covariance.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from tart_gnss.stage5.errors import DynamicsError

#: Geocentric gravitational constant, m^3 s^-2 (IERS Conventions 2010 Table 1.1 = WGS84).
GM_EARTH_M3_S2 = 3.986004418e14
#: Second zonal harmonic (IERS Conventions 2010 Table 1.1).
J2_EARTH = 1.0826359e-3
#: Equatorial radius, m (IERS Conventions 2010 Table 1.1; = skyfield.constants.ERAD).
R_EARTH_M = 6378136.6

DYNAMICS = ("two-body", "constant-velocity", "j2")
RTOL = 1e-10
ATOL = 1e-6
INTEGRATOR = "scipy.integrate.solve_ivp DOP853"
_Z = np.array([0.0, 0.0, 1.0])
_J2_GRADIENT_STEP_M = 100.0

#: Default --process-accel-sigma per model (m/s^2), chosen from physics, never tuned against TLEs.
DEFAULT_PROCESS_ACCEL_SIGMA = {
    # Constant velocity leaves gravity itself unmodelled: GM/r^2 ~ 0.57 m/s^2 at GNSS radius.
    "constant-velocity": 1.0,
    # Two body leaves J2 (~5e-5 m/s^2 at GNSS radius), lunisolar (~5e-6) and radiation pressure (~1e-7).
    "two-body": 1e-4,
    # J2 included: lunisolar perturbations dominate the remainder.
    "j2": 1e-5,
}


def two_body_acceleration(r: np.ndarray) -> np.ndarray:
    """``-GM r / |r|^3`` (m/s^2)."""
    r = np.asarray(r, dtype=float)
    return -GM_EARTH_M3_S2 * r / np.linalg.norm(r) ** 3


def j2_acceleration(r: np.ndarray, pole: np.ndarray = _Z) -> np.ndarray:
    """J2 perturbing acceleration about the unit rotation axis ``pole`` (m/s^2).

    ``a = -3/2 J2 GM a_E^2 / r^5 [(1 - 5 z^2 / r^2) r + 2 z k]`` with ``z = r . k``.
    """
    r = np.asarray(r, dtype=float)
    k = np.asarray(pole, dtype=float)
    rn = np.linalg.norm(r)
    z = float(r @ k)
    c = -1.5 * J2_EARTH * GM_EARTH_M3_S2 * R_EARTH_M ** 2 / rn ** 5
    return c * ((1.0 - 5.0 * z * z / (rn * rn)) * r + 2.0 * z * k)


def acceleration(r: np.ndarray, dynamics: str, pole: Optional[np.ndarray] = None) -> np.ndarray:
    """Model acceleration at position ``r`` (zero for constant velocity)."""
    if dynamics == "constant-velocity":
        return np.zeros(3)
    a = two_body_acceleration(r)
    if dynamics == "j2":
        a = a + j2_acceleration(r, _Z if pole is None else pole)
    return a


def gravity_gradient(r: np.ndarray, dynamics: str, pole: Optional[np.ndarray] = None) -> np.ndarray:
    """``d a / d r`` (s^-2)."""
    r = np.asarray(r, dtype=float)
    if dynamics == "constant-velocity":
        return np.zeros((3, 3))
    rn = np.linalg.norm(r)
    u = r / rn
    g = GM_EARTH_M3_S2 / rn ** 3 * (3.0 * np.outer(u, u) - np.eye(3))
    if dynamics == "j2":
        k = _Z if pole is None else pole
        h = _J2_GRADIENT_STEP_M
        for j in range(3):
            d = np.zeros(3)
            d[j] = h
            g[:, j] += (j2_acceleration(r + d, k) - j2_acceleration(r - d, k)) / (2 * h)
    return g


def constant_velocity_transition(dt_s: float) -> np.ndarray:
    """``F = [[I, dt I], [0, I]]``."""
    f = np.eye(6)
    f[:3, 3:] = dt_s * np.eye(3)
    return f


def propagate_constant_velocity(state: np.ndarray, dt_s: float) -> Tuple[np.ndarray, np.ndarray]:
    """Constant-velocity propagation; return ``(state, F)``."""
    f = constant_velocity_transition(dt_s)
    return f @ np.asarray(state, dtype=float), f


def _orbital(state, dt_s: float, dynamics: str, pole, rtol: float, atol: float):
    from scipy.integrate import solve_ivp

    def rhs(_t, y):
        r, v = y[:3], y[3:6]
        phi = y[6:].reshape(6, 6)
        a_mat = np.zeros((6, 6))
        a_mat[:3, 3:] = np.eye(3)
        a_mat[3:, :3] = gravity_gradient(r, dynamics, pole)
        return np.concatenate([v, acceleration(r, dynamics, pole), (a_mat @ phi).ravel()])

    y0 = np.concatenate([np.asarray(state, dtype=float), np.eye(6).ravel()])
    if not np.all(np.isfinite(y0)) or np.linalg.norm(y0[:3]) < R_EARTH_M * 0.5:
        raise DynamicsError(f"cannot propagate state {np.asarray(state)} (non-finite or inside the Earth)")
    sol = solve_ivp(rhs, (0.0, float(dt_s)), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
        raise DynamicsError(f"ODE integration failed over dt = {dt_s} s: {sol.message}")
    y = sol.y[:, -1]
    return y[:6], y[6:].reshape(6, 6)


def propagate_two_body(state: np.ndarray, dt_s: float, rtol: float = RTOL,
                       atol: float = ATOL) -> Tuple[np.ndarray, np.ndarray]:
    """Propagate position/velocity under two-body gravity and return the state-transition Jacobian."""
    if dt_s == 0:
        return np.asarray(state, dtype=float).copy(), np.eye(6)
    return _orbital(state, dt_s, "two-body", None, rtol, atol)


def propagate(state: np.ndarray, dt_s: float, dynamics: str, pole: Optional[np.ndarray] = None,
              rtol: float = RTOL, atol: float = ATOL) -> Tuple[np.ndarray, np.ndarray]:
    """Propagate under ``dynamics`` (one of :data:`DYNAMICS`); return ``(state, Phi)``."""
    if dynamics == "constant-velocity":
        return propagate_constant_velocity(state, dt_s)
    if dynamics not in DYNAMICS:
        raise DynamicsError(f"unknown dynamics {dynamics!r}; choose from {DYNAMICS}")
    if dt_s == 0:
        return np.asarray(state, dtype=float).copy(), np.eye(6)
    return _orbital(state, dt_s, dynamics, pole, rtol, atol)


def process_noise(dt_s: float, sigma_a_m_s2: float) -> np.ndarray:
    """White-acceleration process noise ``Q(dt)`` (6x6, SI units; see the module docstring)."""
    dt = abs(float(dt_s))
    q = np.zeros((6, 6))
    i3 = np.eye(3)
    q[:3, :3] = dt ** 4 / 4.0 * i3
    q[:3, 3:] = q[3:, :3] = dt ** 3 / 2.0 * i3
    q[3:, 3:] = dt ** 2 * i3
    return sigma_a_m_s2 ** 2 * q
