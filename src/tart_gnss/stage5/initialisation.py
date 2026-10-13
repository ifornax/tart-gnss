"""TLE-independent initial state from the first Stage-4 measurements.

Only Stage-4 positions and covariances enter. With the reference epoch ``t0``
of the first usable measurement, the measurements within
``--initial-velocity-window`` seconds are fitted by generalised least squares
with their full 3x3 covariances::

    y_k = r0 + v0 (t_k - t0) + 1/2 a(r0) (t_k - t0)^2 + noise

where ``a`` is the model acceleration at ``y_0`` (zero for constant velocity;
the gravity term is a known correction, not a fitted parameter). The fit gives
the 6x6 covariance ``N^-1`` of ``(r0, v0)``.

Hierarchy (``method``):
    ``linear_fit``         >= 3 measurements in the window
    ``finite_difference``  2 measurements (in the window, else the first two
                           measurements); the GLS fit through two points is
                           the covariance-weighted finite difference
    ``zero_velocity``      a single measurement: velocity zero with
                           ``LAST_RESORT_VELOCITY_SIGMA_M_S`` per axis

``INITIAL_VELOCITY_WEAK`` is raised for anything but ``linear_fit``, or when
any velocity sigma exceeds ``weak_velocity_sigma_m_s``.

Prior for the filter
    The window measurements are assimilated again by the filter, so using the
    fit covariance ``N^-1`` itself as the prior would count them twice and
    make the filter over-confident. The prior is therefore
    ``P0 = prior_inflation * N^-1`` (default ``1e4``): the fit fixes the
    linearisation point and the *shape* of the initial uncertainty (including
    the position-velocity correlation) but contributes at most
    ``1 / prior_inflation`` of the information of the window data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np

from tart_gnss.stage5.dynamics import acceleration
from tart_gnss.stage5.errors import InitialisationError

LAST_RESORT_VELOCITY_SIGMA_M_S = 1.0e4
DEFAULT_PRIOR_INFLATION = 1.0e4
DEFAULT_WEAK_VELOCITY_SIGMA_M_S = 1000.0
METHODS = ("linear_fit", "finite_difference", "zero_velocity")


@dataclass
class InitialState:
    """Initial GCRS state at ``time_s``, its fit covariance and the inflated filter prior."""

    time_s: float
    state: np.ndarray
    fit_covariance: np.ndarray
    prior_covariance: np.ndarray
    method: str
    n_points: int
    window_s: float
    weak: bool
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """JSON-serialisable summary."""
        sv = np.sqrt(np.clip(np.diag(self.fit_covariance), 0, None))
        return {
            "method": self.method, "n_points": self.n_points, "window_s": self.window_s,
            "reference_time_unix_s": self.time_s,
            "position_gcrs_m": [float(v) for v in self.state[:3]],
            "velocity_gcrs_m_s": [float(v) for v in self.state[3:]],
            "fit_sigma_position_m": [float(v) for v in sv[:3]],
            "fit_sigma_velocity_m_s": [float(v) for v in sv[3:]],
            "weak": self.weak, "notes": list(self.notes),
        }


def _gls(times, positions, covariances, t0, a0) -> tuple:
    """GLS fit of ``(r0, v0)``; return ``(state, covariance)``."""
    n_mat = np.zeros((6, 6))
    b = np.zeros(6)
    for t, y, c in zip(times, positions, covariances):
        dt = t - t0
        h = np.hstack([np.eye(3), dt * np.eye(3)])
        w = np.linalg.inv(c)  # 3x3 weight from a validated, positive-definite covariance
        yc = y - 0.5 * a0 * dt * dt
        n_mat += h.T @ w @ h
        b += h.T @ w @ yc
    try:
        cov = np.linalg.inv(n_mat)
        x = np.linalg.solve(n_mat, b)
    except np.linalg.LinAlgError as exc:
        raise InitialisationError(f"singular initial-state normal matrix: {exc}") from exc
    return x, 0.5 * (cov + cov.T)


def initial_state(
    times_s: Sequence[float],
    positions_gcrs_m: Sequence[np.ndarray],
    covariances_gcrs_m2: Sequence[np.ndarray],
    window_s: float = 30.0,
    dynamics: str = "two-body",
    prior_inflation: float = DEFAULT_PRIOR_INFLATION,
    weak_velocity_sigma_m_s: float = DEFAULT_WEAK_VELOCITY_SIGMA_M_S,
    pole: Optional[np.ndarray] = None,
) -> InitialState:
    """Initial state from usable measurements (time ordered, GCRS); see the module docstring."""
    t = np.asarray(times_s, dtype=float)
    y = np.asarray(positions_gcrs_m, dtype=float)
    c = np.asarray(covariances_gcrs_m2, dtype=float)
    if len(t) == 0:
        raise InitialisationError("no usable Stage-4 measurement to initialise from")
    if prior_inflation < 1:
        raise InitialisationError("prior_inflation must be >= 1")
    t0 = float(t[0])
    a0 = acceleration(y[0], dynamics, pole)
    notes: List[str] = []
    distinct = np.concatenate([[True], np.diff(t) > 0])
    t, y, c = t[distinct], y[distinct], c[distinct]
    in_win = t - t0 <= window_s
    if in_win.sum() >= 3:
        method, sel = "linear_fit", in_win
    elif len(t) >= 2:
        method = "finite_difference"
        sel = in_win if in_win.sum() == 2 else np.arange(len(t)) < 2
        if in_win.sum() < 2:
            notes.append(f"fewer than 2 measurements within {window_s} s; first two used "
                         f"(span {t[1] - t0:.1f} s)")
    else:
        method, sel = "zero_velocity", None
    if sel is not None:
        x, cov = _gls(t[sel], y[sel], c[sel], t0, a0)
        n = int(sel.sum())
    else:
        x = np.concatenate([y[0], np.zeros(3)])
        cov = np.zeros((6, 6))
        cov[:3, :3] = c[0]
        cov[3:, 3:] = LAST_RESORT_VELOCITY_SIGMA_M_S ** 2 * np.eye(3)
        n = 1
        notes.append("single measurement: zero velocity with "
                     f"{LAST_RESORT_VELOCITY_SIGMA_M_S:g} m/s sigma per axis (last resort)")
    sig_v = np.sqrt(np.clip(np.diag(cov)[3:], 0, None))
    weak = method != "linear_fit" or bool(np.max(sig_v) > weak_velocity_sigma_m_s)
    if weak and method == "linear_fit":
        notes.append(f"velocity sigma {np.max(sig_v):.3g} m/s exceeds {weak_velocity_sigma_m_s:g} m/s")
    return InitialState(time_s=t0, state=x, fit_covariance=cov, prior_covariance=prior_inflation * cov,
                        method=method, n_points=n, window_s=float(window_s), weak=weak, notes=notes)
