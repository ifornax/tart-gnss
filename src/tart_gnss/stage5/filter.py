"""Forward extended Kalman filter over Stage-4 measurements (GCRS).

Measurement model ``y_k = H x_k + v_k``, ``H = [I_3 0]``, ``v_k ~ N(0, R_k)``
with ``R_k`` the full Stage-4 covariance rotated into GCRS. For every epoch
(all Stage-4 epochs from the first usable one on, in time order):

    prediction   x_k|k-1 = f(x_k-1|k-1, dt_k),  P_k|k-1 = F_k P F_k^T + Q(dt_k)
    innovation   nu_k = y_k - H x_k|k-1,        S_k = H P_k|k-1 H^T + R_k
    NIS          D_k^2 = nu_k^T S_k^-1 nu_k     (chi^2_3 if the model is right)
    gain         K_k = P_k|k-1 H^T S_k^-1       (Cholesky solve, no inverse)
    update       x_k|k = x_k|k-1 + K_k nu_k
                 P_k|k = (I - K H) P (I - K H)^T + K R K^T   (Joseph form)

``F_k`` is the exact Jacobian of the integrated flow (variational equations).
Only measurements accepted by the quality policy are assimilated; the others
(and policy-rejected epochs with a finite position) still get an innovation
for the diagnostics. A measurement with ``D_k^2`` above the chi^2_3 quantile
of ``--innovation-gate-probability`` is an ``INNOVATION_OUTLIER``: with
``--robust-update none`` (default) it is not assimilated; with ``inflate``
it is assimilated with ``R_k`` multiplied by ``D_k^2 / gate``, which brings
its NIS down to the gate. Gated measurements always stay in the output.

Divergence recovery: with near-deterministic orbital dynamics an early
unlucky gate decision can lock the filter onto a wrong state, after which the
predicted error and its sigma grow together and every later measurement is
gated too. After ``DIVERGENCE_RUN`` (3) consecutive gated measurements the
filter therefore restarts at the first of them, re-initialised from the data
(same TLE-independent initialisation, inflated prior); those rows are marked
``diverged`` (``FILTER_DIVERGENCE``) and the step is a segment boundary
(``reset``) that the smoother does not cross. ``*_first`` arrays keep the innovation and NIS of
each measurement's first (sequential, predictive) evaluation, before any
restart, so the consistency statistics still see the gated values.

After every prediction and update the covariance is symmetrised and checked
(finite, smallest eigenvalue not below ``-1e-9`` times the largest); a failure
raises :class:`FilterDivergenceError` instead of continuing silently.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from tart_gnss.stage5.dynamics import DEFAULT_PROCESS_ACCEL_SIGMA, propagate, process_noise
from tart_gnss.stage5.errors import FilterDivergenceError, InitialisationError
from tart_gnss.stage5.initialisation import (
    DEFAULT_PRIOR_INFLATION,
    DEFAULT_WEAK_VELOCITY_SIGMA_M_S,
    InitialState,
    initial_state,
)

H = np.hstack([np.eye(3), np.zeros((3, 3))])
ROBUST_UPDATES = ("none", "inflate")
EIG_TOL = 1e-9
#: Consecutive gated measurements that trigger a re-initialisation (FILTER_DIVERGENCE).
DIVERGENCE_RUN = 3


@dataclass
class TrajectoryMeasurement:
    """One Stage-4 epoch as a GCRS position measurement.

    ``used`` is the quality-policy decision; ``reason`` says why not.
    """

    time_s: float
    position_gcrs_m: np.ndarray
    covariance_gcrs_m2: np.ndarray
    used: bool = True
    reason: str = ""
    index: int = -1
    stage4_flags: List[str] = field(default_factory=list)
    covariance_problem: str = ""

    @property
    def finite(self) -> bool:
        """Position and covariance are finite."""
        return bool(np.all(np.isfinite(self.position_gcrs_m)) and np.all(np.isfinite(self.covariance_gcrs_m2)))


@dataclass
class FilterConfig:
    """Forward-filter settings (all recorded in the Stage-5 manifest)."""

    dynamics: str = "two-body"
    process_accel_sigma_m_s2: Optional[float] = None
    gate_probability: float = 0.9973
    robust_update: str = "none"
    max_gap_s: float = 120.0
    initial_velocity_window_s: float = 30.0
    prior_inflation: float = DEFAULT_PRIOR_INFLATION
    weak_velocity_sigma_m_s: float = DEFAULT_WEAK_VELOCITY_SIGMA_M_S
    pole_gcrs: Optional[np.ndarray] = None

    @property
    def sigma_a(self) -> float:
        """Process acceleration sigma (the per-model default if unset)."""
        if self.process_accel_sigma_m_s2 is not None:
            return float(self.process_accel_sigma_m_s2)
        return DEFAULT_PROCESS_ACCEL_SIGMA[self.dynamics]


@dataclass
class FilterResult:
    """Per-step forward-filter output; step ``k`` is ``measurements[k]``."""

    measurements: List[TrajectoryMeasurement]
    config: FilterConfig
    initial: InitialState
    gate_threshold: float
    times_s: np.ndarray
    x_pred: np.ndarray
    p_pred: np.ndarray
    x_filt: np.ndarray
    p_filt: np.ndarray
    transition: np.ndarray
    process: np.ndarray
    innovation: np.ndarray
    innovation_cov: np.ndarray
    nis: np.ndarray
    attempted: np.ndarray
    accepted: np.ndarray
    outlier: np.ndarray
    inflation: np.ndarray
    gap_s: np.ndarray
    reset: np.ndarray
    diverged: np.ndarray
    nis_first: np.ndarray
    outlier_first: np.ndarray
    innovation_first: np.ndarray
    innovation_cov_first: np.ndarray
    initials: List[InitialState]


def gate_threshold(probability: float, dof: int = 3) -> float:
    """chi^2 quantile used as the innovation gate."""
    from scipy.stats import chi2

    return float(chi2.ppf(probability, dof))


def check_covariance(p: np.ndarray, what: str) -> np.ndarray:
    """Symmetrise ``p`` and verify it is finite and positive semidefinite (within tolerance)."""
    p = 0.5 * (p + p.T)
    if not np.all(np.isfinite(p)):
        raise FilterDivergenceError(f"{what}: non-finite covariance")
    w = np.linalg.eigvalsh(p)
    if w[-1] <= 0 or w[0] < -EIG_TOL * w[-1]:
        raise FilterDivergenceError(f"{what}: covariance not positive semidefinite (eigenvalues {w[0]:.3g} .. "
                                    f"{w[-1]:.3g})")
    return p


def solve_spd(s: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``S^-1 b`` for symmetric positive-definite ``S`` (Cholesky; LU fallback)."""
    from scipy.linalg import LinAlgError, cho_factor, cho_solve

    try:
        return cho_solve(cho_factor(s, lower=True), b)
    except LinAlgError:
        return np.linalg.solve(s, b)


def kalman_update(x, p, y, r):
    """One Joseph-form position update; return ``(x, P, nu, S, nis)``."""
    nu = y - H @ x
    s = H @ p @ H.T + r
    s = 0.5 * (s + s.T)
    k = solve_spd(s, H @ p).T  # K = P H^T S^-1, S and P symmetric
    nis = float(nu @ solve_spd(s, nu))
    a = np.eye(6) - k @ H
    p_new = a @ p @ a.T + k @ r @ k.T
    return x + k @ nu, p_new, nu, s, nis


def ekf_filter(measurements: List[TrajectoryMeasurement], config: FilterConfig) -> FilterResult:
    """Run forward extended Kalman filtering over Stage-4 measurements.

    ``measurements`` must be time ordered and start with a usable one (see
    :func:`tart_gnss.stage5.cli`); the initial state comes from the usable
    measurements only.
    """
    if config.robust_update not in ROBUST_UPDATES:
        raise ValueError(f"robust_update must be one of {ROBUST_UPDATES}")
    if not measurements or not (measurements[0].used and measurements[0].finite):
        raise InitialisationError("the first measurement passed to the filter must be usable")
    gate = gate_threshold(config.gate_probability)
    n = len(measurements)
    out = {
        "x_pred": np.full((n, 6), np.nan), "p_pred": np.full((n, 6, 6), np.nan),
        "x_filt": np.full((n, 6), np.nan), "p_filt": np.full((n, 6, 6), np.nan),
        "transition": np.full((n, 6, 6), np.nan), "process": np.full((n, 6, 6), np.nan),
        "innovation": np.full((n, 3), np.nan), "innovation_cov": np.full((n, 3, 3), np.nan),
        "nis": np.full(n, np.nan), "attempted": np.zeros(n, bool), "accepted": np.zeros(n, bool),
        "outlier": np.zeros(n, bool), "inflation": np.ones(n), "gap_s": np.zeros(n),
        "reset": np.zeros(n, bool), "diverged": np.zeros(n, bool),
        "nis_first": np.full(n, np.nan), "outlier_first": np.zeros(n, bool),
        "innovation_first": np.full((n, 3), np.nan), "innovation_cov_first": np.full((n, 3, 3), np.nan),
    }
    times = np.array([m.time_s for m in measurements], dtype=float)
    if np.any(np.diff(times) < 0):
        raise ValueError("measurements must be time ordered")
    initials: List[InitialState] = []
    gated_run: List[int] = []
    last_update = times[0]
    x = p = None
    k = 0
    while k < n:
        meas = measurements[k]
        if k == 0 or out["reset"][k]:
            init = _initialise(measurements[k:], config)
            initials.append(init)
            x, p = init.state.copy(), check_covariance(init.prior_covariance, f"prior at step {k}")
            f, q = np.eye(6), np.zeros((6, 6))
            gated_run = []
        else:
            dt = times[k] - times[k - 1]
            x, f = propagate(x, dt, config.dynamics, config.pole_gcrs)
            q = process_noise(dt, config.sigma_a)
            p = check_covariance(f @ p @ f.T + q, f"prediction at step {k}")
        out["x_pred"][k], out["p_pred"][k], out["transition"][k], out["process"][k] = x, p, f, q
        out["gap_s"][k] = times[k] - last_update
        out["outlier"][k] = out["accepted"][k] = False
        out["inflation"][k] = 1.0
        if meas.finite:
            r = meas.covariance_gcrs_m2
            x_new, p_new, nu, s, nis = kalman_update(x, p, meas.position_gcrs_m, r)
            out["innovation"][k], out["innovation_cov"][k], out["nis"][k] = nu, s, nis
            if not np.isfinite(out["nis_first"][k]):  # first (sequential, predictive) evaluation
                out["nis_first"][k] = nis
                out["innovation_first"][k], out["innovation_cov_first"][k] = nu, s
                out["outlier_first"][k] = bool(meas.used and nis > gate)
            if meas.used:
                out["attempted"][k] = True
                if nis > gate:
                    out["outlier"][k] = True
                    if config.robust_update == "inflate":
                        lam = nis / gate
                        out["inflation"][k] = lam
                        x_new, p_new, _, _, _ = kalman_update(x, p, meas.position_gcrs_m, lam * r)
                    else:
                        x_new = None
                        gated_run.append(k)
                if x_new is not None:
                    x, p = x_new, check_covariance(p_new, f"update at step {k}")
                    out["accepted"][k] = True
                    last_update = times[k]
                    gated_run = []
        out["x_filt"][k], out["p_filt"][k] = x, p
        if len(gated_run) >= DIVERGENCE_RUN:
            # The data disagree with the state, not the other way round: restart from the gated run.
            k = gated_run[0]
            out["reset"][k] = True
            out["diverged"][gated_run] = True
            gated_run = []
            continue
        k += 1
    return FilterResult(measurements=measurements, config=config, initial=initials[0], gate_threshold=gate,
                        times_s=times, initials=initials, **out)


def _initialise(measurements: List[TrajectoryMeasurement], config: FilterConfig) -> InitialState:
    usable = [m for m in measurements if m.used and m.finite]
    return initial_state([m.time_s for m in usable], [m.position_gcrs_m for m in usable],
                         [m.covariance_gcrs_m2 for m in usable], config.initial_velocity_window_s,
                         config.dynamics, config.prior_inflation, config.weak_velocity_sigma_m_s,
                         config.pole_gcrs)
