"""Backward Rauch-Tung-Striebel smoother over the forward-filter output.

For ``k = N-2 .. 0``, with the stored nonlinear prediction
``x_k+1|k = f(x_k|k)``, its covariance ``P_k+1|k`` and the flow Jacobian
``F_k+1``::

    C_k   = P_k|k F_k+1^T P_k+1|k^-1           (linear solve, no inverse)
    x_k|N = x_k|k + C_k (x_k+1|N - x_k+1|k)
    P_k|N = P_k|k + C_k (P_k+1|N - P_k+1|k) C_k^T

The smoothed state at every epoch uses all measurements of its filter
segment (the whole run unless the filter was re-initialised after a
divergence, see :mod:`tart_gnss.stage5.filter`); it is the primary Stage-5
product. Covariances are symmetrised and checked at every step.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tart_gnss.stage5.errors import FilterDivergenceError, SmootherError
from tart_gnss.stage5.filter import FilterResult, check_covariance


@dataclass
class SmoothedResult:
    """Smoothed states ``(N, 6)`` and covariances ``(N, 6, 6)`` and the smoother gains."""

    x_smooth: np.ndarray
    p_smooth: np.ndarray
    gains: np.ndarray


def rts_smooth(filtered: FilterResult) -> SmoothedResult:
    """Run backward Rauch-Tung-Striebel smoothing."""
    xf, pf = filtered.x_filt, filtered.p_filt
    n = len(xf)
    xs, ps = xf.copy(), pf.copy()
    gains = np.zeros((n, 6, 6))
    try:
        for k in range(n - 2, -1, -1):
            if filtered.reset[k + 1]:
                continue  # segment boundary (filter re-initialised): the segment ends at its filtered state
            f = filtered.transition[k + 1]
            pp = filtered.p_pred[k + 1]
            # C = P_f F^T P_pred^-1  <=>  P_pred C^T = F P_f
            c = np.linalg.solve(pp, f @ pf[k]).T
            xs[k] = xf[k] + c @ (xs[k + 1] - filtered.x_pred[k + 1])
            ps[k] = check_covariance(pf[k] + c @ (ps[k + 1] - pp) @ c.T, f"smoothed covariance at step {k}")
            gains[k] = c
    except (np.linalg.LinAlgError, FilterDivergenceError) as exc:
        raise SmootherError(str(exc)) from exc
    if not np.all(np.isfinite(xs)):
        raise SmootherError("non-finite smoothed state")
    return SmoothedResult(x_smooth=xs, p_smooth=ps, gains=gains)
