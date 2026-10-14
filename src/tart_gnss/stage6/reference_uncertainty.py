"""Empirical reference dispersion from an ensemble of nearby TLEs.

Each ensemble element set is propagated with SGP4 to every Stage-5 epoch and
transformed to GCRS exactly like the primary TLE. At each epoch the sample
mean and the sample covariance (``ddof = 1``) of the propagated 6-D states
give the **empirical reference dispersion**. This is not a formal SGP4 or
orbit-determination covariance, and it is labelled as such everywhere.

* At least :data:`MIN_ENSEMBLE` (3) members are required before the
  total-covariance statistics are produced; smaller ensembles are flagged
  ``TLE_ENSEMBLE_SMALL``.
* With fewer than 7 members the 6x6 sample covariance is rank deficient; it
  is only ever used as an additive term in ``C_obs + C_TLE,empirical``, which
  stays positive definite.
* The primary TLE remains the reference; the ensemble mean never replaces
  it. The offset ``primary - ensemble mean`` is reported.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence

import numpy as np

from tart_gnss.stage6.errors import ReferenceCovarianceError
from tart_gnss.stage6.propagation import propagate_tle
from tart_gnss.stage6.tle_reference import TLEReference

MIN_ENSEMBLE = 3
LABEL = "empirical reference dispersion"


@dataclass
class ReferenceCovariance:
    """Ensemble mean ``(N, 6)`` and empirical covariance ``(N, 6, 6)`` per epoch."""

    n: int
    mean_state: np.ndarray
    covariance: np.ndarray
    members: List[TLEReference] = field(default_factory=list)
    label: str = LABEL

    @property
    def position_dispersion_m(self) -> np.ndarray:
        """``sqrt(trace)`` of the position block per epoch."""
        return np.sqrt(np.trace(self.covariance[:, :3, :3], axis1=1, axis2=2))

    @property
    def velocity_dispersion_m_s(self) -> np.ndarray:
        """``sqrt(trace)`` of the velocity block per epoch."""
        return np.sqrt(np.trace(self.covariance[:, 3:, 3:], axis1=1, axis2=2))

    def summary(self) -> Dict[str, object]:
        """JSON-serialisable summary."""
        return {
            "label": self.label, "n_tle_ensemble": self.n,
            "member_epochs_utc": [m.epoch_utc.strftime("%Y-%m-%dT%H:%M:%S.%fZ") for m in self.members],
            "position_dispersion_rms_m_median": float(np.median(self.position_dispersion_m)),
            "velocity_dispersion_rms_m_s_median": float(np.median(self.velocity_dispersion_m_s)),
        }


def empirical_tle_covariance(references: Sequence[TLEReference], times: np.ndarray) -> ReferenceCovariance:
    """Estimate propagated-state dispersion from an ensemble of TLEs."""
    if len(references) < 2:
        raise ReferenceCovarianceError(f"need at least 2 ensemble TLEs, got {len(references)}")
    states = np.array([propagate_tle(r, times).state_gcrs for r in references])  # (M, N, 6)
    mean = states.mean(axis=0)
    d = states - mean
    cov = np.einsum("mni,mnj->nij", d, d) / (len(references) - 1)
    return ReferenceCovariance(n=len(references), mean_state=mean, covariance=0.5 * (cov + np.swapaxes(cov, 1, 2)),
                               members=list(references))
