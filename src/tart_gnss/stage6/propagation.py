"""SGP4 propagation of the reference TLE at the exact Stage-5 epochs.

No interpolation: every Stage-5 epoch is propagated directly with Skyfield's
``EarthSatellite`` (``sgp4`` library) and converted to GCRS (see
:mod:`tart_gnss.stage6.frames`). TLE age is ``t_obs - t_TLE epoch``: positive
after the epoch (forward propagation), negative before it (back-propagation).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tart_gnss.stage6.errors import TLEPropagationError
from tart_gnss.stage6.frames import gcrs_to_ecef_positions, itrs_to_gcrs_matrices, skyfield_times, timescale
from tart_gnss.stage6.tle_reference import TLEReference


@dataclass
class TLETrajectory:
    """TLE-predicted states at the requested epochs (SI units)."""

    times_s: np.ndarray
    position_gcrs_m: np.ndarray
    velocity_gcrs_m_s: np.ndarray
    position_ecef_m: np.ndarray
    matrices: np.ndarray
    age_days: np.ndarray
    tle: TLEReference

    @property
    def state_gcrs(self) -> np.ndarray:
        """``(N, 6)`` GCRS states."""
        return np.hstack([self.position_gcrs_m, self.velocity_gcrs_m_s])


def tle_age_days(times_s, tle: TLEReference) -> np.ndarray:
    """``(t_obs - t_epoch)`` in days; positive = propagated forward from the epoch."""
    return (np.asarray(times_s, dtype=float) - tle.epoch_unix_s) / 86400.0


def propagate_tle(tle: TLEReference, times: np.ndarray) -> TLETrajectory:
    """Propagate the TLE with SGP4 at exact Stage-5 epochs (Unix seconds)."""
    from skyfield.api import EarthSatellite

    t = np.atleast_1d(np.asarray(times, dtype=float))
    try:
        sat = EarthSatellite(tle.line1, tle.line2, tle.name, timescale())
        geo = sat.at(skyfield_times(t))
    except Exception as exc:  # noqa: BLE001 - sgp4/skyfield raise various errors
        raise TLEPropagationError(f"SGP4 propagation failed: {exc}") from exc
    r = np.asarray(geo.position.m, dtype=float).reshape(3, -1).T
    v = np.asarray(geo.velocity.m_per_s, dtype=float).reshape(3, -1).T
    msg = getattr(geo, "message", None)
    if not (np.all(np.isfinite(r)) and np.all(np.isfinite(v))):
        raise TLEPropagationError(f"SGP4 returned non-finite states ({msg})")
    mats = itrs_to_gcrs_matrices(t)
    return TLETrajectory(times_s=t, position_gcrs_m=r, velocity_gcrs_m_s=v,
                         position_ecef_m=gcrs_to_ecef_positions(r, mats), matrices=mats,
                         age_days=tle_age_days(t, tle), tle=tle)
