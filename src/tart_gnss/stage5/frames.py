"""ITRS/ECEF <-> GCRS transformations of Stage-4 positions and covariances.

Implementation
    Skyfield's ``skyfield.framelib.itrs.rotation_at(t)`` returns the 3x3
    rotation ``R(t)`` from the GCRS axes to the ITRS (Earth-fixed) axes,
    built from the IAU 2000A/2006 precession-nutation matrix, Greenwich
    apparent sidereal time (UT1) and polar motion. Hence
    ``r_GCRS = R(t)^T r_ITRS`` and ``r_ITRS = R(t) r_GCRS``. Times come from
    ``skyfield.api.load.timescale(builtin=True)``: the leap-second and
    UT1-UTC tables shipped with Skyfield, so no network access is needed.
    That timescale carries no polar-motion table, so polar motion is zero
    (at most ~0.3 arcsec, i.e. ~40 m at GNSS range); the bundled UT1-UTC
    predictions are good to well under a millisecond (~0.5 m at the
    equator for Earth rotation). Both are negligible against the km-scale
    Stage-4 uncertainties. Unix seconds are converted to UTC as whole days
    plus seconds of day, so leap seconds are handled by Skyfield.

Covariance
    The transformation at a fixed epoch is a pure rotation, so a covariance
    transforms as ``C_GCRS = M C_ITRS M^T`` with ``M = R(t)^T``. Position
    and covariance are always transformed together.

Velocity
    Only positions are transformed back to ECEF. An Earth-fixed velocity
    would need the ``omega x r`` term of the rotating frame and is not an
    output of Stage 5.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Tuple

import numpy as np

from tart_gnss.stage5.errors import FrameTransformError

_TIMESCALE = None

FRAME_DESCRIPTION = (
    "GCRS axes via Skyfield itrs.rotation_at(t) (IAU 2000A/2006 precession-nutation, GAST from UT1; "
    "polar motion zero: builtin timescale has no polar-motion table); r_GCRS = R(t)^T r_ITRS, "
    "C_GCRS = R(t)^T C_ITRS R(t)"
)


def timescale():
    """Skyfield timescale built from the bundled data files (no network)."""
    global _TIMESCALE
    if _TIMESCALE is None:
        from skyfield.api import load

        _TIMESCALE = load.timescale(builtin=True)
    return _TIMESCALE


def skyfield_times(t_unix_s):
    """Skyfield ``Time`` for Unix seconds (UTC, leap seconds handled)."""
    t = np.atleast_1d(np.asarray(t_unix_s, dtype=float))
    days = np.floor(t / 86400.0)
    return timescale().utc(1970, 1, 1 + days.astype(np.int64), 0, 0, t - days * 86400.0)


def itrs_to_gcrs_matrices(t_unix_s) -> np.ndarray:
    """Rotation matrices ``M_k`` (shape ``(N, 3, 3)``) with ``r_GCRS = M_k r_ITRS``."""
    from skyfield.framelib import itrs

    t = np.atleast_1d(np.asarray(t_unix_s, dtype=float))
    if not np.all(np.isfinite(t)):
        raise FrameTransformError("non-finite epoch for the ITRS -> GCRS rotation")
    r = np.asarray(itrs.rotation_at(skyfield_times(t)), dtype=float)
    if r.ndim == 2:
        r = r[:, :, None]
    m = np.transpose(r, (2, 1, 0))  # (N, 3, 3), M_k = R_k^T
    if not np.all(np.isfinite(m)):
        raise FrameTransformError("non-finite ITRS -> GCRS rotation")
    return m


def _unix(obstime) -> float:
    if isinstance(obstime, datetime):
        if obstime.tzinfo is None:
            obstime = obstime.replace(tzinfo=timezone.utc)
        return obstime.timestamp()
    return float(obstime)


def ecef_measurement_to_gcrs(
    position_ecef_m: np.ndarray, covariance_ecef_m2: np.ndarray, obstime
) -> Tuple[np.ndarray, np.ndarray]:
    """Transform one ECEF position and covariance into the inertial frame.

    ``obstime`` is an aware ``datetime`` (naive means UTC) or Unix seconds.
    """
    m = itrs_to_gcrs_matrices([_unix(obstime)])[0]
    return rotate_measurement(m, position_ecef_m, covariance_ecef_m2)


def rotate_measurement(m: np.ndarray, position_m, covariance_m2) -> Tuple[np.ndarray, np.ndarray]:
    """``(M r, M C M^T)`` for one rotation ``M``."""
    r = m @ np.asarray(position_m, dtype=float)
    c = m @ np.asarray(covariance_m2, dtype=float) @ m.T
    return r, 0.5 * (c + c.T)


def gcrs_to_ecef_positions(positions_gcrs_m: np.ndarray, matrices: np.ndarray) -> np.ndarray:
    """ECEF positions ``M_k^T r_k`` for GCRS positions ``(N, 3)`` and matrices ``(N, 3, 3)``."""
    return np.einsum("nji,nj->ni", matrices, np.asarray(positions_gcrs_m, dtype=float))


def earth_pole_gcrs(t_unix_s: float) -> np.ndarray:
    """Unit vector of the ITRS z axis (Earth's rotation pole) in GCRS at one epoch."""
    return itrs_to_gcrs_matrices([t_unix_s])[0][:, 2]
