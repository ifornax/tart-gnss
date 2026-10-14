"""Radial / along-track / cross-track (RTN = RIC) frame of the TLE-predicted state.

Convention (RTN = RIC here)::

    R = r / |r|                       radial
    N = (r x v) / |r x v|             cross-track (normal, C)
    T = N x R                         along-track / in-track (I)

``Q = [R^T; T^T; N^T]`` maps a GCRS vector to (R, T, N); it is orthonormal
with ``det Q = +1`` and ``R x T = N``. A covariance rotates as
``C_RTN = Q C Q^T``.

Velocity residuals are the inertial ``dv`` projected on the same axes. This
is *not* the time derivative of the RTN position residual, which also
contains the rotation of the frame (``omega x dr``, ~0.15 m/s per km at GNSS
orbital rates).
"""

from __future__ import annotations

import numpy as np

from tart_gnss.stage6.errors import RTNGeometryError

LABELS = ("radial", "along_track", "cross_track")


def rtn_basis(position: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    """Return the orthonormal R/T/N basis ``Q`` (rows R, T, N) from a reference state."""
    r = np.asarray(position, dtype=float)
    v = np.asarray(velocity, dtype=float)
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    if not (np.isfinite(rn) and rn > 0 and np.isfinite(hn)) or hn <= 1e-12 * rn * np.linalg.norm(v):
        raise RTNGeometryError("RTN basis undefined (zero position or position parallel to velocity)")
    rh = r / rn
    nh = h / hn
    th = np.cross(nh, rh)
    return np.vstack([rh, th, nh])


def rtn_bases(positions: np.ndarray, velocities: np.ndarray) -> np.ndarray:
    """``(N, 3, 3)`` RTN bases."""
    return np.array([rtn_basis(r, v) for r, v in zip(positions, velocities)])


def to_rtn(q: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Components of a GCRS vector on the R, T, N axes."""
    return q @ np.asarray(vector, dtype=float)


def covariance_to_rtn(q: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """``Q C Q^T`` (symmetrised)."""
    c = q @ np.asarray(cov, dtype=float) @ q.T
    return 0.5 * (c + c.T)
