"""Coordinate conversions for Stage 4 (ITRS/ECEF, metres).

Conventions are those of Stage 3 (:mod:`tart_gnss.stage3.geometry`): local
frame East-North-Up; azimuth from North (0 deg) through East (90 deg);
elevation 0 deg at the horizon and 90 deg at the zenith, so

    s_ENU = (cos E sin A, cos E cos A, sin E) = (l, m, n).

An ENU vector is rotated to ECEF with the transpose of the ECEF->ENU matrix at
the station's geodetic latitude/longitude (WGS84). Directions are rotated only;
station positions are converted with :func:`geodetic_to_ecef`.

Angular uncertainty of a sight line is carried as a 3x3 covariance of the unit
vector in ECEF (rank 2, tangent to the sphere), obtained from the Stage-3
covariance of ``(l, m)`` (preferred: finite at the zenith) or of Az/El.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from tart_gnss.stage3.geometry import (  # noqa: F401  (re-exported)
    SPEED_OF_LIGHT_M_S,
    WGS84_A_M,
    WGS84_F,
    azel_to_unit,
    ecef_to_enu_matrix,
    ecef_to_geodetic,
    geodetic_to_ecef,
    lm_to_unit,
    unit_to_azel,
)
from tart_gnss.stage4.errors import CovarianceError, TriangulationError


def enu_to_ecef_vector(v_enu, latitude_deg: float, longitude_deg: float) -> np.ndarray:
    """Rotate ENU vector(s) (no translation) to ECEF."""
    return np.asarray(v_enu, dtype=float) @ ecef_to_enu_matrix(latitude_deg, longitude_deg)


def ecef_to_enu_vector(v_ecef, latitude_deg: float, longitude_deg: float) -> np.ndarray:
    """Rotate ECEF vector(s) (no translation) to ENU."""
    return np.asarray(v_ecef, dtype=float) @ ecef_to_enu_matrix(latitude_deg, longitude_deg).T


def azel_to_ecef_unit_vector(az_deg: float, el_deg: float, lat_deg: float, lon_deg: float) -> np.ndarray:
    """Convert local azimuth/elevation to an ECEF unit line-of-sight vector."""
    if not (np.all(np.isfinite(az_deg)) and np.all(np.isfinite(el_deg))):
        raise TriangulationError("non-finite azimuth/elevation")
    return enu_to_ecef_vector(azel_to_unit(az_deg, el_deg), lat_deg, lon_deg)


def lm_to_ecef_unit_vector(l: float, m: float, lat_deg: float, lon_deg: float) -> np.ndarray:
    """Convert Stage-3 direction cosines ``(l, m)`` (East, North) to an ECEF unit vector."""
    if not (np.isfinite(l) and np.isfinite(m)):
        raise TriangulationError("non-finite direction cosines")
    return enu_to_ecef_vector(lm_to_unit(l, m), lat_deg, lon_deg)


def ecef_to_azel(target_ecef_m, station_ecef_m, lat_deg: float, lon_deg: float) -> Tuple[float, float, float]:
    """Azimuth, elevation (deg) and range (m) of an ECEF point seen from a station (geometric)."""
    d = np.asarray(target_ecef_m, dtype=float) - np.asarray(station_ecef_m, dtype=float)
    az, el = unit_to_azel(ecef_to_enu_vector(d, lat_deg, lon_deg))
    return float(az), float(el), float(np.linalg.norm(d))


def tangent_basis(s) -> Tuple[np.ndarray, np.ndarray]:
    """Orthonormal ``e1, e2`` perpendicular to unit vector ``s`` (any frame).

    The helper axis is the coordinate axis least aligned with ``s``, so the
    construction never degenerates (no special case at poles or zenith).
    """
    s = np.asarray(s, dtype=float)
    s = s / np.linalg.norm(s)
    helper = np.zeros(3)
    helper[int(np.argmin(np.abs(s)))] = 1.0
    e1 = np.cross(helper, s)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(s, e1)
    return e1, e2 / np.linalg.norm(e2)


def lm_jacobian_enu(l: float, m: float) -> np.ndarray:
    """``d s_ENU / d(l, m)`` (3x2); its columns are tangent to the unit sphere."""
    n = np.sqrt(max(1.0 - l * l - m * m, 0.0))
    if n < 1e-9:
        raise CovarianceError("direction on the horizon: (l, m) Jacobian undefined")
    return np.array([[1.0, 0.0], [0.0, 1.0], [-l / n, -m / n]])


def azel_jacobian_enu(az_deg: float, el_deg: float) -> np.ndarray:
    """``d s_ENU / d(Az, El)`` (3x2) per radian."""
    a, e = np.radians(az_deg), np.radians(el_deg)
    return np.array([
        [np.cos(e) * np.cos(a), -np.sin(e) * np.sin(a)],
        [-np.cos(e) * np.sin(a), -np.sin(e) * np.cos(a)],
        [0.0, np.cos(e)],
    ])


def _check_cov2(cov: np.ndarray, what: str) -> np.ndarray:
    cov = np.asarray(cov, dtype=float)
    if cov.shape != (2, 2) or not np.all(np.isfinite(cov)):
        raise CovarianceError(f"{what} covariance is not a finite 2x2 matrix")
    cov = 0.5 * (cov + cov.T)
    if np.min(np.linalg.eigvalsh(cov)) <= 0.0:
        raise CovarianceError(f"{what} covariance is not positive definite")
    return cov


def lm_covariance_to_ecef(l: float, m: float, cov_lm, lat_deg: float, lon_deg: float) -> np.ndarray:
    """3x3 ECEF covariance (rad^2, tangent, rank 2) of a sight line from the ``(l, m)`` covariance."""
    j = enu_to_ecef_vector(lm_jacobian_enu(l, m).T, lat_deg, lon_deg).T
    return j @ _check_cov2(cov_lm, "(l, m)") @ j.T


def azel_covariance_to_ecef(az_deg: float, el_deg: float, cov_azel_deg2, lat_deg: float,
                            lon_deg: float) -> np.ndarray:
    """3x3 ECEF covariance (rad^2) of a sight line from an Az/El covariance in deg^2."""
    cov = _check_cov2(np.asarray(cov_azel_deg2, dtype=float) * np.radians(1.0) ** 2, "Az/El")
    j = enu_to_ecef_vector(azel_jacobian_enu(az_deg, el_deg).T, lat_deg, lon_deg).T
    return j @ cov @ j.T


def tangent_covariance(s, cov3: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Express a 3x3 tangent covariance in the :func:`tangent_basis` of ``s``.

    Returns ``(E, C)`` with ``E = [e1 e2]`` (3x2) and ``C = E^T cov3 E`` (2x2).
    """
    e1, e2 = tangent_basis(s)
    e = np.column_stack([e1, e2])
    c = e.T @ np.asarray(cov3, dtype=float) @ e
    return e, 0.5 * (c + c.T)


def angular_sigma_deg(cov3: np.ndarray) -> float:
    """Major-axis 1-sigma angular uncertainty (deg) of a 3x3 tangent covariance."""
    ev = np.linalg.eigvalsh(np.asarray(cov3, dtype=float))
    return float(np.degrees(np.sqrt(max(ev[-1], 0.0))))


def pairwise_angles_deg(sight_vectors) -> np.ndarray:
    """All pairwise angles ``arccos(s_i . s_j)`` (deg) between unit vectors, i < j."""
    s = np.asarray(sight_vectors, dtype=float)
    out = []
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            out.append(np.degrees(np.arctan2(np.linalg.norm(np.cross(s[i], s[j])), float(np.dot(s[i], s[j])))))
    return np.asarray(out, dtype=float)
