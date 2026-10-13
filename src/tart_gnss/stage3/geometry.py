"""Coordinate conventions and conversions for Stage 3.

Local frame
-----------
Every direction is expressed in the station's local topocentric frame

* x = East, y = North, z = Up (ENU),

so a unit direction vector is ``s = (l, m, n)`` with ``n = sqrt(1 - l^2 - m^2)``
for a direction above the horizon. ``l`` and ``m`` are therefore the East and
North direction cosines of a zenith-centred tangent plane.

Azimuth is measured from North (0 deg) through East (90 deg); elevation is
0 deg at the horizon and 90 deg at the zenith:

    l = sin(Az) cos(El),   m = cos(Az) cos(El),   n = sin(El)
    Az = atan2(l, m),      El = arccos(sqrt(l^2 + m^2)) = arcsin(n)

Earth-fixed coordinates (ITRS/ECEF, metres) are rotated to ENU at geodetic
latitude ``phi`` and longitude ``lambda`` with

    [E]   [ -sin(lambda)           cos(lambda)          0       ] [dX]
    [N] = [ -sin(phi)cos(lambda)  -sin(phi)sin(lambda)  cos(phi)] [dY]
    [U]   [  cos(phi)cos(lambda)   cos(phi)sin(lambda)  sin(phi)] [dZ]

where ``(dX, dY, dZ)`` is a position relative to a reference point.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from tart_gnss.stage3.errors import CoordinateTransformError

#: WGS84 ellipsoid.
WGS84_A_M = 6378137.0
WGS84_F = 1.0 / 298.257223563
SPEED_OF_LIGHT_M_S = 299792458.0


def ecef_to_enu_matrix(latitude_deg: float, longitude_deg: float) -> np.ndarray:
    """Rotation matrix ``R`` such that ``enu = R @ d_ecef`` (see module docstring)."""
    phi, lam = np.radians(latitude_deg), np.radians(longitude_deg)
    return np.array(
        [
            [-np.sin(lam), np.cos(lam), 0.0],
            [-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)],
            [np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)],
        ]
    )


def ecef_to_enu(positions_m: np.ndarray, latitude_deg: float, longitude_deg: float,
                reference_m: np.ndarray = None) -> np.ndarray:
    """Convert ECEF positions (N, 3) to ENU (N, 3) relative to ``reference_m``.

    The reference defaults to the mean of ``positions_m`` (the array centre).
    """
    pos = np.atleast_2d(np.asarray(positions_m, dtype=float))
    ref = pos.mean(axis=0) if reference_m is None else np.asarray(reference_m, dtype=float)
    return (pos - ref) @ ecef_to_enu_matrix(latitude_deg, longitude_deg).T


def geodetic_to_ecef(latitude_deg: float, longitude_deg: float, height_m: float = 0.0) -> np.ndarray:
    """WGS84 ECEF position (m) of a geodetic point."""
    phi, lam = np.radians(latitude_deg), np.radians(longitude_deg)
    e2 = WGS84_F * (2.0 - WGS84_F)
    n_rad = WGS84_A_M / np.sqrt(1.0 - e2 * np.sin(phi) ** 2)
    return np.array([
        (n_rad + height_m) * np.cos(phi) * np.cos(lam),
        (n_rad + height_m) * np.cos(phi) * np.sin(lam),
        (n_rad * (1.0 - e2) + height_m) * np.sin(phi),
    ])


def enu_to_ecef(enu_m: np.ndarray, latitude_deg: float, longitude_deg: float, height_m: float = 0.0) -> np.ndarray:
    """ECEF positions (N, 3) of ENU offsets (N, 3) from a geodetic reference point."""
    ref = geodetic_to_ecef(latitude_deg, longitude_deg, height_m)
    return ref + np.atleast_2d(enu_m) @ ecef_to_enu_matrix(latitude_deg, longitude_deg)


def ecef_to_geodetic(position_m: np.ndarray) -> Tuple[float, float, float]:
    """WGS84 geodetic latitude, longitude (deg) and height (m) of one ECEF point.

    Uses Bowring's iteration, converged far below a millimetre after a few steps.
    """
    x, y, z = (float(v) for v in position_m)
    e2 = WGS84_F * (2.0 - WGS84_F)
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1.0 - e2))
    for _ in range(10):
        n_rad = WGS84_A_M / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n_rad
        lat = np.arctan2(z, p * (1.0 - e2 * n_rad / (n_rad + h)))
    n_rad = WGS84_A_M / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - n_rad
    return float(np.degrees(lat)), float(np.degrees(lon)), float(h)


def lm_to_unit(l, m) -> np.ndarray:
    """ENU unit vector(s) ``(l, m, n)`` for direction cosines above the horizon.

    Raises
    ------
    CoordinateTransformError
        If ``l^2 + m^2 > 1`` (no real direction).
    """
    l, m = np.asarray(l, dtype=float), np.asarray(m, dtype=float)
    r2 = l * l + m * m
    if np.any(r2 > 1.0 + 1e-12):
        raise CoordinateTransformError("l^2 + m^2 > 1: not a direction on the sky")
    n = np.sqrt(np.clip(1.0 - r2, 0.0, None))
    return np.stack([l, m, n], axis=-1)


def azel_to_unit(az_deg, el_deg) -> np.ndarray:
    """ENU unit vector(s) for azimuth/elevation in degrees."""
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.stack([np.sin(az) * np.cos(el), np.cos(az) * np.cos(el), np.sin(el)], axis=-1)


def unit_to_azel(s) -> Tuple[np.ndarray, np.ndarray]:
    """Azimuth [0, 360) and elevation (deg) of ENU vector(s) (need not be normalised)."""
    s = np.asarray(s, dtype=float)
    e, nth, u = s[..., 0], s[..., 1], s[..., 2]
    az = np.degrees(np.arctan2(e, nth)) % 360.0
    el = np.degrees(np.arctan2(u, np.hypot(e, nth)))
    return az, el


def lm_to_azel(l, m) -> Tuple[np.ndarray, np.ndarray]:
    """Azimuth/elevation (deg) of direction cosines ``(l, m)`` above the horizon."""
    return unit_to_azel(lm_to_unit(l, m))


def azel_to_lm(az_deg, el_deg) -> Tuple[np.ndarray, np.ndarray]:
    """Direction cosines ``(l, m)`` of azimuth/elevation in degrees."""
    s = azel_to_unit(az_deg, el_deg)
    return s[..., 0], s[..., 1]


def angular_separation_deg(s1, s2) -> np.ndarray:
    """Great-circle angle (deg) between ENU vectors, robust at small angles."""
    a = np.asarray(s1, dtype=float)
    b = np.asarray(s2, dtype=float)
    a = a / np.linalg.norm(a, axis=-1, keepdims=True)
    b = b / np.linalg.norm(b, axis=-1, keepdims=True)
    cross = np.linalg.norm(np.cross(a, b), axis=-1)
    dot = np.sum(a * b, axis=-1)
    return np.degrees(np.arctan2(cross, dot))


def tangent_basis(s0) -> Tuple[np.ndarray, np.ndarray]:
    """Two unit vectors orthogonal to ``s0`` (and to each other).

    ``e1`` points along increasing azimuth and ``e2`` along increasing elevation
    (at the zenith, East and North).
    """
    s0 = np.asarray(s0, dtype=float)
    s0 = s0 / np.linalg.norm(s0)
    up = np.array([0.0, 0.0, 1.0])
    e1 = np.cross(up, s0)
    if np.linalg.norm(e1) < 1e-9:
        e1 = np.array([1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(s0, e1)
    return e1, e2 / np.linalg.norm(e2)


def offset_directions(s0, x_rad, y_rad) -> np.ndarray:
    """Directions displaced from ``s0`` by tangent-plane offsets (radians).

    The displacement is an exact great-circle rotation by angle
    ``rho = hypot(x, y)`` towards ``x e1 + y e2``, so ``rho`` is the true
    angular distance from ``s0``.
    """
    s0 = np.asarray(s0, dtype=float)
    s0 = s0 / np.linalg.norm(s0)
    e1, e2 = tangent_basis(s0)
    x, y = np.asarray(x_rad, dtype=float), np.asarray(y_rad, dtype=float)
    rho = np.hypot(x, y)
    safe = np.where(rho > 0, rho, 1.0)
    dirn = (x[..., None] * e1 + y[..., None] * e2) / safe[..., None]
    return np.cos(rho)[..., None] * s0 + np.sin(rho)[..., None] * dirn


def azel_jacobian(l: float, m: float) -> np.ndarray:
    """Jacobian ``d(Az, El)/d(l, m)`` in radians per unit direction cosine.

    With ``r = sqrt(l^2 + m^2)``:

        dAz/dl =  m / r^2,                  dAz/dm = -l / r^2
        dEl/dl = -l / (r sqrt(1 - r^2)),    dEl/dm = -m / (r sqrt(1 - r^2))

    Singular at the zenith (``r = 0``, azimuth undefined) and the horizon
    (``r = 1``).

    Raises
    ------
    CoordinateTransformError
        At either singularity.
    """
    r2 = l * l + m * m
    r = np.sqrt(r2)
    if r < 1e-12 or r2 >= 1.0:
        raise CoordinateTransformError("Az/El Jacobian is singular at the zenith or horizon")
    q = r * np.sqrt(1.0 - r2)
    return np.array([[m / r2, -l / r2], [-l / q, -m / q]])


def wrap_azimuth_near(az_deg, reference_deg: float) -> np.ndarray:
    """Shift azimuths by multiples of 360 deg to lie within 180 deg of ``reference_deg``."""
    az = np.asarray(az_deg, dtype=float)
    return reference_deg + (az - reference_deg + 180.0) % 360.0 - 180.0
