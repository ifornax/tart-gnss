"""Frames for the TLE comparison.

SGP4 source frame
    SGP4 (``sgp4`` library, WGS72 constants as required by the TLE theory)
    returns position and velocity in TEME (true equator, mean equinox of
    date).
Transformation library
    Skyfield ``EarthSatellite.at(t)``: TEME -> ITRS with GMST (IAU 1982),
    then ITRS -> GCRS with ``itrs.rotation_at(t)``'s inverse, i.e. the IAU
    2000A/2006 precession-nutation and GAST from UT1. Velocities are carried
    through the same chain, including the Earth-rotation term, so the result
    is the *inertial* GCRS velocity.
Comparison frame
    GCRS, identical to the Stage-5 state frame (Stage 5 uses the same
    ``itrs.rotation_at``). ITRS/ECEF positions are ``R(t) r_GCRS`` with the
    Stage-5 matrices, so observed and predicted ECEF positions share one
    rotation. Earth-fixed *velocities* are not compared: they differ from
    inertial ones by ``omega x r``.
Time scale
    UTC epochs from the Stage-5 table, converted with
    :func:`tart_gnss.stage5.frames.skyfield_times` (whole days plus seconds of
    day on Skyfield's builtin timescale: leap seconds and UT1 from the bundled
    tables, no network).
Earth orientation
    Polar motion zero (builtin timescale), UT1 from Skyfield's bundled
    predictions: ~40 m and well below 1 m respectively at GNSS range.

TEME is never treated as GCRS/J2000: the two differ by precession since
J2000 plus the equation of the equinoxes (about 0.36 deg in 2026, ~170 km
at GNSS radius).
"""

from __future__ import annotations

import numpy as np

from tart_gnss.stage5.frames import (  # noqa: F401  (re-exported: one rotation for Stages 5 and 6)
    FRAME_DESCRIPTION,
    gcrs_to_ecef_positions,
    itrs_to_gcrs_matrices,
    skyfield_times,
    timescale,
)

OMEGA_EARTH_RAD_S = 7.292115e-5

FRAMES = {
    "sgp4_source_frame": "TEME (true equator, mean equinox of date), sgp4 library, WGS72",
    "transformation_library": "Skyfield EarthSatellite.at(t): TEME -> ITRS (GMST 1982) -> GCRS "
                              "(itrs.rotation_at, IAU 2000A/2006, GAST from UT1), velocity included",
    "comparison_frame": "GCRS (same rotation as Stage 5); ITRS/ECEF positions as R(t) r_GCRS",
    "time_scale": "UTC from the Stage-5 table via tart_gnss.stage5.frames.skyfield_times (builtin timescale, "
                  "leap seconds handled)",
    "earth_orientation": "polar motion zero, UT1 from Skyfield's bundled tables; no network",
    "velocity": "inertial GCRS velocities compared; Earth-fixed velocities never compared with inertial ones",
}


def ecef_velocity_from_gcrs(r_gcrs: np.ndarray, v_gcrs: np.ndarray, matrices: np.ndarray) -> np.ndarray:
    """Approximate Earth-fixed velocity ``R v - omega x (R r)`` (diagnostic only; never compared)."""
    r_e = np.einsum("nji,nj->ni", matrices, r_gcrs)
    v_e = np.einsum("nji,nj->ni", matrices, v_gcrs)
    w = np.array([0.0, 0.0, OMEGA_EARTH_RAD_S])
    return v_e - np.cross(w, r_e)
