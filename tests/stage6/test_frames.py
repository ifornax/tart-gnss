"""Frame consistency: Skyfield's TEME -> GCRS chain against an independent route (spec section 35)."""

import numpy as np
import pytest

pytest.importorskip("skyfield")

from tart_gnss.stage6 import synthetic as syn  # noqa: E402
from tart_gnss.stage6.frames import ecef_velocity_from_gcrs  # noqa: E402
from tart_gnss.stage6.propagation import propagate_tle  # noqa: E402

from .helpers import L1, L2, reference, times  # noqa: E402

T = times(12, step=300.0)


def test_skyfield_gcrs_matches_independent_teme_route():
    # Raw sgp4 TEME -> explicit GMST82 -> Stage-5 ITRS->GCRS matrices, versus Skyfield's EarthSatellite.
    tr = propagate_tle(reference(), T)
    ind = syn.independent_gcrs_states(L1, L2, T)
    assert np.max(np.abs(tr.position_gcrs_m - ind[:, :3])) < 0.1
    assert np.max(np.abs(tr.velocity_gcrs_m_s - ind[:, 3:])) < 1e-3


def test_teme_is_not_gcrs():
    r_teme, _ = syn.teme_states(L1, L2, T)
    tr = propagate_tle(reference(), T)
    # Precession since J2000 (~0.36 deg) displaces the position by ~170 km at GNSS radius.
    assert np.min(np.linalg.norm(r_teme - tr.position_gcrs_m, axis=1)) > 1e5


def test_velocity_is_the_inertial_derivative_of_position():
    h = 0.5
    tr = propagate_tle(reference(), T)
    fd = (propagate_tle(reference(), T + h).position_gcrs_m - propagate_tle(reference(), T - h).position_gcrs_m) / (2 * h)
    # SGP4's analytic velocity is not the exact derivative of its position: agreement to a few cm/s.
    assert np.max(np.abs(fd - tr.velocity_gcrs_m_s)) < 0.05


def test_ecef_position_and_velocity_frames():
    tr = propagate_tle(reference(), T)
    assert np.allclose(np.einsum("nij,nj->ni", tr.matrices, tr.position_ecef_m), tr.position_gcrs_m, atol=1e-6)
    from skyfield.api import EarthSatellite
    from skyfield.framelib import itrs

    from tart_gnss.stage5.frames import skyfield_times, timescale

    pos, vel = EarthSatellite(L1, L2, "x", timescale()).at(skyfield_times(T)).frame_xyz_and_velocity(itrs)
    assert np.allclose(pos.m.T, tr.position_ecef_m, atol=1e-3)
    v_ecef = ecef_velocity_from_gcrs(tr.position_gcrs_m, tr.velocity_gcrs_m_s, tr.matrices)
    assert np.allclose(vel.m_per_s.T, v_ecef, atol=0.05)
    # Earth-fixed and inertial velocities differ by |omega x r| ~ 1.9 km/s: never compare them.
    v_inertial_in_ecef = np.einsum("nji,nj->ni", tr.matrices, tr.velocity_gcrs_m_s)
    assert np.min(np.linalg.norm(v_inertial_in_ecef - v_ecef, axis=1)) > 1e3


def test_velocity_units_are_m_per_s():
    tr = propagate_tle(reference(), T)
    speed = np.linalg.norm(tr.velocity_gcrs_m_s, axis=1)
    radius = np.linalg.norm(tr.position_gcrs_m, axis=1)
    assert np.all((3.5e3 < speed) & (speed < 4.2e3)) and np.all((2.6e7 < radius) & (radius < 2.7e7))
