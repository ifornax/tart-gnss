"""ITRS/ECEF <-> GCRS rotation of positions and covariances."""

from datetime import datetime, timezone

import numpy as np
import pytest

pytest.importorskip("skyfield")

from tart_gnss.stage4.geometry import geodetic_to_ecef  # noqa: E402
from tart_gnss.stage5 import frames  # noqa: E402
from tart_gnss.stage5.errors import FrameTransformError  # noqa: E402
from tart_gnss.stage5.synthetic import T0_UNIX_S  # noqa: E402

T = T0_UNIX_S + np.array([0.0, 37.5, 3600.0])


def test_rotation_matches_independent_skyfield_topos_route():
    # Skyfield's own GCRS position of a ground point (not via our matrices) must equal M r_ITRS.
    from skyfield.api import wgs84

    lat, lon, alt = -33.319413, 26.507688, 642.26
    m = frames.itrs_to_gcrs_matrices(T)
    g = wgs84.latlon(lat, lon, alt).at(frames.skyfield_times(T)).position.m
    e = geodetic_to_ecef(lat, lon, alt)
    for k in range(len(T)):
        assert np.allclose(m[k] @ e, g[:, k], atol=1e-3)


def test_unix_seconds_handle_leap_seconds():
    from skyfield.api import load

    ts = load.timescale(builtin=True)
    ref = ts.from_datetimes([datetime.fromtimestamp(t, tz=timezone.utc) for t in T])
    assert np.allclose((frames.skyfield_times(T).tt - ref.tt) * 86400.0, 0.0, atol=1e-5)


def test_matrices_are_rotations():
    for m in frames.itrs_to_gcrs_matrices(T):
        assert np.allclose(m @ m.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(m) == pytest.approx(1.0, abs=1e-12)


def test_measurement_round_trip_and_covariance_invariants():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(3, 3))
    c = a @ a.T * 1e8 + np.eye(3) * 1e6
    r = np.array([1.2e7, 9e6, -2e7])
    rg, cg = frames.ecef_measurement_to_gcrs(r, c, datetime.fromtimestamp(T[0], tz=timezone.utc))
    m = frames.itrs_to_gcrs_matrices([T[0]])[0]
    assert np.allclose(frames.gcrs_to_ecef_positions(rg[None], m[None])[0], r, atol=1e-6)
    assert np.allclose(m.T @ cg @ m, c, rtol=1e-10)
    assert np.allclose(np.linalg.eigvalsh(cg), np.linalg.eigvalsh(c), rtol=1e-10)
    assert np.linalg.norm(rg) == pytest.approx(np.linalg.norm(r), rel=1e-14)
    # Off-diagonal terms are carried through (the covariance is not diagonalised).
    assert np.allclose(cg, cg.T)


def test_naive_datetime_is_utc_and_unix_equivalent():
    r, c = np.array([1e7, 2e7, 3e6]), np.eye(3)
    naive = datetime.fromtimestamp(T[1], tz=timezone.utc).replace(tzinfo=None)
    a = frames.ecef_measurement_to_gcrs(r, c, naive)[0]
    b = frames.ecef_measurement_to_gcrs(r, c, T[1])[0]
    assert np.allclose(a, b, atol=1e-6)


def test_covariance_transform_is_statistically_consistent():
    # Rotating samples drawn from C_ITRS gives samples whose covariance is M C M^T.
    rng = np.random.default_rng(1)
    c = np.diag([9e8, 1e6, 4e6])
    m = frames.itrs_to_gcrs_matrices([T[0]])[0]
    x = rng.multivariate_normal(np.zeros(3), c, size=20000) @ m.T
    _, cg = frames.rotate_measurement(m, np.zeros(3), c)
    assert np.allclose(np.cov(x.T), cg, rtol=0.05, atol=0.03 * 9e8)


def test_pole_is_close_to_gcrs_z():
    k = frames.earth_pole_gcrs(T[0])
    assert np.degrees(np.arccos(k[2])) < 0.5  # precession since J2000 is ~0.35 deg


def test_non_finite_time_raises():
    with pytest.raises(FrameTransformError):
        frames.itrs_to_gcrs_matrices([np.nan])
