"""Coordinate conventions: cardinal directions, frames, tangent bases, covariance mapping."""

import numpy as np
import pytest

from tart_gnss.stage4 import geometry as g
from tart_gnss.stage4.errors import CovarianceError, TriangulationError


def test_enu_formula_matches_spec():
    a, e = np.radians(37.0), np.radians(21.0)
    want = [np.cos(e) * np.sin(a), np.cos(e) * np.cos(a), np.sin(e)]
    assert np.allclose(g.azel_to_unit(37.0, 21.0), want)


@pytest.mark.parametrize("az, el, want", [
    (0.0, 0.0, (0, 0, 1)),     # North at (0, 0) is +Z
    (90.0, 0.0, (0, 1, 0)),    # East is +Y
    (180.0, 0.0, (0, 0, -1)),  # South is -Z
    (270.0, 0.0, (0, -1, 0)),  # West is -Y
    (0.0, 90.0, (1, 0, 0)),    # zenith is +X
])
def test_cardinal_directions_at_equator_greenwich(az, el, want):
    assert np.allclose(g.azel_to_ecef_unit_vector(az, el, 0.0, 0.0), want, atol=1e-12)


def test_cardinal_directions_at_north_pole_and_lon90():
    assert np.allclose(g.azel_to_ecef_unit_vector(0.0, 90.0, 90.0, 0.0), (0, 0, 1), atol=1e-12)
    # At (0, 90E): East is -X, North +Z, Up +Y.
    assert np.allclose(g.azel_to_ecef_unit_vector(90.0, 0.0, 0.0, 90.0), (-1, 0, 0), atol=1e-12)
    assert np.allclose(g.azel_to_ecef_unit_vector(0.0, 90.0, 0.0, 90.0), (0, 1, 0), atol=1e-12)


def test_southern_station_sees_up_away_from_earth_centre():
    lat, lon, alt = -33.319413, 26.507688, 642.26
    up = g.azel_to_ecef_unit_vector(0.0, 90.0, lat, lon)
    r = g.geodetic_to_ecef(lat, lon, alt)
    assert np.dot(up, r / np.linalg.norm(r)) > 0.99
    north = g.azel_to_ecef_unit_vector(0.0, 0.0, lat, lon)
    assert north[2] > 0.8  # North points towards +Z (the pole) from the southern hemisphere


def test_lm_and_azel_paths_agree():
    l, m = 0.31, -0.42
    az, el = g.unit_to_azel(g.lm_to_unit(l, m))
    assert np.allclose(g.lm_to_ecef_unit_vector(l, m, -20.0, 57.0), g.azel_to_ecef_unit_vector(az, el, -20.0, 57.0))


def test_geodetic_round_trip():
    for lat, lon, h in [(-33.3, 26.5, 642.0), (89.9, -120.0, 10.0), (0.0, 179.0, 2.0e7)]:
        back = g.ecef_to_geodetic(g.geodetic_to_ecef(lat, lon, h))
        assert np.allclose(back, (lat, lon, h), atol=1e-6)


def test_ecef_to_azel_against_skyfield():
    """Non-circular check: skyfield's own WGS84 and alt/az against ours."""
    api = pytest.importorskip("skyfield.api")
    ts = api.load.timescale()
    t = ts.utc(2026, 6, 25, 9, 14, 0)
    station_geo = (-20.2587508, 57.7591989, 20.0)
    target_geo = (-5.0, 40.0, 2.0e7)
    station = api.wgs84.latlon(station_geo[0], station_geo[1], elevation_m=station_geo[2])
    target = api.wgs84.latlon(target_geo[0], target_geo[1], elevation_m=target_geo[2])
    alt, az, dist = (target - station).at(t).altaz()
    x_t = g.geodetic_to_ecef(*target_geo)
    r_s = g.geodetic_to_ecef(*station_geo)
    assert np.allclose(x_t, target.itrs_xyz.m, atol=1e-3)
    az2, el2, rng = g.ecef_to_azel(x_t, r_s, station_geo[0], station_geo[1])
    assert az2 == pytest.approx(az.degrees, abs=1e-6)
    assert el2 == pytest.approx(alt.degrees, abs=1e-6)
    assert rng == pytest.approx(dist.m, rel=1e-9)
    # And the Az/El -> ECEF sight line passes through the target.
    s = g.azel_to_ecef_unit_vector(az.degrees, alt.degrees, station_geo[0], station_geo[1])
    assert np.linalg.norm(np.cross(s, x_t - r_s)) / rng < 1e-9


@pytest.mark.parametrize("s", [(0, 0, 1), (1, 0, 0), (0, -1, 0), (0.6, 0.0, 0.8), (1e-12, 1e-12, 1.0)])
def test_tangent_basis_orthonormal_everywhere(s):
    s = np.asarray(s, float) / np.linalg.norm(s)
    e1, e2 = g.tangent_basis(s)
    m = np.column_stack([e1, e2, s])
    assert np.allclose(m.T @ m, np.eye(3), atol=1e-12)


def test_lm_covariance_is_tangent_and_matches_isotropic_sigma():
    from tart_gnss.stage4.synthetic import isotropic_lm_covariance

    l, m = 0.2, 0.5
    sig = np.radians(0.1)
    cov3 = g.lm_covariance_to_ecef(l, m, isotropic_lm_covariance(l, m, sig), -20.0, 57.0)
    s = g.lm_to_ecef_unit_vector(l, m, -20.0, 57.0)
    assert np.allclose(cov3 @ s, 0.0, atol=1e-18)
    _, c2 = g.tangent_covariance(s, cov3)
    assert np.allclose(np.linalg.eigvalsh(c2), sig ** 2, rtol=1e-9)
    assert g.angular_sigma_deg(cov3) == pytest.approx(0.1, rel=1e-9)


def test_lm_and_azel_covariances_agree():
    from tart_gnss.stage3.geometry import azel_jacobian

    l, m = -0.3, 0.4
    cov_lm = np.array([[2e-6, 3e-7], [3e-7, 1e-6]])
    jac = azel_jacobian(l, m)
    cov_ae = np.degrees(1.0) ** 2 * jac @ cov_lm @ jac.T
    az, el = g.unit_to_azel(g.lm_to_unit(l, m))
    a = g.lm_covariance_to_ecef(l, m, cov_lm, 10.0, 20.0)
    b = g.azel_covariance_to_ecef(az, el, cov_ae, 10.0, 20.0)
    assert np.allclose(a, b, rtol=1e-8, atol=1e-18)


def test_zenith_lm_covariance_is_finite():
    cov3 = g.lm_covariance_to_ecef(0.0, 0.0, np.eye(2) * 1e-6, -33.0, 26.0)
    assert np.all(np.isfinite(cov3))


def test_invalid_inputs_raise():
    with pytest.raises(TriangulationError):
        g.azel_to_ecef_unit_vector(np.nan, 10.0, 0.0, 0.0)
    with pytest.raises(CovarianceError):
        g.lm_covariance_to_ecef(0.1, 0.1, [[1e-6, 0.0], [0.0, -1e-6]], 0.0, 0.0)
    with pytest.raises(CovarianceError):
        g.lm_covariance_to_ecef(0.1, 0.1, [[np.nan, 0.0], [0.0, 1e-6]], 0.0, 0.0)


def test_pairwise_angles():
    s = np.array([[1, 0, 0], [0, 1, 0], [np.cos(0.1), np.sin(0.1), 0]])
    ang = g.pairwise_angles_deg(s)
    assert np.allclose(sorted(ang), sorted([90.0, np.degrees(0.1), 90.0 - np.degrees(0.1)]))
