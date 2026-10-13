"""Coordinate conventions: ENU, (l, m), azimuth from North through East."""

import numpy as np
import pytest

from tart_gnss.stage3 import geometry as g
from tart_gnss.stage3.errors import CoordinateTransformError


@pytest.mark.parametrize("az,el,expected", [
    (0.0, 0.0, (0, 1, 0)),     # North horizon
    (90.0, 0.0, (1, 0, 0)),    # East horizon
    (180.0, 0.0, (0, -1, 0)),  # South
    (270.0, 0.0, (-1, 0, 0)),  # West
    (123.0, 90.0, (0, 0, 1)),  # zenith
])
def test_cardinal_directions(az, el, expected):
    np.testing.assert_allclose(g.azel_to_unit(az, el), expected, atol=1e-12)


def test_lm_azel_round_trip():
    rng = np.random.default_rng(0)
    az = rng.uniform(0, 360, 200)
    el = rng.uniform(1, 89, 200)
    l, m = g.azel_to_lm(az, el)
    az2, el2 = g.lm_to_azel(l, m)
    np.testing.assert_allclose(g.wrap_azimuth_near(az2, 0) % 360, az, atol=1e-9)
    np.testing.assert_allclose(el2, el, atol=1e-9)


def test_east_is_positive_l():
    l, m = g.azel_to_lm(90.0, 30.0)
    assert l > 0 and abs(m) < 1e-12
    l, m = g.azel_to_lm(0.0, 30.0)
    assert m > 0 and abs(l) < 1e-12


def test_invalid_lm_raises():
    with pytest.raises(CoordinateTransformError):
        g.lm_to_unit(0.9, 0.9)


def test_ecef_enu_round_trip_and_up_axis():
    lat, lon, h = -33.319413, 26.507688, 642.26
    ref = g.geodetic_to_ecef(lat, lon, h)
    lat2, lon2, h2 = g.ecef_to_geodetic(ref)
    assert abs(lat2 - lat) < 1e-9 and abs(lon2 - lon) < 1e-9 and abs(h2 - h) < 1e-4
    enu = np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, 0.0]])
    back = g.ecef_to_enu(g.enu_to_ecef(enu, lat, lon, h), lat, lon, reference_m=ref)
    np.testing.assert_allclose(back, enu, atol=1e-6)
    up = g.ecef_to_enu(g.geodetic_to_ecef(lat, lon, h + 10.0), lat, lon, reference_m=ref)[0]
    np.testing.assert_allclose(up, [0, 0, 10.0], atol=1e-6)


def test_angular_separation_and_offsets():
    s0 = g.azel_to_unit(200.0, 45.0)
    for rho_deg, psi in ((0.0, 0.0), (5.0, 0.3), (10.0, 2.0)):
        x, y = np.radians(rho_deg) * np.cos(psi), np.radians(rho_deg) * np.sin(psi)
        s = g.offset_directions(s0, np.array(x), np.array(y))
        assert abs(np.linalg.norm(s) - 1) < 1e-12
        assert abs(g.angular_separation_deg(s, s0) - rho_deg) < 1e-9


def test_azel_jacobian_matches_finite_differences():
    l, m = 0.3, -0.4
    jac = g.azel_jacobian(l, m)
    h = 1e-7
    az0, el0 = g.lm_to_azel(l, m)
    for j, (dl, dm) in enumerate(((h, 0), (0, h))):
        az1, el1 = g.lm_to_azel(l + dl, m + dm)
        np.testing.assert_allclose(np.radians([az1 - az0, el1 - el0]) / h, jac[:, j], rtol=1e-5)
    with pytest.raises(CoordinateTransformError):
        g.azel_jacobian(0.0, 0.0)


def test_wrap_azimuth_near():
    np.testing.assert_allclose(g.wrap_azimuth_near([359.0, 1.0], 0.0), [-1.0, 1.0])
