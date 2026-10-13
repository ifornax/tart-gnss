"""Monte Carlo validation of the predicted position covariance (spec section 26)."""

import numpy as np
import pytest

from tart_gnss.stage4.synthetic import EXAMPLE_POSITION_M
from tart_gnss.stage4.triangulation import triangulate_lines, triangulate_weighted
from tart_gnss.stage4.uncertainty import bootstrap_position, monte_carlo_check

from .helpers import lines


@pytest.mark.parametrize("method", ["weighted", "unweighted"])
@pytest.mark.parametrize("sigma_deg", [0.05, 0.15])
def test_empirical_scatter_matches_predicted_covariance(method, sigma_deg):
    r, s, c = lines(sigma_deg=sigma_deg)
    mc = monte_carlo_check(EXAMPLE_POSITION_M, r, c, 1500, np.random.default_rng(42), method=method)
    # sigma ratio sampling error ~ 1/sqrt(2 n) ~ 2 %.
    assert np.all(np.abs(mc.sigma_ratio - 1.0) < 0.08), mc.sigma_ratio
    assert mc.los_axis_empirical_sigma_m / mc.los_axis_predicted_sigma_m == pytest.approx(1.0, abs=0.08)
    assert 2.7 < mc.mean_mahalanobis2 < 3.3
    assert np.all(np.abs(mc.mean_error_m / mc.predicted_sigma_m) < 0.15)


def test_unequal_station_errors():
    r, s, c = lines(sigmas=[0.03, 0.1, 0.2])
    mc = monte_carlo_check(EXAMPLE_POSITION_M, r, c, 1500, np.random.default_rng(1))
    assert np.all(np.abs(mc.sigma_ratio - 1.0) < 0.08)


def test_angular_uncertainty_propagates_linearly():
    r, s, c1 = lines(sigma_deg=0.05)
    _, _, c2 = lines(sigma_deg=0.10)
    a = triangulate_weighted(r, s, c1).sigma_ecef_m
    b = triangulate_weighted(r, s, c2).sigma_ecef_m
    assert np.allclose(b / a, 2.0, rtol=1e-6)
    assert np.allclose(triangulate_lines(r, s, c2).sigma_ecef_m / triangulate_lines(r, s, c1).sigma_ecef_m, 2.0)


def test_along_line_of_sight_is_least_constrained():
    r, s, c = lines()
    cov = triangulate_weighted(r, s, c).covariance_ecef_m2
    w, v = np.linalg.eigh(cov)
    los = s.mean(axis=0) / np.linalg.norm(s.mean(axis=0))
    assert abs(v[:, -1] @ los) > 0.99
    assert np.sqrt(w[-1]) > 5 * np.sqrt(w[-2])


def test_nonlinear_regime_detected_by_bootstrap():
    # At 2 deg the angular errors approach the crossing angles (3-9 deg): the estimator becomes biased.
    r, s, c = lines(sigma_deg=2.0)
    res = triangulate_weighted(r, s, c)
    bs = bootstrap_position(r, s, c, 400, np.random.default_rng(3), linear_sigma_m=res.sigma_ecef_m,
                            point_estimate_m=res.position_ecef_m)
    assert bs.max_median_shift_sigma > 0.25
    r, s, c = lines(sigma_deg=0.1)
    res = triangulate_weighted(r, s, c)
    bs = bootstrap_position(r, s, c, 400, np.random.default_rng(3), linear_sigma_m=res.sigma_ecef_m,
                            point_estimate_m=res.position_ecef_m)
    assert bs.max_median_shift_sigma < 0.2 and bs.max_sigma_ratio < 1.3
