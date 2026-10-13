"""Sight-line intersection: exact recovery, degenerate geometry, ranges, weighting, errors."""

import numpy as np
import pytest

from tart_gnss.stage4.errors import CovarianceError, InsufficientStationsError, TriangulationError
from tart_gnss.stage4.synthetic import EXAMPLE_POSITION_M
from tart_gnss.stage4.triangulation import geometry_metrics, triangulate_lines, triangulate_weighted

from .helpers import lines


@pytest.mark.parametrize("method", ["unweighted", "weighted"])
def test_noise_free_exact_recovery(method):
    r, s, c = lines()
    res = triangulate_lines(r, s, c) if method == "unweighted" else triangulate_weighted(r, s, c)
    assert res.success and res.converged
    assert np.linalg.norm(res.position_ecef_m - EXAMPLE_POSITION_M) < 1e-3
    assert np.allclose(res.line_residuals_m, 0.0, atol=1e-3)
    assert np.allclose(res.station_ranges_m, np.linalg.norm(EXAMPLE_POSITION_M - r, axis=1), rtol=1e-12)
    assert res.chi2 == pytest.approx(0.0, abs=1e-12)
    assert res.dof == 3
    assert res.covariance_ok and np.all(res.sigma_ecef_m > 0)
    assert res.solver == "solve" and res.matrix_rank == 3


def test_two_skew_lines_midpoint():
    r = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    s = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    res = triangulate_lines(r, s)
    assert np.allclose(res.position_ecef_m, [0.0, 0.0, 0.5])
    assert np.allclose(res.line_residuals_m, 0.5)
    assert res.min_crossing_angle_deg == pytest.approx(90.0)
    assert res.chi2 is None and not res.covariance_ok


def test_parallel_lines_are_rank_deficient():
    r = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    s = np.tile([0.0, 0.0, 1.0], (3, 1))
    res = triangulate_lines(r, s)
    assert res.matrix_rank == 2 and not res.success
    assert res.solver == "lstsq"
    assert res.min_crossing_angle_deg == pytest.approx(0.0, abs=1e-9)
    assert not np.isfinite(res.condition_number) or res.condition_number > 1e12


def test_reversed_lines_give_negative_ranges():
    r, s, _ = lines()
    res = triangulate_lines(r, -s)
    # Same lines, so the same point, but it lies behind every station.
    assert np.linalg.norm(res.position_ecef_m - EXAMPLE_POSITION_M) < 1e-3
    assert np.all(res.station_ranges_m < 0)


def test_condition_number_grows_as_lines_become_parallel():
    conds = []
    for ang in (10.0, 1.0, 0.1):
        a = np.radians(ang)
        s = np.array([[0, 0, 1.0], [np.sin(a), 0, np.cos(a)], [0, np.sin(a), np.cos(a)]])
        conds.append(geometry_metrics(s)[0])
    assert conds[0] < conds[1] < conds[2]


def test_weighting_downweights_noisy_station():
    rng = np.random.default_rng(5)
    r, s, c = lines(sigmas=[0.05, 0.05, 1.0])
    err_w, err_u = [], []
    from .helpers import perturb

    for _ in range(200):
        sp = perturb(s, c, rng)
        err_w.append(np.linalg.norm(triangulate_weighted(r, sp, c).position_ecef_m - EXAMPLE_POSITION_M))
        err_u.append(np.linalg.norm(triangulate_lines(r, sp, c).position_ecef_m - EXAMPLE_POSITION_M))
    assert np.median(err_w) < 0.7 * np.median(err_u)


def test_unweighted_and_weighted_covariances_agree_for_equal_errors():
    r, s, c = lines()
    cu = triangulate_lines(r, s, c).covariance_ecef_m2
    cw = triangulate_weighted(r, s, c).covariance_ecef_m2
    # Equal angular errors and similar ranges: the two estimators are nearly equivalent.
    assert np.allclose(np.sqrt(np.diag(cu)), np.sqrt(np.diag(cw)), rtol=0.1)


def test_invalid_inputs():
    r, s, c = lines()
    with pytest.raises(InsufficientStationsError):
        triangulate_lines(r[:1], s[:1])
    bad = s.copy()
    bad[0, 0] = np.nan
    with pytest.raises(TriangulationError):
        triangulate_lines(r, bad)
    with pytest.raises(TriangulationError):
        triangulate_lines(r, 2 * s)
    with pytest.raises(TriangulationError):
        triangulate_lines(r, s[:2])
    with pytest.raises(CovarianceError):
        triangulate_weighted(r, s, None)
    zero = c.copy()
    zero[1] = 0.0
    with pytest.raises(CovarianceError):
        triangulate_weighted(r, s, zero)
