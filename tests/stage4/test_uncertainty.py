"""Covariance summaries, direction perturbation and the parametric bootstrap."""

import numpy as np
import pytest

from tart_gnss.stage4.synthetic import EXAMPLE_POSITION_M
from tart_gnss.stage4.triangulation import triangulate_weighted
from tart_gnss.stage4.uncertainty import (
    bootstrap_position,
    covariance_components,
    perturb_directions,
    triangulator,
)

from .helpers import lines


def test_covariance_components():
    c = np.array([[4.0, 1.0, 0.5], [1.0, 9.0, 0.2], [0.5, 0.2, 16.0]])
    out = covariance_components(c)
    assert out["cov_xy"] == 1.0 and out["cov_yz"] == 0.2 and out["cov_zz"] == 16.0
    assert (out["sigma_x_m"], out["sigma_y_m"], out["sigma_z_m"]) == (2.0, 3.0, 4.0)
    assert np.isnan(covariance_components(np.full((3, 3), np.nan))["sigma_x_m"])


def test_perturbation_statistics():
    r, s, c = lines(sigma_deg=0.2)
    rng = np.random.default_rng(1)
    draws = np.array([perturb_directions(s, c, rng)[0] for _ in range(4000)])
    assert np.allclose(np.linalg.norm(draws, axis=1), 1.0)
    ang = np.degrees(np.arccos(np.clip(draws @ s[0], -1, 1)))
    # Isotropic 2-D Gaussian: RMS angle = sqrt(2) sigma.
    assert np.sqrt(np.mean(ang ** 2)) == pytest.approx(0.2 * np.sqrt(2), rel=0.05)


def test_bootstrap_brackets_truth_and_matches_linear_sigma():
    r, s, c = lines()
    res = triangulate_weighted(r, s, c)
    bs = bootstrap_position(r, s, c, 400, np.random.default_rng(2), linear_sigma_m=res.sigma_ecef_m)
    assert bs.n_success == 400 and bs.failure_fraction == 0.0
    for i, ax in enumerate("xyz"):
        assert bs.percentiles[f"{ax}_p16"] < EXAMPLE_POSITION_M[i] < bs.percentiles[f"{ax}_p84"]
    assert bs.max_sigma_ratio < 1.3


def test_bootstrap_is_reproducible_with_seed():
    r, s, c = lines()
    a = bootstrap_position(r, s, c, 20, np.random.default_rng(7))
    b = bootstrap_position(r, s, c, 20, np.random.default_rng(7))
    assert a.percentiles == b.percentiles


def test_triangulator_names():
    assert triangulator("weighted") is triangulate_weighted
    with pytest.raises(ValueError):
        triangulator("magic")
