"""Covariance propagation and the baseline bootstrap."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage3.estimator import estimate_direction  # noqa: E402
from tart_gnss.stage3.geometry import azel_jacobian  # noqa: E402
from tart_gnss.stage3.uncertainty import (  # noqa: E402
    bootstrap_direction,
    estimate_direction_uncertainty,
    propagate_to_azel,
)

from .helpers import make_samples  # noqa: E402


def test_propagation_matches_jacobian():
    l, m = 0.4, -0.3
    cov = np.array([[1e-6, 2e-7], [2e-7, 4e-6]])
    out = propagate_to_azel(l, m, cov)
    g = azel_jacobian(l, m)
    expected = np.degrees(1) ** 2 * g @ cov @ g.T
    assert np.isclose(out.sigma_az_deg ** 2, expected[0, 0])
    assert np.isclose(out.sigma_el_deg ** 2, expected[1, 1])
    assert np.isclose(out.cov_az_el, expected[0, 1])
    assert not out.zenith_az_uncertain


def test_near_zenith_azimuth_flagged():
    out = propagate_to_azel(1e-4, 0.0, np.eye(2) * 1e-6)
    assert out.zenith_az_uncertain and np.isnan(out.sigma_az_deg)


def test_covariance_shrinks_with_noise(layout):
    sig = []
    for noise in (0.1, 0.4):
        sel, _ = make_samples(layout, 150.0, 40.0, amp=1.0, noise=noise, seed=1)
        e = estimate_direction(sel.samples, 150.0, 40.0)
        sig.append(estimate_direction_uncertainty(e, sel.samples).sigma_el_deg)
    assert 2.5 < sig[1] / sig[0] < 6


def test_bootstrap_brackets_estimate_and_agrees_with_hessian(layout):
    sel, _ = make_samples(layout, 330.0, 55.0, amp=1.0, noise=0.5, seed=4)
    e = estimate_direction(sel.samples, 330.0, 55.0)
    cov = estimate_direction_uncertainty(e, sel.samples)
    boot = bootstrap_direction(e, sel.samples, 60, np.random.default_rng(0))
    assert boot.n_success >= 55
    assert boot.az_p16 <= e.az_deg + 0.5 * cov.sigma_az_deg and boot.az_p84 >= e.az_deg - 0.5 * cov.sigma_az_deg
    width = (boot.el_p84 - boot.el_p16) / 2
    assert 0.4 < width / cov.sigma_el_deg < 2.5
    assert bootstrap_direction(e, sel.samples, 0, np.random.default_rng(0)) is None


def test_bootstrap_unwraps_azimuth_near_north(layout):
    sel, _ = make_samples(layout, 0.2, 50.0, amp=1.0, noise=0.5, seed=5)
    e = estimate_direction(sel.samples, 0.0, 50.0)
    boot = bootstrap_direction(e, sel.samples, 30, np.random.default_rng(1))
    assert boot.az_p84 - boot.az_p16 < 20  # no 0/360 wrap blow-up
