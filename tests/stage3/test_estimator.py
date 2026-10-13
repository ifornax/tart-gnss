"""Estimator mechanics: search grid, options, failure modes."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage3 import estimator as est_mod  # noqa: E402
from tart_gnss.stage3.errors import InsufficientBaselinesError  # noqa: E402

from .helpers import make_samples  # noqa: E402


def test_search_grid_stays_inside_cap_and_above_horizon():
    region = est_mod.SearchRegion.around_azel(100.0, 5.0, 10.0)
    grid = est_mod.search_grid(region, np.radians(0.5))
    assert len(grid) > 100
    assert np.all(region.separation_deg(grid) <= 10.0 + 1e-9)
    assert np.all(grid[:, 2] > 0)
    np.testing.assert_allclose(np.linalg.norm(grid, axis=1), 1.0)


def test_options_validation():
    with pytest.raises(ValueError):
        est_mod.EstimatorOptions(robust_loss="cauchy-ish")
    with pytest.raises(ValueError):
        est_mod.EstimatorOptions(search_radius_deg=0)


def test_too_few_samples(layout):
    sel, _ = make_samples(layout, 200.0, 45.0)
    with pytest.raises(InsufficientBaselinesError):
        est_mod.estimate_direction(sel.samples.subset(np.arange(2)), 200.0, 45.0)


def test_statistics_definitions(layout):
    sel, _ = make_samples(layout, 200.0, 45.0, amp=1.0, noise=0.5, seed=2)
    e = est_mod.estimate_direction(sel.samples, 200.0, 45.0, robust_loss="linear")
    assert e.n_samples == 276 and e.n_baselines_used == 276
    assert np.isclose(e.reduced_chi2, e.chi2 / (2 * 276 - 4))
    assert abs(e.sigma_vis - 0.5) < 0.1  # robust noise recovers the injected per-component sigma
    assert abs(e.snr - abs(e.complex_amplitude) * np.sqrt(276 / 2) / e.sigma_vis) < 1e-9
    assert 0.8 < e.residual_tail_ratio < 1.2
    assert e.max_abs_residual >= e.median_abs_residual


def test_robust_f_scale_tracks_noise_level(layout):
    sel, _ = make_samples(layout, 200.0, 45.0, amp=1.0, noise=0.01, seed=3)
    e = est_mod.estimate_direction(sel.samples, 200.0, 45.0)
    assert e.robust_loss == "huber" and 0.005 < e.f_scale < 0.03
