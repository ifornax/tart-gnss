"""Correlation statistics: time-trend control, effective sample size, Holm, runs."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage7.statistics import correlation, holm, runs  # noqa: E402


def test_shared_trend_is_not_a_correlation():
    rng = np.random.default_rng(0)
    t = np.arange(60.0)
    x = t + rng.normal(0, 3, 60)
    y = 2 * t + rng.normal(0, 3, 60)
    c = correlation(x, y, t)
    assert c["coefficient"] > 0.9  # the naive rank correlation is fooled by the trend
    assert abs(c["partial_coefficient"]) < 0.4


def test_smooth_series_have_small_effective_samples():
    rng = np.random.default_rng(1)
    t = np.arange(80.0)
    x = np.cumsum(rng.normal(size=80))
    y = np.cumsum(rng.normal(size=80))
    c = correlation(x, y, t)
    assert c["n_effective"] < 0.5 * c["n_samples"]
    assert c["p_value"] > c["p_value_naive"] or not np.isfinite(c["p_value_naive"])


def test_real_relation_survives():
    rng = np.random.default_rng(2)
    t = np.arange(50.0)
    x = rng.normal(size=50)
    y = x + 0.3 * rng.normal(size=50)
    c = correlation(x, y, t)
    assert c["partial_coefficient"] > 0.8 and c["p_value"] < 1e-6


def test_degenerate_inputs():
    assert np.isnan(correlation(np.ones(10), np.arange(10.0))["coefficient"])
    assert correlation([1, 2], [3, 4])["n_samples"] == 2
    with pytest.raises(ValueError):
        correlation(np.arange(10.0), np.arange(10.0), method="magic")


def test_holm_and_runs():
    p = holm([0.01, 0.04, np.nan, 0.03])
    assert p[0] == pytest.approx(0.03) and p[3] == pytest.approx(0.06) and p[1] == pytest.approx(0.06)
    assert np.isnan(p[2])
    r = runs([False, True, True, False, True], np.arange(5.0) * 10)
    assert r == [(1, 2, 2, 10.0), (4, 4, 1, 0.0)]
