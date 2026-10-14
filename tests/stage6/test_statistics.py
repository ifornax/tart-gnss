"""Descriptive statistics, trends, autocorrelation."""

import numpy as np
import pytest

from tart_gnss.stage6.statistics import describe, effective_sample_size, lag1_autocorrelation, trend


def test_describe():
    d = describe([1.0, 2.0, 3.0, 4.0, np.nan, 100.0])
    assert d["n"] == 5 and d["median"] == 3.0 and d["min"] == 1.0 and d["max"] == 100.0
    assert d["mad_sigma"] == pytest.approx(1.4826 * 1.0)
    assert d["rms"] == pytest.approx(np.sqrt(np.mean(np.array([1, 4, 9, 16, 10000.0]))))
    assert describe([])["n"] == 0


def test_weighted_trend():
    t = np.arange(20.0) * 10
    y = 5.0 + 0.3 * t
    tr = trend(t, y, np.ones(20))
    assert tr["slope"] == pytest.approx(0.3) and tr["intercept"] == pytest.approx(5.0)
    assert tr["change_over_span"] == pytest.approx(0.3 * 190)
    # Weights matter: a point with a huge sigma barely moves the fit.
    y2 = y.copy()
    y2[5] += 1000
    s = np.ones(20)
    s[5] = 1e6
    assert trend(t, y2, s)["slope"] == pytest.approx(0.3, abs=1e-6)


def test_autocorrelation_and_effective_size():
    rng = np.random.default_rng(0)
    assert abs(lag1_autocorrelation(rng.normal(size=5000))) < 0.05
    smooth = np.sin(np.linspace(0, 3, 100))
    rho = lag1_autocorrelation(smooth)
    assert rho > 0.9
    assert effective_sample_size(100, rho) < 5
    assert effective_sample_size(100, -0.2) == 100
