"""Statistics helpers."""

import numpy as np
import pytest

from tart_gnss.campaign.statistics import compare_groups, dependence, describe, sample_tier, wilson_interval


def test_describe_and_tiers():
    d = describe([1.0, 2, 3, np.nan, 4])
    assert d["n"] == 4 and d["median"] == 2.5 and d["p95"] == pytest.approx(3.85)
    assert [sample_tier(n) for n in (2, 5, 12, 40)] == ["descriptive_only", "preliminary", "basic_inference",
                                                       "grouped_comparisons"]


def test_kruskal_and_effect():
    rng = np.random.default_rng(0)
    g = {"a": rng.normal(0, 1, 10), "b": rng.normal(3, 1, 10), "c": rng.normal(6, 1, 10)}
    c = compare_groups(g)
    assert c["test"] == "Kruskal-Wallis" and c["p_value"] < 1e-4 and c["effect_size"] > 0.5


def test_small_sample_p_values_are_exact():
    # 3 vs 13 runs, perfectly separated, tied within groups: of the C(16, 3) = 560 label assignments only the
    # observed one (all three low values in group a) is as extreme (the mirror extreme is impossible with these
    # ties), so the exact p is 1/560 (asymptotic Mann-Whitney gives ~1e-4).
    c = compare_groups({"a": [12.0] * 3, "b": [20.0] * 13})
    assert c["p_method"] == "exact permutation" and c["p_value"] == pytest.approx(1 / 560)
    c = compare_groups({"a": [1.0, 2, 3], "b": [4.0, 5, 6]})  # untied: both extremes count, 2 / C(6, 3)
    assert c["p_value"] == pytest.approx(2 / 20)
    # three points with a perfect rank correlation: exact two-sided p = 2 / 3! = 1/3
    d = dependence([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])
    assert d["spearman_r"] == pytest.approx(-1.0) and d["spearman_p"] == pytest.approx(1 / 3)
    assert d["note"] is not None


def test_dependence_and_wilson():
    x = np.arange(10.0)
    d = dependence(x, 2 * x + 1)
    assert d["spearman_r"] == pytest.approx(1.0) and d["theil_sen_slope"] == pytest.approx(2.0)
    assert dependence([1, 2], [3, 4])["spearman_r"] is None
    lo, hi = wilson_interval(5, 10)
    assert lo < 0.5 < hi and wilson_interval(0, 0) == (None, None)
