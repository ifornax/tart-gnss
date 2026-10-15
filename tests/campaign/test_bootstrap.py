"""Run-level bootstrap: coverage, reproducibility, clustering (specification section 48)."""

import numpy as np
import pandas as pd
import pytest

from tart_gnss.campaign.bootstrap import bootstrap_run_level, median_of
from tart_gnss.campaign.errors import BootstrapError


def test_coverage_of_the_median():
    rng = np.random.default_rng(5)
    hits = 0
    reps = 200
    for i in range(reps):
        runs = pd.DataFrame({"run_id": range(25), "v": rng.lognormal(0, 1, 25)})  # true median 1
        b = bootstrap_run_level(runs, median_of("v"), 400, seed=i)
        hits += b.ci_low <= 1.0 <= b.ci_high
    assert 0.88 <= hits / reps <= 0.99


def test_reproducible_with_seed():
    runs = pd.DataFrame({"run_id": range(10), "v": np.arange(10.0)})
    a = bootstrap_run_level(runs, median_of("v"), 200, seed=3)
    b = bootstrap_run_level(runs, median_of("v"), 200, seed=3)
    assert (a.ci_low, a.ci_high) == (b.ci_low, b.ci_high)


def test_cluster_resampling_widens_correlated_epochs():
    rng = np.random.default_rng(1)
    run_means = rng.normal(0, 1, 8)
    epochs = pd.DataFrame([{"run_id": r, "v": m + rng.normal(0, 0.05)} for r, m in enumerate(run_means)
                           for _ in range(200)])
    runs = pd.DataFrame({"run_id": range(8)})
    clustered = bootstrap_run_level(runs, median_of("v"), 300, seed=2, rows=epochs)
    naive = bootstrap_run_level(epochs.assign(run_id=range(len(epochs))), median_of("v"), 300, seed=2)
    assert (clustered.ci_high - clustered.ci_low) > 5 * (naive.ci_high - naive.ci_low)


def test_errors_and_small_samples():
    runs = pd.DataFrame({"run_id": [0], "v": [1.0]})
    assert np.isnan(bootstrap_run_level(runs, median_of("v"), 100).ci_low)
    with pytest.raises(BootstrapError):
        bootstrap_run_level(runs, median_of("v"), -1)
