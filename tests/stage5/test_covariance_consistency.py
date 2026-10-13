"""Monte Carlo coverage and NEES (specification sections 33-34; synthetic truth only)."""

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")
from scipy.stats import chi2  # noqa: E402

from tart_gnss.stage5 import synthetic as syn  # noqa: E402

from .helpers import run_track, times, truth  # noqa: E402

N_RUNS = 60


@pytest.fixture(scope="module")
def mc():
    t = times(20, step=15.0)
    x = truth(t)
    rng = np.random.default_rng(2026)
    nees, dpos, dvel, nees_f = [], [], [], []
    for _ in range(N_RUNS):
        res = run_track(syn.stage4_rows(t, x, 3e4, 2e3, rng))
        for k in (0, 10, 19):
            e = res.smoothed.x_smooth[k] - x[k]
            p = res.smoothed.p_smooth[k]
            nees.append(e @ np.linalg.solve(p, e))
            dpos.append(e[:3] @ np.linalg.solve(p[:3, :3], e[:3]))
            dvel.append(e[3:] @ np.linalg.solve(p[3:, 3:], e[3:]))
        ef = res.filtered.x_filt[5] - x[5]
        nees_f.append(ef @ np.linalg.solve(res.filtered.p_filt[5], ef))
    return np.array(nees), np.array(dpos), np.array(dvel), np.array(nees_f)


def test_nees_is_chi2_6(mc):
    nees, _, _, nees_f = mc
    # Epochs within one run are correlated, so judge by the per-run mean (sd of the mean ~ sqrt(12/N)).
    lim = 3.5 * np.sqrt(12.0 / N_RUNS)
    assert abs(nees.mean() - 6.0) < lim
    # Early filtered states: the initialisation window must not be counted twice.
    assert abs(nees_f.mean() - 6.0) < lim


@pytest.mark.parametrize("k", [1, 2, 3])
def test_position_and_velocity_coverage(mc, k):
    _, dpos, dvel, _ = mc
    expected = chi2.cdf(k * k, 3)  # 0.199, 0.739, 0.971 for the 1, 2, 3-sigma ellipsoids
    tol = 0.12 if k < 3 else 0.05
    assert abs(np.mean(dpos <= k * k) - expected) < tol
    assert abs(np.mean(dvel <= k * k) - expected) < tol
