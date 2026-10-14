"""SGP4 at exact epochs; age convention."""

import numpy as np
import pytest

pytest.importorskip("skyfield")

from tart_gnss.stage6.propagation import propagate_tle, tle_age_days  # noqa: E402

from .helpers import TLE_EPOCH_S, reference, times  # noqa: E402


def test_exact_epochs_no_interpolation():
    t = np.sort(times(8) + np.random.default_rng(0).uniform(-3, 3, 8))
    tr = propagate_tle(reference(), t)
    for k in (0, 3, 7):
        one = propagate_tle(reference(), t[k:k + 1])
        assert np.allclose(one.state_gcrs[0], tr.state_gcrs[k], atol=1e-6)
    assert np.array_equal(tr.times_s, t)


def test_age_sign_convention():
    ref = reference()
    assert tle_age_days([TLE_EPOCH_S + 86400], ref)[0] == pytest.approx(1.0)
    assert tle_age_days([TLE_EPOCH_S - 2 * 86400], ref)[0] == pytest.approx(-2.0)
