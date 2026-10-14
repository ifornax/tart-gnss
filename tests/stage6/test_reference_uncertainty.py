"""Empirical TLE reference dispersion (spec section 16)."""

import numpy as np
import pytest

pytest.importorskip("skyfield")

from tart_gnss.stage6 import synthetic as syn  # noqa: E402
from tart_gnss.stage6.errors import ReferenceCovarianceError  # noqa: E402
from tart_gnss.stage6.propagation import propagate_tle  # noqa: E402
from tart_gnss.stage6.reference_uncertainty import LABEL, empirical_tle_covariance  # noqa: E402

from .helpers import L1, L2, reference, times  # noqa: E402

T = times(5)


def test_identical_sets_have_zero_dispersion():
    rc = empirical_tle_covariance([reference(), reference(), reference()], T)
    assert np.allclose(rc.covariance, 0.0) and rc.n == 3 and rc.label == LABEL == "empirical reference dispersion"


def test_dispersion_matches_sample_covariance():
    sets = [syn.perturbed_tle(L1, L2, d_mean_anomaly_deg=d) for d in (-0.002, 0.0, 0.001, 0.003)]
    refs = [reference(*s) for s in sets]
    rc = empirical_tle_covariance(refs, T)
    st = np.array([propagate_tle(r, T).state_gcrs for r in refs])
    assert np.allclose(rc.mean_state, st.mean(axis=0))
    assert np.allclose(rc.covariance[2], np.cov(st[:, 2, :].T), rtol=1e-9, atol=1e-6)
    # A mean-anomaly spread is mostly along-track: ~0.002 deg x 26 600 km ~ 1 km.
    assert 300.0 < np.median(rc.position_dispersion_m) < 3000.0


def test_needs_two_members():
    with pytest.raises(ReferenceCovarianceError):
        empirical_tle_covariance([reference()], T)
