"""Shared Stage-3 test helpers: synthetic TART-like data only (no real data)."""

import importlib.util

import numpy as np
import pytest

from tart_gnss.stage3 import synthetic as syn

HAS_CASACORE = importlib.util.find_spec("casacore") is not None
HAS_SCIPY = importlib.util.find_spec("scipy") is not None
needs_casacore = pytest.mark.skipif(not HAS_CASACORE, reason="python-casacore not installed")

LAT, LON, ALT = -20.26, 57.76, 20.0
T0 = syn.utc_to_mjd_seconds(2026, 6, 25, 9, 14, 0.0)


def make_samples(layout, az, el, amp=0.5 + 0.2j, noise=0.0, seed=0, flag_fraction=0.0, n_times=1, **kw):
    """Valid samples of a synthetic snapshot (or a few times)."""
    from tart_gnss.stage3 import loader

    times = T0 + np.arange(n_times, dtype=float)
    obs = syn.simulate_observation(layout, LAT, LON, ALT, times, az, el, amplitude=amp, noise_sigma=noise,
                                   flag_fraction=flag_fraction, rng=np.random.default_rng(seed), **kw)
    ds = syn.to_dataset(obs)
    return loader.select_samples(ds, np.arange(ds.n_rows)), ds
