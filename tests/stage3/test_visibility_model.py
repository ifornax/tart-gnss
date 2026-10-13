"""Forward model, sign convention and complex-amplitude elimination."""

import numpy as np

from tart_gnss.stage3 import visibility_model as vm
from tart_gnss.stage3.geometry import azel_to_lm, azel_to_unit

F = 1575.42e6


def test_sign_convention_east_source_gives_negative_phase_on_east_baseline():
    """b = r_q - r_p pointing East and a source East of zenith: phase of V is -2 pi b.s/lambda."""
    b = np.array([[1.0, 0.0, 0.0]])
    l, m = azel_to_lm(90.0, 60.0)
    v = vm.predict_visibilities(b, np.array([F]), float(l), float(m), 1.0)
    expected = -2 * np.pi * F / 299792458.0 * 1.0 * np.cos(np.radians(60.0))
    assert np.isclose(np.angle(v[0]), np.angle(np.exp(1j * expected)))
    # The tart2ms UVW for this pair is -b, so V = exp(+2 pi i u l / lambda).
    uvw = -b
    assert np.isclose(np.angle(v[0]), np.angle(np.exp(2j * np.pi * F / 299792458.0 * uvw[0, 0] * l)))


def test_baselines_from_uvw():
    np.testing.assert_allclose(vm.baselines_from_uvw([[1.0, -2.0, 0.5]]), [[-1.0, 2.0, -0.5]])


def test_profiled_amplitude_and_chi2(layout):
    p, q = np.triu_indices(len(layout), 1)
    b = layout[q] - layout[p]
    freq = np.full(b.shape[0], F)
    l, m = azel_to_lm(150.0, 50.0)
    c = 0.3 - 0.7j
    vis = vm.predict_visibilities(b, freq, float(l), float(m), c)
    s = vm.VisibilitySamples(b, freq, vis, np.ones(len(vis)), p, q, np.zeros(len(vis)))
    chi2, c_hat = vm.profiled_chi2(s, azel_to_unit(150.0, 50.0))
    assert np.isclose(c_hat, c) and chi2 < 1e-20
    chi2_off, _ = vm.profiled_chi2(s, azel_to_unit(160.0, 50.0))
    assert chi2_off > 1e-3


def test_per_channel_frequency_is_used():
    b = np.array([[2.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    freqs = np.array([1.5e9, 1.6e9])
    v = vm.predict_visibilities(b, freqs, 0.5, 0.0, 1.0)
    assert not np.isclose(np.angle(v[0]), np.angle(v[1]))


def test_jacobian_matches_finite_differences(layout):
    p, q = np.triu_indices(len(layout), 1)
    b = layout[q] - layout[p]
    freq = np.full(b.shape[0], F)
    rng = np.random.default_rng(3)
    vis = rng.standard_normal(len(b)) + 1j * rng.standard_normal(len(b))
    s = vm.VisibilitySamples(b, freq, vis, rng.uniform(0.5, 2, len(b)), p, q, np.zeros(len(b)))
    x0 = np.array([0.2, -0.3, 0.4, 0.1])

    def f(x):
        r = vm.residuals(s, x[0], x[1], complex(x[2], x[3]))
        return np.concatenate([r.real, r.imag])

    jac = vm.residual_jacobian(s, x0[0], x0[1], complex(x0[2], x0[3]))
    for j in range(4):
        h = np.zeros(4)
        h[j] = 1e-7
        np.testing.assert_allclose((f(x0 + h) - f(x0 - h)) / 2e-7, jac[:, j], rtol=1e-4, atol=1e-6)
