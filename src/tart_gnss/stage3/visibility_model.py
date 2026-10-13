"""Point-source visibility forward model.

Measurement equation (the single place where the sign convention lives)
-----------------------------------------------------------------------
For antennas ``p = ANTENNA1`` and ``q = ANTENNA2`` with ENU positions
``r_p``, ``r_q`` (metres), define the baseline

    b_pq = r_q - r_p.

A far-field point source in ENU direction ``s`` contributes

    V_pq(nu) = C exp(-i phi_pq),     phi_pq = 2 pi (nu / c) b_pq . s,

where ``C = a + i b`` is a complex nuisance amplitude.

The Measurement Sets written by ``tart2ms`` store ``UVW = r_p - r_q`` (i.e.
``UVW = -b_pq``) for a zenith phase centre, so in terms of the stored UVW the
same equation reads ``V = C exp(+2 pi i (u l + v m + w n) / lambda)``. This
was verified on real data (``docs/stage3_direction_estimation.md``); the loader
checks the ``UVW = -b`` relation for every Measurement Set it reads.

Complex-amplitude elimination
-----------------------------
For a fixed direction the weighted least-squares amplitude is

    C_hat = sum_k w_k V_k exp(+i phi_k) / sum_k w_k,

and the profiled objective is

    chi2(l, m) = sum_k w_k |V_k|^2 - |sum_k w_k V_k exp(+i phi_k)|^2 / sum_k w_k,

so minimising ``chi2`` over ``(l, m)`` is a two-dimensional problem.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from tart_gnss.stage3.geometry import SPEED_OF_LIGHT_M_S, lm_to_unit


@dataclass
class VisibilitySamples:
    """Flattened, valid complex visibility samples for one fit (one time bin).

    Each element ``k`` is one (row, channel) sample.

    Attributes
    ----------
    baseline_enu_m
        ``b_k = r_q - r_p`` in metres, shape (N, 3).
    frequency_hz
        Channel frequency of each sample, shape (N,).
    vis
        Complex visibilities, shape (N,).
    weight
        Non-negative statistical weights, shape (N,).
    antenna1, antenna2
        Antenna indices of each sample (used to resample baselines).
    time_mjd_seconds
        MS ``TIME`` of each sample (MJD seconds).
    """

    baseline_enu_m: np.ndarray
    frequency_hz: np.ndarray
    vis: np.ndarray
    weight: np.ndarray
    antenna1: np.ndarray
    antenna2: np.ndarray
    time_mjd_seconds: np.ndarray

    @property
    def n(self) -> int:
        """Number of samples."""
        return int(self.vis.size)

    def subset(self, index: np.ndarray) -> "VisibilitySamples":
        """Samples selected by an index or boolean mask (indices may repeat)."""
        return VisibilitySamples(
            self.baseline_enu_m[index], self.frequency_hz[index], self.vis[index],
            self.weight[index], self.antenna1[index], self.antenna2[index], self.time_mjd_seconds[index],
        )

    def baseline_ids(self) -> np.ndarray:
        """Integer id of each sample's (antenna1, antenna2) pair."""
        return self.antenna1.astype(np.int64) * 100000 + self.antenna2.astype(np.int64)

    def max_baseline_m(self) -> float:
        """Longest baseline among the samples (m)."""
        return float(np.linalg.norm(self.baseline_enu_m, axis=1).max()) if self.n else 0.0


def baselines_from_uvw(uvw_m: np.ndarray) -> np.ndarray:
    """Convert ``tart2ms`` UVW (``r_p - r_q``, zenith phase centre) to ``b = r_q - r_p``."""
    return -np.asarray(uvw_m, dtype=float)


def wavenumber(frequency_hz) -> np.ndarray:
    """``2 pi nu / c`` in radians per metre."""
    return 2.0 * np.pi * np.asarray(frequency_hz, dtype=float) / SPEED_OF_LIGHT_M_S


def phases(baseline_enu_m: np.ndarray, frequency_hz: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Geometric phases ``phi = 2 pi nu/c b . s``.

    ``s`` may be a single direction (3,) giving shape (N,), or a stack (G, 3)
    giving shape (G, N).
    """
    s = np.asarray(s, dtype=float)
    proj = s @ np.asarray(baseline_enu_m, dtype=float).T
    return proj * wavenumber(frequency_hz)


def predict_visibilities(
    baseline_enu_m: np.ndarray,
    frequency_hz: np.ndarray,
    l: float,
    m: float,
    complex_amplitude: complex,
) -> np.ndarray:
    """Predict far-field complex visibilities ``C exp(-i phi)`` of a point source at ``(l, m)``.

    ``frequency_hz`` gives each sample's own channel frequency, so no single
    frequency is assumed.
    """
    s = lm_to_unit(l, m)
    return complex(complex_amplitude) * np.exp(-1j * phases(baseline_enu_m, frequency_hz, s))


def best_amplitude(vis: np.ndarray, weight: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Weighted least-squares ``C_hat`` for fixed phases (vectorised over leading axes of ``phi``)."""
    return (np.exp(1j * phi) @ (weight * vis)) / weight.sum()


def profiled_chi2(samples: VisibilitySamples, s: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Profiled ``chi2`` and ``C_hat`` for direction(s) ``s`` (see module docstring)."""
    phi = phases(samples.baseline_enu_m, samples.frequency_hz, s)
    w = samples.weight
    c_hat = best_amplitude(samples.vis, w, phi)
    chi2 = float(np.sum(w * np.abs(samples.vis) ** 2)) - np.abs(c_hat) ** 2 * w.sum()
    return np.maximum(chi2, 0.0), c_hat


def residuals(samples: VisibilitySamples, l: float, m: float, c: complex) -> np.ndarray:
    """Complex weighted residuals ``sqrt(w) (V - C exp(-i phi))``."""
    model = predict_visibilities(samples.baseline_enu_m, samples.frequency_hz, l, m, c)
    return np.sqrt(samples.weight) * (samples.vis - model)


def residual_jacobian(samples: VisibilitySamples, l: float, m: float, c: complex) -> np.ndarray:
    """Jacobian of the stacked real residual vector ``[Re r, Im r]`` w.r.t. ``(l, m, Re C, Im C)``.

    With ``n = sqrt(1 - l^2 - m^2)``, ``k = 2 pi nu / c`` and model ``M = C e^{-i phi}``:

        dphi/dl = k (b_E - b_U l / n),     dphi/dm = k (b_N - b_U m / n)
        dr/dl = i sqrt(w) M dphi/dl,       dr/dm = i sqrt(w) M dphi/dm
        dr/da = -sqrt(w) e^{-i phi},       dr/db = -i sqrt(w) e^{-i phi}

    Returns shape (2N, 4).
    """
    b = samples.baseline_enu_m
    k = wavenumber(samples.frequency_hz)
    n = np.sqrt(max(1.0 - l * l - m * m, 1e-12))
    s = np.array([l, m, n])
    e = np.exp(-1j * k * (b @ s))
    sw = np.sqrt(samples.weight)
    model = c * e
    dphi_dl = k * (b[:, 0] - b[:, 2] * l / n)
    dphi_dm = k * (b[:, 1] - b[:, 2] * m / n)
    cols = [1j * sw * model * dphi_dl, 1j * sw * model * dphi_dm, -sw * e, -1j * sw * e]
    jc = np.stack(cols, axis=1)
    return np.concatenate([jc.real, jc.imag], axis=0)


def simulate_visibilities(
    baseline_enu_m: np.ndarray,
    frequency_hz: np.ndarray,
    l: float,
    m: float,
    complex_amplitude: complex,
    noise_sigma: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Synthetic visibilities from the same forward model plus circular complex Gaussian noise.

    ``noise_sigma`` is the standard deviation of each of the real and imaginary parts.
    """
    v = predict_visibilities(baseline_enu_m, frequency_hz, l, m, complex_amplitude)
    if noise_sigma > 0:
        rng = rng or np.random.default_rng()
        v = v + noise_sigma * (rng.standard_normal(v.shape) + 1j * rng.standard_normal(v.shape))
    return v
