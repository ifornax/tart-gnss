"""Image-domain direction estimator (validation path).

A dirty image is formed on a regular ``(l, m)`` pixel grid,

    I(l, m) = Re[ sum_k w_k V_k exp(+i phi_k(l, m)) ] / sum_k w_k,

using the same sign convention as the forward model (a point source at
``s0`` gives a peak ``Re C`` at ``s0``). No deconvolution is performed.

The estimate is the brightest pixel inside the search cap, refined to
sub-pixel precision by independent parabolic interpolation along ``l`` and
``m`` through the peak and its neighbours.

``image_snr = peak / (1.4826 MAD)`` of a coarser dirty image of the whole
visible hemisphere (one pixel per resolution element ``lambda / B_max``). It is a dynamic-range style
measure: the "noise" includes the sidelobes of every source in the sky.

Relation to the visibility estimator: for the linear loss, minimising the
profiled ``chi2`` is equivalent to maximising ``|sum w V e^{+i phi}|``, i.e.
the dirty-image *amplitude*. This estimator differs by using the real part
on a fixed pixel grid with interpolation instead of a continuous complex fit
with a free phase and robust loss, so it is an independent implementation,
not a statistically independent measurement.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from tart_gnss.stage3.estimator import MAD_TO_SIGMA, SearchRegion, resolution_rad, search_grid
from tart_gnss.stage3.geometry import lm_to_azel
from tart_gnss.stage3.visibility_model import VisibilitySamples, phases

#: Pixels per resolution element of the search-region image and of the noise image.
IMAGE_OVERSAMPLE = 6.0
NOISE_OVERSAMPLE = 1.0
PIXEL_CHUNK = 4096


@dataclass
class ImageDirectionEstimate:
    """Result of the image-domain estimator."""

    l: float
    m: float
    az_deg: float
    el_deg: float
    peak_value: float
    image_noise: float
    image_snr: float
    pixel_deg: float
    boundary: bool
    success: bool


def dirty_image_at(samples: VisibilitySamples, s: np.ndarray) -> np.ndarray:
    """Real dirty-image values at ENU directions ``s`` (G, 3)."""
    w = samples.weight
    wv = w * samples.vis
    out = np.empty(len(s))
    for i in range(0, len(s), PIXEL_CHUNK):
        ph = phases(samples.baseline_enu_m, samples.frequency_hz, s[i:i + PIXEL_CHUNK])
        out[i:i + PIXEL_CHUNK] = np.real(np.exp(1j * ph) @ wv) / w.sum()
    return out


def dirty_image(samples: VisibilitySamples, l_axis: np.ndarray, m_axis: np.ndarray) -> np.ndarray:
    """Dirty image on the grid ``m_axis x l_axis`` (rows = m). NaN where ``l^2 + m^2 >= 1``."""
    ll, mm = np.meshgrid(l_axis, m_axis)
    r2 = ll ** 2 + mm ** 2
    img = np.full(ll.shape, np.nan)
    inside = r2 < 1.0
    s = np.stack([ll[inside], mm[inside], np.sqrt(1.0 - r2[inside])], axis=1)
    img[inside] = dirty_image_at(samples, s)
    return img


def _parabolic_offset(ym: float, y0: float, yp: float) -> float:
    denom = ym - 2.0 * y0 + yp
    if not np.isfinite(denom) or denom >= 0:
        return 0.0
    return float(np.clip(0.5 * (ym - yp) / denom, -0.5, 0.5))


def image_noise(samples: VisibilitySamples) -> Tuple[float, float]:
    """Robust (MAD) rms of a hemisphere dirty image and its pixel size (rad)."""
    pix = resolution_rad(samples) / NOISE_OVERSAMPLE
    axis = np.arange(-1.0, 1.0 + pix, pix)
    img = dirty_image(samples, axis, axis)
    v = img[np.isfinite(img)]
    return float(MAD_TO_SIGMA * np.median(np.abs(v - np.median(v)))), pix


def image_direction(
    samples: VisibilitySamples,
    search_region: Optional[SearchRegion] = None,
) -> ImageDirectionEstimate:
    """Estimate source direction from a dirty-image peak inside ``search_region``.

    ``None`` searches the whole visible hemisphere.
    """
    region = search_region or SearchRegion.hemisphere()
    pix = resolution_rad(samples) / IMAGE_OVERSAMPLE
    # Bounding box of the cap in (l, m), from a coarse sampling of the cap itself.
    cap = search_grid(region, max(pix * 4, 1e-3))
    if len(cap) == 0:
        nan = float("nan")
        return ImageDirectionEstimate(nan, nan, nan, nan, nan, nan, nan, float(np.degrees(pix)), False, False)
    pad = 4 * pix
    l_axis = np.arange(max(cap[:, 0].min() - pad, -1.0), min(cap[:, 0].max() + pad, 1.0) + pix, pix)
    m_axis = np.arange(max(cap[:, 1].min() - pad, -1.0), min(cap[:, 1].max() + pad, 1.0) + pix, pix)
    img = dirty_image(samples, l_axis, m_axis)
    ll, mm = np.meshgrid(l_axis, m_axis)
    nn = np.sqrt(np.clip(1.0 - ll ** 2 - mm ** 2, 0.0, None))
    sep = region.separation_deg(np.stack([ll, mm, nn], axis=-1))
    masked = np.where(np.isfinite(img) & (sep <= region.radius_deg), img, -np.inf)
    j, i = np.unravel_index(int(np.argmax(masked)), masked.shape)
    peak = float(img[j, i])
    dl = dm = 0.0
    if 0 < i < len(l_axis) - 1:
        dl = _parabolic_offset(img[j, i - 1], peak, img[j, i + 1])
    if 0 < j < len(m_axis) - 1:
        dm = _parabolic_offset(img[j - 1, i], peak, img[j + 1, i])
    l, m = l_axis[i] + dl * pix, m_axis[j] + dm * pix
    if l * l + m * m >= 1.0:
        l, m = l_axis[i], m_axis[j]
    az, el = lm_to_azel(l, m)
    noise, _ = image_noise(samples)
    boundary = bool(region.radius_deg < 90 and sep[j, i] > region.radius_deg - np.degrees(pix))
    return ImageDirectionEstimate(
        l=float(l), m=float(m), az_deg=float(az), el_deg=float(el), peak_value=peak, image_noise=noise,
        image_snr=float(peak / noise) if noise > 0 else float("inf"), pixel_deg=float(np.degrees(pix)),
        boundary=boundary, success=bool(np.isfinite(peak)),
    )
