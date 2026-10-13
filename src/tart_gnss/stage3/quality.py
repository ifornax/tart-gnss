"""Deterministic quality flags for Stage-3 estimates.

Flags (several may apply; an estimate with no flag is ``GOOD``):

``LOW_BASELINE_COUNT``      fewer than ``min_baselines`` valid baselines: no fit is attempted
``HIGH_FLAGGED_FRACTION``   more than ``max_flagged_fraction`` of the bin's cross-correlation samples excluded
``FIT_FAILED``              the optimiser did not converge, or the covariance is singular/non-finite
``BOUNDARY_SOLUTION``       the solution (or coarse optimum) is on the search-cap edge or the horizon
``LOW_SNR``                 ``snr < min_snr``
``HIGH_RESIDUAL``           ``residual_tail_ratio > max_residual_tail_ratio`` (heavy-tailed residuals:
                            outliers or unmodelled sources); with calibrated weights, also
                            ``reduced_chi2 > max_reduced_chi2``
``IMAGE_DISAGREEMENT``      visibility/image separation larger than
                            ``max(image_disagreement_beams x lambda/B_max, 3 sigma_position)``
``ZENITH_AZ_UNCERTAIN``     zenith distance below 3 sigma: azimuth undefined, ``sigma_az`` is NaN
``INVALID_GEOMETRY``        solution at/below the horizon (or the minimum elevation)
``CALIBRATION_WARNING``     Stage-2 calibration issue for the station (see the Stage-3 manifest)

An estimate is *usable* for Stage 4 when it carries none of the
:data:`BLOCKING_FLAGS`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

GOOD = "GOOD"
LOW_BASELINE_COUNT = "LOW_BASELINE_COUNT"
HIGH_FLAGGED_FRACTION = "HIGH_FLAGGED_FRACTION"
FIT_FAILED = "FIT_FAILED"
BOUNDARY_SOLUTION = "BOUNDARY_SOLUTION"
LOW_SNR = "LOW_SNR"
HIGH_RESIDUAL = "HIGH_RESIDUAL"
IMAGE_DISAGREEMENT = "IMAGE_DISAGREEMENT"
ZENITH_AZ_UNCERTAIN = "ZENITH_AZ_UNCERTAIN"
INVALID_GEOMETRY = "INVALID_GEOMETRY"
CALIBRATION_WARNING = "CALIBRATION_WARNING"

ALL_FLAGS = (LOW_BASELINE_COUNT, HIGH_FLAGGED_FRACTION, FIT_FAILED, BOUNDARY_SOLUTION, LOW_SNR,
             HIGH_RESIDUAL, IMAGE_DISAGREEMENT, ZENITH_AZ_UNCERTAIN, INVALID_GEOMETRY, CALIBRATION_WARNING)

#: Flags that make an estimate unusable as a measurement. The others are
#: warnings: the estimate is a measurement but should be down-weighted or
#: cross-checked in Stage 4. ``LOW_SNR`` blocks because below SNR ~ 5 the fit
#: can lock onto noise or sidelobe peaks (heavy-tailed errors in the
#: synthetic pull tests).
BLOCKING_FLAGS = frozenset({LOW_BASELINE_COUNT, FIT_FAILED, BOUNDARY_SOLUTION, INVALID_GEOMETRY, LOW_SNR})


@dataclass
class QualityThresholds:
    """Thresholds behind the quality flags (recorded in the Stage-3 manifest)."""

    min_baselines: int = 10
    max_flagged_fraction: float = 0.5
    min_snr: float = 5.0
    max_residual_tail_ratio: float = 1.5
    max_reduced_chi2: float = 3.0
    weights_calibrated: bool = False
    image_disagreement_beams: float = 0.5

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        return dict(self.__dict__)


def assess(
    *,
    thresholds: QualityThresholds,
    n_baselines_used: int,
    flagged_fraction: float,
    fitted: bool,
    fit_success: bool = False,
    boundary: bool = False,
    invalid_geometry: bool = False,
    snr: float = float("nan"),
    residual_tail_ratio: float = float("nan"),
    reduced_chi2: float = float("nan"),
    covariance_finite: bool = True,
    zenith_az_uncertain: bool = False,
    vis_image_separation_deg: Optional[float] = None,
    resolution_deg: Optional[float] = None,
    sigma_position_deg: Optional[float] = None,
    calibration_warning: bool = False,
) -> List[str]:
    """Return the sorted-by-definition list of quality flags for one estimate."""
    flags = []
    if n_baselines_used < thresholds.min_baselines:
        flags.append(LOW_BASELINE_COUNT)
    if flagged_fraction > thresholds.max_flagged_fraction:
        flags.append(HIGH_FLAGGED_FRACTION)
    if fitted:
        if not fit_success or not covariance_finite:
            flags.append(FIT_FAILED)
        if boundary:
            flags.append(BOUNDARY_SOLUTION)
        if not (snr >= thresholds.min_snr):
            flags.append(LOW_SNR)
        heavy = residual_tail_ratio > thresholds.max_residual_tail_ratio
        if heavy or (thresholds.weights_calibrated and reduced_chi2 > thresholds.max_reduced_chi2):
            flags.append(HIGH_RESIDUAL)
        if vis_image_separation_deg is not None and np.isfinite(vis_image_separation_deg):
            limit = thresholds.image_disagreement_beams * (resolution_deg or 0.0)
            if sigma_position_deg is not None and np.isfinite(sigma_position_deg):
                limit = max(limit, 3.0 * sigma_position_deg)
            if vis_image_separation_deg > limit:
                flags.append(IMAGE_DISAGREEMENT)
        if zenith_az_uncertain:
            flags.append(ZENITH_AZ_UNCERTAIN)
        if invalid_geometry:
            flags.append(INVALID_GEOMETRY)
    if calibration_warning:
        flags.append(CALIBRATION_WARNING)
    return [f for f in ALL_FLAGS if f in flags]


def is_usable(flags: List[str]) -> bool:
    """True if no blocking flag is present."""
    return not any(f in BLOCKING_FLAGS for f in flags)
