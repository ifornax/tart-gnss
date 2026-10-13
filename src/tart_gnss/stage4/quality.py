"""Stage-3 input policy and Stage-4 quality flags.

Stage-3 flags (input policy, ``--quality-policy``)
--------------------------------------------------
Always fatal (the row carries no usable direction), under every policy:
``FIT_FAILED``, ``INVALID_GEOMETRY``, ``LOW_BASELINE_COUNT``.

``strict`` (default) also rejects the other Stage-3 *blocking* flags, i.e.
uses only rows that Stage 3 marked ``usable``: ``BOUNDARY_SOLUTION`` (search
cap edge or sidelobe), ``LOW_SNR`` (fit may sit on a noise peak) and
``TLE_UNASSOCIATED`` (the measured source is > 5 deg from the target's TLE
prediction, i.e. some other source). The latter is a post-fit *association*
decision taken in Stage 3; it selects which measurements describe the target
but never enters the Stage-4 optimisation.

``permissive`` keeps rows with those three flags and marks the epoch
``STAGE3_QUALITY_WARNING``.

Advisory under both policies (row kept, epoch marked
``INPUT_DIRECTION_WARNING``): ``HIGH_FLAGGED_FRACTION``, ``HIGH_RESIDUAL``,
``IMAGE_DISAGREEMENT``, ``ZENITH_AZ_UNCERTAIN``. A Stage-3
``CALIBRATION_WARNING`` becomes the Stage-4 ``CALIBRATION_WARNING``.

Stage-4 flags (per epoch; an epoch with none is ``GOOD``)
----------------------------------------------------------
``INSUFFICIENT_STATIONS``    fewer than ``min_stations`` matched (two-station
                             solutions are degraded diagnostics only)
``TIME_MISMATCH``            target motion during a station's time offset
                             exceeds ``time_mismatch_sigma`` x its angular sigma
``POOR_GEOMETRY``            condition number of ``sum P_i`` above
                             ``max_condition_number`` (or rank < 3)
``NEGATIVE_RANGE``           a station range ``rho_i <= 0`` (point behind a station)
``HIGH_LINE_RESIDUAL``       a station's normalised line residual above
                             ``max_normalised_residual`` (or, without
                             covariances, its line residual above ``max_line_residual_m``)
``HIGH_CHI2``                ``reduced_chi2 > max_reduced_chi2``
``INPUT_DIRECTION_WARNING``  a contributing Stage-3 row has an advisory flag
``CALIBRATION_WARNING``      a contributing station has a Stage-2/3 calibration warning
``STAGE3_QUALITY_WARNING``   a contributing row has a Stage-3 blocking flag (permissive policy)
``COVARIANCE_SINGULAR``      directional or position covariance singular / not positive definite
``BOOTSTRAP_UNSTABLE``       bootstrap failures, bootstrap spread inconsistent with the
                             linear covariance, or bootstrap median shifted from the point
                             estimate (non-linear regime)
``TRIANGULATION_FAILED``     no finite solution (or the weighted iteration did not converge)

An epoch is *valid* (a 3-D measurement for Stage 5) when it has none of the
:data:`BLOCKING_FLAGS`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List

GOOD = "GOOD"
INSUFFICIENT_STATIONS = "INSUFFICIENT_STATIONS"
TIME_MISMATCH = "TIME_MISMATCH"
POOR_GEOMETRY = "POOR_GEOMETRY"
NEGATIVE_RANGE = "NEGATIVE_RANGE"
HIGH_LINE_RESIDUAL = "HIGH_LINE_RESIDUAL"
HIGH_CHI2 = "HIGH_CHI2"
INPUT_DIRECTION_WARNING = "INPUT_DIRECTION_WARNING"
CALIBRATION_WARNING = "CALIBRATION_WARNING"
STAGE3_QUALITY_WARNING = "STAGE3_QUALITY_WARNING"
COVARIANCE_SINGULAR = "COVARIANCE_SINGULAR"
BOOTSTRAP_UNSTABLE = "BOOTSTRAP_UNSTABLE"
TRIANGULATION_FAILED = "TRIANGULATION_FAILED"

ALL_FLAGS = (INSUFFICIENT_STATIONS, TIME_MISMATCH, POOR_GEOMETRY, NEGATIVE_RANGE, HIGH_LINE_RESIDUAL, HIGH_CHI2,
             INPUT_DIRECTION_WARNING, CALIBRATION_WARNING, STAGE3_QUALITY_WARNING, COVARIANCE_SINGULAR,
             BOOTSTRAP_UNSTABLE, TRIANGULATION_FAILED)

#: Flags that make an epoch unusable as a 3-D measurement.
BLOCKING_FLAGS = frozenset({INSUFFICIENT_STATIONS, POOR_GEOMETRY, NEGATIVE_RANGE, COVARIANCE_SINGULAR,
                            TRIANGULATION_FAILED})

#: Stage-3 flags whose rows are never used.
STAGE3_FATAL_FLAGS = frozenset({"FIT_FAILED", "INVALID_GEOMETRY", "LOW_BASELINE_COUNT"})
#: Further Stage-3 blocking flags rejected by the strict policy only.
STAGE3_STRICT_FLAGS = frozenset({"BOUNDARY_SOLUTION", "LOW_SNR", "TLE_UNASSOCIATED"})
#: Stage-3 advisory flags (row kept, epoch marked INPUT_DIRECTION_WARNING).
STAGE3_ADVISORY_FLAGS = frozenset({"HIGH_FLAGGED_FRACTION", "HIGH_RESIDUAL", "IMAGE_DISAGREEMENT",
                                   "ZENITH_AZ_UNCERTAIN"})
STAGE3_CALIBRATION_FLAG = "CALIBRATION_WARNING"

POLICIES = ("strict", "permissive")


@dataclass
class Stage4Thresholds:
    """Thresholds behind the Stage-4 flags (recorded in the manifest)."""

    min_stations: int = 3
    max_condition_number: float = 1e8
    max_normalised_residual: float = 5.0
    max_line_residual_m: float = 1.0e6
    max_reduced_chi2: float = 5.0
    time_mismatch_sigma: float = 1.0
    bootstrap_max_failure_fraction: float = 0.1
    bootstrap_max_sigma_ratio: float = 2.0
    bootstrap_max_median_shift_sigma: float = 0.25

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        return dict(self.__dict__)


def stage3_rejection_reason(flags: Iterable[str], policy: str) -> str:
    """Why a Stage-3 row is rejected under ``policy`` ('' if it is kept)."""
    if policy not in POLICIES:
        raise ValueError(f"unknown quality policy {policy!r}")
    flags = set(flags)
    fatal = sorted(flags & STAGE3_FATAL_FLAGS)
    if fatal:
        return fatal[0]
    if policy == "strict":
        strict = sorted(flags & STAGE3_STRICT_FLAGS)
        if strict:
            return strict[0]
    return ""


def input_flags(stage3_flags: Iterable[str]) -> List[str]:
    """Stage-4 epoch flags implied by one kept Stage-3 row."""
    f = set(stage3_flags)
    out = []
    if f & STAGE3_ADVISORY_FLAGS:
        out.append(INPUT_DIRECTION_WARNING)
    if STAGE3_CALIBRATION_FLAG in f:
        out.append(CALIBRATION_WARNING)
    if f & STAGE3_STRICT_FLAGS:
        out.append(STAGE3_QUALITY_WARNING)
    return out


def ordered(flags: Iterable[str]) -> List[str]:
    """Unique flags in definition order."""
    f = set(flags)
    return [x for x in ALL_FLAGS if x in f]


def is_valid(flags: Iterable[str]) -> bool:
    """True if no blocking flag is present."""
    return not any(f in BLOCKING_FLAGS for f in flags)
