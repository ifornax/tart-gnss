"""Stage-5 input policy, Stage-6 flags and precision-validation status.

Stage-5 input policy (``--quality-policy``)
-------------------------------------------
``strict`` (default): only Stage-5 rows that are ``valid`` (none of
``NOT_INITIALISED``, ``SMOOTHER_FAILED``, ``COVARIANCE_INVALID``,
``INSUFFICIENT_MEASUREMENTS``) with a finite smoothed state.

``permissive``: also rows carrying those flags if they still have a finite
state; where the smoothed state is missing (e.g. ``SMOOTHER_FAILED``) the
filtered state and covariance are used (``stage5_state_source = filtered``).

Advisory Stage-5 flags (``FILTER_DIVERGENCE``, ``INNOVATION_OUTLIER``,
``LARGE_TIME_GAP``, ``EXTRAPOLATED``, ``HIGH_NIS``, ``DYNAMICS_MISMATCH``,
``STAGE4_QUALITY_WARNING``, ``INITIAL_VELOCITY_WEAK``) are carried into the
row as ``STAGE5_QUALITY_WARNING`` (the Stage-5 flags are kept verbatim in
``stage5_quality_flags``).

Stage-6 flags
-------------
Per row:
``INVALID_STAGE5_STATE``          row not compared: no finite Stage-5 state, or (strict) Stage-5
                                  blocking flags; ``stage5_rejection_reason`` says which
``INVALID_STAGE5_COVARIANCE``     covariance unusable: Mahalanobis distances NaN
``TLE_AGE_WARNING``               |age| > ``--tle-age-warning-days``
``TLE_STALE``                     |age| > ``--tle-age-stale-days``
``TLE_BACKPROPAGATED``            age < 0 (propagated before the TLE epoch)
``TLE_REFERENCE_DISPERSION_HIGH`` ensemble position dispersion exceeds the Stage-5 position sigma
``HIGH_POSITION_RESIDUAL``        |dr| above ``--position-residual-warning-km`` (only if given)
``HIGH_VELOCITY_RESIDUAL``        |dv| above ``--velocity-residual-warning-m-s`` (only if given)
``HIGH_CONDITIONAL_MAHALANOBIS``  position D2 (Stage-5 covariance only) above chi2_3(0.9973) = 14.16
``STAGE5_QUALITY_WARNING``        advisory Stage-5 flags present
Run level (on every compared row):
``REFERENCE_UNCERTAINTY_UNKNOWN`` no ensemble of >= 3 TLEs: reference uncertainty unknown
``TLE_ENSEMBLE_SMALL``            an ensemble was requested but has fewer than 3 members
``STRONG_RESIDUAL_TREND``         descriptive trend over the arc above 3 median sigmas
``RESIDUAL_AUTOCORRELATED``       lag-1 autocorrelation above the white-noise bound
``NOT_PRECISION_VALIDATION``      the run does not meet the precision-validation conditions
Errors (run aborted, named in the message): ``TLE_NOT_FOUND``,
``TLE_NORAD_MISMATCH`` (exit 2), ``TLE_FRAME_FAILURE`` (propagation/frame error, exit 1).

``valid`` is true for a compared row with a finite residual.
"""

from __future__ import annotations

from typing import Dict, Iterable, List

POLICIES = ("strict", "permissive")
STAGE5_BLOCKING = frozenset({"NOT_INITIALISED", "SMOOTHER_FAILED", "COVARIANCE_INVALID", "INSUFFICIENT_MEASUREMENTS"})
STAGE5_ADVISORY = frozenset({"FILTER_DIVERGENCE", "INNOVATION_OUTLIER", "LARGE_TIME_GAP", "EXTRAPOLATED", "HIGH_NIS",
                             "DYNAMICS_MISMATCH", "STAGE4_QUALITY_WARNING", "INITIAL_VELOCITY_WEAK"})
CONDITIONAL_D2_LIMIT = 14.156  # chi2_3 at 0.9973

GOOD = "GOOD"
INVALID_STAGE5_STATE = "INVALID_STAGE5_STATE"
INVALID_STAGE5_COVARIANCE = "INVALID_STAGE5_COVARIANCE"
TLE_AGE_WARNING = "TLE_AGE_WARNING"
TLE_STALE = "TLE_STALE"
TLE_BACKPROPAGATED = "TLE_BACKPROPAGATED"
TLE_ENSEMBLE_SMALL = "TLE_ENSEMBLE_SMALL"
TLE_REFERENCE_DISPERSION_HIGH = "TLE_REFERENCE_DISPERSION_HIGH"
REFERENCE_UNCERTAINTY_UNKNOWN = "REFERENCE_UNCERTAINTY_UNKNOWN"
HIGH_POSITION_RESIDUAL = "HIGH_POSITION_RESIDUAL"
HIGH_VELOCITY_RESIDUAL = "HIGH_VELOCITY_RESIDUAL"
HIGH_CONDITIONAL_MAHALANOBIS = "HIGH_CONDITIONAL_MAHALANOBIS"
STRONG_RESIDUAL_TREND = "STRONG_RESIDUAL_TREND"
RESIDUAL_AUTOCORRELATED = "RESIDUAL_AUTOCORRELATED"
STAGE5_QUALITY_WARNING = "STAGE5_QUALITY_WARNING"
NOT_PRECISION_VALIDATION = "NOT_PRECISION_VALIDATION"
ERROR_FLAGS = ("TLE_NOT_FOUND", "TLE_NORAD_MISMATCH", "TLE_FRAME_FAILURE")

FLAG_ORDER = (
    INVALID_STAGE5_STATE, INVALID_STAGE5_COVARIANCE, NOT_PRECISION_VALIDATION, TLE_STALE, TLE_AGE_WARNING,
    TLE_BACKPROPAGATED, REFERENCE_UNCERTAINTY_UNKNOWN, TLE_ENSEMBLE_SMALL, TLE_REFERENCE_DISPERSION_HIGH,
    HIGH_POSITION_RESIDUAL, HIGH_VELOCITY_RESIDUAL, HIGH_CONDITIONAL_MAHALANOBIS, STRONG_RESIDUAL_TREND,
    RESIDUAL_AUTOCORRELATED, STAGE5_QUALITY_WARNING,
)


def stage5_row_decision(flags: Iterable[str], valid: bool, has_smoothed: bool, has_filtered: bool,
                        policy: str):
    """``(use, state_source, reason)`` for one Stage-5 row under ``policy``."""
    if policy not in POLICIES:
        raise ValueError(f"unknown quality policy {policy!r}")
    flags = set(flags)
    if policy == "strict":
        if not valid or flags & STAGE5_BLOCKING:
            return False, "", "STAGE5_BLOCKING:" + ",".join(sorted(flags & STAGE5_BLOCKING) or ["not valid"])
        if not has_smoothed:
            return False, "", INVALID_STAGE5_STATE
        return True, "smoothed", ""
    if has_smoothed:
        return True, "smoothed", ""
    if has_filtered:
        return True, "filtered", ""
    return False, "", INVALID_STAGE5_STATE


def age_flags(age_days: float, warning_days: float, stale_days: float) -> List[str]:
    """TLE-age flags of one epoch."""
    out = []
    if abs(age_days) > stale_days:
        out.append(TLE_STALE)
    elif abs(age_days) > warning_days:
        out.append(TLE_AGE_WARNING)
    if age_days < 0:
        out.append(TLE_BACKPROPAGATED)
    return out


def ordered(flags: Iterable[str]) -> List[str]:
    """Unique flags in canonical order; ``[GOOD]`` if none."""
    s = set(flags) - {GOOD}
    out = [f for f in FLAG_ORDER if f in s] + sorted(s - set(FLAG_ORDER))
    return out or [GOOD]


def precision_validation(conditions: Dict[str, bool]) -> bool:
    """True only if every condition holds."""
    return all(bool(v) for v in conditions.values())
