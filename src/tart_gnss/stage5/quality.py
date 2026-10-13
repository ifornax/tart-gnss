"""Stage-4 input policy, covariance integrity checks and Stage-5 quality flags.

Stage-4 flags (input policy, ``--quality-policy``)
--------------------------------------------------
Always rejected (no usable position), under every policy, together with rows
whose ``fit_success`` is false or whose position/covariance is not finite:
``TRIANGULATION_FAILED``, ``INSUFFICIENT_STATIONS``, ``NEGATIVE_RANGE``,
``COVARIANCE_SINGULAR``.

``strict`` (default) also rejects ``POOR_GEOMETRY`` (a Stage-4 blocking flag:
strict uses exactly the Stage-4 ``valid`` epochs) and
``STAGE3_QUALITY_WARNING`` (the epoch rests on Stage-3 directions that Stage 3
itself marked unusable; on the 2026-06-25 data such epochs come from another
source). ``permissive`` accepts both.

Advisory under both policies (measurement used, row marked
``STAGE4_QUALITY_WARNING``): ``TIME_MISMATCH``, ``HIGH_LINE_RESIDUAL``,
``HIGH_CHI2``, ``INPUT_DIRECTION_WARNING``, ``CALIBRATION_WARNING``,
``BOOTSTRAP_UNSTABLE`` and, under ``permissive``, the two flags above.

Exactly repeated timestamps: the first row is used, later ones are rejected
(``DUPLICATE_TIMESTAMP``).

Covariance integrity
--------------------
Every Stage-4 covariance must be finite, symmetric to a relative tolerance of
``1e-6`` (smaller asymmetry is removed by ``(R + R^T)/2``), positive definite
and with condition number at most ``1e12``. A failing row is rejected and the
row flagged ``INVALID_STAGE4_COVARIANCE``; covariances are never repaired.

Stage-5 flags (per trajectory row; a row with none is ``GOOD``)
---------------------------------------------------------------
``INVALID_STAGE4_COVARIANCE``  the Stage-4 covariance failed the checks
``MEASUREMENT_REJECTED``       the Stage-4 row was not used (quality policy or duplicate)
``NO_MEASUREMENT``             the Stage-4 row has no finite position
``INITIAL_VELOCITY_WEAK``      initial velocity not from a >= 3-point fit, or poorly determined
``LARGE_TIME_GAP``             more than ``--max-gap-s`` since the last assimilated measurement
``INNOVATION_OUTLIER``         NIS above the gate (not assimilated unless ``--robust-update inflate``)
``HIGH_NIS``                   run-level: mean NIS too high for chi^2_3
``DYNAMICS_MISMATCH``          run-level: innovations significantly autocorrelated
``FILTER_DIVERGENCE``          one of >= 3 consecutive gated measurements after which the filter
                               was re-initialised from the data (a smoother segment boundary)
``EXTRAPOLATED``               after the last assimilated measurement (prediction only)
``STAGE4_QUALITY_WARNING``     the Stage-4 row carried advisory flags
``COVARIANCE_INVALID``         the state covariance failed its checks
``SMOOTHER_FAILED``            the RTS smoother failed (filtered state reported)
``NOT_INITIALISED``            before the first usable measurement (no state)
``INSUFFICIENT_MEASUREMENTS``  fewer than ``MIN_MEASUREMENTS`` usable measurements

Blocking (row not ``valid``): ``COVARIANCE_INVALID``, ``SMOOTHER_FAILED``,
``NOT_INITIALISED``, ``INSUFFICIENT_MEASUREMENTS``. ``FILTER_DIVERGENCE`` is
advisory: the re-initialised filter assimilates those rows again.
"""

from __future__ import annotations

from typing import Iterable, List, Tuple

import numpy as np

POLICIES = ("strict", "permissive")
MIN_MEASUREMENTS = 3
COVARIANCE_ASYMMETRY_TOL = 1e-6
MAX_COVARIANCE_CONDITION = 1e12

STAGE4_FATAL_FLAGS = frozenset({"TRIANGULATION_FAILED", "INSUFFICIENT_STATIONS", "NEGATIVE_RANGE",
                                "COVARIANCE_SINGULAR"})
STAGE4_STRICT_FLAGS = frozenset({"POOR_GEOMETRY", "STAGE3_QUALITY_WARNING"})
STAGE4_ADVISORY_FLAGS = frozenset({"TIME_MISMATCH", "HIGH_LINE_RESIDUAL", "HIGH_CHI2", "INPUT_DIRECTION_WARNING",
                                   "CALIBRATION_WARNING", "BOOTSTRAP_UNSTABLE"})

GOOD = "GOOD"
INVALID_STAGE4_COVARIANCE = "INVALID_STAGE4_COVARIANCE"
MEASUREMENT_REJECTED = "MEASUREMENT_REJECTED"
NO_MEASUREMENT = "NO_MEASUREMENT"
INITIAL_VELOCITY_WEAK = "INITIAL_VELOCITY_WEAK"
LARGE_TIME_GAP = "LARGE_TIME_GAP"
INNOVATION_OUTLIER = "INNOVATION_OUTLIER"
HIGH_NIS = "HIGH_NIS"
DYNAMICS_MISMATCH = "DYNAMICS_MISMATCH"
FILTER_DIVERGENCE = "FILTER_DIVERGENCE"
EXTRAPOLATED = "EXTRAPOLATED"
STAGE4_QUALITY_WARNING = "STAGE4_QUALITY_WARNING"
COVARIANCE_INVALID = "COVARIANCE_INVALID"
SMOOTHER_FAILED = "SMOOTHER_FAILED"
NOT_INITIALISED = "NOT_INITIALISED"
INSUFFICIENT_MEASUREMENTS = "INSUFFICIENT_MEASUREMENTS"

FLAG_ORDER = (
    INSUFFICIENT_MEASUREMENTS, NOT_INITIALISED, SMOOTHER_FAILED, COVARIANCE_INVALID, FILTER_DIVERGENCE,
    INVALID_STAGE4_COVARIANCE, MEASUREMENT_REJECTED, NO_MEASUREMENT, INNOVATION_OUTLIER, LARGE_TIME_GAP,
    EXTRAPOLATED, INITIAL_VELOCITY_WEAK, HIGH_NIS, DYNAMICS_MISMATCH, STAGE4_QUALITY_WARNING,
)
BLOCKING_FLAGS = frozenset({COVARIANCE_INVALID, SMOOTHER_FAILED, NOT_INITIALISED, INSUFFICIENT_MEASUREMENTS})


def stage4_rejection_reason(flags: Iterable[str], fit_success: bool, policy: str) -> str:
    """Why a Stage-4 row is not used under ``policy`` (empty string if it is used)."""
    if policy not in POLICIES:
        raise ValueError(f"unknown quality policy {policy!r}")
    flags = set(flags)
    for f in sorted(flags & STAGE4_FATAL_FLAGS):
        return f
    if not fit_success:
        return "FIT_FAILED"
    if policy == "strict":
        for f in sorted(flags & STAGE4_STRICT_FLAGS):
            return f
    return ""


def has_advisory(flags: Iterable[str], policy: str) -> bool:
    """The row carries flags that warrant ``STAGE4_QUALITY_WARNING``."""
    flags = set(flags)
    adv = STAGE4_ADVISORY_FLAGS | (STAGE4_STRICT_FLAGS if policy == "permissive" else frozenset())
    return bool(flags & adv)


def validate_covariance(r: np.ndarray, max_condition: float = MAX_COVARIANCE_CONDITION) -> Tuple[np.ndarray, str]:
    """Check one 3x3 covariance; return ``(symmetrised R, "")`` or ``(R, reason)``."""
    r = np.asarray(r, dtype=float)
    if r.shape != (3, 3) or not np.all(np.isfinite(r)):
        return r, "non-finite covariance"
    scale = float(np.max(np.abs(r)))
    if scale == 0:
        return r, "zero covariance"
    if float(np.max(np.abs(r - r.T))) > COVARIANCE_ASYMMETRY_TOL * scale:
        return r, "asymmetric covariance"
    r = 0.5 * (r + r.T)
    w = np.linalg.eigvalsh(r)
    if w[0] <= 0:
        return r, f"covariance not positive definite (smallest eigenvalue {w[0]:.3g})"
    if w[-1] / w[0] > max_condition:
        return r, f"covariance condition number {w[-1] / w[0]:.3g} above {max_condition:.3g}"
    return r, ""


def ordered(flags: Iterable[str]) -> List[str]:
    """Unique flags in canonical order; ``[GOOD]`` if none."""
    s = set(flags) - {GOOD}
    out = [f for f in FLAG_ORDER if f in s] + sorted(s - set(FLAG_ORDER))
    return out or [GOOD]


def is_valid(flags: Iterable[str]) -> bool:
    """No blocking flag."""
    return not (set(flags) & BLOCKING_FLAGS)
