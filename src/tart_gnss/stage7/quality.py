"""Stage-7 input policy and data-quality flags.

``strict`` (default): an epoch is diagnosed only if its Stage-6 row is
``valid`` and it links to Stage-5 and Stage-4 rows; upstream-invalid rows are
listed but never diagnosed (no diagnosis is manufactured from an invalid
chain). ``permissive``: also epochs with a finite Stage-6 residual whose
Stage-6 row is not ``valid``, and provenance hash mismatches are tolerated
(recorded) instead of fatal. Advisory upstream flags stay visible in
``stage4_flags``/``stage5_flags``/``stage6_flags``.

Quality flags (never evidence tags; missing evidence is not evidence of absence):

``CALIBRATION_METADATA_MISSING`` no Stage-2 per-file calibration: lowers confidence, blocks Level 5
``STAGE3_LINK_MISSING``          no Stage-3 row for a contributing station (no re-triangulation for it)
``PROVENANCE_UNVERIFIED``        an upstream link without a recorded hash, or (permissive) a mismatch
``SUBSET_TEST_UNAVAILABLE``      fewer than 3 stations: no leave-one-out/station test at this epoch
``NOT_DIAGNOSED``                excluded by the policy
"""

from __future__ import annotations

POLICIES = ("strict", "permissive")
CALIBRATION_METADATA_MISSING = "CALIBRATION_METADATA_MISSING"
STAGE3_LINK_MISSING = "STAGE3_LINK_MISSING"
PROVENANCE_UNVERIFIED = "PROVENANCE_UNVERIFIED"
SUBSET_TEST_UNAVAILABLE = "SUBSET_TEST_UNAVAILABLE"
NOT_DIAGNOSED = "NOT_DIAGNOSED"
MIN_DIAGNOSED_EPOCHS = 3


def epoch_used(row, policy: str) -> bool:
    """Whether an epoch-table row is diagnosed under ``policy``."""
    if policy not in POLICIES:
        raise ValueError(f"unknown quality policy {policy!r}")
    linked = bool(row.get("link_stage5")) and bool(row.get("link_stage4"))
    if policy == "strict":
        return bool(row.get("stage6_valid")) and linked
    import numpy as np

    return linked and np.isfinite(float(row.get("position_residual_3d_m", np.nan)))
