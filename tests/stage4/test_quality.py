"""Stage-3 quality policy and Stage-4 flags."""

import pytest

from tart_gnss.stage4 import quality as q


@pytest.mark.parametrize("flags, strict, permissive", [
    ([], "", ""),
    (["FIT_FAILED"], "FIT_FAILED", "FIT_FAILED"),
    (["INVALID_GEOMETRY", "LOW_SNR"], "INVALID_GEOMETRY", "INVALID_GEOMETRY"),
    (["LOW_BASELINE_COUNT"], "LOW_BASELINE_COUNT", "LOW_BASELINE_COUNT"),
    (["LOW_SNR"], "LOW_SNR", ""),
    (["TLE_UNASSOCIATED"], "TLE_UNASSOCIATED", ""),
    (["BOUNDARY_SOLUTION"], "BOUNDARY_SOLUTION", ""),
    (["HIGH_RESIDUAL", "IMAGE_DISAGREEMENT", "CALIBRATION_WARNING"], "", ""),
])
def test_rejection_policy(flags, strict, permissive):
    assert q.stage3_rejection_reason(flags, "strict") == strict
    assert q.stage3_rejection_reason(flags, "permissive") == permissive


def test_unknown_policy():
    with pytest.raises(ValueError):
        q.stage3_rejection_reason([], "lenient")


def test_input_flags():
    assert q.input_flags([]) == []
    assert q.input_flags(["IMAGE_DISAGREEMENT"]) == [q.INPUT_DIRECTION_WARNING]
    assert q.input_flags(["CALIBRATION_WARNING"]) == [q.CALIBRATION_WARNING]
    assert q.input_flags(["LOW_SNR"]) == [q.STAGE3_QUALITY_WARNING]


def test_validity_and_ordering():
    assert q.is_valid([])
    assert q.is_valid([q.HIGH_CHI2, q.TIME_MISMATCH, q.CALIBRATION_WARNING])
    for f in (q.INSUFFICIENT_STATIONS, q.POOR_GEOMETRY, q.NEGATIVE_RANGE, q.COVARIANCE_SINGULAR,
              q.TRIANGULATION_FAILED):
        assert not q.is_valid([f])
    assert q.ordered([q.TRIANGULATION_FAILED, q.HIGH_CHI2, q.HIGH_CHI2]) == [q.HIGH_CHI2, q.TRIANGULATION_FAILED]
    spec = {"INSUFFICIENT_STATIONS", "TIME_MISMATCH", "POOR_GEOMETRY", "NEGATIVE_RANGE", "HIGH_LINE_RESIDUAL",
            "HIGH_CHI2", "INPUT_DIRECTION_WARNING", "CALIBRATION_WARNING", "STAGE3_QUALITY_WARNING",
            "COVARIANCE_SINGULAR", "BOOTSTRAP_UNSTABLE", "TRIANGULATION_FAILED"}
    assert set(q.ALL_FLAGS) == spec


def test_policy_sets_cover_stage3_flags():
    from tart_gnss.stage3.quality import ALL_FLAGS as S3, BLOCKING_FLAGS as S3_BLOCK

    assert q.STAGE3_FATAL_FLAGS | q.STAGE3_STRICT_FLAGS == set(S3_BLOCK)
    covered = q.STAGE3_FATAL_FLAGS | q.STAGE3_STRICT_FLAGS | q.STAGE3_ADVISORY_FLAGS | {q.STAGE3_CALIBRATION_FLAG}
    assert covered == set(S3)
