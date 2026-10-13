"""Stage-4 input policy, covariance integrity and Stage-5 flag bookkeeping."""

import numpy as np
import pytest

from tart_gnss.stage5 import quality as q


@pytest.mark.parametrize("flag", sorted(q.STAGE4_FATAL_FLAGS))
def test_fatal_flags_rejected_under_both_policies(flag):
    for policy in q.POLICIES:
        assert q.stage4_rejection_reason([flag], True, policy) == flag


@pytest.mark.parametrize("flag", sorted(q.STAGE4_STRICT_FLAGS))
def test_strict_only_flags(flag):
    assert q.stage4_rejection_reason([flag], True, "strict") == flag
    assert q.stage4_rejection_reason([flag], True, "permissive") == ""
    assert q.has_advisory([flag], "permissive")


@pytest.mark.parametrize("flag", sorted(q.STAGE4_ADVISORY_FLAGS))
def test_advisory_flags_never_reject(flag):
    for policy in q.POLICIES:
        assert q.stage4_rejection_reason([flag], True, policy) == ""
        assert q.has_advisory([flag], policy)


def test_fit_failure_and_good():
    assert q.stage4_rejection_reason(["GOOD"], False, "permissive") == "FIT_FAILED"
    assert q.stage4_rejection_reason(["GOOD"], True, "strict") == ""
    assert not q.has_advisory(["GOOD"], "strict")
    with pytest.raises(ValueError):
        q.stage4_rejection_reason([], True, "lenient")


def test_validate_covariance():
    good = np.diag([1e10, 1e8, 1e8])
    good[0, 1] = good[1, 0] = 1e7
    r, why = q.validate_covariance(good)
    assert why == "" and np.array_equal(r, good)
    tiny = good.copy()
    tiny[0, 1] += 1e-3  # numerical asymmetry far below 1e-6 relative: symmetrised
    r, why = q.validate_covariance(tiny)
    assert why == "" and np.allclose(r, r.T)
    asym = good.copy()
    asym[0, 1] += 1e6
    assert "asymmetric" in q.validate_covariance(asym)[1]
    assert "positive definite" in q.validate_covariance(np.diag([1.0, -1.0, 1.0]))[1]
    assert "non-finite" in q.validate_covariance(np.full((3, 3), np.nan))[1]
    assert "condition" in q.validate_covariance(np.diag([1e14, 1.0, 1.0]))[1]
    assert "zero" in q.validate_covariance(np.zeros((3, 3)))[1]


def test_ordering_and_validity():
    assert q.ordered([]) == ["GOOD"]
    assert q.ordered(["HIGH_NIS", "INNOVATION_OUTLIER", "HIGH_NIS"]) == ["INNOVATION_OUTLIER", "HIGH_NIS"]
    assert q.is_valid(["HIGH_NIS", "INNOVATION_OUTLIER", "LARGE_TIME_GAP"])
    for f in q.BLOCKING_FLAGS:
        assert not q.is_valid([f])
