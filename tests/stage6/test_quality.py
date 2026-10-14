"""Stage-5 policy, TLE-age flags, precision validation."""

import pytest

from tart_gnss.stage6 import quality as q


def test_strict_and_permissive():
    assert q.stage5_row_decision(["GOOD"], True, True, True, "strict") == (True, "smoothed", "")
    ok, _, why = q.stage5_row_decision(["SMOOTHER_FAILED"], False, False, True, "strict")
    assert not ok and "SMOOTHER_FAILED" in why
    assert q.stage5_row_decision(["SMOOTHER_FAILED"], False, False, True, "permissive") == (True, "filtered", "")
    assert q.stage5_row_decision(["NOT_INITIALISED"], False, False, False, "permissive")[0] is False
    assert q.stage5_row_decision(["HIGH_NIS"], True, True, True, "strict")[0]
    with pytest.raises(ValueError):
        q.stage5_row_decision([], True, True, True, "lenient")


@pytest.mark.parametrize("age,expected", [
    (0.5, []), (-0.5, ["TLE_BACKPROPAGATED"]), (5.0, ["TLE_AGE_WARNING"]), (20.0, ["TLE_STALE"]),
    (-103.0, ["TLE_STALE", "TLE_BACKPROPAGATED"]),
])
def test_age_flags(age, expected):  # 34.6
    assert q.age_flags(age, 3.0, 14.0) == expected


def test_precision_and_ordering():
    assert q.precision_validation({"a": True, "b": True})
    assert not q.precision_validation({"a": True, "b": False})
    assert q.ordered([]) == ["GOOD"]
    assert q.ordered(["STAGE5_QUALITY_WARNING", "TLE_STALE", "NOT_PRECISION_VALIDATION"])[:2] == [
        "NOT_PRECISION_VALIDATION", "TLE_STALE"]
