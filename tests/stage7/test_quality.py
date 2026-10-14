"""Stage-7 input policy."""

import pytest

from tart_gnss.stage7 import quality as q


def test_epoch_policy():
    row = {"stage6_valid": True, "link_stage5": True, "link_stage4": True, "position_residual_3d_m": 5.0}
    assert q.epoch_used(row, "strict") and q.epoch_used(row, "permissive")
    assert not q.epoch_used({**row, "stage6_valid": False}, "strict")
    assert q.epoch_used({**row, "stage6_valid": False}, "permissive")
    assert not q.epoch_used({**row, "link_stage4": False}, "permissive")
    assert not q.epoch_used({**row, "stage6_valid": False, "position_residual_3d_m": float("nan")}, "permissive")
    with pytest.raises(ValueError):
        q.epoch_used(row, "lenient")
