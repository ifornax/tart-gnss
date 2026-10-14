"""Confidence rules."""

from tart_gnss.stage7.scoring import COMPONENTS, confidence


def test_levels():
    allc = {c: True for c in COMPONENTS}
    assert confidence(allc) == "HIGH"
    assert confidence({**allc, "persistence": False}) == "MODERATE"
    assert confidence({**allc, "reference_quality": False, "data_completeness": False}) == "LOW"
    assert confidence({**allc, "reference_quality": False}) == "MODERATE"
    assert confidence({}) == "LOW"
