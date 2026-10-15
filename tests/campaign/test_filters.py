"""Filters, policy and exclusion reasons."""

import pytest

from tart_gnss.campaign.filters import RunFilters, exclusion_reason, split_runs

ROW = {"run_id": "r", "manifest_sha256": "a", "n_valid_epochs": 10, "precision_validation": True, "norad_id": 1,
       "constellation": "GPS", "stations": ["x", "y"], "start_utc": "2026-06-25 09:00:00+00:00",
       "median_elevation_deg": 40.0, "max_abs_tle_age_days": 2.0, "overall_evidence_level": 3, "quality_flags": []}


@pytest.mark.parametrize("change,filters,reason", [
    ({"quality_flags": ["BROKEN_PROVENANCE"]}, RunFilters(), "BROKEN_PROVENANCE"),
    ({"n_valid_epochs": 2}, RunFilters(), "INSUFFICIENT_EPOCHS"),
    ({"precision_validation": False}, RunFilters(precision_validation_only=True), "NO_PRECISION_VALIDATION"),
    ({}, RunFilters(norad_ids=[2]), "FILTERED_NORAD_ID"),
    ({}, RunFilters(constellations=["galileo"]), "FILTERED_CONSTELLATION"),
    ({}, RunFilters(stations=["z"]), "FILTERED_STATION"),
    ({}, RunFilters(start="2026-07-01T00:00:00Z"), "FILTERED_TIME_RANGE"),
    ({}, RunFilters(min_elevation_deg=50), "FILTERED_ELEVATION"),
    ({}, RunFilters(max_tle_age_days=1), "FILTERED_TLE_AGE"),
    ({}, RunFilters(max_evidence_level=2), "FILTERED_EVIDENCE_LEVEL"),
    ({}, RunFilters(constellations=["GPS"], stations=["y"]), ""),
])
def test_reasons(change, filters, reason):
    assert exclusion_reason({**ROW, **change}, filters, "strict") == reason


def test_permissive_and_duplicates():
    broken = {**ROW, "quality_flags": ["BROKEN_PROVENANCE"]}
    assert exclusion_reason(broken, RunFilters(), "permissive") == ""
    inc, exc = split_runs([ROW, dict(ROW), {**ROW, "manifest_sha256": "b", "n_valid_epochs": 0}], RunFilters(), "strict")
    assert len(inc) == 1 and [e["reason"] for e in exc] == ["DUPLICATE_RUN", "INSUFFICIENT_EPOCHS"]
    with pytest.raises(ValueError):
        exclusion_reason(ROW, RunFilters(), "lenient")
