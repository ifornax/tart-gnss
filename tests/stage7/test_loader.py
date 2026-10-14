"""Provenance-chain loading."""

import json
import shutil

import pytest

pytest.importorskip("pandas")

from tart_gnss.stage7.errors import InvalidStage6ManifestError, ProvenanceChainError  # noqa: E402
from tart_gnss.stage7.loader import calibration_at, load_diagnostic_chain, stage2_calibration_table  # noqa: E402


def test_chain_verified(runs):
    inp = load_diagnostic_chain(runs("clean").chain["stage6"])
    st = {lk.stage: lk.status for lk in inp.links}
    assert st == {6: "input", 5: "verified", 4: "verified", 3: "verified", 2: "verified"}
    assert inp.calibration_metadata_available and set(inp.stations) == {"madagascar", "mauritius", "south_africa"}
    assert len(inp.residuals) == len(inp.trajectory) == len(inp.positions) == 36
    assert inp.directions is not None and len(inp.contributions) == 108


def test_missing_stage2_recorded(runs):
    inp = load_diagnostic_chain(runs("nostage2").chain["stage6"])
    assert inp.link(2).status == "missing" and not inp.calibration_metadata_available


def test_hash_mismatch_strict_vs_permissive(runs, tmp_path):
    chain = runs("clean").chain
    work = tmp_path / "copy"
    shutil.copytree(chain["stage6"].parents[3], work)
    m6 = next((work / "stage6").glob("*/*/stage6_manifest.json"))
    doc = json.loads(m6.read_text())
    m5 = next((work / "stage5").glob("*/*/stage5_manifest.json"))
    doc["input_stage5_manifest"] = str(m5)
    m6.write_text(json.dumps(doc))
    m5.write_text(m5.read_text() + "\n")  # modified after Stage 6 recorded its hash
    with pytest.raises(ProvenanceChainError):
        load_diagnostic_chain(m6, "strict")
    inp = load_diagnostic_chain(m6, "permissive")
    assert inp.link(5).status == "mismatch"


def test_not_a_stage6_manifest(tmp_path):
    p = tmp_path / "m.json"
    p.write_text(json.dumps({"stage": 5, "schema_version": 1}))
    with pytest.raises(InvalidStage6ManifestError):
        load_diagnostic_chain(p)
    with pytest.raises(InvalidStage6ManifestError):
        load_diagnostic_chain(tmp_path / "none.json")


def test_stage2_calibration_lookup():
    doc = {"stations": {"a": {"observation": {"files": [
        {"first_timestamp_utc": "2026-06-25T09:00:00Z", "last_timestamp_utc": "2026-06-25T09:01:00Z",
         "zero_gain_antennas": [1, 2], "calibration_fingerprint": "x"},
        {"first_timestamp_utc": "2026-06-25T09:01:01Z", "last_timestamp_utc": "2026-06-25T09:02:00Z",
         "zero_gain_antennas": [], "calibration_fingerprint": "y"}]}}}}
    tab = stage2_calibration_table(doc)
    from datetime import datetime, timezone

    t = datetime(2026, 6, 25, 9, 0, 30, tzinfo=timezone.utc).timestamp()
    assert calibration_at(tab, "a", t)["zero_gain_count"] == 2
    assert calibration_at(tab, "a", t + 60)["solution_index"] == 1
    assert calibration_at(tab, "a", t + 3600) is None and calibration_at(tab, "b", t) is None
