"""Reference-TLE selection and provenance (spec section 5)."""

import json

import pytest

from tart_gnss.stage6 import synthetic as syn
from tart_gnss.stage6.errors import TLENoradMismatchError, TLEReferenceError
from tart_gnss.stage6.tle_reference import (
    download_tle_reference,
    load_ensemble,
    load_tle_reference,
    stage1_tle_reference,
)

from .helpers import L1, L2, TLE_EPOCH_S

OTHER = ("1 25544U 98067A   26279.50000000  .00010000  00000+0  18000-3 0  9990",
         "2 25544  51.6400 100.0000 0005000  90.0000 270.0000 15.50000000400000")


def _fix(l1, l2):
    from tart_gnss.tle import tle_checksum

    return l1[:68] + str(tle_checksum(l1)), l2[:68] + str(tle_checksum(l2))


def test_single_file(tmp_path):
    p = syn.write_tle_file(tmp_path / "a.tle", [(L1, L2)])
    ref = load_tle_reference(40534, p)
    assert ref.line1 == L1 and ref.sha256 and ref.source_kind == "tle-file"
    assert ref.epoch_unix_s == pytest.approx(TLE_EPOCH_S)


def test_history_file_selects_nearest_epoch(tmp_path):
    sets = [syn.perturbed_tle(L1, L2, d_epoch_days=d) for d in (-2.0, 0.0, 3.0)]
    p = syn.write_tle_file(tmp_path / "h.tle", sets)
    ref = load_tle_reference(40534, p, reference_time_s=TLE_EPOCH_S + 2.6 * 86400)
    assert ref.line1 == sets[2][0]
    assert "nearest" in ref.selection_rule and len(ref.provenance["candidates"]) == 3


def test_norad_mismatch_and_not_found(tmp_path):
    other = _fix(*OTHER)
    with pytest.raises(TLENoradMismatchError):
        load_tle_reference(40534, syn.write_tle_file(tmp_path / "o.tle", [other]))
    p = syn.write_tle_file(tmp_path / "two.tle", [other, _fix(*OTHER)])
    with pytest.raises(TLEReferenceError, match="TLE_NOT_FOUND"):
        load_tle_reference(40534, p)
    with pytest.raises(TLEReferenceError):
        load_tle_reference(40534, tmp_path / "missing.tle")


def test_modified_file_rejected(tmp_path):
    p = syn.write_tle_file(tmp_path / "a.tle", [(L1, L2)])
    (tmp_path / "a.tle.meta.json").write_text(json.dumps({"sha256": "0" * 64}))
    with pytest.raises(TLEReferenceError, match="modified"):
        load_tle_reference(40534, p)


def test_stage1_chain(tmp_path):
    m4 = tmp_path / "s4" / "stage4_manifest.json"
    m4.parent.mkdir()
    m4.write_text(json.dumps({"stage": 4, "satellite_reference": {"satellite": {
        "satellite_name": "GPS BIIF-9", "tle_line1": L1, "tle_line2": L2, "tle_source": "CelesTrak"}}}))
    m5 = tmp_path / "s5" / "stage5_manifest.json"
    doc = {"input_stage4_manifest": "s4/stage4_manifest.json", "provenance": {"working_directory": str(tmp_path)}}
    ref = stage1_tle_reference(m5, doc, 40534)
    assert ref.line1 == L1 and ref.source_kind == "stage1" and str(m4) in ref.provenance["chain"]
    with pytest.raises(TLENoradMismatchError):
        stage1_tle_reference(m5, doc, 99999)
    with pytest.raises(TLEReferenceError):
        stage1_tle_reference(m5, {"input_stage4_manifest": None}, 40534)


def test_download_is_explicit_and_archived(tmp_path):
    def fetch(url, timeout):
        return f"GPS BIIF-9\n{L1}\n{L2}\n"

    ref = download_tle_reference(40534, tmp_path / "tle", None, fetch=fetch)
    assert ref.source_kind == "download" and ref.provenance["downloaded"]
    assert list((tmp_path / "tle").glob("*.meta.json"))


def test_ensemble_window_and_dedup(tmp_path):
    sets = [syn.perturbed_tle(L1, L2, d_epoch_days=d) for d in (-1.0, 0.5, 2.0, 30.0)]
    syn.write_tle_file(tmp_path / "e" / "a.tle", sets[:3])
    syn.write_tle_file(tmp_path / "e" / "b.tle", [sets[0], sets[3], _fix(*OTHER)])
    members, info = load_ensemble(40534, dirs=[tmp_path / "e"], centre_time_s=TLE_EPOCH_S, window_days=7)
    assert len(members) == 3 and info["n_outside_window"] == 1 and info["n_sets_for_norad"] == 5
