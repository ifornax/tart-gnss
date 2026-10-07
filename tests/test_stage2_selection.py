import json
from datetime import timedelta

import pytest

from stage2_helpers import stage1_doc, utc
from tart_gnss.stage2.errors import Stage1SelectionError
from tart_gnss.stage2.selection import load_stage1_results, query_interval, select_window, signal_relevance


@pytest.fixture
def results(tmp_path):
    path = tmp_path / "common_gnss_test.json"
    path.write_text(json.dumps(stage1_doc()))
    return path


def test_valid_rank(results):
    doc, path = load_stage1_results(results)
    w = select_window(doc, path, rank=1)
    assert (w.norad_id, w.rank, w.window_index) == (40534, 1, 1)
    assert w.start_utc == utc(2026, 6, 25, 10, 0)
    assert len(w.stations) == 3 and len(w.results_sha256) == 64


def test_invalid_rank(results):
    doc, path = load_stage1_results(results)
    with pytest.raises(Stage1SelectionError, match="Rank 99 not found"):
        select_window(doc, path, rank=99)


def test_valid_norad_single_window(results):
    doc, path = load_stage1_results(results)
    assert select_window(doc, path, norad_id=40534).rank == 1


def test_ambiguous_norad(results):
    doc, path = load_stage1_results(results)
    with pytest.raises(Stage1SelectionError, match="2 windows; add --window-index"):
        select_window(doc, path, norad_id=41549)


def test_norad_with_window_index(results):
    doc, path = load_stage1_results(results)
    w = select_window(doc, path, norad_id=41549, window_index=2)
    assert w.rank == 3 and w.window_index == 2 and w.n_windows_for_satellite == 2
    with pytest.raises(Stage1SelectionError, match="out of range"):
        select_window(doc, path, norad_id=41549, window_index=3)


def test_inconsistent_rank_and_norad(results):
    doc, path = load_stage1_results(results)
    with pytest.raises(Stage1SelectionError, match="is NORAD 40534"):
        select_window(doc, path, rank=1, norad_id=41549)


def test_no_selection(results):
    doc, path = load_stage1_results(results)
    with pytest.raises(Stage1SelectionError):
        select_window(doc, path)


def test_csv_resolves_to_json(results):
    csv = results.with_suffix(".csv")
    csv.write_text("rank\n1\n")
    _, path = load_stage1_results(csv)
    assert path == results


def test_missing_results_file(tmp_path):
    with pytest.raises(Stage1SelectionError, match="not found"):
        load_stage1_results(tmp_path / "missing.json")


def test_padding_arithmetic():
    a, b = utc(2026, 6, 25, 10, 0), utc(2026, 6, 25, 10, 20)
    q0, q1 = query_interval(a, b, 120)
    assert q0 == a - timedelta(seconds=120) and q1 == b + timedelta(seconds=120)
    assert query_interval(a, b, 0) == (a, b)
    with pytest.raises(ValueError):
        query_interval(a, b, -1)


def test_max_duration_centres_on_midpoint():
    a, b = utc(2026, 6, 25, 10, 0), utc(2026, 6, 25, 12, 0)
    q0, q1 = query_interval(a, b, 120, max_duration_s=600)
    assert (q1 - q0).total_seconds() == 600
    assert q0 + (q1 - q0) / 2 == utc(2026, 6, 25, 11, 0)


@pytest.mark.parametrize("const,name,expected", [
    ("GPS", "GPS BIIF-9  (PRN 26)", "GPS_L1"),
    ("Galileo", "GSAT0211 (GALILEO 14)", "GALILEO_E1"),
    ("BeiDou", "BEIDOU-3 M21 (C43)", "BEIDOU_B1C"),
    ("BeiDou", "BEIDOU-3 IGSO-1 (C38)", "BEIDOU_B1C"),
    ("BeiDou", "BEIDOU-2 G6 (C02)", "OUT_OF_BAND"),
    ("BeiDou", "BEIDOU-3 G2 (C60)", "UNKNOWN"),
    ("GLONASS", "COSMOS 2433 (720)", "UNKNOWN"),
    ("UNKNOWN", "X", "UNKNOWN"),
])
def test_signal_relevance(const, name, expected):
    assert signal_relevance(const, name)[0] == expected
