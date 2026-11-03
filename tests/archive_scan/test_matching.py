"""Stage-1 matching and ranking through the full scanner (mocked archive, synthetic Stage 1)."""

from datetime import timedelta

import pytest

from tart_gnss.archive_scan.matching import load_windows
from tart_gnss.stage2 import selection as sel

from .helpers import TRIAD, MockArchive, read_csv, run_of_files, scan, t, window, write_stage1

START, END = t("2026-10-09T00:00:00"), t("2026-10-11T00:00:00")
EPOCH = t("2026-10-10T00:00:00")


def _all(first_end, n):
    return {s: run_of_files(first_end, n) for s in TRIAD}


def _scan(tmp_path, files, windows, *argv, archive=None):
    s1 = write_stage1(tmp_path / "s1.json", windows, START, END)
    return scan(tmp_path, archive or MockArchive(files), "--visibility-results", str(s1), *argv)


def test_no_overlap_despite_geometric_visibility(tmp_path):
    files = {"mg-tana": run_of_files(t("2026-10-09T10:01:00"), 30),
             "za-rhodes": run_of_files(t("2026-10-09T12:01:00"), 30),
             "mu-udm": run_of_files(t("2026-10-09T14:01:00"), 30)}
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T09:00:00"), t("2026-10-09T16:00:00"), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    assert code == 0 and m["exit_classification"] == "NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS"
    assert read_csv(p["candidates_csv"]) == [] and m["summary"]["n_segments"] == 0


def test_all_station_common_overlap(tmp_path):
    files = _all(t("2026-10-09T12:01:00"), 20)
    files["mu-udm"] = run_of_files(t("2026-10-09T12:06:00"), 20)
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00"), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    rows = read_csv(p["candidates_csv"])
    assert code == 0 and m["exit_classification"] == "CANDIDATES_FOUND" and len(rows) == 1
    r = rows[0]
    assert r["overlap_start_utc"] == "2026-10-09T12:04:00.000000Z"
    assert r["overlap_end_utc"] == "2026-10-09T12:20:00.000000Z"
    assert float(r["possible_overlap_duration_s"]) == 960.0 and r["verified_overlap_duration_s"] == ""
    assert r["coverage_status"] == "FILENAME_CANDIDATE" and r["evidence_class"] == "ARCHIVE_COVERAGE_CANDIDATE"
    assert {"files_mg_tana", "files_za_rhodes", "files_mu_udm"} <= set(r)
    long_form = read_csv(p["candidate_stations_csv"])
    assert sorted(x["station_id"] for x in long_form) == sorted(TRIAD)


def test_two_station_overlap_is_not_three_station_valid(tmp_path):
    files = _all(t("2026-10-09T12:01:00"), 20)
    files["mu-udm"] = run_of_files(t("2026-10-09T18:01:00"), 20)
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00"), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    assert m["exit_classification"] == "NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS"
    code, m, p = _scan(tmp_path, files, w, "--min-stations", "2")
    rows = read_csv(p["candidates_csv"])
    assert rows and all(r["coverage_status"] == "TWO_STATION_DIAGNOSTIC" for r in rows)
    assert all("PARTIAL_NETWORK_DIAGNOSTIC" in r["flags"] for r in rows)
    assert m["exit_classification"] == "NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS"  # all-station claim only
    assert m["summary"]["n_candidates"] == 0 and m["summary"]["n_diagnostic_candidates"] == len(rows)


def test_separated_segments_stay_separate(tmp_path):
    files = {s: run_of_files(t("2026-10-09T12:01:00"), 10) + run_of_files(t("2026-10-09T13:01:00"), 10) for s in TRIAD}
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T15:00:00"), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    rows = read_csv(p["candidates_csv"])
    assert len(rows) == 2 and len({r["segment_id"] for r in rows}) == 2
    assert sum(float(r["possible_overlap_duration_s"]) for r in rows) == 2 * 660.0  # the 50-min gap is not covered
    assert m["summary"]["shared_possible_s"] == 1320.0


@pytest.mark.parametrize("seconds,kept", [(300.0, True), (299.9, False)])
def test_min_common_duration_boundary(tmp_path, seconds, kept):
    files = _all(t("2026-10-09T12:01:00"), 30)
    a = t("2026-10-09T12:10:00")
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", a, a + timedelta(seconds=seconds), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    assert (len(read_csv(p["candidates_csv"])) == 1) == kept
    if not kept:
        assert m["filters"]["n_intersections_below_min_duration"] == 1


def test_partial_overlap_keeps_subsecond_bounds(tmp_path):
    files = _all(t("2026-10-09T12:01:00"), 20)  # candidate 11:59:00 .. 12:20:00
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T12:10:00.3"), t("2026-10-09T13:00:00"), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    r = read_csv(p["candidates_csv"])[0]
    assert r["overlap_start_utc"] == "2026-10-09T12:10:00.300000Z"
    assert float(r["possible_overlap_duration_s"]) == pytest.approx(599.7, abs=1e-9)
    assert r["predicted_start_utc"] == "2026-10-09T12:10:00.300Z"


def test_tle_age_from_satellite_epoch_not_catalogue_time(tmp_path):
    files = _all(t("2026-10-09T12:01:00"), 20)
    epoch = t("2026-10-07T06:00:00")
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00"), epoch),
         window(2, 2, "GPS B (PRN 2)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00"), None)]
    code, m, p = _scan(tmp_path, files, w)
    rows = {r["rank_stage1"]: r for r in read_csv(p["candidates_csv"])}
    mid = t("2026-10-09T12:09:30")  # overlap 11:59:00 .. 12:20:00
    assert float(rows["1"]["tle_age_days"]) == pytest.approx((mid - epoch).total_seconds() / 86400, abs=1e-6)
    assert rows["2"]["tle_age_days"] == "" and "TLE_AGE_MISSING" in rows["2"]["flags"]


def test_signal_relevance_filters(tmp_path):
    files = _all(t("2026-10-09T12:01:00"), 20)
    a, b = t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00")
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", a, b, EPOCH),
         window(2, 2, "BEIDOU-2 G5 (C05)", "BeiDou", a, b, EPOCH),
         window(3, 3, "BEIDOU-3 G1 (C59)", "BeiDou", a, b, EPOCH),
         window(4, 4, "BEIDOU-3 M1 (C19)", "BeiDou", a, b, EPOCH)]
    code, m, p = _scan(tmp_path, files, w, "--signal-relevance", "GPS_L1", "BEIDOU_B1C")
    rows = {r["norad_id"]: r for r in read_csv(p["candidates_csv"])}
    assert set(rows) == {"1", "4"} and rows["4"]["signal_relevance"] == "BEIDOU_B1C"
    assert m["filters"]["window_exclusions"] == {"signal_relevance:OUT_OF_BAND": 1, "signal_relevance:UNKNOWN": 1}
    code, m, p = _scan(tmp_path, files, w)  # default with Stage-1 input: GPS_L1 only
    assert [r["norad_id"] for r in read_csv(p["candidates_csv"])] == ["1"]
    code, m, p = _scan(tmp_path, files, w, "--signal-relevance", "UNKNOWN")
    r = read_csv(p["candidates_csv"])[0]
    assert r["norad_id"] == "3" and r["signal_relevance"] == "UNKNOWN" and "SIGNAL_UNKNOWN" in r["flags"]


def test_rank_and_window_index_preserved_from_stage1(tmp_path):
    files = {s: run_of_files(t("2026-10-09T06:01:00"), 20) + run_of_files(t("2026-10-09T18:01:00"), 20) for s in TRIAD}
    w = [window(7, 9, "GPS A (PRN 1)", "GPS", t("2026-10-09T17:00:00"), t("2026-10-09T19:00:00"), EPOCH),
         window(3, 9, "GPS A (PRN 1)", "GPS", t("2026-10-09T05:00:00"), t("2026-10-09T07:00:00"), EPOCH)]
    s1 = write_stage1(tmp_path / "s1.json", w, START, END)
    code, m, p = scan(tmp_path, MockArchive(files), "--visibility-results", str(s1))
    rows = {r["rank_stage1"]: r for r in read_csv(p["candidates_csv"])}
    doc, path = sel.load_stage1_results(s1)
    for rank in (3, 7):
        ref = sel.select_window(doc, path, rank=rank)
        assert int(rows[str(rank)]["window_index"]) == ref.window_index
        assert int(rows[str(rank)]["n_windows_for_satellite"]) == 2
    assert {v.rank: v.window_index for v in load_windows(doc)} == {3: 1, 7: 2}


def test_long_geo_window_does_not_dominate_by_length(tmp_path):
    files = _all(t("2026-10-09T08:01:00"), 300)  # 08:00 .. 13:00
    w = [window(1, 1, "BEIDOU-3 IGSO-1 (C38)", "BeiDou", t("2026-10-09T00:00:00"), t("2026-10-10T00:00:00"),
                t("2026-10-01T00:00:00"), elevation=60.0),
         window(2, 2, "BEIDOU-3 G2 (C60)", "BeiDou", t("2026-10-09T00:00:00"), t("2026-10-10T00:00:00"), EPOCH),
         window(50, 3, "GPS A (PRN 1)", "GPS", t("2026-10-09T10:00:00"), t("2026-10-09T10:15:00"), EPOCH,
                elevation=25.0)]
    code, m, p = _scan(tmp_path, files, w, "--signal-relevance", "GPS_L1", "BEIDOU_B1C", "UNKNOWN")
    rows = read_csv(p["candidates_csv"])
    order = [r["norad_id"] for r in sorted(rows, key=lambda r: int(r["ranking_order"]))]
    # the 15-min GPS overlap beats the ~5-h IGSO overlap (both capped at 600 s; fresher TLE), and the
    # in-band IGSO beats the UNKNOWN GEO despite equal length.
    assert order == ["3", "1", "2"]
    assert m["ranking"]["duration_cap_s"] == 600.0 and len(m["ranking"]["rule"]) == 7


def test_no_double_counting_across_neighbouring_windows(tmp_path):
    files = _all(t("2026-10-09T12:01:00"), 30)  # one 31-min segment
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00"), EPOCH),
         window(2, 2, "GPS B (PRN 2)", "GPS", t("2026-10-09T11:30:00"), t("2026-10-09T14:00:00"), EPOCH)]
    code, m, p = _scan(tmp_path, files, w)
    assert len(read_csv(p["candidates_csv"])) == 2
    assert m["summary"]["shared_possible_s"] == 1860.0
    assert m["summary"]["in_band_candidate_union_s"] == 1860.0
    seg = read_csv(p["common_segments_csv"])
    assert len(seg) == 1 and seg[0]["n_matched_candidates"] == "2"
