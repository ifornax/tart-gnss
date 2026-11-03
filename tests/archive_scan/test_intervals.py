"""Candidate segments: gaps kept, clipping, intersections, k-of-n sweep, file counts."""

from datetime import timedelta

from tart_gnss.archive_scan.intervals import (
    count_files,
    gap_groups,
    intersect_all,
    k_station_segments,
    merge_intervals,
    station_candidate_segments,
)
from tart_gnss.stage2.synchronise import total_seconds

from .helpers import run_of_files, t

W = (t("2026-10-09T00:00:00"), t("2026-10-10T00:00:00"))


def test_gaps_are_not_filled():
    files = run_of_files(t("2026-10-09T12:01:00"), 5) + run_of_files(t("2026-10-09T12:30:00"), 5)
    segs = station_candidate_segments(files, 120, W)
    assert segs == [(t("2026-10-09T11:59:00"), t("2026-10-09T12:05:00")),
                    (t("2026-10-09T12:28:00"), t("2026-10-09T12:34:00"))]
    assert total_seconds(segs) == 720.0


def test_merge_only_overlapping_or_touching():
    a = [(t("2026-10-09T01:00:00"), t("2026-10-09T01:02:00")), (t("2026-10-09T01:02:00"), t("2026-10-09T01:03:00")),
         (t("2026-10-09T01:03:00.000001"), t("2026-10-09T01:04:00"))]
    assert merge_intervals(a) == [(a[0][0], a[1][1]), a[2]]


def test_clipped_to_half_open_scan_interval():
    files = [t("2026-10-09T00:00:30"), t("2026-10-10T00:01:00")]  # crosses both ends
    segs = station_candidate_segments(files, 120, W)
    assert segs == [(W[0], t("2026-10-09T00:00:30")), (t("2026-10-09T23:59:00"), W[1])]


def test_zero_duration_boundary_intersection_discarded():
    a = [(t("2026-10-09T01:00:00"), t("2026-10-09T02:00:00"))]
    b = [(t("2026-10-09T02:00:00"), t("2026-10-09T03:00:00"))]
    assert intersect_all([a, b]) == []
    assert k_station_segments({"a": a, "b": b}, 2) == []


def test_k_of_n_sweep_matches_intersection_and_labels_station_sets():
    a = [(t("2026-10-09T01:00:00"), t("2026-10-09T03:00:00"))]
    b = [(t("2026-10-09T02:00:00"), t("2026-10-09T04:00:00"))]
    c = [(t("2026-10-09T02:30:00"), t("2026-10-09T02:40:00"))]
    full = k_station_segments({"a": a, "b": b, "c": c}, 3)
    assert [s for s, _ in full] == intersect_all([a, b, c]) and full[0][1] == ("a", "b", "c")
    two = k_station_segments({"a": a, "b": b, "c": c}, 2)
    assert [ids for _, ids in two] == [("a", "b"), ("a", "b", "c"), ("a", "b")]
    assert total_seconds([s for s, _ in two]) == 3600.0  # disjoint pieces, no double counting


def test_count_files_uses_stage2_envelope():
    times = run_of_files(t("2026-10-09T12:00:00"), 10)  # 12:00 .. 12:09
    assert count_files(times, t("2026-10-09T12:03:00"), t("2026-10-09T12:04:00"), 120) == 4  # 12:03..12:06


def test_gap_groups_report_covered_time_separately():
    segs = [(t("2026-10-09T01:00:00"), t("2026-10-09T01:10:00")), (t("2026-10-09T01:12:00"), t("2026-10-09T01:20:00"))]
    g = gap_groups(segs, 300)
    assert len(g) == 1 and g[0]["covered_s"] == 1080.0 and g[0]["span_s"] == 1200.0
    assert len(gap_groups(segs, 0)) == 2
