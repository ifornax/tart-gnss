"""``tart-archive-scan`` end to end: inputs, stations, handoff, honesty of reporting, determinism."""

import shlex
import subprocess

import pytest

from tart_gnss.archive_scan import cli
from tart_gnss.archive_scan.inventory import day_prefix
from tart_gnss.stage2 import cli as stage2_cli

from .helpers import (
    TRIAD,
    MockArchive,
    read_csv,
    registry_station,
    run_of_files,
    scan,
    t,
    window,
    write_stage1,
)

START, END = t("2026-10-09T00:00:00"), t("2026-10-11T00:00:00")
EPOCH = t("2026-10-10T00:00:00")


def _files():
    return {s: run_of_files(t("2026-10-09T12:01:00"), 20) for s in TRIAD}


def _s1(tmp_path, stations=None):
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T11:00:00"), t("2026-10-09T13:00:00"), EPOCH)]
    return write_stage1(tmp_path / "s1.json", w, START, END, stations)


def test_help_and_public_name(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["--help"])
    assert exc.value.code == 0
    text = capsys.readouterr().out
    assert text.startswith("usage: tart-archive-scan")
    for opt in ("--stations", "--visibility-results", "--min-common-duration", "--signal-relevance",
                "--cache-policy", "--file-span-assumption-s", "--min-stations", "--top"):
        assert opt in text


def test_aliases_resolve_to_frozen_stage1_geometry(tmp_path):
    stations = [registry_station(s) for s in TRIAD]
    stations[0]["latitude_deg"] += 1e-4  # frozen Stage-1 value differs from the registry
    s1 = _s1(tmp_path, stations)
    code, m, _ = scan(tmp_path, MockArchive(_files()), "--visibility-results", str(s1),
                      "--stations", "madagascar", "south_africa", "mauritius")
    assert code == 0
    mg = next(s for s in m["stations"] if s["id"] == "mg-tana")
    assert mg["latitude_deg"] == stations[0]["latitude_deg"] and mg["station_metadata_changed"] is True
    assert mg["archive_name"] == "mg-tana"


def test_station_set_mismatch_is_rejected(tmp_path, capsys):
    s1 = _s1(tmp_path)
    argv = ["--visibility-results", str(s1), "--stations", "mg-tana", "za-rhodes", "na-unam",
            "--output", str(tmp_path / "o.csv"), "--cache-dir", str(tmp_path / "c"), "--cache-policy", "offline"]
    assert cli.main(argv) == cli.EXIT_USAGE
    assert "do not match the frozen Stage-1 stations" in capsys.readouterr().err


@pytest.mark.parametrize("extra,msg", [
    (["--start", "2026-10-08T00:00:00Z"], "may only restrict"),
    (["--start", "2026-10-10T00:00:00Z", "--end", "2026-10-09T12:00:00Z"], "must be after"),
    (["--start", "2026-10-09T05:00:00"], "no timezone"),
    (["--min-stations", "1"], "between 2 and 3"),
    (["--file-span-assumption-s", "0"], "> 0"),
])
def test_invalid_inputs_exit_2(tmp_path, capsys, extra, msg):
    argv = ["--visibility-results", str(_s1(tmp_path)), "--output", str(tmp_path / "o.csv"),
            "--cache-dir", str(tmp_path / "c"), "--cache-policy", "offline", *extra]
    assert cli.main(argv) == cli.EXIT_USAGE
    assert msg in capsys.readouterr().err


def test_missing_stage1_and_inventory_only_requirements(tmp_path):
    base = ["--output", str(tmp_path / "o.csv"), "--cache-dir", str(tmp_path / "c"), "--cache-policy", "offline"]
    assert cli.main(["--visibility-results", str(tmp_path / "missing.json"), *base]) == cli.EXIT_USAGE
    assert cli.main(["--start", "2026-10-09T00:00:00Z", "--end", "2026-10-10T00:00:00Z", *base]) == cli.EXIT_USAGE
    assert cli.main(["--stations", *TRIAD, *base]) == cli.EXIT_USAGE


def test_inventory_only_scan(tmp_path):
    code, m, p = scan(tmp_path, MockArchive(_files()), "--stations", "madagascar", "south_africa", "mauritius",
                      "--start", "2026-10-09T00:00:00Z", "--end", "2026-10-10T00:00:00Z")
    assert code == 0 and m["summary"]["mode"] == "INVENTORY_ONLY"
    assert m["exit_classification"] == "COMMON_CANDIDATE_SEGMENTS_FOUND" and m["stage1"] is None
    assert read_csv(p["candidates_csv"]) == []
    segs = read_csv(p["common_segments_csv"])
    assert len(segs) == 1 and segs[0]["coverage_confidence"] == "FILENAME_CANDIDATE"
    days = read_csv(p["station_days_csv"])
    assert len(days) == 3 * 3 and {d["status"] for d in days} == {"FOUND", "EMPTY"}
    summary = {r["station_id"]: r for r in read_csv(p["station_summary_csv"])}
    assert summary["mg-tana"]["files"] == "20" and float(summary["ALL_STATIONS"]["possible_coverage_fraction"]) <= 1


def test_retrieval_command_is_a_valid_safe_stage2_dry_run(tmp_path):
    s1 = _s1(tmp_path)
    code, m, p = scan(tmp_path, MockArchive(_files()), "--visibility-results", str(s1))
    cmd = read_csv(p["candidates_csv"])[0]["retrieval_command"]
    words = shlex.split(cmd)
    assert words[0] == "tart-retrieve" and "--dry-run" in words and "--max-duration" not in words
    args = stage2_cli.build_parser().parse_args(words[1:])
    assert args.rank == 1 and args.dry_run and not args.download and args.visibility_results == str(s1)


def test_scanner_never_runs_stage2_or_downloads(tmp_path, monkeypatch):
    calls = []
    real = subprocess.run

    def spy(cmd, *a, **k):
        calls.append(cmd)
        return real(cmd, *a, **k)

    monkeypatch.setattr(subprocess, "run", spy)
    monkeypatch.setattr(stage2_cli, "run", lambda *a, **k: pytest.fail("tart-retrieve must not run"))
    monkeypatch.setattr("tart_gnss.stage2.download.download_observation",
                        lambda *a, **k: pytest.fail("no download allowed"))
    code, m, p = scan(tmp_path, MockArchive(_files()), "--visibility-results", str(_s1(tmp_path)))
    assert code == 0  # MockArchive.download would raise if called
    assert all(c[0] == "git" for c in calls)


def test_empty_complete_archive_reports_honestly(tmp_path, capsys):
    code, m, p = scan(tmp_path, MockArchive(), "--visibility-results", str(_s1(tmp_path)))
    out = capsys.readouterr().out
    assert code == 0 and m["exit_classification"] == "NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS"
    assert "Inventory completeness: COMPLETE" in out and "no common data" not in out.lower()
    assert {d["status"] for d in read_csv(p["station_days_csv"])} == {"EMPTY"}


def test_partial_failed_day_flags_incompleteness_and_exits_3(tmp_path, capsys):
    archive = MockArchive(_files(), errors={day_prefix("mu-udm", t("2026-10-09T00:00:00").date()): -1})
    code, m, p = scan(tmp_path, archive, "--visibility-results", str(_s1(tmp_path)))
    out = capsys.readouterr().out
    assert code == cli.EXIT_INCOMPLETE
    assert m["exit_classification"] == "NO_CANDIDATE_OVERLAP_BUT_INCOMPLETE_INVENTORY"
    assert "Inventory completeness: INCOMPLETE" in out and m["listing"]["status_counts"]["ERROR"] == 1
    errors = [d for d in read_csv(p["station_days_csv"]) if d["status"] == "ERROR"]
    assert len(errors) == 1 and errors[0]["n_files"] == "0" and "mock failure" in errors[0]["error"]


def test_incomplete_neighbouring_day_marks_candidate(tmp_path):
    files = {s: run_of_files(t("2026-10-09T23:51:00"), 5) for s in TRIAD}  # ends 23:55, envelope to 23:57
    files2 = {s: run_of_files(t("2026-10-09T00:01:00"), 10) for s in TRIAD}  # from the scan start
    for s in TRIAD:
        files[s] += files2[s]
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T00:00:00"), t("2026-10-10T00:00:00"), EPOCH)]
    s1 = write_stage1(tmp_path / "s1.json", w, START, END)
    archive = MockArchive(files, errors={day_prefix("za-rhodes", t("2026-10-08T00:00:00").date()): -1})
    code, m, p = scan(tmp_path, archive, "--visibility-results", str(s1))
    rows = {r["overlap_start_utc"][:16]: r for r in read_csv(p["candidates_csv"])}
    assert rows["2026-10-09T00:00"]["coverage_status"] == "INVENTORY_INCOMPLETE"  # 10-08 could hold files for it
    assert rows["2026-10-09T23:49"]["coverage_status"] == "FILENAME_CANDIDATE"
    assert code == cli.EXIT_INCOMPLETE and m["exit_classification"] == "CANDIDATES_FOUND_INVENTORY_INCOMPLETE"
    first = read_csv(p["candidates_csv"])[0]
    assert first["inventory_complete"] == "True"  # complete-inventory candidates rank first


def test_all_listings_failed_is_a_listing_failure(tmp_path):
    days = [t("2026-10-08T00:00:00"), t("2026-10-09T00:00:00"), t("2026-10-10T00:00:00"), t("2026-10-11T00:00:00")]
    errors = {day_prefix(s, d.date()): -1 for s in TRIAD for d in days}
    code, m, p = scan(tmp_path, MockArchive(_files(), errors=errors), "--visibility-results", str(_s1(tmp_path)))
    assert code == cli.EXIT_FAILURE and m["exit_classification"] == "ARCHIVE_LISTING_FAILURE"


def test_offline_without_cache_is_unknown_not_empty(tmp_path):
    code, m, p = scan(tmp_path, None, "--visibility-results", str(_s1(tmp_path)), "--cache-policy", "offline")
    assert code == cli.EXIT_INCOMPLETE and m["exit_classification"] == "ARCHIVE_LISTING_FAILURE"
    assert {d["status"] for d in read_csv(p["station_days_csv"])} == {"UNKNOWN"}


def test_output_independent_of_listing_order(tmp_path):
    files = {s: run_of_files(t("2026-10-09T12:01:00"), 15) + run_of_files(t("2026-10-10T03:01:00"), 15) for s in TRIAD}
    w = [window(1, 1, "GPS A (PRN 1)", "GPS", t("2026-10-09T00:00:00"), t("2026-10-11T00:00:00"), EPOCH),
         window(2, 2, "GPS B (PRN 2)", "GPS", t("2026-10-10T02:00:00"), t("2026-10-10T05:00:00"), EPOCH)]
    s1 = write_stage1(tmp_path / "s1.json", w, START, END)
    texts = []
    for seed in (1, 2):
        d = tmp_path / f"run{seed}"
        code, m, p = scan(d, MockArchive(files, shuffle_seed=seed), "--visibility-results", str(s1))
        days = [{k: v for k, v in r.items() if k != "cache_path"} for r in read_csv(p["station_days_csv"])]
        texts.append((p["candidates_csv"].read_text(), p["common_segments_csv"].read_text(), days))
    assert texts[0] == texts[1]


def test_previous_manifest_is_kept(tmp_path):
    s1 = _s1(tmp_path)
    scan(tmp_path, MockArchive(_files()), "--visibility-results", str(s1))
    scan(tmp_path, MockArchive(_files()), "--visibility-results", str(s1), "--cache-policy", "reuse",
         now=t("2026-11-10T13:00:00").replace(tzinfo=t("2026-01-01T00:00:00").tzinfo))
    assert len(list((tmp_path / "out").glob("archive_scan_manifest*.json"))) == 2
