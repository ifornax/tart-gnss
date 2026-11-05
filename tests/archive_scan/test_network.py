"""``tart-archive-scan --all-stations``: eligibility, physical sites, baselines, sweep, triads, outputs.

Synthetic five-station demo (day 2026-10-09, 120 s file span, one-minute files):

=========  ==========  ======================================  =========================
station    position    candidate coverage (UTC)                 note
=========  ==========  ======================================  =========================
aa         (0, 0)      10:00-12:00
bb         (0, 10)     10:30-11:30
cc         (10, 0)     11:00-11:45
cc-2       (10, .001)  10:00-12:00                              co-located with cc
dd         (0, 2)      10:40-10:44, 11:10-11:50                 outage; 222 km from aa
ee         (5, 5)      -                                        no archive_name
ff         lat 95      -                                        invalid coordinates
=========  ==========  ======================================  =========================

Expected (hand-computed): network intervals with >= 3 sites cover 10:30-11:50
(4800 s, 7 constant-membership pieces); triad aa+bb+cc is ONE 1800 s segment
11:00-11:30 although dd joins at 11:10; aa+bb+cc-2 is 3600 s; every triad
containing dd loses the 240 s 10:40-10:44 piece to the 300 s minimum and is
removed by a 500 km baseline cut.
"""

import json
from datetime import timedelta
from itertools import combinations

import pytest
import yaml

from tart_gnss.archive_scan import cli
from tart_gnss.archive_scan import network as net
from tart_gnss.archive_scan.inventory import day_prefix
from tart_gnss.archive_scan.network_run import network_output_paths
from tart_gnss.stations import default_station_file, load_stations

from .helpers import NOW, MockArchive, read_csv, run_of_files, t

DAY = "2026-10-09"
START, END = f"{DAY}T00:00:00Z", "2026-10-10T00:00:00Z"
POSITIONS = {"aa": (0.0, 0.0), "bb": (0.0, 10.0), "cc": (10.0, 0.0), "cc-2": (10.0, 0.001), "dd": (0.0, 2.0)}


def _span(first_end: str, last_end: str):
    """One-minute files whose 120 s candidate envelopes cover [first_end - 120 s, last_end]."""
    a, b = t(f"{DAY}T{first_end}"), t(f"{DAY}T{last_end}")
    return run_of_files(a, int((b - a).total_seconds() // 60) + 1)


def demo_files():
    return {
        "aa": _span("10:02:00", "12:00:00"),
        "bb": _span("10:32:00", "11:30:00"),
        "cc": _span("11:02:00", "11:45:00"),
        "cc-2": _span("10:02:00", "12:00:00"),
        "dd": _span("10:42:00", "10:44:00") + _span("11:12:00", "11:50:00"),
    }


def write_registry(path, positions=POSITIONS, extra=True):
    stations = {}
    for sid, (lat, lon) in positions.items():
        stations[sid] = {"name": sid.upper(), "tart_id": sid, "latitude_deg": lat, "longitude_deg": lon,
                         "altitude_m": 10.0, "archive_name": sid, "archive_verified": True, "source": "test"}
    if extra:
        stations["ee"] = {"name": "EE", "tart_id": "ee", "latitude_deg": 5.0, "longitude_deg": 5.0,
                          "altitude_m": 0.0, "archive_name": None, "source": "test"}
        stations["ff"] = {"name": "FF", "tart_id": "ff", "latitude_deg": 95.0, "longitude_deg": 0.0,
                          "altitude_m": 0.0, "archive_name": "ff", "source": "test"}
    path.write_text(yaml.safe_dump({"metadata": {"source": "test"}, "stations": stations}), encoding="utf-8")
    return path


def nscan(tmp_path, archive, *argv, registry=None, now=NOW, cache_dir=None, out_name="network_overlap.csv"):
    """Run a network scan in-process; returns (exit code, manifest, output paths)."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    reg = registry or write_registry(tmp_path / "registry.yaml")
    out = tmp_path / "out" / out_name
    args = ["--all-stations", "--station-file", str(reg), "--start", START, "--end", END,
            "--output", str(out), "--cache-dir", str(cache_dir or tmp_path / "cache"), *argv]
    ns = cli.build_parser().parse_args(args)
    code = cli.run(ns, ["tart-archive-scan", *args], archive=archive, now=lambda: now, sleep=lambda s: None)
    paths = network_output_paths(out)
    manifest = json.loads(paths["manifest"].read_text()) if paths["manifest"].exists() else None
    return code, manifest, paths


def triads_by_id(paths):
    return {r["triad_id"]: r for r in read_csv(paths["triads_csv"])}


# --------------------------------------------------------------------------- geodesy

def test_geodesic_reference_values():
    dms = lambda d, m, s: d + m / 60 + s / 3600  # noqa: E731
    # Vincenty (1975) / Geoscience Australia test line Flinders Peak -> Buninyong: 54 972.271 m.
    d = net.geodesic_distance_km(-dms(37, 57, 3.72030), dms(144, 25, 29.52440),
                                 -dms(37, 39, 10.15610), dms(143, 55, 35.38390))
    assert d * 1000 == pytest.approx(54972.271, abs=0.001)
    assert net.geodesic_distance_km(0, 0, 0, 1) * 1000 == pytest.approx(111319.4908, abs=1e-3)  # a * pi / 180
    assert net.geodesic_distance_km(0, 0, 90, 0) * 1000 == pytest.approx(10001965.729, abs=1e-2)  # WGS84 quadrant
    assert net.geodesic_distance_km(1, 2, 1, 2) == 0.0


def test_every_registry_pair_converges_and_chord_is_shorter():
    stations = load_stations()
    b = net.Baselines(stations)
    for x, y in combinations(sorted(stations), 2):
        g, c = b.geodesic_km(x, y), b.chord_km(x, y)
        dh_km = abs(stations[x].altitude_m - stations[y].altitude_m) / 1000.0
        assert g >= 0 and c <= g + dh_km + 1e-3  # the chord includes station heights


def test_nearly_antipodal_points_raise():
    with pytest.raises(net.GeodesicError):
        net.geodesic_distance_km(0, 0, 0.5, 179.7)


# --------------------------------------------------------------------------- registry and sites

def test_packaged_registry_eligibility_and_sites():
    catalogue, invalid = net.load_registry_tolerant(default_station_file())
    eligible, excluded = net.select_eligible(catalogue)
    reasons = {e.station_id: e.reason for e in invalid + excluded}
    assert {sid for sid, r in reasons.items() if r == net.NO_ARCHIVE_NAME} == {s.id for s in catalogue.values()
                                                                               if not s.archive_name}
    assert reasons.get("bw-biust") == net.INVALID_COORDINATES  # skipped by tart-stations-update
    ids = {s.id for s in eligible}
    assert all(catalogue[i].archive_name for i in ids)
    site = net.site_groups(eligible, net.Baselines({s.id: s for s in eligible}), 1.0)
    assert site["nz-nelson"] == site["nz-nelson-32"]  # identical coordinates: one site
    assert site["za-rhodes"] != site["za-hammanskraal"]  # two distinct South African stations
    assert site["mg-tana"] != site["mu-udm"]


def test_colocated_stations_are_one_site_and_tolerance_matters():
    from tart_gnss.stations import Station

    mk = lambda sid, lat, lon: Station(sid, sid, lat, lon, 0.0, "test", archive_name=sid)  # noqa: E731
    sts = [mk("a", 0, 0), mk("a2", 0, 0), mk("a3", 0, 0.005), mk("b", 0, 1)]  # a3 ~0.56 km from a
    b = net.Baselines({s.id: s for s in sts})
    assert set(net.site_groups(sts, b, 1.0).values()) == {"SITE_a", "SITE_b"}
    assert len(set(net.site_groups(sts, b, 0.0).values())) == 3  # tolerance 0: identical coordinates only


def test_duplicate_archive_name_listed_once():
    from tart_gnss.stations import Station

    cat = {s: Station(s, s, i, 0, 0, "test", archive_name="shared") for i, s in enumerate(("p", "q"))}
    eligible, excluded = net.select_eligible(cat)
    assert [s.id for s in eligible] == ["p"] and excluded[0].reason == net.DUPLICATE_ARCHIVE_NAME


# --------------------------------------------------------------------------- sweep

def test_half_open_touch_gives_no_overlap():
    segs = {"x": [(t(f"{DAY}T10:00:00"), t(f"{DAY}T10:30:00"))],
            "y": [(t(f"{DAY}T10:30:00"), t(f"{DAY}T11:00:00"))],
            "z": [(t(f"{DAY}T09:00:00"), t(f"{DAY}T12:00:00"))]}
    site = {s: f"SITE_{s}" for s in segs}
    assert net.network_active_segments(segs, site, 3) == []


def test_demo_end_to_end_matches_ground_truth(tmp_path):
    arch = MockArchive(demo_files())
    code, m, paths = nscan(tmp_path, arch)
    assert code == 0 and m["exit_classification"] == "TRIADS_FOUND"
    s = m["summary"]
    assert s["n_eligible_stations"] == 5 and s["n_distinct_sites"] == 4
    assert s["active_possible_s"] == 4800.0 and s["n_active_segments"] == 7
    active = read_csv(paths["active_segments_csv"])
    assert [(r["start_utc"][11:19], r["end_utc"][11:19], r["n_distinct_sites"]) for r in active] == [
        ("10:30:00", "10:40:00", "3"), ("10:40:00", "10:44:00", "4"), ("10:44:00", "11:00:00", "3"),
        ("11:00:00", "11:10:00", "3"), ("11:10:00", "11:30:00", "4"), ("11:30:00", "11:45:00", "3"),
        ("11:45:00", "11:50:00", "3")]
    assert json.loads(active[3]["active_station_ids_json"]) == ["aa", "bb", "cc", "cc-2"]

    tri = triads_by_id(paths)
    expected = {"aa+bb+cc-2": 3600, "aa+cc-2+dd": 2400, "aa+cc+dd": 2100, "aa+bb+cc": 1800,
                "aa+bb+dd": 1200, "bb+cc+dd": 1200, "bb+cc-2+dd": 1200}
    assert {k: float(v["possible_shared_duration_s"]) for k, v in tri.items()} == expected
    assert not any("cc+cc-2" in k for k in tri)  # co-located stations never form a triad
    # dd joins at 11:10: aa+bb+cc stays ONE contiguous segment.
    abc = tri["aa+bb+cc"]
    assert abc["n_contiguous_segments"] == "1" and abc["longest_contiguous_duration_s"] == "1800.0"
    segs = [r for r in read_csv(paths["triad_segments_csv"]) if r["triad_id"] == "aa+bb+cc"]
    assert len(segs) == 1 and segs[0]["co_active_additional_stations"] == "cc-2;dd"
    assert segs[0]["verified_duration_s"] == "" and segs[0]["coverage_confidence"] == "FILENAME_CANDIDATE"
    # The 240 s piece is below the 300 s minimum: excluded and counted, never bridged.
    assert tri["aa+bb+dd"]["n_short_segments_excluded"] == "1"
    assert tri["aa+bb+dd"]["first_shared_start_utc"].startswith(f"{DAY}T11:10:00")
    # Ranking: total duration, then longest, then min baseline (bb+cc+dd ~890 km beats aa+bb+dd ~222 km).
    order = [r["triad_id"] for r in read_csv(paths["triads_csv"])]
    assert order[:4] == ["aa+bb+cc-2", "aa+cc-2+dd", "aa+cc+dd", "aa+bb+cc"]
    assert order.index("bb+cc+dd") < order.index("aa+bb+dd")
    # Triads sharing a time interval do not double count the network total.
    assert s["triad_union_possible_s"] <= s["active_possible_s"]
    assert tri["aa+bb+cc"]["coverage_class"] == "FILENAME_CANDIDATE"
    assert tri["aa+bb+cc"]["verified_shared_duration_s"] == ""


def test_baseline_filter_requires_all_pairs(tmp_path):
    code, m, paths = nscan(tmp_path, MockArchive(demo_files()), "--min-baseline-km", "500")
    assert code == 0
    tri = triads_by_id(paths)
    # aa-dd = 222 km removes every triad with both aa and dd; bb+cc+dd (bb-dd ~890 km) passes.
    assert set(tri) == {"aa+bb+cc-2", "aa+bb+cc", "bb+cc+dd", "bb+cc-2+dd"}
    assert m["summary"]["n_triads_before_baseline"] == 7 and m["summary"]["n_triads_rejected_by_baseline"] == 3
    assert all(float(r["min_pairwise_baseline_km"]) >= 500 for r in tri.values())
    doc = json.loads(paths["triads_json"].read_text())
    assert {r["triad_id"] for r in doc["baseline_rejected"]} == {"aa+bb+dd", "aa+cc+dd", "aa+cc-2+dd"}
    # A 1000 km cut also removes the triads whose single short pair is bb-dd.
    code, m, paths = nscan(tmp_path, MockArchive(demo_files()), "--min-baseline-km", "1000",
                           cache_dir=tmp_path / "c2")
    assert set(triads_by_id(paths)) == {"aa+bb+cc-2", "aa+bb+cc"}


def test_baseline_filter_excludes_all_keeps_active_segments(tmp_path):
    code, m, paths = nscan(tmp_path, MockArchive(demo_files()), "--min-baseline-km", "5000")
    assert code == 0 and m["exit_classification"] == "THREE_SITE_OVERLAP_BUT_BASELINE_FILTER_EXCLUDES_ALL"
    assert read_csv(paths["triads_csv"]) == [] and len(read_csv(paths["active_segments_csv"])) == 7


def test_only_short_overlap_is_no_overlap_with_complete_inventory(tmp_path):
    files = {"aa": _span("10:02:00", "10:06:00"), "bb": _span("10:02:00", "10:06:00"),
             "cc": _span("10:02:00", "10:06:00")}  # 10:00-10:06 = 360 s
    code, m, paths = nscan(tmp_path, MockArchive(files), "--min-common-duration", "400")
    assert code == 0 and m["exit_classification"] == "NO_THREE_SITE_CANDIDATE_OVERLAP_WITH_COMPLETE_INVENTORY"
    assert m["summary"]["n_active_segments"] == 1 and m["summary"]["n_triads"] == 0


def test_listing_calls_once_per_station_day_not_per_triad(tmp_path):
    arch = MockArchive(demo_files())
    code, m, _ = nscan(tmp_path, arch)
    days = m["interval"]["days_listed"]
    assert days == ["2026-10-08", "2026-10-09", "2026-10-10"]
    assert len(arch.calls) == len(set(arch.calls)) == 5 * len(days)
    assert not any(c.startswith(("ee/", "ff/")) for c in arch.calls)
    # Cache reuse: no live calls on the second run; offline without a cache: UNKNOWN, exit 3.
    arch2 = MockArchive(demo_files())
    code, m, _ = nscan(tmp_path, arch2, "--cache-policy", "reuse")
    assert code == 0 and arch2.calls == [] and m["summary"]["cache_hits"] == 15
    code, m, _ = nscan(tmp_path, MockArchive(demo_files()), "--cache-policy", "offline",
                       cache_dir=tmp_path / "empty")
    assert code == 3 and m["exit_classification"] == "ARCHIVE_LISTING_FAILURE"  # nothing known at all


def test_no_overlap_with_incomplete_inventory_is_undetermined(tmp_path):
    files = {"aa": _span("10:02:00", "11:00:00"), "bb": _span("10:02:00", "11:00:00")}
    bad = day_prefix("cc", t(f"{DAY}T00:00:00").date())
    code, m, _ = nscan(tmp_path, MockArchive(files, errors={bad: -1}))
    assert code == 3 and m["exit_classification"] == "OVERLAP_UNDETERMINED_DUE_TO_INCOMPLETE_INVENTORY"
    code, m, _ = nscan(tmp_path, MockArchive(files), cache_dir=tmp_path / "c2")
    assert code == 0 and m["exit_classification"] == "NO_THREE_SITE_CANDIDATE_OVERLAP_WITH_COMPLETE_INVENTORY"


def test_listing_error_is_not_an_empty_day(tmp_path):
    bad = day_prefix("bb", t(f"{DAY}T00:00:00").date())
    code, m, paths = nscan(tmp_path, MockArchive(demo_files(), errors={bad: -1}))
    assert code == 3 and m["exit_classification"] == "TRIADS_FOUND_INVENTORY_INCOMPLETE"
    day = next(r for r in read_csv(paths["station_days_csv"]) if r["prefix"] == bad)
    assert day["status"] == "ERROR" and day["n_files"] == "0"
    summary = {r["station_id"]: r for r in read_csv(paths["station_summary_csv"])}
    assert summary["bb"]["days_error"] == "1" and summary["bb"]["inventory_complete"] == "False"
    tri = triads_by_id(paths)
    assert all(not k.startswith("bb") and "+bb+" not in k for k in tri)  # bb has no files left at all
    assert tri["aa+cc+dd"]["inventory_complete"] == "True"  # bb's error does not affect sets without bb
    active = read_csv(paths["active_segments_csv"])
    assert all(r["inventory_complete"] == "False" for r in active)  # a missing station could add coverage


def test_no_eligible_stations(tmp_path):
    code, m, paths = nscan(tmp_path, MockArchive({}), "--include-stations", "ee")
    assert code == 0 and m["exit_classification"] == "NO_ELIGIBLE_ARCHIVE_STATIONS"
    assert read_csv(paths["triads_csv"]) == []
    reasons = {r["station_id"]: r["reason"] for r in read_csv(paths["excluded_stations_csv"])}
    assert reasons["ee"] == "NO_ARCHIVE_NAME" and reasons["ff"] == "INVALID_COORDINATES"
    assert reasons["aa"] == "NOT_IN_INCLUDE_LIST"


def test_verified_archive_name_does_not_imply_files(tmp_path):
    files = demo_files()
    del files["dd"]
    code, m, paths = nscan(tmp_path, MockArchive(files))
    dd = next(r for r in read_csv(paths["station_summary_csv"]) if r["station_id"] == "dd")
    assert dd["archive_verified"] == "True" and dd["files"] == "0" and dd["days_empty"] == "3"
    assert dd["inventory_complete"] == "True"


def test_filename_time_not_upload_time(tmp_path):
    code, m, paths = nscan(tmp_path, MockArchive(demo_files()))
    aa = next(r for r in read_csv(paths["station_summary_csv"]) if r["station_id"] == "aa")
    assert aa["earliest_filename_utc"].startswith(f"{DAY}T10:02:00")  # mock last_modified is 3 min later


def test_station_summary_counts_only_files_reaching_the_interval(tmp_path):
    files = demo_files()
    files["aa"] = files["aa"] + [t("2026-10-08T23:59:00"), t("2026-10-10T00:01:00"), t("2026-10-10T00:03:00")]
    code, m, paths = nscan(tmp_path, MockArchive(files))
    aa = next(r for r in read_csv(paths["station_summary_csv"]) if r["station_id"] == "aa")
    # 23:59 (ends before start) and 00:03 (starts after end) are outside; 00:01 reaches back into the interval.
    assert aa["files_listed"] == "122" and aa["files"] == "120"
    assert aa["earliest_filename_utc"].startswith(f"{DAY}T10:02:00")
    assert aa["latest_filename_utc"].startswith("2026-10-10T00:01:00")


def test_colocated_only_network_is_rejected(tmp_path, capsys):
    reg = write_registry(tmp_path / "co.yaml", {"aa": (0, 0), "aa-2": (0, 0), "aa-3": (0, 0.0001), "bb": (0, 9)},
                         extra=False)
    files = {s: _span("10:02:00", "11:00:00") for s in ("aa", "aa-2", "aa-3", "bb")}
    args = ["--all-stations", "--station-file", str(reg), "--start", START, "--end", END,
            "--output", str(tmp_path / "o.csv"), "--cache-dir", str(tmp_path / "c")]
    ns = cli.build_parser().parse_args(args)
    with pytest.raises(cli.ArchiveScanInputError, match="2 distinct physical sites"):
        cli.run(ns, args, archive=MockArchive(files), now=lambda: NOW)
    # With --min-stations 2 the pair aa-family + bb is found, but never aa+aa-2.
    code, m, paths = nscan(tmp_path, MockArchive(files), "--min-stations", "2", registry=reg)
    assert code == 0 and set(triads_by_id(paths)) == {"aa+bb", "aa-2+bb", "aa-3+bb"}


def test_deterministic_regardless_of_listing_order(tmp_path):
    outputs = []
    reg = write_registry(tmp_path / "registry.yaml")
    for seed in (1, 2):
        _, _, paths = nscan(tmp_path / str(seed), MockArchive(demo_files(), shuffle_seed=seed), registry=reg)
        outputs.append([paths[k].read_bytes() for k in ("triads_csv", "triad_segments_csv", "active_segments_csv")])
    assert outputs[0] == outputs[1]


def test_max_triads_limits_export_only(tmp_path):
    code, m, paths = nscan(tmp_path, MockArchive(demo_files()), "--max-triads", "2")
    assert len(read_csv(paths["triads_csv"])) == 2
    assert m["summary"]["n_triads"] == 7 and m["summary"]["n_triads_exported"] == 2


def test_outputs_hashes_and_handoff_commands(tmp_path):
    code, m, paths = nscan(tmp_path, MockArchive(demo_files()))
    import hashlib

    for name, digest in m["output_sha256"].items():
        assert hashlib.sha256(paths[name].read_bytes()).hexdigest() == digest
    assert set(m["outputs"]) == set(paths) and m["ranking"]["version"] == net.NETWORK_RANKING_VERSION
    assert m["site_grouping"]["groups"]["SITE_cc"] == ["cc", "cc-2"]
    assert len(m["station_file_sha256"]) == 64 and m["listing"]["snapshots"][0]["checksum"]
    top = read_csv(paths["triads_csv"])[0]
    assert top["handoff_status"] == "PROPOSED_NOT_RUN"
    cmd = top["stage1_command"]
    assert cmd.startswith("tart-common-gnss --stations aa bb cc-2 --start 2026-10-09T10:30:00Z "
                          "--end 2026-10-09T11:30:00Z --min-elevation 20 --step 10 --tle-file")
    assert "--station-file" in cmd
    assert top["archive_scan_command"].startswith("tart-archive-scan --stations aa bb cc-2 --visibility-results")
    doc = json.loads(paths["triads_json"].read_text())
    assert doc["triads"][0]["segments"][0]["file_counts_per_station"] == {"aa": 63, "bb": 59, "cc-2": 63}  # Stage-2 rule: [start, end + span]
    assert {r["triad_id"] for r in doc["baseline_rejected"]} == set()
    # Re-running keeps the previous manifest under its run tag (never silently overwritten).
    nscan(tmp_path, MockArchive(demo_files()), now=NOW + timedelta(hours=1))
    assert any(p.name.startswith("network_archive_scan_manifest.2026") for p in paths["manifest"].parent.iterdir())


# --------------------------------------------------------------------------- CLI errors

@pytest.mark.parametrize("extra, message", [
    (["--include-stations", "nope"], "Unknown station"),
    (["--include-stations", "aa", "bb", "--exclude-stations", "bb"], "both included and excluded"),
    (["--visibility-results", "x.json"], "cannot be used with --all-stations"),
    (["--min-stations", "1"], ">= 2"),
    (["--min-stations", "5"], "exceeds the 4 distinct physical sites"),
    (["--min-baseline-km", "-1"], "--min-baseline-km must be >= 0"),
])
def test_invalid_network_inputs(tmp_path, extra, message):
    with pytest.raises(cli.ArchiveScanInputError, match=message):
        nscan(tmp_path, MockArchive(demo_files()), *extra)


def test_cli_conflicts_and_explicit_mode_rejects_network_options(tmp_path, capsys):
    with pytest.raises(SystemExit) as exc:
        cli.build_parser().parse_args(["--all-stations", "--stations", "mg-tana"])
    assert exc.value.code == 2
    rc = cli.main(["--stations", "mg-tana", "za-rhodes", "--start", START, "--end", END,
                   "--min-baseline-km", "500", "--output", str(tmp_path / "x.csv")])
    assert rc == 2 and "--min-baseline-km needs --all-stations" in capsys.readouterr().err
    rc = cli.main(["--all-stations", "--output", str(tmp_path / "x.csv")])
    assert rc == 2 and "--start and --end are required" in capsys.readouterr().err
