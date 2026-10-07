import csv
import json

import pytest

from tart_gnss.cli import main, parse_utc, CLIError


def base_args(station_yaml, tle_path, tmp_path, **kw):
    args = [
        "--stations", "alpha", "beta", "gamma",
        "--start", kw.get("start", "2026-10-01T00:00:00Z"),
        "--end", kw.get("end", "2026-10-01T12:00:00Z"),
        "--station-file", str(station_yaml),
        "--tle-file", str(tle_path),
        "--output", str(tmp_path / "out" / "common_gnss_test.csv"),
        "--step", "30",
    ]
    return args


def test_start_after_end_fails(station_yaml, sample_tle_path, tmp_path, capsys):
    args = base_args(station_yaml, sample_tle_path, tmp_path,
                     start="2026-10-02T00:00:00Z", end="2026-10-01T00:00:00Z")
    assert main(args) == 2
    assert "must be later" in capsys.readouterr().err


def test_invalid_date_fails(station_yaml, sample_tle_path, tmp_path, capsys):
    args = base_args(station_yaml, sample_tle_path, tmp_path, start="yesterday")
    assert main(args) == 2
    assert "ISO-8601" in capsys.readouterr().err


def test_unknown_station_fails(station_yaml, sample_tle_path, tmp_path, capsys):
    args = base_args(station_yaml, sample_tle_path, tmp_path)
    args[args.index("gamma")] = "atlantis"
    assert main(args) == 2
    assert "Unknown station 'atlantis'" in capsys.readouterr().err


def test_wrong_number_of_stations_fails(station_yaml, sample_tle_path, tmp_path, capsys):
    args = base_args(station_yaml, sample_tle_path, tmp_path)
    args.remove("gamma")
    assert main(args) == 2
    assert "exactly 3" in capsys.readouterr().err


def test_missing_tle_file_fails(station_yaml, tmp_path):
    args = base_args(station_yaml, tmp_path / "missing.tle", tmp_path)
    assert main(args) == 1


def test_parse_utc_conventions():
    assert parse_utc("2026-10-01T00:00:00Z") == parse_utc("2026-10-01T00:00:00")
    assert parse_utc("2026-10-01T02:00:00+02:00") == parse_utc("2026-10-01T00:00:00Z")
    with pytest.raises(CLIError):
        parse_utc("2026-13-01T00:00:00")


def test_end_to_end(station_yaml, sample_tle_path, tmp_path, capsys):
    assert main(base_args(station_yaml, sample_tle_path, tmp_path)) == 0
    out = tmp_path / "out"
    rows = list(csv.DictReader(open(out / "common_gnss_test.csv")))
    doc = json.loads((out / "common_gnss_test.json").read_text())
    meta = json.loads((out / "run_metadata_test.json").read_text())
    assert rows, "the sample includes a geostationary BeiDou satellite visible from all test sites"
    durations = [float(r["duration_seconds"]) for r in rows]
    assert durations == sorted(durations, reverse=True)
    assert [int(r["rank"]) for r in rows] == list(range(1, len(rows) + 1))
    assert all(r["coordinate_frame"] == "ITRS" for r in rows)

    for w in doc["windows"]:
        assert w["tle"]["line1"].startswith("1 ") and w["tle"]["line2"].startswith("2 ")
        assert w["satellite_position"]["coordinate_frame"] == "ITRS"
        # every station sees the satellite at >= 20 deg (to refinement tolerance) at all epochs
        for geom in w["stations"].values():
            for epoch in ("start", "mid", "end"):
                assert geom[epoch]["elevation_deg"] >= 20.0 - 1e-3
        # a refined (non-clipped) boundary sits on the threshold at one station
        if not w["start_clipped_by_search_interval"]:
            assert min(g["start"]["elevation_deg"] for g in w["stations"].values()) == pytest.approx(20.0, abs=0.01)
    assert meta["tle_sha256"] and meta["software"]["skyfield"]
    assert meta["station_ids"] == ["alpha", "beta", "gamma"]
    assert "COMMON GNSS VISIBILITY" in capsys.readouterr().out


def test_no_windows_reported_honestly(tmp_path, sample_tle_path, capsys):
    stations = tmp_path / "far.yaml"
    stations.write_text(
        "stations:\n"
        "  north: {name: n, latitude_deg: 80, longitude_deg: 0, altitude_m: 0, source: test}\n"
        "  south: {name: s, latitude_deg: -80, longitude_deg: 180, altitude_m: 0, source: test}\n"
        "  east: {name: e, latitude_deg: 0, longitude_deg: 90, altitude_m: 0, source: test}\n"
    )
    args = [
        "--stations", "north", "south", "east", "--start", "2026-10-01T00:00:00Z",
        "--end", "2026-10-01T01:00:00Z", "--station-file", str(stations),
        "--tle-file", str(sample_tle_path), "--output", str(tmp_path / "r.csv"),
    ]
    assert main(args) == 0
    assert "No common GNSS visibility windows found." in capsys.readouterr().out
