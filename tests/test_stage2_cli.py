import json
from datetime import timedelta
from pathlib import Path

import pytest

from stage2_helpers import STATIONS_YAML, MockArchive, populate, stage1_doc, utc
from tart_gnss.stage2 import cli
from tart_gnss.stage2.manifest import validate_manifest

WINDOW_START = utc(2026, 6, 25, 10, 0)  # rank 1 in stage1_doc: 10:00-10:20


@pytest.fixture
def env(tmp_path):
    results = tmp_path / "common_gnss_test.json"
    results.write_text(json.dumps(stage1_doc()))
    stations = tmp_path / "stations.yaml"
    stations.write_text(STATIONS_YAML)
    return tmp_path, results, stations


def fake_converter(files, output_ms, station, **kwargs):
    Path(output_ms).mkdir(parents=True)
    return {"path": str(output_ms), "valid": True, "number_of_rows": 276 * 60 * len(files),
            "number_of_antennas": 24, "calibration_applied_by_tart2ms": kwargs["apply_calibration"]}


def run(env, archive, *extra, converter=fake_converter):
    tmp, results, stations = env
    argv = ["--visibility-results", str(results), "--rank", "1", "--padding", "120",
            "--output-dir", str(tmp / "obs"), "--station-file", str(stations), *extra]
    args = cli.build_parser().parse_args(argv)
    code = cli.run(args, argv, archive=archive, converter=converter)
    manifest = tmp / "obs" / "20260625" / "NORAD_40534" / "manifest.json"
    return code, (json.loads(manifest.read_text()) if manifest.exists() else None)


def full_archive(tmp):
    arch = MockArchive()
    start = WINDOW_START - timedelta(minutes=3)
    populate(arch, tmp / "src", "alpha", start, 26)
    populate(arch, tmp / "src", "beta", start, 26, offset_s=0.25)
    populate(arch, tmp / "src", "gamma", start, 26, offset_s=0.5)
    return arch


def test_full_three_station_run(env):
    code, m = run(env, full_archive(env[0]), "--download", "--convert-ms", "--keep-native")
    assert code == cli.EXIT_VALID
    assert m["three_station_valid"] is True and m["validity_reasons"] == []
    assert validate_manifest(m) == []
    s = m["synchronisation"]
    assert s["actual_common_duration_seconds"] > 20 * 60
    assert s["three_station_coverage_fraction"] == pytest.approx(1.0)
    assert s["station_pair_time_offset_statistics"]["alpha__beta"]["median_offset_s"] == pytest.approx(0.25, abs=1e-3)
    for sid, st in m["stations"].items():
        assert st["archive_query_status"] == "FOUND" and st["measurement_set_status"] == "VALID"
        assert all(len(f["sha256"]) == 64 and f["retained"] for f in st["raw_files"])
        assert st["observation"]["coverage"]["frequencies_hz"] == [1575.42e6]
    assert m["satellite"]["signal_relevance"] == "GPS_L1"
    assert m["prediction"]["query_start_utc"] == "2026-06-25T09:58:00.000Z"
    run_dir = env[0] / "obs" / "20260625" / "NORAD_40534"
    assert (run_dir / "stage1_window.json").is_file()
    assert (run_dir / "alpha" / "ms" / "alpha.ms").is_dir()


def test_missing_station_is_not_three_station_valid(env):
    arch = full_archive(env[0])
    arch.objects = {k: v for k, v in arch.objects.items() if not k.startswith("bb-beta")}
    code, m = run(env, arch, "--download", "--convert-ms")
    assert code == cli.EXIT_NOT_VALID
    assert m["three_station_valid"] is False and m["synchronisation"]["three_station_valid"] is False
    assert m["stations"]["beta"]["archive_query_status"] == "NOT_FOUND"
    assert any("beta" in r for r in m["validity_reasons"])
    assert arch.downloads == [], "nothing is downloaded without --allow-partial"


def test_allow_partial_still_invalid(env):
    arch = full_archive(env[0])
    arch.objects = {k: v for k, v in arch.objects.items() if not k.startswith("bb-beta")}
    code, m = run(env, arch, "--download", "--allow-partial")
    assert code == cli.EXIT_NOT_VALID and m["three_station_valid"] is False
    assert arch.downloads and m["stations"]["alpha"]["raw_files"]


def test_no_common_time(env):
    tmp = env[0]
    arch = MockArchive()
    populate(arch, tmp / "src", "alpha", WINDOW_START - timedelta(minutes=2), 8)
    populate(arch, tmp / "src", "beta", WINDOW_START + timedelta(minutes=10), 8)
    populate(arch, tmp / "src", "gamma", WINDOW_START - timedelta(minutes=2), 24)
    code, m = run(env, arch, "--download")
    assert code == cli.EXIT_NOT_VALID and m["three_station_valid"] is False
    assert m["synchronisation"]["actual_common_start_utc"] is None
    assert any("no common overlap" in r for r in m["validity_reasons"])


def test_query_only_is_not_valid(env):
    code, m = run(env, full_archive(env[0]))
    assert code == cli.EXIT_NOT_VALID
    assert any("not retrieved" in r for r in m["validity_reasons"])


def test_dry_run_writes_nothing(env, capsys):
    arch = full_archive(env[0])
    code, m = run(env, arch, "--dry-run", "--download", "--convert-ms")
    assert m is None and not (env[0] / "obs").exists() and arch.downloads == []
    assert "DRY RUN" in capsys.readouterr().out


def test_conversion_failure_marks_invalid(env):
    from tart_gnss.stage2.errors import MeasurementSetConversionError

    def failing(files, output_ms, station, **kw):
        if station.id == "gamma":
            raise MeasurementSetConversionError("tart2ms failed for gamma (exit 1)")
        return fake_converter(files, output_ms, station, **kw)

    code, m = run(env, full_archive(env[0]), "--download", "--convert-ms", converter=failing)
    assert code == cli.EXIT_NOT_VALID
    assert m["stations"]["gamma"]["measurement_set_status"] == "FAILED"
    # native files are kept for the failed station but removed for converted ones (no --keep-native)
    assert all(f["retained"] for f in m["stations"]["gamma"]["raw_files"])
    assert not any(f["retained"] for f in m["stations"]["alpha"]["raw_files"])


def test_rerun_reuses_files_and_keeps_old_manifest(env):
    arch = full_archive(env[0])
    run(env, arch, "--download")
    n = len(arch.downloads)
    code, m = run(env, arch, "--download")
    assert len(arch.downloads) == n, "verified local files are reused"
    run_dir = env[0] / "obs" / "20260625" / "NORAD_40534"
    assert len(list(run_dir.glob("manifest_*.json"))) == 1


def test_main_reports_selection_errors(env, capsys):
    tmp, results, stations = env
    code = cli.main(["--visibility-results", str(results), "--rank", "99", "--station-file", str(stations),
                     "--output-dir", str(tmp / "obs"), "--dry-run"])
    assert code == cli.EXIT_USAGE and "Rank 99 not found" in capsys.readouterr().err


def test_convert_requires_download(env, capsys):
    tmp, results, stations = env
    code = cli.main(["--visibility-results", str(results), "--rank", "1", "--station-file", str(stations),
                     "--output-dir", str(tmp / "obs"), "--convert-ms"])
    assert code == cli.EXIT_USAGE


def test_manifest_validation_catches_problems():
    assert validate_manifest({}) and "missing 'stage'" in validate_manifest({})[0]
    good = {"stage": 2, "manifest_version": "1", "three_station_valid": True, "validity_reasons": [],
            "satellite": {}, "prediction": {}, "stations": {"a": {}, "b": {}}, "software": {},
            "provenance": {}, "synchronisation": {"three_station_valid": False}}
    problems = validate_manifest(good)
    assert any("expected 3 stations" in p for p in problems)
    assert any("differs between top level" in p for p in problems)
    assert any("satellite.norad_id" in p for p in problems)
