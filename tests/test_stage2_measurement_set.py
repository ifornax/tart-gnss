import subprocess
from pathlib import Path

import pytest

from stage2_helpers import make_tart_hdf5, utc
from tart_gnss.stage2 import measurement_set as ms
from tart_gnss.stage2.errors import MeasurementSetConversionError
from tart_gnss.stations import Station

ALPHA = Station("alpha", "a", -20.0, 57.0, 0.0, "test", archive_name="aa-alpha")

GOOD_SUMMARY = {"path": "x.ms", "number_of_rows": 16560, "has_time": True, "has_uvw": True, "has_data": True,
                "has_flag": True, "has_antenna_table": True, "has_spectral_window_table": True,
                "number_of_antennas": 24, "number_of_channels": [1], "channel_frequencies_hz": [1575.42e6],
                "first_time_unix": 1782381600.0, "last_time_unix": 1782381659.0, "flagged_fraction": 0.0}


def fake_runner(returncode=0, create=True, stderr=""):
    calls = []

    def runner(args, **kwargs):
        calls.append(args)
        if create:
            Path(args[args.index("--ms") + 1]).mkdir(parents=True)
        return subprocess.CompletedProcess(args, returncode, stdout="ok", stderr=stderr)

    runner.calls = calls
    return runner


def test_successful_conversion(tmp_path):
    runner = fake_runner()
    meta = ms.convert_to_ms(["a.hdf", "b.hdf"], tmp_path / "alpha.ms", ALPHA, tart2ms_command=["tart2ms"],
                            expected_antennas=24, log_path=tmp_path / "log.txt", runner=runner,
                            inspector=lambda path, python=None: dict(GOOD_SUMMARY, path=str(path)))
    args = runner.calls[0]
    assert args[:4] == ["tart2ms", "--hdf", "a.hdf", "b.hdf"] and "--uncalibrated" not in args
    assert "--skip-online-source-catalog" in args
    assert meta["valid"] and meta["number_of_rows"] == 16560 and meta["first_time_utc"].endswith("Z")
    assert meta["calibration_applied_by_tart2ms"] is True
    assert (tmp_path / "log.txt").read_text().startswith("$ tart2ms")


def test_uncalibrated_flag(tmp_path):
    runner = fake_runner()
    ms.convert_to_ms(["a.hdf"], tmp_path / "alpha.ms", ALPHA, tart2ms_command=["tart2ms"],
                     apply_calibration=False, runner=runner,
                     inspector=lambda path, python=None: dict(GOOD_SUMMARY, path=str(path)))
    assert runner.calls[0][-1] == "--uncalibrated"


def test_failed_conversion(tmp_path):
    with pytest.raises(MeasurementSetConversionError, match="exit 1"):
        ms.convert_to_ms(["a.hdf"], tmp_path / "alpha.ms", ALPHA, tart2ms_command=["tart2ms"],
                         runner=fake_runner(returncode=1, create=False, stderr="boom"))


def test_success_without_output_is_failure(tmp_path):
    with pytest.raises(MeasurementSetConversionError, match="not created"):
        ms.convert_to_ms(["a.hdf"], tmp_path / "alpha.ms", ALPHA, tart2ms_command=["tart2ms"],
                         runner=fake_runner(create=False))


def test_existing_output_never_overwritten(tmp_path):
    (tmp_path / "alpha.ms").mkdir()
    with pytest.raises(MeasurementSetConversionError, match="already exists"):
        ms.convert_to_ms(["a.hdf"], tmp_path / "alpha.ms", ALPHA, tart2ms_command=["tart2ms"],
                         runner=fake_runner())


def test_invalid_ms_rejected(tmp_path):
    bad = dict(GOOD_SUMMARY, number_of_rows=0, has_uvw=False)
    with pytest.raises(MeasurementSetConversionError, match="no rows; missing UVW column"):
        ms.convert_to_ms(["a.hdf"], tmp_path / "alpha.ms", ALPHA, tart2ms_command=["tart2ms"],
                         runner=fake_runner(), inspector=lambda path, python=None: dict(bad, path=str(path)))


def test_antenna_count_checked():
    with pytest.raises(MeasurementSetConversionError, match="expected 24"):
        ms.validate_ms(dict(GOOD_SUMMARY, number_of_antennas=12), expected_antennas=24)


def _tart2ms_available():
    try:
        ms.find_tart2ms()
        import casacore.tables  # noqa: F401
        return True
    except (MeasurementSetConversionError, ImportError):
        return False


@pytest.mark.skipif(not _tart2ms_available(), reason="tart2ms / python-casacore not installed")
def test_real_tart2ms_conversion(tmp_path):
    """End-to-end conversion of a synthetic native file with the real tart2ms."""
    src = make_tart_hdf5(tmp_path / "vis.hdf", "alpha", utc(2026, 6, 25, 10, 1), n=5)
    meta = ms.convert_to_ms([str(src)], tmp_path / "alpha.ms", ALPHA, expected_antennas=24,
                            log_path=tmp_path / "tart2ms.log")
    assert meta["number_of_rows"] == 5 * 276
    assert meta["number_of_antennas"] == 24
    assert meta["frequency_hz"] == [1575.42e6]
    from tart_gnss.stage2.selection import parse_utc

    # 5 samples, 1 s apart, ending 10:01:00 -> first at 10:00:56 (MS TIME is an MJD double: ~us precision)
    assert abs((parse_utc(meta["first_time_utc"]) - utc(2026, 6, 25, 10, 0, 56)).total_seconds()) < 1e-3
    assert abs((parse_utc(meta["last_time_utc"]) - utc(2026, 6, 25, 10, 1, 0)).total_seconds()) < 1e-3
