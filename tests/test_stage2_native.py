from datetime import timedelta

import h5py
import numpy as np
import pytest

from stage2_helpers import make_tart_hdf5, utc
from tart_gnss.stage2.errors import ObservationValidationError
from tart_gnss.stage2.native import inspect_hdf5, summarise_coverage
from tart_gnss.stations import Station

ALPHA = Station("alpha", "a", -20.0, 57.0, 0.0, "test", archive_name="aa-alpha")
T = utc(2026, 6, 25, 10, 1)


def test_inspect_valid_file(tmp_path):
    obs = inspect_hdf5(make_tart_hdf5(tmp_path / "a.hdf", "alpha", T, zero_gains=(3,)), ALPHA)
    assert obs.n_timestamps == 60 and obs.last_timestamp_utc == T
    assert obs.n_antennas == 24 and obs.n_baselines == 276 and obs.vis_shape == (60, 276)
    assert obs.frequency_hz == 1575.42e6 and obs.frequency_consistent is True
    assert obs.median_spacing_s == pytest.approx(1.0)
    assert obs.zero_gain_antennas == [3] and obs.calibration_fingerprint
    assert obs.position_offset_deg == 0.0 and not obs.warnings


def test_frequency_outside_band_flagged(tmp_path):
    obs = inspect_hdf5(make_tart_hdf5(tmp_path / "a.hdf", "alpha", T, frequency=1602e6), ALPHA)
    assert obs.frequency_consistent is False and obs.warnings


def test_wrong_station_rejected(tmp_path):
    with pytest.raises(ObservationValidationError, match="differs from station"):
        inspect_hdf5(make_tart_hdf5(tmp_path / "a.hdf", "alpha", T, lat=-25.0), ALPHA)


def test_non_monotonic_timestamps_rejected(tmp_path):
    ts = [T, T - timedelta(seconds=1)]
    with pytest.raises(ObservationValidationError, match="increasing"):
        inspect_hdf5(make_tart_hdf5(tmp_path / "a.hdf", "alpha", T, timestamps=ts), ALPHA)


def test_missing_vis_rejected(tmp_path):
    p = make_tart_hdf5(tmp_path / "a.hdf", "alpha", T)
    with h5py.File(p, "a") as f:
        del f["vis"]
    with pytest.raises(ObservationValidationError, match="missing 'vis'"):
        inspect_hdf5(p, ALPHA)


def test_all_zero_vis_rejected(tmp_path):
    p = make_tart_hdf5(tmp_path / "a.hdf", "alpha", T)
    with h5py.File(p, "a") as f:
        f["vis"][...] = 0
    with pytest.raises(ObservationValidationError, match="all zero"):
        inspect_hdf5(p, ALPHA)


def test_not_hdf5(tmp_path):
    p = tmp_path / "x.hdf"
    p.write_text("not hdf5")
    with pytest.raises(ObservationValidationError, match="not a readable HDF5"):
        inspect_hdf5(p, ALPHA)


def test_coverage_detects_gap_and_merges_files(tmp_path):
    a = inspect_hdf5(make_tart_hdf5(tmp_path / "a.hdf", "alpha", T), ALPHA)
    b = inspect_hdf5(make_tart_hdf5(tmp_path / "b.hdf", "alpha", T + timedelta(seconds=60)), ALPHA)
    c = inspect_hdf5(make_tart_hdf5(tmp_path / "c.hdf", "alpha", T + timedelta(seconds=300)), ALPHA)
    cov = summarise_coverage("alpha", [c, a, b])
    assert cov.number_of_timestamps == 180
    assert cov.median_time_spacing_seconds == pytest.approx(1.0)
    assert len(cov.gaps) == 1 and len(cov.segments) == 2
    assert cov.segments[0] == (a.first_timestamp_utc, b.last_timestamp_utc)
    assert cov.frequencies_hz == [1575.42e6] and cov.n_antennas == [24]


def test_coverage_requires_data():
    with pytest.raises(ObservationValidationError):
        summarise_coverage("alpha", [])
