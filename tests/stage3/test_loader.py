"""Stage-2 manifest validation, MS loading, UVW convention check, time binning."""

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from tart_gnss.stage3 import loader
from tart_gnss.stage3 import synthetic as syn
from tart_gnss.stage3.errors import InvalidManifestError, InvalidVisibilityDataError, MissingMeasurementSetError

from .helpers import ALT, LAT, LON, T0, needs_casacore

EXAMPLES = Path(__file__).resolve().parents[2] / "docs" / "examples" / "stage2"
VALID = EXAMPLES / "valid_20260625_NORAD_40534" / "manifest.json"
INVALID = EXAMPLES / "invalid_20261001_NORAD_38091" / "manifest.json"


def test_load_real_stage2_manifest_and_hash():
    doc, sha = loader.load_stage2_manifest(VALID)
    assert doc["three_station_valid"] is True and len(sha) == 64


def test_not_a_stage2_manifest(tmp_path):
    p = tmp_path / "m.json"
    p.write_text(json.dumps({"stage": 1}))
    with pytest.raises(InvalidManifestError):
        loader.load_stage2_manifest(p)
    p.write_text("{not json")
    with pytest.raises(InvalidManifestError):
        loader.load_stage2_manifest(p)


def test_station_aliases_resolve_tart_id_and_station_id():
    doc, _ = loader.load_stage2_manifest(VALID)
    assert loader.resolve_station_keys(doc, ["mu-udm", "SOUTH_AFRICA"]) == ["mauritius", "south_africa"]
    assert loader.resolve_station_keys(doc, None) == ["mauritius", "madagascar", "south_africa"]
    with pytest.raises(InvalidManifestError):
        loader.resolve_station_keys(doc, ["nowhere"])


def test_invalid_three_station_dataset_refused_without_opt_in():
    doc, _ = loader.load_stage2_manifest(INVALID)
    with pytest.raises(InvalidManifestError, match="not three-station valid"):
        loader.station_input(doc, "south_africa", INVALID)
    # Opting in still refuses a station whose own MS is not VALID.
    with pytest.raises(InvalidManifestError, match="not VALID"):
        loader.station_input(doc, "south_africa", INVALID, allow_invalid_stage2=True)


def test_valid_manifest_without_local_ms_raises_missing(tmp_path):
    copy = tmp_path / "manifest.json"
    shutil.copy(VALID, copy)
    doc, _ = loader.load_stage2_manifest(copy)
    with pytest.raises(MissingMeasurementSetError):
        loader.station_input(doc, "mauritius", copy, data_root=tmp_path)


def test_calibration_warnings_from_stage2_metadata():
    doc, _ = loader.load_stage2_manifest(VALID)
    doc = json.loads(json.dumps(doc))
    for st in doc["stations"].values():  # pretend the MS exist (path resolution tested separately)
        st["measurement_set"]["path"] = str(Path(__file__).parent)
    za = loader.station_input(doc, "south_africa", VALID)
    assert any("zero gain" in w for w in za.calibration_warnings) and len(za.zero_gain_antennas) == 9
    mu = loader.station_input(doc, "mauritius", VALID)
    assert any("solutions changed" in w for w in mu.calibration_warnings)
    assert mu.stage2_three_station_valid is True


def test_ms_path_resolution_relative_to_manifest(tmp_path):
    run = tmp_path / "obs" / "20260625" / "NORAD_1"
    ms = run / "alpha" / "ms" / "alpha.ms"
    ms.mkdir(parents=True)
    manifest = run / "manifest.json"
    found = loader.resolve_ms_path("data/observations/20260625/NORAD_1/alpha/ms/alpha.ms", manifest,
                                   "data/observations/20260625/NORAD_1")
    assert found == run / "alpha" / "ms" / "alpha.ms"
    assert loader.resolve_ms_path("alpha/ms/alpha.ms", manifest, None, data_root=run) == run / "alpha/ms/alpha.ms"


def _obs(layout, **kw):
    return syn.simulate_observation(layout, LAT, LON, ALT, T0 + np.arange(5.0), 200.0, 45.0, **kw)


def test_uvw_sign_mismatch_is_rejected(layout):
    obs = _obs(layout)
    kwargs = dict(station_id="x", time_mjd_seconds=obs.time_mjd_seconds, antenna1=obs.antenna1,
                  antenna2=obs.antenna2, data=obs.data, flag=obs.flag, weight=obs.weight,
                  frequency_hz=obs.frequency_hz, antenna_positions_ecef_m=obs.antenna_positions_ecef_m,
                  latitude_deg=LAT, longitude_deg=LON)
    ds = loader.dataset_from_arrays(uvw_m=obs.uvw_m, **kwargs)
    assert ds.uvw_check["consistent"]
    with pytest.raises(InvalidVisibilityDataError, match="UVW"):
        loader.dataset_from_arrays(uvw_m=-obs.uvw_m, **kwargs)  # opposite convention
    rotated = obs.uvw_m.copy()
    rotated[:, 2] = rotated[:, 0]  # not a zenith phase centre (w != 0)
    with pytest.raises(InvalidVisibilityDataError):
        loader.dataset_from_arrays(uvw_m=rotated, **kwargs)


def test_time_bins_never_fabricate_timestamps(layout):
    times = T0 + np.array([0.0, 1.03, 2.05, 7.2, 8.24])  # gap: no data between 2.05 and 7.2 s
    obs = syn.simulate_observation(layout, LAT, LON, ALT, times, 200.0, 45.0)
    ds = syn.to_dataset(obs)
    bins = list(loader.time_bins(ds, 1.0))
    assert len(bins) == 5
    bins3 = list(loader.time_bins(ds, 3.0))
    assert [np.unique(ds.time_mjd_seconds[b.rows]).size for b in bins3] == [3, 2]
    for b in bins3:
        t = ds.time_mjd_seconds[b.rows]
        assert b.start_mjd_seconds <= t.min() and t.max() < b.end_mjd_seconds


def test_select_samples_excludes_bad_data(layout):
    obs = _obs(layout, flagged_antennas=[2])
    obs.data[0, 0, 0] = np.nan
    obs.weight[3, 0] = 0.0  # baseline (0, 4)
    ds = syn.to_dataset(obs)
    rows = np.flatnonzero(ds.time_mjd_seconds == ds.time_mjd_seconds.min())
    sel = loader.select_samples(ds, rows, invalid_antennas=[5])
    s = sel.samples
    assert not np.isin(s.antenna1, [2, 5]).any() and not np.isin(s.antenna2, [2, 5]).any()
    assert np.all(np.isfinite(s.vis)) and np.all(s.weight > 0)
    assert sel.n_baselines_total == 276
    assert sel.n_baselines_used == 276 - 23 - 22 - 2  # antennas 2 and 5 removed, plus 2 bad samples
    assert np.isclose(sel.flagged_fraction, 1 - sel.n_baselines_used / 276)


@needs_casacore
def test_load_synthetic_measurement_set(tmp_path, layout):
    obs = _obs(layout, noise_sigma=0.01, flag_fraction=0.1, rng=np.random.default_rng(1))
    path = syn.write_ms(obs, tmp_path / "alpha.ms")
    ds = loader.load_station_ms(path, "alpha", LAT, LON)
    assert ds.n_rows == obs.time_mjd_seconds.size
    np.testing.assert_allclose(ds.vis[:, 0], obs.data[:, 0, 0])
    np.testing.assert_array_equal(ds.flag[:, 0], obs.flag[:, 0, 0])
    np.testing.assert_allclose(ds.frequency_hz, obs.frequency_hz)
    assert ds.uvw_check["consistent"]
    with pytest.raises(MissingMeasurementSetError):
        loader.load_station_ms(tmp_path / "nope.ms", "alpha")
