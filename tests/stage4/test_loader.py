"""Loading Stage-3 outputs, frozen coordinates, quality policy and sight lines."""

import json

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage4 import synthetic as syn  # noqa: E402
from tart_gnss.stage4.errors import InvalidStage3ManifestError, MissingDirectionDataError  # noqa: E402
from tart_gnss.stage4.geometry import geodetic_to_ecef  # noqa: E402
from tart_gnss.stage4.loader import (  # noqa: E402
    apply_quality_policy,
    load_stage3_manifest,
    load_stage3_run,
    sha256_file,
    sight_line,
)


@pytest.fixture
def run_dir(tmp_path):
    times = {sid: [syn.T0_UNIX_S + k for k in range(5)] for sid in syn.EXAMPLE_STATIONS}
    rows = syn.make_directions(syn.linear_track(), times, 0.1, None)
    rows[0]["quality_flags"] = ["LOW_SNR"]
    rows[1]["quality_flags"] = ["FIT_FAILED"]
    rows[1]["fit_success"] = False
    rows[2]["quality_flags"] = ["IMAGE_DISAGREEMENT"]
    return syn.write_stage3_run(tmp_path / "s3", rows)


def test_load_run_uses_manifest_coordinates(run_dir):
    run = load_stage3_run(run_dir)
    assert set(run.stations) == set(syn.EXAMPLE_STATIONS)
    st = run.stations["south_africa"]
    assert (st.latitude_deg, st.longitude_deg, st.altitude_m) == syn.EXAMPLE_STATIONS["south_africa"]
    assert np.allclose(st.ecef_m, geodetic_to_ecef(*syn.EXAMPLE_STATIONS["south_africa"]))
    assert "stations.south_africa" in st.source
    assert len(run.directions) == 15 and run.norad_id == 40534
    assert run.sha256 == sha256_file(run_dir)


def test_station_subset_and_unknown(run_dir):
    run = load_stage3_run(run_dir, ["mauritius", "madagascar"])
    assert set(run.directions["station_id"]) == {"mauritius", "madagascar"}
    with pytest.raises(InvalidStage3ManifestError):
        load_stage3_run(run_dir, ["atlantis"])


def test_csv_fallback(tmp_path):
    rows = syn.make_directions(syn.linear_track(), {"mauritius": [syn.T0_UNIX_S]}, 0.1, None)
    m = syn.write_stage3_run(tmp_path, rows, csv_only=True)
    run = load_stage3_run(m)
    assert run.directions_path.suffix == ".csv"


def test_quality_policy(run_dir):
    run = load_stage3_run(run_dir)
    kept, rej = apply_quality_policy(run.directions, "strict")
    assert len(kept) == 13 and rej == {"LOW_SNR": 1, "FIT_FAILED": 1}
    kept, rej = apply_quality_policy(run.directions, "permissive")
    assert len(kept) == 14 and rej == {"FIT_FAILED": 1}


def test_invalid_manifests(tmp_path, run_dir):
    with pytest.raises(InvalidStage3ManifestError):
        load_stage3_manifest(tmp_path / "missing.json")
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"stage": 2, "schema_version": 1, "stations": {"a": {}}}))
    with pytest.raises(InvalidStage3ManifestError):
        load_stage3_manifest(bad)
    doc = json.loads(run_dir.read_text())
    doc["stations"]["mauritius"].pop("latitude_deg")
    run_dir.write_text(json.dumps(doc))
    with pytest.raises(InvalidStage3ManifestError):
        load_stage3_run(run_dir)


def test_missing_columns(tmp_path, run_dir):
    df = pd.read_csv(run_dir.parent / "directions.csv").drop(columns=["cov_lm"])
    (run_dir.parent / "directions.parquet").unlink(missing_ok=True)
    df.to_csv(run_dir.parent / "directions.csv", index=False)
    with pytest.raises(MissingDirectionDataError):
        load_stage3_run(run_dir)


def test_stage2_cross_check(tmp_path, run_dir):
    s2 = tmp_path / "stage2_manifest.json"
    stations = {sid: {"latitude_deg": lat, "longitude_deg": lon, "altitude_m": alt}
                for sid, (lat, lon, alt) in syn.EXAMPLE_STATIONS.items()}
    s2.write_text(json.dumps({"stations": stations}))
    doc = json.loads(run_dir.read_text())
    doc["input_manifest"] = str(s2)
    doc["input_manifest_sha256"] = sha256_file(s2)
    run_dir.write_text(json.dumps(doc))
    assert load_stage3_run(run_dir).stations["mauritius"].cross_check.startswith("verified")
    # The frozen Stage-2 manifest disagrees with the Stage-3 copy: refuse.
    stations["mauritius"]["latitude_deg"] += 0.01
    s2.write_text(json.dumps({"stations": stations}))
    doc["input_manifest_sha256"] = sha256_file(s2)
    run_dir.write_text(json.dumps(doc))
    with pytest.raises(InvalidStage3ManifestError):
        load_stage3_run(run_dir)
    # A Stage-2 manifest that changed after Stage 3 ran is not trusted for the check.
    doc["input_manifest_sha256"] = "0" * 64
    run_dir.write_text(json.dumps(doc))
    assert "SHA-256 mismatch" in load_stage3_run(run_dir).stations["mauritius"].cross_check


def test_sight_line_lm_and_azel_and_zenith(run_dir):
    run = load_stage3_run(run_dir)
    row = run.directions.iloc[3].copy()
    st = run.stations[row["station_id"]]
    s1, c1, src1 = sight_line(row, st)
    row2 = row.copy()
    row2["l"] = row2["m"] = np.nan
    s2, c2, src2 = sight_line(row2, st)
    assert (src1, src2) == ("lm", "azel")
    assert np.allclose(s1, s2, atol=1e-9)
    assert np.allclose(c1, c2, rtol=1e-6, atol=1e-16)
    zen = row.copy()
    zen["l"], zen["m"], zen["sigma_az_deg"], zen["cov_az_el"] = 0.0, 0.0, np.nan, np.nan
    s3, c3, src3 = sight_line(zen, st)
    assert src3 == "lm" and np.all(np.isfinite(c3))
