"""``tart-localise`` end to end on a synthetic Measurement Set, plus exit codes."""

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("pandas")

from tart_gnss.stage3 import cli  # noqa: E402
from tart_gnss.stage3 import synthetic as syn  # noqa: E402
from tart_gnss.stage3.geometry import angular_separation_deg, azel_to_unit  # noqa: E402
from tart_gnss.stage3.output import read_directions  # noqa: E402

from .helpers import needs_casacore  # noqa: E402

EXAMPLES = Path(__file__).resolve().parents[2] / "docs" / "examples" / "stage2"
VALID = EXAMPLES / "valid_20260625_NORAD_40534" / "manifest.json"
INVALID = EXAMPLES / "invalid_20261001_NORAD_38091" / "manifest.json"
TLE = (json.loads(VALID.read_text())["satellite"]["tle_line1"], json.loads(VALID.read_text())["satellite"]["tle_line2"])
SITE = {"latitude_deg": -20.2587508, "longitude_deg": 57.7591989, "altitude_m": 20.0}


def _synthetic_run(tmp_path, offset_az=5.0, offset_el=2.0, n_t=6, three_valid=True, ms_status="VALID"):
    """Synthetic station observing a source offset from the real TLE track of NORAD 40534."""
    t = syn.utc_to_mjd_seconds(2026, 6, 25, 9, 14, 0.0) + np.arange(n_t) * 1.0
    tle = cli.TleReference(*TLE, SITE["latitude_deg"], SITE["longitude_deg"], SITE["altitude_m"])
    az_t, el_t, _ = tle.predict(t - 3506716800.0)
    obs = syn.simulate_observation(syn.tart_like_layout(), SITE["latitude_deg"], SITE["longitude_deg"],
                                   SITE["altitude_m"], t, az_t + offset_az, el_t + offset_el, amplitude=0.5 + 0.3j,
                                   noise_sigma=0.3, rng=np.random.default_rng(1))
    run = tmp_path / "obs"
    (run / "alpha" / "ms").mkdir(parents=True)
    syn.write_ms(obs, run / "alpha" / "ms" / "alpha.ms")
    manifest = syn.write_manifest(
        run / "manifest.json",
        {"alpha": {**SITE, "ms_path": "alpha/ms/alpha.ms", "n_antennas": 24, "tart_id": "zz-alpha",
                   "first_time_utc": "2026-06-25T09:14:00Z", "last_time_utc": "2026-06-25T09:14:05Z",
                   "measurement_set_status": ms_status}},
        norad_id=40534, tle_lines=TLE, three_station_valid=three_valid)
    return manifest, run, obs


@needs_casacore
def test_end_to_end_synthetic(tmp_path):
    manifest, run, obs = _synthetic_run(tmp_path)
    out = tmp_path / "results"
    code = cli.main(["--manifest", str(manifest), "--station", "zz-alpha", "--output-dir", str(out),
                     "--image-validation", "--bootstrap", "5", "--data-root", str(run)])
    assert code == cli.EXIT_OK
    d = out / "20260625" / "NORAD_40534"
    m = json.loads((d / "stage3_manifest.json").read_text())
    assert m["stage"] == 3 and m["schema_version"] == 1 and len(m["input_manifest_sha256"]) == 64
    assert m["estimator"]["robust_loss"] == "huber" and m["estimator"]["bootstrap"]["random_seed"] == 20260625
    assert m["stations"]["alpha"]["near_field"]["far_field_valid"] is True
    assert m["stations"]["alpha"]["uvw_check"]["consistent"] is True
    assert "git_commit" in m["software"] and m["software"]["scipy"]
    assert (d / "logs" / "stage3.log").read_text().count("alpha") > 2
    assert len(list((d / "diagnostics").glob("alpha_*.png"))) == 4
    df = read_directions(d / "directions.csv")
    assert len(df) == 6
    for _, r in df.iterrows():
        k = int(round((r["time_mjd"] * 86400 - obs.time_mjd_seconds[0])))
        true = azel_to_unit(obs.true_az_deg[k], obs.true_el_deg[k])
        err = angular_separation_deg(azel_to_unit(r["az_obs_deg"], r["el_obs_deg"]), true)
        tle_sep = r["tle_separation_deg"]
        assert err < 0.3 and tle_sep > 4.0  # measurement-driven, not the TLE direction
        assert r["vis_image_separation_deg"] < 1.0
        assert r["usable"] and np.isfinite(r["sigma_az_deg"]) and np.isfinite(r["az_boot_p50"])


@needs_casacore
def test_invalid_stage2_requires_opt_in_and_is_marked(tmp_path):
    manifest, run, _ = _synthetic_run(tmp_path, n_t=2, three_valid=False)
    args = ["--manifest", str(manifest), "--all-stations", "--output-dir", str(tmp_path / "r"), "--no-plots",
            "--data-root", str(run)]
    assert cli.main(args) == cli.EXIT_USAGE
    assert cli.main(args + ["--allow-invalid-stage2"]) == cli.EXIT_OK
    df = read_directions(tmp_path / "r" / "20260625" / "NORAD_40534" / "directions.csv")
    assert not df["stage2_three_station_valid"].any()


@needs_casacore
def test_no_usable_estimates_exit_3(tmp_path):
    manifest, run, _ = _synthetic_run(tmp_path, n_t=2)
    code = cli.main(["--manifest", str(manifest), "--all-stations", "--output-dir", str(tmp_path / "r"),
                     "--no-plots", "--min-baselines", "1000", "--data-root", str(run)])
    assert code == cli.EXIT_NO_ESTIMATES


def test_real_invalid_example_manifest_exits_2(tmp_path):
    assert cli.main(["--manifest", str(INVALID), "--all-stations", "--output-dir", str(tmp_path)]) == 2


def test_valid_manifest_without_ms_exits_2(tmp_path):
    copy = tmp_path / "manifest.json"
    shutil.copy(VALID, copy)
    assert cli.main(["--manifest", str(copy), "--station", "mu-udm", "--output-dir", str(tmp_path / "r"),
                     "--data-root", str(tmp_path)]) == 2


def test_usage_errors_exit_2(tmp_path):
    assert cli.main(["--manifest", str(VALID)]) == 2  # neither --station nor --all-stations
    assert cli.main(["--manifest", str(tmp_path / "missing.json"), "--all-stations"]) == 2
    assert cli.main(["--manifest", str(VALID), "--station", "nowhere"]) == 2


def test_near_field_diagnostic_values():
    d = cli.near_field_diagnostic(4.9, 1575.42e6, 20000.0, "test")
    assert 5e-6 < d["eta_B2_over_lambda_R"] < 7e-6 and d["far_field_valid"]
    assert not cli.near_field_diagnostic(4.9, 1575.42e6, 0.01, "test")["far_field_valid"]


@needs_casacore
def test_unrelated_source_is_measured_but_flagged_unassociated(tmp_path):
    """A source ~8 deg from the TLE is fitted where it is, then flagged TLE_UNASSOCIATED (not usable)."""
    manifest, run, obs = _synthetic_run(tmp_path, offset_az=10.0, offset_el=5.0, n_t=2)
    out = tmp_path / "r"
    code = cli.main(["--manifest", str(manifest), "--all-stations", "--output-dir", str(out), "--no-plots",
                     "--data-root", str(run)])
    assert code == cli.EXIT_NO_ESTIMATES
    df = read_directions(out / "20260625" / "NORAD_40534" / "directions.csv")
    for _, r in df.iterrows():
        assert "TLE_UNASSOCIATED" in r["quality_flags"] and not r["usable"]
        err = angular_separation_deg(azel_to_unit(r["az_obs_deg"], r["el_obs_deg"]),
                                     azel_to_unit(obs.true_az_deg[0], obs.true_el_deg[0]))
        assert err < 0.5 and r["tle_separation_deg"] > 7.0
    # disabling the association check makes the same measurements usable
    assert cli.main(["--manifest", str(manifest), "--all-stations", "--output-dir", str(tmp_path / "r2"),
                     "--no-plots", "--data-root", str(run), "--association-radius-deg", "0"]) == cli.EXIT_OK
