"""``tart-validate`` end to end: exit codes, sources, flags, provenance, no modification of Stage 5."""

import json

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("skyfield")

from tart_gnss.stage4.loader import sha256_file  # noqa: E402
from tart_gnss.stage6 import cli  # noqa: E402

from .helpers import L1, L2, TLE_EPOCH_S, inject_rtn, run_cli, stage5_run, times, tle_file, truth  # noqa: E402


def _fresh(tmp_path, n=12, **kw):
    t = times(n)
    return stage5_run(tmp_path, t, inject_rtn(truth(t), (500.0, 0, 0)), **kw), t


def test_fresh_reference_exit_0(tmp_path):
    m, _ = _fresh(tmp_path)
    code, res, summ, man = run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path)), "--no-plots")
    assert code == cli.EXIT_OK
    assert summ["precision_validation"] is True and man["precision_validation"]["value"] is True
    assert res["precision_validation"].all() and res["valid"].all()
    assert np.allclose(res["dr_radial_m"], 500.0, atol=0.5)
    assert man["residual_definition"] == "observed_minus_tle" == summ["residual_definition"]
    assert man["tle_reference"]["line1"] == L1 and len(man["tle_reference"]["sha256"]) == 64
    assert "sgp4" in man["tle_reference"]["sgp4_implementation"] and man["frames"]["sgp4_source_frame"].startswith("TEME")
    assert man["stage5_inputs_unchanged"] is True
    # Single TLE: no reference covariance pretended.
    assert res["position_mahalanobis_total_cov"].isna().all() and res["cov_tle_dispersion_x_x"].isna().all()
    assert all("REFERENCE_UNCERTAINTY_UNKNOWN" in f for f in res["quality_flags"])
    # No residual thresholds given -> no HIGH_*_RESIDUAL flags.
    assert not any(f.startswith("HIGH_POSITION") for fl in res["quality_flags"] for f in fl)


def test_stale_reference_exit_3_and_allow(tmp_path):
    t = times(6, t0=TLE_EPOCH_S - 103 * 86400)
    m = stage5_run(tmp_path, t, truth(t))
    f = str(tle_file(tmp_path))
    code, res, summ, _ = run_cli(tmp_path, m, "--tle-file", f, "--no-plots")
    assert code == cli.EXIT_NOT_PRECISION and summ["precision_validation"] is False
    for flags in res["quality_flags"]:
        assert {"TLE_STALE", "TLE_BACKPROPAGATED", "NOT_PRECISION_VALIDATION"} <= set(flags)
    assert (res["tle_age_days"] < -102.9).all() and (res["tle_propagation_direction"] == "backward").all()
    code2, _, summ2, _ = run_cli(tmp_path, m, "--tle-file", f, "--no-plots", "--allow-stale-reference")
    assert code2 == cli.EXIT_OK and summ2["precision_validation"] is False


def test_thresholds_set_flags(tmp_path):
    m, _ = _fresh(tmp_path)
    _, res, _, _ = run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path)), "--no-plots",
                           "--position-residual-warning-km", "0.1", "--velocity-residual-warning-m-s", "1000")
    assert all("HIGH_POSITION_RESIDUAL" in f for f in res["quality_flags"])
    assert not any("HIGH_VELOCITY_RESIDUAL" in f for f in res["quality_flags"])


def test_source_arguments_are_exclusive_and_required(tmp_path):
    m, _ = _fresh(tmp_path, 3)
    base = ["--stage5-manifest", str(m), "--output-dir", str(tmp_path / "o")]
    assert cli.main(base) == cli.EXIT_USAGE
    assert cli.main(base + ["--tle-file", str(tle_file(tmp_path)), "--tle-source", "stage1"]) == cli.EXIT_USAGE


def test_missing_and_mismatched_tle_exit_2(tmp_path):
    m, _ = _fresh(tmp_path, 3)
    assert run_cli(tmp_path, m, "--tle-file", str(tmp_path / "none.tle"))[0] == cli.EXIT_USAGE
    other = ("1 25544U 98067A   26279.50000000  .00010000  00000+0  18000-3 0  9990",
             "2 25544  51.6400 100.0000 0005000  90.0000 270.0000 15.50000000400000")
    from tart_gnss.tle import tle_checksum

    other = tuple(x[:68] + str(tle_checksum(x)) for x in other)
    assert run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path, [other], "o.tle")))[0] == cli.EXIT_USAGE


def test_stage1_source(tmp_path):
    m4 = tmp_path / "stage4_manifest.json"
    m4.write_text(json.dumps({"stage": 4, "satellite_reference": {"satellite": {
        "satellite_name": "GPS BIIF-9", "tle_line1": L1, "tle_line2": L2}}}))
    t = times(5)
    m = stage5_run(tmp_path, t, truth(t), stage4_manifest=str(m4))
    code, res, _, man = run_cli(tmp_path, m, "--tle-source", "stage1", "--no-plots")
    assert code == cli.EXIT_OK and man["tle_reference"]["source_kind"] == "stage1"
    assert np.all(res["position_residual_3d_m"] < 1.0)


def test_stage5_not_modified_and_plots(tmp_path):
    from tart_gnss.stage6.plotting import PLOT_NAMES

    pytest.importorskip("matplotlib")
    m, _ = _fresh(tmp_path, 8)
    before = {p.name: sha256_file(p) for p in m.parent.iterdir() if p.is_file()}
    code, _, _, man = run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path)))
    after = {p.name: sha256_file(p) for p in m.parent.iterdir() if p.is_file()}
    assert code == 0 and before == after
    assert sorted(p.split("/")[-1][:-4] for p in man["outputs"]["diagnostics"]) == sorted(PLOT_NAMES)


def test_strict_policy_skips_invalid_rows(tmp_path):
    t = times(6)
    m = stage5_run(tmp_path, t, truth(t), flags={0: ["NOT_INITIALISED"]}, valid=[False] + [True] * 5)
    _, res, summ, _ = run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path)), "--no-plots")
    assert not res["valid"].iloc[0] and "INVALID_STAGE5_STATE" in res["quality_flags"].iloc[0]
    assert "NOT_INITIALISED" in res["stage5_rejection_reason"].iloc[0]
    assert summ["n_compared"] == 5


def test_ensemble_small_flag(tmp_path):
    from tart_gnss.stage6 import synthetic as syn

    m, _ = _fresh(tmp_path, 4)
    ens = syn.write_tle_file(tmp_path / "ens.tle", [(L1, L2), syn.perturbed_tle(L1, L2, d_epoch_days=1)])
    _, res, summ, _ = run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path)), "--no-plots",
                              "--tle-ensemble-files", str(ens))
    assert summ["ensemble"]["n_distinct_in_window"] == 2
    assert all({"TLE_ENSEMBLE_SMALL", "REFERENCE_UNCERTAINTY_UNKNOWN"} <= set(f) for f in res["quality_flags"])
    assert res["position_mahalanobis_total_cov"].isna().all()


def test_prog_name():
    assert cli.build_parser().prog == "tart-validate"
