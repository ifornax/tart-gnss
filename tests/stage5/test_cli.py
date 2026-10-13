"""``tart-track`` end to end: exit codes, outputs, manifest, flags."""

import json

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")

from tart_gnss.stage5 import cli  # noqa: E402
from tart_gnss.stage5 import synthetic as syn  # noqa: E402
from tart_gnss.stage5.output import state_covariance  # noqa: E402

from .helpers import run_cli, times, truth  # noqa: E402


def _rows(n=20, seed=0, **kw):
    t = times(n)
    return syn.stage4_rows(t, truth(t), 2e4, 2e3, np.random.default_rng(seed), **kw)


def test_successful_run_and_manifest(tmp_path):
    code, traj, inn, man = run_cli(tmp_path, _rows(), "--no-plots")
    assert code == cli.EXIT_OK
    assert len(traj) == 20 and len(inn) == 20
    assert traj["valid"].all()
    assert man["stage"] == 5 and man["state_frame"] == "GCRS" and man["measurement_frame"] == "ITRS/ECEF"
    assert man["tle_used"] is False and man["filter"]["initialisation"]["tle_used"] is False
    assert man["state_model"]["dynamics"] == "two-body"
    assert man["state_model"]["rtol"] == 1e-10
    assert man["process_noise"]["process_accel_sigma_m_s2"] == 1e-4
    assert man["innovation_gate"]["threshold"] == pytest.approx(14.156, abs=1e-3)
    assert len(man["input_sha256"]) == 64
    assert man["software"]["git_commit"] is None or len(man["software"]["git_commit"]) == 40
    assert man["provenance"]["command_line"][0] == "tart-track"
    p = state_covariance(traj.iloc[5], "smoothed")
    assert np.allclose(p, p.T) and np.linalg.eigvalsh(p).min() > 0
    # ECEF output is the GCRS state rotated back (radius preserved).
    r_g = np.linalg.norm(traj[["x_gcrs_smoothed_m", "y_gcrs_smoothed_m", "z_gcrs_smoothed_m"]].to_numpy(), axis=1)
    r_e = np.linalg.norm(traj[["x_ecef_smoothed_m", "y_ecef_smoothed_m", "z_ecef_smoothed_m"]].to_numpy(), axis=1)
    assert np.allclose(r_g, r_e, rtol=1e-12)
    assert np.allclose(traj["acceleration_m_s2"], 3.986004418e14 / r_g ** 2)


def test_plots_written(tmp_path):
    from tart_gnss.stage5.plotting import PLOT_NAMES

    pytest.importorskip("matplotlib")
    code, _, _, man = run_cli(tmp_path, _rows(12))
    assert code == cli.EXIT_OK
    names = sorted(p.split("/")[-1][:-4] for p in man["outputs"]["diagnostics"])
    assert names == sorted(PLOT_NAMES)


def test_insufficient_measurements_exit_3(tmp_path):
    rows = _rows(5, flags={0: ["INSUFFICIENT_STATIONS"], 1: ["INSUFFICIENT_STATIONS"], 2: ["NEGATIVE_RANGE"]})
    code, traj, _, man = run_cli(tmp_path, rows, "--no-plots")
    assert code == cli.EXIT_NO_VALID
    assert all("INSUFFICIENT_MEASUREMENTS" in f for f in traj["quality_flags"])
    assert not traj["valid"].any()
    assert man["filter"]["success"] is False


def test_leading_rejected_rows_are_not_initialised(tmp_path):
    rows = _rows(10, flags={0: ["POOR_GEOMETRY"], 1: ["POOR_GEOMETRY"]})
    code, traj, inn, _ = run_cli(tmp_path, rows, "--no-plots")
    assert code == cli.EXIT_OK
    assert traj["quality_flags"].iloc[0][:2] == ["NOT_INITIALISED", "MEASUREMENT_REJECTED"]
    assert not traj["valid"].iloc[:2].any() and traj["valid"].iloc[2:].all()
    assert len(inn) == 8


def test_no_smoother(tmp_path):
    code, traj, _, man = run_cli(tmp_path, _rows(8), "--no-plots", "--no-smoother")
    assert code == cli.EXIT_OK
    assert traj["x_gcrs_smoothed_m"].isna().all() and traj["x_gcrs_filtered_m"].notna().all()
    assert man["smoother"]["enabled"] is False


def test_constant_velocity_option(tmp_path):
    code, traj, _, man = run_cli(tmp_path, _rows(8), "--no-plots", "--dynamics", "constant-velocity",
                                 "--process-accel-sigma", "0.5")
    assert code == cli.EXIT_OK
    assert man["state_model"]["integrator"] == "closed form"
    assert man["process_noise"]["process_accel_sigma_m_s2"] == 0.5
    assert (traj["acceleration_m_s2"] == 0).all()


@pytest.mark.parametrize("extra", [["--innovation-gate-probability", "1.5"], ["--max-gap-s", "0"],
                                   ["--dynamics", "kepler"], ["--process-accel-sigma", "-1"]])
def test_bad_arguments_exit_2(tmp_path, extra):
    m = syn.write_stage4_run(tmp_path / "s4", _rows(5))
    assert cli.main(["--stage4-manifest", str(m), "--output-dir", str(tmp_path / "o"), *extra]) == cli.EXIT_USAGE


def test_bad_manifest_exit_2(tmp_path):
    p = tmp_path / "m.json"
    p.write_text(json.dumps({"stage": 3, "schema_version": 1}))
    assert cli.main(["--stage4-manifest", str(p), "--output-dir", str(tmp_path / "o")]) == cli.EXIT_USAGE
    assert cli.main(["--stage4-manifest", str(tmp_path / "none.json")]) == cli.EXIT_USAGE


def test_prog_name():
    assert cli.build_parser().prog == "tart-track"
