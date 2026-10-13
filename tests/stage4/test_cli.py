"""``tart-triangulate`` end to end: exit codes, outputs, manifest, policy, TLE independence."""

import json

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage4 import cli  # noqa: E402
from tart_gnss.stage4 import synthetic as syn  # noqa: E402

from .helpers import run_cli  # noqa: E402

TRUTH = syn.linear_track()
T = {sid: [syn.T0_UNIX_S + k + 0.1 * i for k in range(10)] for i, sid in enumerate(syn.EXAMPLE_STATIONS)}


def rows(rng=None):
    return syn.make_directions(TRUTH, T, 0.1, rng)


def test_outputs_and_manifest(tmp_path):
    m = syn.write_stage3_run(tmp_path / "s3", rows(np.random.default_rng(0)))
    out = tmp_path / "s4"
    code = cli.main(["--stage3-manifest", str(m), "--output-dir", str(out), "--bootstrap", "30", "--seed", "5"])
    assert code == 0
    d = out / "20260625" / "NORAD_40534"
    for name in ("positions.csv", "station_contributions.csv", "stage4_manifest.json", "logs/stage4.log"):
        assert (d / name).is_file(), name
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pass
    else:
        assert (d / "positions.parquet").is_file() and (d / "station_contributions.parquet").is_file()
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        pass
    else:
        for p in ("position_track_ecef", "geometry_condition", "station_residuals", "position_uncertainty",
                  "timing_offsets"):
            assert (d / "diagnostics" / f"{p}.png").is_file()
    man = json.loads((d / "stage4_manifest.json").read_text())
    for key in ("schema_version", "stage", "input_stage3_manifest", "input_sha256", "norad_id", "coordinate_frame",
                "time_matching", "triangulation", "uncertainty", "quality_thresholds", "stations", "software",
                "provenance", "outputs"):
        assert key in man, key
    assert man["stage"] == 4 and man["coordinate_frame"] == "ITRS/ECEF"
    assert man["uncertainty"] == {**man["uncertainty"], "bootstrap": 30, "random_seed": 5}
    assert man["time_matching"]["tolerance_s"] == 0.75
    assert man["quality_policy"]["policy"] == "strict"
    assert man["triangulation"]["tle_used"] is False and man["triangulation"]["method"] == "weighted"
    assert man["provenance"]["command_line"][0] == "tart-triangulate"
    assert "git_commit" in man["software"] and man["software"]["tart_gnss"]
    st = man["stations"]["mauritius"]
    assert st["latitude_deg"] == syn.EXAMPLE_STATIONS["mauritius"][0] and "stations.mauritius" in st["source"]
    log = (d / "logs" / "stage4.log").read_text()
    for text in ("input Stage-3 manifest", "matched epochs", "median line residual", "outputs:"):
        assert text in log
    pos = pd.read_csv(d / "positions.csv")
    assert pos["x_p16"].notna().all() and (pos["bootstrap_n_success"] == 30).all()


def test_tle_never_used(tmp_path):
    """Identical positions with and without any TLE/satellite information or TLE columns."""
    r = rows(np.random.default_rng(1))
    for row in r:
        row.update({"az_tle_deg": 10.0, "el_tle_deg": 20.0, "tle_separation_deg": 30.0})
    _, a, _, man_a = run_cli(tmp_path / "a", r)
    stripped = [{k: v for k, v in row.items() if "tle" not in k} for row in r]
    m = syn.write_stage3_run(tmp_path / "b" / "stage3", stripped)
    doc = json.loads(m.read_text())
    doc.pop("satellite")
    m.write_text(json.dumps(doc))
    cli.main(["--stage3-manifest", str(m), "--output-dir", str(tmp_path / "b" / "stage4"), "--no-plots"])
    b = pd.read_csv(next((tmp_path / "b" / "stage4").glob("*/*/positions.csv")))
    assert np.array_equal(a[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy(), b[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy())


def test_quality_policy_strict_vs_permissive(tmp_path):
    r = rows()
    for row in r:
        if row["station_id"] == "madagascar":
            row["quality_flags"] = ["TLE_UNASSOCIATED"]
    code, pos, _, man = run_cli(tmp_path / "s", r)
    assert code == 3 and (pos["n_stations"] == 2).all()
    assert man["quality_policy"]["rejected"] == {"TLE_UNASSOCIATED": 10}
    code, pos, _, man = run_cli(tmp_path / "p", r, "--quality-policy", "permissive")
    assert code == 0 and (pos["n_stations"] == 3).all()
    assert all("STAGE3_QUALITY_WARNING" in f for f in pos["quality_flags"])


def test_advisory_and_calibration_flags(tmp_path):
    r = rows()
    r[0]["quality_flags"] = ["IMAGE_DISAGREEMENT"]
    code, pos, con, _ = run_cli(tmp_path, r, calibration_warnings={"south_africa": ["9/24 antennas zero gain"]})
    assert code == 0 and pos["valid"].all()
    assert all("CALIBRATION_WARNING" in f for f in pos["quality_flags"])
    assert sum("INPUT_DIRECTION_WARNING" in f for f in pos["quality_flags"]) == 1


def test_station_selection_reference_time_and_unweighted(tmp_path):
    code, pos, _, man = run_cli(tmp_path / "a", rows(), "--stations", "mauritius", "madagascar",
                                "--min-stations", "2")
    assert code == 0 and set(man["stations"]) == {"mauritius", "madagascar"}
    code, pos, _, _ = run_cli(tmp_path / "b", rows(), "--reference-time", "2026-06-25T09:14:05.2Z")
    assert len(pos) == 1 and pos.loc[0, "time_utc"].second == 5
    code, pos, _, man = run_cli(tmp_path / "c", rows(), "--method", "unweighted")
    assert code == 0 and (pos["method"] == "unweighted").all() and pos["sigma_x_m"].notna().all()


@pytest.mark.parametrize("extra", [
    ["--time-match-tolerance", "-1"], ["--min-stations", "1"], ["--bootstrap", "-2"], ["--method", "magic"],
    ["--stations", "atlantis"], ["--reference-time", "not-a-time"],
])
def test_invalid_configuration_exit_2(tmp_path, extra):
    code, *_ = run_cli(tmp_path, rows(), *extra)
    assert code == 2


def test_missing_inputs_exit_2(tmp_path):
    assert cli.main(["--stage3-manifest", str(tmp_path / "nope.json")]) == 2
    m = syn.write_stage3_run(tmp_path, rows())
    for f in ("directions.csv", "directions.parquet"):
        (tmp_path / f).unlink(missing_ok=True)
    assert cli.main(["--stage3-manifest", str(m), "--output-dir", str(tmp_path / "o")]) == 2


def test_no_valid_estimates_exit_3(tmp_path):
    r = [row for row in rows() if row["station_id"] == "mauritius"]
    r += [dict(row, station_id="madagascar") for row in rows() if row["station_id"] == "madagascar"][:0]
    code, pos, _, man = run_cli(tmp_path, r)
    assert code == 3 and man["summary"]["n_valid"] == 0


def test_covariance_fallback(tmp_path):
    r = rows()
    for row in r:
        if row["station_id"] == "south_africa":
            row["sigma_l"] = row["sigma_m"] = row["cov_lm"] = row["sigma_az_deg"] = float("nan")
    code, pos, con, _ = run_cli(tmp_path, r)
    assert code == 3
    assert all("COVARIANCE_SINGULAR" in f for f in pos["quality_flags"])
    assert np.isfinite(pos["x_ecef_m"]).all()  # unweighted fallback still gives the point
