"""``tart-diagnose`` end to end: exit codes, manifest, no physical-cause labels, inputs untouched."""

import json

import pytest

pytest.importorskip("pandas")

from tart_gnss.stage4.loader import sha256_file  # noqa: E402
from tart_gnss.stage7 import cli  # noqa: E402

FORBIDDEN = ("SATELLITE_MANOEUVRE", "CLOCK_FAILURE", "RFI_EVENT", "HARDWARE_FAILURE")


def test_manifest_and_summary(runs):
    r = runs("network")
    m, s = r.manifest, r.summary
    assert m["stage"] == 7 and m["schema_version"] == 1
    assert m["alignment"]["row_number_joins"] is False and m["alignment"]["keys"] == ["time_utc", "station_id",
                                                                                      "norad_id"]
    assert m["evidence_rules"]["anomaly_threshold_used"] > 14 and m["evidence_rules"]["min_persistent_epochs"] == 5
    assert {lk["status"] for lk in m["provenance_chain"]} <= {"input", "verified"}
    assert len(m["input_hashes"]["stage5_manifest"]) == 64
    for key in ("overall_evidence_level", "overall_diagnostic_tags", "persistent_anomaly_intervals",
                "most_influential_station", "station_summaries", "strongest_calibration_correlation",
                "strongest_timing_correlation", "strongest_geometry_correlation",
                "strongest_trajectory_model_correlation", "reference_quality", "precision_validation",
                "unexplained_intervals"):
        assert key in s
    st = s["station_summaries"]["mauritius"]
    for key in ("usable_fraction", "median_direction_residual_deg", "median_line_residual_m",
                "median_leave_one_out_shift_m", "median_leave_one_out_mahalanobis", "zero_gain_fraction",
                "median_flagged_fraction", "calibration_change_count", "median_snr",
                "median_image_visibility_offset_deg", "timing_offset_median_s", "timing_offset_mad_s",
                "evidence_tags", "health_flags"):
        assert key in st


def test_no_physical_cause_labels(runs):
    for name in ("network", "bias", "dynamics"):
        text = (runs(name).dir / "summary.json").read_text() + (runs(name).dir / "diagnosis.csv").read_text()
        assert not any(f in text for f in FORBIDDEN)


def test_inputs_untouched_and_plots(runs, tmp_path):
    from tart_gnss.stage7.plotting import PLOT_NAMES

    pytest.importorskip("matplotlib")
    chain = runs("bias").chain
    files = [p for k in ("stage3", "stage4", "stage5", "stage6") for p in chain[k].parent.iterdir() if p.is_file()]
    before = {p: sha256_file(p) for p in files}
    assert cli.main(["--stage6-manifest", str(chain["stage6"]), "--output-dir", str(tmp_path / "o")]) == 0
    assert before == {p: sha256_file(p) for p in files}
    m = json.loads(next((tmp_path / "o").glob("*/*/stage7_manifest.json")).read_text())
    assert sorted(p.split("/")[-1][:-4] for p in m["outputs"]["diagnostics"]) == sorted(PLOT_NAMES)


def test_invalid_upstream_chain_exit_3(runs, tmp_path):
    # Strict policy with every Stage-6 row invalid: no diagnosis is manufactured.
    from tart_gnss.stage7.synthetic import Scenario, build_chain

    chain = build_chain(Scenario("insufficient", n_epochs=2), tmp_path / "c")
    code = cli.main(["--stage6-manifest", str(chain["stage6"]), "--output-dir", str(tmp_path / "o")])
    assert code == cli.EXIT_INSUFFICIENT


def test_bad_arguments(runs, tmp_path):
    m = str(runs("clean").chain["stage6"])
    assert cli.main(["--stage6-manifest", m, "--min-samples", "1"]) == cli.EXIT_USAGE
    assert cli.main(["--stage6-manifest", m, "--correlation-method", "magic"]) == cli.EXIT_USAGE
    assert cli.main(["--stage6-manifest", str(tmp_path / "none.json")]) == cli.EXIT_USAGE
    assert cli.main(["--stage6-manifest", m, "--output-dir", str(tmp_path / "x"), "--no-leave-one-out",
                     "--no-pairwise", "--no-plots"]) == cli.EXIT_OK


def test_prog_name():
    assert cli.build_parser().prog == "tart-diagnose"
