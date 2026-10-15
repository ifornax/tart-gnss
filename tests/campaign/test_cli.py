"""``tart-campaign`` end to end."""

import json

import pandas as pd
import pytest

from tart_gnss.campaign import cli


def test_outputs_and_manifest(campaign):
    code, d, s, m = campaign
    for name in ("campaign_runs.parquet", "campaign_epochs.parquet", "campaign_stations.parquet", "excluded_runs.csv",
                 "summary.json", "campaign_manifest.json", "logs/campaign.log"):
        assert (d / name).is_file(), name
    for t in ("station_summary", "constellation_summary", "satellite_summary", "evidence_summary", "tle_age_summary",
              "geometry_summary"):
        assert (d / "tables" / f"{t}.csv").is_file()
    assert m["analysis"] == "campaign" and m["bootstrap"]["unit"] == "run" and len(m["input_hashes"]) == 17
    assert len(m["included_runs"]) == 17 and m["software"]["git_commit"] is not None
    from tart_gnss.campaign.plotting import PLOT_NAMES

    assert sorted(p.split("/")[-1][:-4] for p in m["outputs"]["diagnostics"]) == sorted(PLOT_NAMES)


def test_no_reprocessing(pool, tmp_path, monkeypatch):
    import tart_gnss.stage3.cli as c3
    import tart_gnss.stage4.cli as c4
    import tart_gnss.stage5.cli as c5
    import tart_gnss.stage6.cli as c6
    import tart_gnss.stage7.cli as c7

    def boom(*a, **k):
        raise AssertionError("a pipeline stage was run")

    for mod in (c3, c4, c5, c6, c7):
        monkeypatch.setattr(mod, "main", boom)
        monkeypatch.setattr(mod, "run", boom)
    args = []
    for k in ("gps_clean_d0", "gps_clean_d1", "gps_clean_d2"):
        args += ["--stage7-manifest", str(pool[k])]
    assert cli.main([*args, "--output-dir", str(tmp_path), "--no-plots", "--bootstrap", "50"]) == cli.EXIT_OK


def test_insufficient_and_excluded(pool, tmp_path):
    code = cli.main(["--stage7-manifest", str(pool["gps_aged_d20"]), "--stage7-manifest", str(pool["gps_clean_d0"]),
                     "--output-dir", str(tmp_path), "--campaign-name", "small", "--no-plots", "--bootstrap", "20",
                     "--precision-validation-only"])
    assert code == cli.EXIT_INSUFFICIENT
    s = json.loads((tmp_path / "small" / "summary.json").read_text())
    assert "INSUFFICIENT_CAMPAIGN_SIZE" in s["campaign_quality_flags"] and s["n_included_runs"] == 1
    exc = pd.read_csv(tmp_path / "small" / "excluded_runs.csv")
    assert list(exc["reason"]) == ["NO_PRECISION_VALIDATION"]


def test_filters_and_publication_plots(pool, tmp_path):
    root = next(iter(pool.values())).parents[4]
    code = cli.main(["--input-root", str(root), "--output-dir", str(tmp_path), "--campaign-name", "gps",
                     "--constellation", "GPS", "--bootstrap", "50", "--publication-plots", "--no-precision-only",
                     "--group-by", "station_set", "norad_id"])
    assert code == cli.EXIT_OK
    s = json.loads((tmp_path / "gps" / "summary.json").read_text())
    assert s["constellations"] == ["GPS"] and s["analysis_cohort"] == "all_included"
    assert (tmp_path / "gps" / "diagnostics" / "residual_distribution.pdf").is_file()
    assert (tmp_path / "gps" / "tables" / "grouped_by_station_set.csv").is_file()


@pytest.mark.parametrize("argv", [[], ["--input-root", "/nonexistent"], ["--stage7-manifest", "x.json",
                                                                      "--bootstrap", "-1"]])
def test_invalid_input(argv, tmp_path):
    assert cli.main([*argv, "--output-dir", str(tmp_path)]) == cli.EXIT_USAGE


def test_prog_name():
    assert cli.build_parser().prog == "tart-campaign"
