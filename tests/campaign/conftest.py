"""Session-wide synthetic campaign: the default 17-run pool, built once through Stages 3-7."""

import json

import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")


@pytest.fixture(scope="session")
def pool(tmp_path_factory):
    """``name -> stage7_manifest.json`` for the default synthetic pool."""
    from tart_gnss.campaign.synthetic import build_runs, default_pool

    return build_runs(tmp_path_factory.mktemp("campaign_pool"), default_pool())


@pytest.fixture(scope="session")
def campaign(pool, tmp_path_factory):
    """The campaign analysis of the whole pool: ``(exit_code, out_dir, summary, manifest)``."""
    from tart_gnss.campaign import cli

    out = tmp_path_factory.mktemp("campaign_out")
    root = next(iter(pool.values())).parents[4]
    code = cli.main(["--input-root", str(root), "--output-dir", str(out), "--campaign-name", "syn",
                     "--bootstrap", "200"])
    d = out / "syn"
    return code, d, json.loads((d / "summary.json").read_text()), json.loads((d / "campaign_manifest.json").read_text())
