"""Run discovery."""

import pytest

from tart_gnss.campaign.discovery import collect_manifests, discover_stage7_runs
from tart_gnss.campaign.errors import ManifestDiscoveryError


def test_discovery(tmp_path):
    (tmp_path / "a" / "b").mkdir(parents=True)
    (tmp_path / "a" / "b" / "stage7_manifest.json").write_text("{}")
    (tmp_path / "a" / "positions.parquet").write_text("x")
    (tmp_path / "a" / "stage6_manifest.json").write_text("{}")
    found = discover_stage7_runs(tmp_path)
    assert [p.name for p in found] == ["stage7_manifest.json"]
    both = collect_manifests([tmp_path], [found[0]])
    assert len(both) == 1  # duplicates listed once


def test_errors(tmp_path):
    with pytest.raises(ManifestDiscoveryError):
        discover_stage7_runs(tmp_path / "none")
    with pytest.raises(ManifestDiscoveryError):
        collect_manifests([tmp_path], [])
    with pytest.raises(ManifestDiscoveryError):
        collect_manifests([], [tmp_path / "missing.json"])
