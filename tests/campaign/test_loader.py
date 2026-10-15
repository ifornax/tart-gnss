"""Loading completed Stage-7 runs."""

import json
import shutil

import pytest

from tart_gnss.campaign.errors import CampaignSchemaError, InvalidStage7ManifestError
from tart_gnss.campaign.loader import load_campaign_run


def test_load(pool):
    r = load_campaign_run(pool["gal_repeat_0"])
    assert r.norad_id == 90001 and r.constellation == "Galileo" and r.satellite["satellite_name"] == "SYN-GAL"
    assert not r.broken and set(r.stage_manifests) >= {2, 3, 4, 5, 6}
    assert r.residuals is not None and r.trajectory is not None and r.contributions is not None


def test_schema_and_type_errors(pool, tmp_path):
    src = pool["gps_clean_d0"]
    doc = json.loads(src.read_text())
    bad = tmp_path / "s7"
    shutil.copytree(src.parent, bad)
    (bad / "stage7_manifest.json").write_text(json.dumps({**doc, "schema_version": 2}))
    with pytest.raises(CampaignSchemaError):
        load_campaign_run(bad / "stage7_manifest.json")
    (bad / "stage7_manifest.json").write_text(json.dumps({**doc, "stage": 6}))
    with pytest.raises(InvalidStage7ManifestError):
        load_campaign_run(bad / "stage7_manifest.json")


def test_broken_provenance(pool, tmp_path):
    src = pool["gps_clean_d1"]
    work = tmp_path / "copy"
    shutil.copytree(src.parents[3], work)  # one run: stage2 ... stage7
    m7 = next(work.glob("stage7/*/*/stage7_manifest.json"))
    doc = json.loads(m7.read_text())
    for link in doc["provenance_chain"]:
        if link["stage"] == 5:
            p = work / "stage5" / next((work / "stage5").glob("*/*/stage5_manifest.json")).relative_to(work / "stage5")
            link["path"] = str(p)
            p.write_text(p.read_text() + "\n")  # changed after Stage 7 hashed it
    m7.write_text(json.dumps(doc))
    r = load_campaign_run(m7)
    assert r.broken and any("stage 5" in x for x in r.provenance_problems)
