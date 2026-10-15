"""Evidence-level and tag statistics with separate denominators."""

import pandas as pd

from tart_gnss.campaign.anomaly_metrics import evidence_summary, level_fractions


def test_denominators_are_separate():
    runs = pd.DataFrame({"run_id": ["a", "b"], "overall_evidence_level": [0, 4],
                         "overall_diagnostic_tags": [["NORMAL"], ["STATISTICAL_ANOMALY", "NETWORK_CONSISTENT"]]})
    epochs = pd.DataFrame({"evidence_level": [0, 0, 0, 4], "diagnostic_tags": [["NORMAL"]] * 3 + [
        ["STATISTICAL_ANOMALY"]]})
    stations = pd.DataFrame({"overall_evidence_level": [0, 0, 4], "station_tags": [[], ["STATION_SPECIFIC"], []]})
    e = evidence_summary(runs, epochs, stations)
    get = lambda u, q: e[(e["unit"] == u) & (e["quantity"] == q)].iloc[0]  # noqa: E731
    assert get("run", "evidence_level_4")["fraction"] == 0.5 and get("epoch", "evidence_level_4")["fraction"] == 0.25
    assert get("station", "STATION_SPECIFIC")["n"] == 3
    r = get("run", "NETWORK_CONSISTENT")
    assert 0 < r["ci95_low"] < 0.5 < r["ci95_high"] < 1
    assert pd.isna(get("epoch", "STATISTICAL_ANOMALY")["ci95_low"])


def test_level_fractions():
    f = level_fractions([0, 0, 5, None])
    assert f["n"] == 3 and f["level_0"] == 2 / 3 and f["level_5"] == 1 / 3


def test_campaign_rates(campaign):
    _, d, _, _ = campaign
    r = pd.read_csv(d / "tables" / "anomaly_rates.csv")
    assert {"constellation", "station", "abs_tle_age_days"} <= set(r["condition"])
    ref = r[(r["condition"] == "abs_tle_age_days") & (r["tag"] == "REFERENCE_LIMITED")]
    assert ref.loc[ref["value"] == ">14", "fraction"].iloc[0] == 1.0
