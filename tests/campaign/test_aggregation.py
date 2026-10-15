"""Run, epoch and station tables from completed products."""

import numpy as np

from tart_gnss.campaign.aggregation import build_campaign_tables, run_id_for
from tart_gnss.campaign.loader import load_campaign_run


def test_tables(pool):
    runs = [load_campaign_run(pool[k]) for k in ("gps_clean_d0", "gal_repeat_0", "gps_wide_d9")]
    ids = [run_id_for(r, str(r.diagnosis["time_utc"].min())) for r in runs]
    assert ids[0] == "20261007_NORAD40534_021211"
    t = build_campaign_tables(runs, ids)
    r = t.runs.set_index("run_id")
    assert list(r["n_valid_epochs"]) == [24, 24, 24] and (r["usable_fraction"] == 1).all()
    assert r.loc[ids[2], "n_stations"] == 4 and r.loc[ids[1], "constellation"] == "Galileo"
    assert (r["median_sigma_major_m"] > r["median_sigma_intermediate_m"]).all()
    assert (r["median_sigma_intermediate_m"] >= r["median_sigma_minor_m"]).all()
    assert np.all(np.isfinite(r["median_velocity_residual_m_s"])) and np.all(np.isfinite(r["median_elevation_deg"]))
    assert np.all(np.isfinite(r["stage5_rmse_vs_stage4_m"])) and (r["cohort"] == "precision_validation").all()
    assert len(t.epochs) == 72 and t.epochs["velocity_residual_3d_m_s"].notna().all()
    assert t.epochs["stage4_sigma_3d_m"].notna().all() and t.epochs["stage5_sigma_3d_m"].notna().all()
    assert len(t.stations) == 3 + 3 + 4 and t.stations["median_elevation_deg"].notna().all()


def test_gcrs_components(pool):
    run = load_campaign_run(pool["gps_clean_d0"])
    a = build_campaign_tables([run], ["x"], "rtn").runs.iloc[0]
    b = build_campaign_tables([run], ["x"], "gcrs").runs.iloc[0]
    assert a["median_radial_m"] != b["median_radial_m"]
