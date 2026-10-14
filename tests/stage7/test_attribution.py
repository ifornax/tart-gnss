"""Evidence rules on hand-built inputs."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage7.attribution import EvidenceRules, assign_evidence, loo_limit, station_pick  # noqa: E402


def _epochs(n=10, d2=1.0, d4=1.0, flags6=(), flags5=()):
    return pd.DataFrame({
        "time_s": np.arange(n) * 10.0, "time_utc": pd.date_range("2026-01-01", periods=n, freq="10s", tz="UTC"),
        "diagnosed": True, "position_mahalanobis_obs_cov": d2, "position_mahalanobis_total_cov": np.nan,
        "stage4_d2_tle": d4, "stage6_flags": [list(flags6)] * n, "stage5_flags": [list(flags5)] * n,
        "stage4_flags": [[] for _ in range(n)],
    })


def _infl(n, d=1.0, stations=("a", "b", "c")):
    return pd.DataFrame([{"epoch_index": e, "station_id": s, "loo_mahalanobis": d, "loo_dof": 2,
                          "loo_d2_tle": 100.0, "full_d2_tle": 100.0} for e in range(n) for s in stations])


def _pairs(n, z=(0.5, 0.5, 0.5)):
    return pd.DataFrame([{"epoch_index": e, "pair_station_1": a, "pair_station_2": b, "pair_separation_z": zz}
                         for e in range(n) for (a, b), zz in zip((("a", "b"), ("a", "c"), ("b", "c")), z)])


def _st(n):
    return pd.DataFrame([{"epoch_index": e, "station_id": s, "time_mismatch": False} for e in range(n)
                         for s in "abc"])


def _run(ep, infl, pairs, **kw):
    args = dict(corr_tags={}, stage5_run_flags=[], stage5_resets=0, precision_validation=True,
                calibration_available=True, data_complete=True, rules=EvidenceRules(),
                reference_uncertainty_known=True)
    args.update(kw)
    return assign_evidence(ep, _st(len(ep)), infl, pairs, **args)


def test_normal():
    rows, info = _run(_epochs(), _infl(10), _pairs(10))
    assert all(r["evidence_level"] == 0 and r["diagnostic_tags"] == ["NORMAL"] for r in rows)
    assert all(r["diagnostic_confidence"] == "HIGH" for r in rows)


def test_network_consistent_to_level_5_needs_persistence():
    rows, _ = _run(_epochs(d2=1e4, d4=500.0), _infl(10), _pairs(10))
    assert all(r["evidence_level"] == 5 and "UNEXPLAINED" in r["diagnostic_tags"] for r in rows)
    rows, _ = _run(_epochs(d2=1e4, d4=1.0), _infl(10), _pairs(10))  # not independently strong per epoch
    assert all(r["evidence_level"] == 4 for r in rows)


def test_caps_and_blocks():
    ep = _epochs(d2=1e4, d4=500.0)
    assert {r["evidence_level"] for r in _run(ep, _infl(10), _pairs(10), precision_validation=False)[0]} == {1}
    assert {r["evidence_level"] for r in _run(ep, _infl(10), _pairs(10), calibration_available=False)[0]} == {4}
    rows = _run(ep, _infl(10), _pairs(10), reference_uncertainty_known=False)[0]
    assert {r["evidence_level"] for r in rows} == {4} and all(r["diagnostic_confidence"] != "HIGH" for r in rows)
    rows, _ = _run(ep, _infl(10), _pairs(10), stage5_run_flags=["DYNAMICS_MISMATCH"])
    assert all("TRAJECTORY_MODEL_CORRELATED" in r["diagnostic_tags"] and r["evidence_level"] == 1 for r in rows)
    rows, _ = _run(ep, _infl(10, d=50.0), _pairs(10))  # a significant leave-one-out shift
    assert all("NETWORK_CONSISTENT" not in r["diagnostic_tags"] for r in rows)


def test_station_pick_methods():
    r = EvidenceRules()
    pairs = _pairs(1, z=(9.0, 8.0, 0.4))  # a-b and a-c inconsistent
    assert station_pick(pairs, None, False, r) == ("a", "pairwise", False)
    pairs1 = _pairs(1, z=(9.0, 0.4, 0.4))  # only a-b: ambiguous
    assert station_pick(pairs1, None, False, r)[0] == "AMBIGUOUS"
    loo = pd.DataFrame({"station_id": ["a", "b", "c"], "loo_d2_tle": [2.0, 300.0, 400.0],
                        "full_d2_tle": [300.0] * 3})
    assert station_pick(pairs1, loo, True, r) == ("a", "reference", False)
    assert station_pick(pairs, loo, True, r) == ("a", "pairwise+reference", True)


def test_station_specific_level_2():
    ep = _epochs(d2=1e4, d4=500.0)
    rows, info = _run(ep, _infl(10, d=80.0), _pairs(10, z=(9.0, 8.0, 0.4)))
    assert info["station_specific"] == ["a"]
    assert all("STATION_SPECIFIC" in r["diagnostic_tags"] and r["evidence_level"] == 2 for r in rows)
    assert all("LOCALISATION_CORRELATED" not in r["diagnostic_tags"] for r in rows)  # not counted twice
    rows, _ = _run(ep, _infl(10, d=80.0), _pairs(10, z=(9.0, 0.4, 0.4)))  # inconsistency, no unique station
    assert all("LOCALISATION_CORRELATED" in r["diagnostic_tags"] for r in rows)


def test_threshold_and_limits():
    r = EvidenceRules()
    assert r.threshold_for(1) == pytest.approx(14.156, abs=0.01)
    assert r.threshold_for(100) > r.threshold_for(10)
    assert EvidenceRules(anomaly_threshold=7.0).threshold_for(50) == 7.0
    assert loo_limit(2) == pytest.approx(11.83, abs=0.01) and loo_limit(0) == np.inf
