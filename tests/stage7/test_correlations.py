"""Correlation analyses and false-positive control (specification section 37)."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage7.correlations import compute_correlations, strongest, supported_tags  # noqa: E402


def _tables(rng, n=40, link=None):
    t = np.arange(n) * 10.0
    rows = []
    for sid in ("a", "b", "c"):
        ff = rng.random(n)
        resid = rng.normal(0, 0.01, n) + (0.05 * ff if sid == link else 0)
        for k in range(n):
            rows.append({"station_id": sid, "time_s": t[k], "flagged_fraction": ff[k], "snr": rng.normal(20, 2),
                         "direction_residual_angle_deg": abs(resid[k]), "line_residual_m": abs(rng.normal(0, 1e3))})
    ep = pd.DataFrame({"time_s": t, "stage4_residual_3d_m": abs(rng.normal(0, 1e4, n)),
                       "predicted_condition_number": 200 + rng.normal(0, 1, n)})
    return pd.DataFrame(rows), ep


def test_injected_calibration_link_is_found():
    rng = np.random.default_rng(0)
    st, ep = _tables(rng, link="b")
    c = compute_correlations(st, ep)
    tags = supported_tags(c)
    assert "CALIBRATION_CORRELATED" in tags
    assert {x["station_id"] for x in tags["CALIBRATION_CORRELATED"]} == {"b"}
    assert strongest(c, "calibration")["station_id"] == "b"


@pytest.mark.parametrize("seed", range(10))
def test_random_metadata_rarely_supports_a_tag(seed):
    rng = np.random.default_rng(100 + seed)
    st, ep = _tables(rng)
    # random-walk metadata (smooth, autocorrelated) as well as white noise
    st["flagged_fraction"] = st.groupby("station_id")["flagged_fraction"].transform(lambda v: np.cumsum(v - 0.5))
    assert not supported_tags(compute_correlations(st, ep))


def test_constant_and_small_samples():
    rng = np.random.default_rng(1)
    st, ep = _tables(rng, n=5)
    c = compute_correlations(st, ep)
    assert (c["interpretation_tag"].isin(["INSUFFICIENT_SAMPLES", "INSUFFICIENT_VARIATION", "EXPECTED"])).all()
    st, ep = _tables(rng)
    st["flagged_fraction"] = 0.2
    c = compute_correlations(st, ep)
    assert (c.loc[c["x_variable"] == "flagged_fraction", "interpretation_tag"] == "INSUFFICIENT_VARIATION").all()
    # geometry that barely changes over the pass is not a geometry test
    assert (c.loc[c["family"] == "geometry", "interpretation_tag"].isin(["INSUFFICIENT_VARIATION"])).all()


def test_binary_indicator_is_allowed():
    rng = np.random.default_rng(4)
    st, ep = _tables(rng)
    flag = rng.random(len(st)) < 0.4
    st["flagged_fraction"] = np.where(flag, 0.6, 0.1)
    st.loc[(st["station_id"] == "a") & flag, "direction_residual_angle_deg"] += 0.1
    assert "CALIBRATION_CORRELATED" in supported_tags(compute_correlations(st, ep))
