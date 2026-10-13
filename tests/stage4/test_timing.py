"""Time matching of independent station timestamps and timing diagnostics."""

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage4.errors import TimeMatchError  # noqa: E402
from tart_gnss.stage4.timing import (  # noqa: E402
    angular_rates_deg_s,
    choose_reference_station,
    match_station_epochs,
    motion_during_offset_deg,
    time_seconds,
)

T0 = pd.Timestamp("2026-06-25T09:14:00Z")


def table(offsets_s, shuffle=None):
    t = [T0 + pd.Timedelta(seconds=float(o)) for o in offsets_s]
    df = pd.DataFrame({"time_utc": t, "v": np.arange(len(t))})
    if shuffle is not None:
        df = df.iloc[shuffle].reset_index(drop=True)
    return df


def test_matches_nearest_not_row_number():
    a = table([0, 1, 2, 3, 4, 5])
    # b is shuffled and starts later: row numbers do not line up with a.
    b = table([2.3, 0.2, 4.1, 1.4, 3.3], shuffle=None)
    eps = match_station_epochs({"a": a, "b": b}, tolerance_s=0.75, reference_station="a")
    got = {round(e.reference_time_s - time_seconds([T0])[0]): (b.loc[e.rows["b"], "v"] if "b" in e.rows else None)
           for e in eps}
    assert got == {0: 1, 1: 3, 2: 0, 3: 4, 4: 2, 5: None}
    off = {round(e.reference_time_s - time_seconds([T0])[0]): e.offsets_s.get("b") for e in eps}
    assert off[1] == pytest.approx(0.4) and off[3] == pytest.approx(0.3)


def test_tolerance_excludes_and_each_sample_used_once():
    a = table([0.0, 0.5, 1.0])
    b = table([0.45])
    eps = match_station_epochs({"a": a, "b": b}, tolerance_s=0.75, reference_station="a")
    users = [e for e in eps if "b" in e.rows]
    assert len(users) == 1 and users[0].offsets_s["b"] == pytest.approx(-0.05)
    eps = match_station_epochs({"a": table([0.0]), "b": table([1.0])}, tolerance_s=0.75, reference_station="a")
    assert "b" not in eps[0].rows


def test_reference_station_and_time():
    tabs = {"x": table([0, 1]), "y": table([0, 1, 2]), "z": table([0, 1, 2])}
    assert choose_reference_station(tabs) == "y"
    eps = match_station_epochs(tabs, 0.5, reference_time_s=time_seconds([T0])[0] + 1.2)
    assert len(eps) == 1 and eps[0].reference_station == "y"
    assert sorted(eps[0].rows) == ["x", "y", "z"]
    with pytest.raises(TimeMatchError):
        match_station_epochs(tabs, -1.0)
    with pytest.raises(TimeMatchError):
        match_station_epochs(tabs, 0.5, reference_station="nope")


def test_max_offset_is_spread():
    eps = match_station_epochs({"a": table([0.0]), "b": table([0.3]), "c": table([-0.4])}, 0.75,
                               reference_station="a")
    assert eps[0].max_offset_s == pytest.approx(0.7)


def test_angular_rate_of_uniform_motion():
    w = np.radians(0.01)  # 0.01 deg/s
    t = np.arange(0.0, 200.0, 1.0)
    s = np.stack([np.cos(w * t), np.sin(w * t), np.zeros_like(t)], axis=1)
    rate = angular_rates_deg_s(t, s)
    assert np.allclose(rate, 0.01, rtol=1e-4)
    assert motion_during_offset_deg(rate[50], -0.5) == pytest.approx(0.005, rel=1e-4)


def test_angular_rate_is_robust_to_measurement_noise():
    rng = np.random.default_rng(0)
    w = np.radians(0.01)
    t = np.arange(0.0, 300.0, 1.0)
    s = np.stack([np.cos(w * t), np.sin(w * t), np.zeros_like(t)], axis=1)
    s += rng.normal(0.0, np.radians(0.1), s.shape)
    s /= np.linalg.norm(s, axis=1)[:, None]
    rate = angular_rates_deg_s(t, s)
    assert np.median(rate) == pytest.approx(0.01, rel=0.15)


def test_angular_rate_needs_samples():
    assert np.all(np.isnan(angular_rates_deg_s(np.array([0.0, 1.0]), np.eye(3)[:2])))
