"""Time alignment across stages (never by row number)."""

import numpy as np
import pytest

pytest.importorskip("pandas")

from tart_gnss.stage7.alignment import build_epoch_table, build_station_table, nearest  # noqa: E402
from tart_gnss.stage7.loader import load_diagnostic_chain  # noqa: E402


def test_nearest():
    t = np.array([0.0, 10.0, 20.0])
    assert nearest(t, 10.004, 0.01) == 1 and nearest(t, 15.0, 0.01) is None and nearest(np.zeros(0), 1, 1) is None


def test_tables_link_every_stage(runs):
    inp = load_diagnostic_chain(runs("clean").chain["stage6"])
    ep = build_epoch_table(inp)
    assert ep["link_stage5"].all() and ep["link_stage4"].all()
    assert np.all(np.isfinite(np.stack(ep["tle_ecef"])))
    st = build_station_table(inp, ep)
    assert len(st) == 108 and st["link_stage3"].all() and st["link_stage2"].all()
    assert (st.groupby("epoch_index")["station_id"].nunique() == 3).all()


def test_shuffled_rows_give_the_same_alignment(runs):
    inp = load_diagnostic_chain(runs("clean").chain["stage6"])
    a = build_epoch_table(inp)
    inp.positions = inp.positions.sample(frac=1.0, random_state=1).reset_index(drop=True)
    inp.trajectory = inp.trajectory.sample(frac=1.0, random_state=2).reset_index(drop=True)
    b = build_epoch_table(inp)
    assert np.allclose(np.stack(a["stage4_position_ecef"]), np.stack(b["stage4_position_ecef"]))
    assert np.allclose(a["stage5_nis"], b["stage5_nis"], equal_nan=True)


def test_timing_offsets_are_carried(runs):
    inp = load_diagnostic_chain(runs("timing").chain["stage6"])
    st = build_station_table(inp, build_epoch_table(inp))
    mu = st[st["station_id"] == "mauritius"]
    assert mu["station_time_offset_s"].between(0.9, 5.1).all() and mu["link_stage3"].all()
