"""Loading Stage-5 runs."""

import json

import numpy as np
import pytest

pytest.importorskip("pandas")

from tart_gnss.stage6 import synthetic as syn  # noqa: E402
from tart_gnss.stage6.errors import InvalidStage5ManifestError, MissingTrajectoryError  # noqa: E402
from tart_gnss.stage6.loader import load_stage5_manifest, load_stage5_trajectory, row_state  # noqa: E402

from .helpers import times, truth  # noqa: E402


def test_manifest_validation(tmp_path):
    with pytest.raises(InvalidStage5ManifestError):
        load_stage5_manifest(tmp_path / "none.json")
    p = tmp_path / "m.json"
    for doc in ({"stage": 4, "schema_version": 1}, {"stage": 5, "schema_version": 2},
                {"stage": 5, "schema_version": 1, "state_frame": "ITRS"}):
        p.write_text(json.dumps(doc))
        with pytest.raises(InvalidStage5ManifestError):
            load_stage5_manifest(p)


def test_full_covariance_rebuilt(tmp_path):
    t = times(4)
    rng = np.random.default_rng(0)
    a = rng.normal(size=(6, 6))
    p = a @ a.T * 1e4
    m = syn.write_stage5_run(tmp_path, t, truth(t), np.repeat(p[None], 4, axis=0))
    s5 = load_stage5_trajectory(m)
    x, pc = row_state(s5.table.iloc[2])
    assert np.allclose(pc, p, rtol=1e-12) and pc[0, 4] != 0  # off-diagonal terms kept
    assert np.allclose(x, truth(t)[2])
    assert s5.norad_id == 40534 and len(s5.sha256) == 64 and np.allclose(s5.times_s, t, atol=1e-5)


def test_moved_manifest_and_missing_columns(tmp_path):
    t = times(3)
    m = syn.write_stage5_run(tmp_path, t, truth(t), np.repeat(np.eye(6)[None], 3, axis=0))
    doc = json.loads(m.read_text())
    doc["outputs"] = {"trajectory_csv": "/nowhere/trajectory.csv"}
    m.write_text(json.dumps(doc))
    assert load_stage5_trajectory(m).trajectory_path.parent == m.parent
    import pandas as pd

    (m.parent / "trajectory.parquet").unlink(missing_ok=True)
    csv = m.parent / "trajectory.csv"
    pd.read_csv(csv).drop(columns=["cov_smoothed_x_vy"]).to_csv(csv, index=False)
    with pytest.raises(MissingTrajectoryError):
        load_stage5_trajectory(m)
