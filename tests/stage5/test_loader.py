"""Loading Stage-4 runs and building GCRS measurements."""

import json

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")

from tart_gnss.stage5 import synthetic as syn  # noqa: E402
from tart_gnss.stage5.errors import InvalidStage4ManifestError, MissingPositionDataError  # noqa: E402
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices  # noqa: E402
from tart_gnss.stage5.loader import (  # noqa: E402
    build_measurements,
    covariance_matrix,
    load_stage4_manifest,
    load_stage4_positions,
)

from .helpers import frame, times, truth  # noqa: E402


def _rows(n=6, **kw):
    t = times(n)
    return syn.stage4_rows(t, truth(t), 3e4, 2e3, np.random.default_rng(0), **kw)


def test_manifest_validation(tmp_path):
    with pytest.raises(InvalidStage4ManifestError):
        load_stage4_manifest(tmp_path / "missing.json")
    p = tmp_path / "m.json"
    p.write_text("{not json")
    with pytest.raises(InvalidStage4ManifestError):
        load_stage4_manifest(p)
    for doc in ({"stage": 3, "schema_version": 1}, {"stage": 4, "schema_version": 2},
                {"stage": 4, "schema_version": 1, "coordinate_frame": "GCRS"}):
        p.write_text(json.dumps(doc))
        with pytest.raises(InvalidStage4ManifestError):
            load_stage4_manifest(p)


def test_positions_found_next_to_moved_manifest(tmp_path):
    m = syn.write_stage4_run(tmp_path / "a", _rows())
    doc = json.loads(m.read_text())
    doc["outputs"] = {"positions_csv": "/nonexistent/positions.csv"}
    syn.write_json(m, doc)
    inp = load_stage4_positions(m)
    assert inp.positions_path.parent == m.parent
    assert len(inp.positions) == 6 and inp.norad_id == 40534
    assert len(inp.sha256) == 64


def test_csv_only_and_missing_table(tmp_path):
    m = syn.write_stage4_run(tmp_path / "a", _rows(), csv_only=True)
    assert load_stage4_positions(m).positions_path.suffix == ".csv"
    (m.parent / "positions.csv").unlink()
    with pytest.raises(MissingPositionDataError):
        load_stage4_positions(m)


def test_missing_columns(tmp_path):
    import pandas as pd

    m = syn.write_stage4_run(tmp_path / "a", _rows(), csv_only=True)
    csv = m.parent / "positions.csv"
    pd.read_csv(csv).drop(columns=["cov_xy"]).to_csv(csv, index=False)
    with pytest.raises(MissingPositionDataError):
        load_stage4_positions(m)


def test_full_covariance_rotated_with_position():
    rows = _rows()
    df = frame(rows)
    c = covariance_matrix(df.iloc[0])
    assert c[0, 1] != 0 and c[1, 2] != 0  # off-diagonal terms kept
    ms = build_measurements(df, "strict")
    m = itrs_to_gcrs_matrices([ms.measurements[0].time_s])[0]
    y = df[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy(float)[0]
    assert np.allclose(ms.measurements[0].position_gcrs_m, m @ y)
    assert np.allclose(ms.measurements[0].covariance_gcrs_m2, m @ c @ m.T)


def test_policy_duplicates_and_sorting():
    rows = _rows(6, flags={1: ["POOR_GEOMETRY"], 2: ["HIGH_CHI2"], 3: ["NEGATIVE_RANGE"]})
    rows[5]["time_utc"] = rows[4]["time_utc"]  # exact duplicate timestamp
    rows = [rows[0], rows[2], rows[1], rows[3], rows[4], rows[5]]  # out of order
    ms = build_measurements(frame(rows), "strict")
    assert not ms.sorted_input
    t = [m.time_s for m in ms.measurements]
    assert t == sorted(t)
    reasons = [m.reason for m in ms.measurements]
    assert reasons == ["", "POOR_GEOMETRY", "", "NEGATIVE_RANGE", "", "DUPLICATE_TIMESTAMP"]
    assert ms.rejected == {"POOR_GEOMETRY": 1, "NEGATIVE_RANGE": 1, "DUPLICATE_TIMESTAMP": 1}
    perm = build_measurements(frame(rows), "permissive")
    assert [m.reason for m in perm.measurements][1] == ""


def test_invalid_covariance_is_rejected_not_repaired():
    rows = _rows()
    rows[2]["cov_xx"] = -1.0
    ms = build_measurements(frame(rows), "permissive")
    m = ms.measurements[2]
    assert not m.used and m.reason == "INVALID_STAGE4_COVARIANCE" and m.covariance_problem
    assert not m.finite
