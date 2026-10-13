"""Stage-5 tables: 6x6 covariance persistence and round trips."""

import numpy as np
import pytest

pytest.importorskip("pandas")

from tart_gnss.stage5.output import (  # noqa: E402
    INNOVATION_COLUMNS,
    TRAJECTORY_COLUMNS,
    covariance_columns,
    covariance_fields,
    read_table,
    state_covariance,
    write_tables,
)


def test_covariance_columns():
    cols = covariance_columns("smoothed")
    assert len(cols) == 21 and cols[0] == "cov_smoothed_x_x" and cols[-1] == "cov_smoothed_vz_vz"
    assert set(cols) <= set(TRAJECTORY_COLUMNS) and set(covariance_columns("filtered")) <= set(TRAJECTORY_COLUMNS)


def test_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    a = rng.normal(size=(6, 6))
    p = a @ a.T * 1e6
    row = {"time_utc": "2026-06-25T09:14:00.000000Z", "valid": True, "quality_flags": ["HIGH_NIS", "EXTRAPOLATED"],
           "stage4_quality_flags": ["GOOD"], **covariance_fields(p, "smoothed"), **covariance_fields(p, "filtered")}
    inn = {"time_utc": row["time_utc"], "nis": 1.0, "stage4_quality_flags": [], "stage5_quality_flags": ["GOOD"]}
    paths = write_tables([row], [inn], tmp_path)
    for key in ("trajectory_csv", "trajectory_parquet"):
        if paths[key] is None:
            continue
        df = read_table(paths[key])
        assert list(df.columns) == TRAJECTORY_COLUMNS
        assert np.allclose(state_covariance(df.iloc[0], "smoothed"), p, rtol=1e-12)
        assert df["quality_flags"].iloc[0] == ["HIGH_NIS", "EXTRAPOLATED"]
        assert str(df["time_utc"].dt.tz) == "UTC"
    assert list(read_table(paths["innovations_csv"]).columns) == INNOVATION_COLUMNS
