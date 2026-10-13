"""Stage-5 tables (CSV/Parquet) and ``stage5_manifest.json``.

``trajectory`` (schema :data:`STAGE5_SCHEMA_VERSION`): one row per Stage-4
epoch. States are GCRS (m, m/s); ``*_ecef_*`` columns are the GCRS position
rotated back to ITRS/ECEF at the same epoch. ``sigma_*`` columns are GCRS
1-sigma values. The full 6x6 state covariances are stored as their 21
upper-triangle elements ``cov_filtered_<a>_<b>`` / ``cov_smoothed_<a>_<b>``
with ``a, b`` in ``x y z vx vy vz`` (m^2, m^2/s, m^2/s^2); use
:func:`state_covariance` to rebuild them. ``cov_ecef_smoothed_*`` is the
smoothed 3x3 position covariance rotated to ECEF. ``quality_flags`` is a
``;``-separated string in CSV and a list in Parquet.

``innovations``: one row per Stage-4 epoch with a finite position, from the
final filter pass; ``nis_first_pass`` is the NIS of the first (sequential)
evaluation, which differs only for rows re-processed after a divergence restart.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from tart_gnss.stage4.output import write_manifest  # noqa: F401  (same JSON conventions as Stage 4)

STAGE5_SCHEMA_VERSION = 1
STATE_LABELS = ("x", "y", "z", "vx", "vy", "vz")
_IJ = [(i, j) for i in range(6) for j in range(i, 6)]


def covariance_columns(kind: str) -> List[str]:
    """The 21 upper-triangle column names of a 6x6 state covariance (``kind`` = filtered/smoothed)."""
    return [f"cov_{kind}_{STATE_LABELS[i]}_{STATE_LABELS[j]}" for i, j in _IJ]


def covariance_fields(p: np.ndarray, kind: str) -> Dict[str, float]:
    """Upper-triangle fields of one 6x6 covariance."""
    return {name: float(p[i, j]) for name, (i, j) in zip(covariance_columns(kind), _IJ)}


def state_covariance(row, kind: str = "smoothed") -> np.ndarray:
    """Rebuild the 6x6 state covariance from a trajectory row."""
    p = np.empty((6, 6))
    for name, (i, j) in zip(covariance_columns(kind), _IJ):
        p[i, j] = p[j, i] = float(row[name])
    return p


ECEF_COV_COLUMNS = ["cov_ecef_smoothed_xx", "cov_ecef_smoothed_xy", "cov_ecef_smoothed_xz",
                    "cov_ecef_smoothed_yy", "cov_ecef_smoothed_yz", "cov_ecef_smoothed_zz"]

TRAJECTORY_COLUMNS: List[str] = [
    "time_utc", "time_mjd", "norad_id", "dynamics",
    "x_gcrs_filtered_m", "y_gcrs_filtered_m", "z_gcrs_filtered_m",
    "vx_gcrs_filtered_m_s", "vy_gcrs_filtered_m_s", "vz_gcrs_filtered_m_s",
    "x_gcrs_smoothed_m", "y_gcrs_smoothed_m", "z_gcrs_smoothed_m",
    "vx_gcrs_smoothed_m_s", "vy_gcrs_smoothed_m_s", "vz_gcrs_smoothed_m_s",
    "x_ecef_smoothed_m", "y_ecef_smoothed_m", "z_ecef_smoothed_m", "speed_smoothed_m_s",
    "sigma_x_filtered_m", "sigma_y_filtered_m", "sigma_z_filtered_m",
    "sigma_x_smoothed_m", "sigma_y_smoothed_m", "sigma_z_smoothed_m",
    "sigma_vx_smoothed_m_s", "sigma_vy_smoothed_m_s", "sigma_vz_smoothed_m_s",
    "acceleration_m_s2", "radius_m",
    "stage4_x_ecef_m", "stage4_y_ecef_m", "stage4_z_ecef_m",
    "stage4_x_gcrs_m", "stage4_y_gcrs_m", "stage4_z_gcrs_m",
    "stage4_sigma_x_gcrs_m", "stage4_sigma_y_gcrs_m", "stage4_sigma_z_gcrs_m",
    "measurement_used", "measurement_accepted", "rejection_reason", "time_since_update_s",
    "nis", "valid", "quality_flags", "stage4_quality_flags",
    *ECEF_COV_COLUMNS, *covariance_columns("filtered"), *covariance_columns("smoothed"),
]

INNOVATION_COLUMNS: List[str] = [
    "time_utc", "innovation_x_m", "innovation_y_m", "innovation_z_m",
    "innovation_cov_xx", "innovation_cov_xy", "innovation_cov_xz",
    "innovation_cov_yy", "innovation_cov_yz", "innovation_cov_zz",
    "whitened_x", "whitened_y", "whitened_z",
    "nis", "nis_first_pass", "gate_threshold", "accepted", "attempted", "inflation",
    "stage4_quality_flags", "stage5_quality_flags",
]

FLAG_COLUMNS = ("quality_flags", "stage4_quality_flags", "stage5_quality_flags")


def _frame(rows: Sequence[Dict[str, object]], columns: List[str]):
    import pandas as pd

    df = pd.DataFrame(list(rows))
    for c in columns:
        if c not in df.columns:
            df[c] = None
    return df[columns]


def _write(df, out_dir: Path, stem: str) -> Dict[str, Optional[str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    c = df.copy()
    for col in FLAG_COLUMNS:
        if col in c.columns:
            c[col] = c[col].map(lambda f: ";".join(f) if isinstance(f, list) else (f or ""))
    p = out_dir / f"{stem}.csv"
    c.to_csv(p, index=False, float_format="%.15g")
    paths: Dict[str, Optional[str]] = {"csv": str(p), "parquet": None}
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return paths
    p = out_dir / f"{stem}.parquet"
    df.to_parquet(p, index=False)
    paths["parquet"] = str(p)
    return paths


def write_tables(trajectory: Sequence[Dict[str, object]], innovations: Sequence[Dict[str, object]],
                 out_dir: Path) -> Dict[str, Optional[str]]:
    """Write ``trajectory.{csv,parquet}`` and ``innovations.{csv,parquet}``."""
    t = _write(_frame(trajectory, TRAJECTORY_COLUMNS), out_dir, "trajectory")
    i = _write(_frame(innovations, INNOVATION_COLUMNS), out_dir, "innovations")
    return {"trajectory_csv": t["csv"], "trajectory_parquet": t["parquet"],
            "innovations_csv": i["csv"], "innovations_parquet": i["parquet"]}


def read_table(path):
    """Read a Stage-5 table (CSV or Parquet); flag columns become lists, times UTC datetimes."""
    import pandas as pd

    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
        for col in FLAG_COLUMNS:
            if col in df.columns:
                df[col] = df[col].map(lambda f: list(f) if f is not None else [])
    else:
        df = pd.read_csv(path)
        for col in FLAG_COLUMNS:
            if col in df.columns:
                df[col] = df[col].map(lambda f: [x for x in str(f).split(";") if x] if isinstance(f, str) else [])
    if "time_utc" in df.columns:
        df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    return df
