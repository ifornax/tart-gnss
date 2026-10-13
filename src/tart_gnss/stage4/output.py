"""Stage-4 tables (CSV/Parquet) and ``stage4_manifest.json``.

``positions`` (schema :data:`STAGE4_SCHEMA_VERSION`): one row per matched
reference epoch. Positions and covariances in ITRS/ECEF metres (m^2);
``latitude_deg``/``longitude_deg``/``altitude_m`` are WGS84 geodetic
conveniences derived from the ECEF solution. ``quality_flags`` is a
``;``-separated string in CSV and a list in Parquet. ``valid`` is true when no
blocking flag is set.

``station_contributions``: one row per matched station per epoch.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence

STAGE4_SCHEMA_VERSION = 1

POSITION_COLUMNS: List[str] = [
    "time_utc", "time_mjd", "norad_id", "reference_station", "method",
    "x_ecef_m", "y_ecef_m", "z_ecef_m",
    "cov_xx", "cov_xy", "cov_xz", "cov_yy", "cov_yz", "cov_zz",
    "sigma_x_m", "sigma_y_m", "sigma_z_m",
    "n_stations", "stations",
    "rms_line_residual_m", "max_line_residual_m", "max_normalised_residual",
    "condition_number", "matrix_rank", "min_crossing_angle_deg", "max_crossing_angle_deg",
    "chi2", "dof", "reduced_chi2",
    "max_time_offset_s", "max_motion_during_offset_deg", "mean_light_time_s",
    "x_p16", "x_p50", "x_p84", "y_p16", "y_p50", "y_p84", "z_p16", "z_p50", "z_p84",
    "bootstrap_n_success", "bootstrap_failure_fraction",
    "solver", "n_iterations", "converged", "fit_success", "valid", "quality_flags",
    "latitude_deg", "longitude_deg", "altitude_m",
]

CONTRIBUTION_COLUMNS: List[str] = [
    "time_utc", "reference_time_utc", "station_id", "station_time_utc", "station_time_offset_s",
    "station_x_ecef_m", "station_y_ecef_m", "station_z_ecef_m",
    "az_obs_deg", "el_obs_deg", "sigma_az_deg", "sigma_el_deg", "cov_az_el",
    "l", "m", "sigma_theta_deg", "direction_source",
    "sight_x", "sight_y", "sight_z",
    "range_m", "line_residual_m", "normalised_residual",
    "angular_rate_deg_s", "motion_during_offset_deg", "light_time_s",
    "stage3_quality_flags", "stage4_station_flags",
]

FLAG_COLUMNS = ("quality_flags", "stage3_quality_flags", "stage4_station_flags")


def _frame(rows: Sequence[Dict[str, object]], columns: List[str]):
    import pandas as pd

    df = pd.DataFrame(list(rows))
    for c in columns:
        if c not in df.columns:
            df[c] = None
    return df[columns]


def _write(df, out_dir: Path, stem: str, csv: bool = True) -> Dict[str, Optional[str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Optional[str]] = {"csv": None, "parquet": None}
    if csv:
        c = df.copy()
        for col in FLAG_COLUMNS:
            if col in c.columns:
                c[col] = c[col].map(lambda f: ";".join(f) if isinstance(f, list) else (f or ""))
        p = out_dir / f"{stem}.csv"
        c.to_csv(p, index=False, float_format="%.12g")
        paths["csv"] = str(p)
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return paths
    p = out_dir / f"{stem}.parquet"
    df.to_parquet(p, index=False)
    paths["parquet"] = str(p)
    return paths


def write_tables(positions: Sequence[Dict[str, object]], contributions: Sequence[Dict[str, object]],
                 out_dir: Path) -> Dict[str, Optional[str]]:
    """Write ``positions.{csv,parquet}`` and ``station_contributions.{csv,parquet}``."""
    p = _write(_frame(positions, POSITION_COLUMNS), out_dir, "positions")
    c = _write(_frame(contributions, CONTRIBUTION_COLUMNS), out_dir, "station_contributions")
    return {"positions_csv": p["csv"], "positions_parquet": p["parquet"],
            "station_contributions_csv": c["csv"], "station_contributions_parquet": c["parquet"]}


def read_table(path):
    """Read a Stage-4 table (CSV or Parquet); flag columns become lists."""
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
    for col in ("time_utc", "reference_time_utc", "station_time_utc"):
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], utc=True)
    return df


def write_manifest(doc: Dict[str, object], path: Path) -> str:
    """Write the Stage-4 manifest as indented JSON (NaN/Inf written as null)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(doc), indent=2) + "\n", encoding="utf-8")
    return str(path)


def _clean(obj):
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if hasattr(obj, "tolist") and not isinstance(obj, (str, bytes)):
        return _clean(obj.tolist())
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj
