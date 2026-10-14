"""Stage-6 tables (CSV/Parquet), ``summary.json`` and ``stage6_manifest.json``.

``residuals`` (schema :data:`STAGE6_SCHEMA_VERSION`): one row per Stage-5
epoch. All residuals are ``observed - TLE``; Cartesian components are GCRS.
``position_/state_mahalanobis_*`` are *squared* Mahalanobis distances
(``D^2``): ``*_obs_cov`` with the Stage-5 covariance only (conditional on the
TLE being fixed), ``*_total_cov`` with the empirical TLE reference dispersion
added (NaN without an ensemble of >= 3 TLEs). ``cov_tle_dispersion_<a>_<b>``
are the 21 upper-triangle terms of the per-epoch 6x6 empirical TLE dispersion
(position and velocity blocks and their cross terms; GCRS; NaN without an
ensemble); rebuild with ``tart_gnss.stage5.output.state_covariance(row,
"tle_dispersion")``. Rows not compared (Stage-5
policy) keep their time and flags with NaN residuals.

``tle_trajectory``: the TLE-predicted state at the same epochs (GCRS position
and velocity, ECEF position) and the RTN basis vectors, so later stages need
not re-propagate.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence

from tart_gnss.stage4.output import write_manifest  # noqa: F401  (same JSON conventions)
from tart_gnss.stage5.output import covariance_columns

STAGE6_SCHEMA_VERSION = 1

RESIDUAL_COLUMNS = [
    "time_utc", "time_mjd", "norad_id", "tle_epoch_utc", "tle_age_days", "tle_propagation_direction",
    "obs_x_gcrs_m", "obs_y_gcrs_m", "obs_z_gcrs_m", "obs_vx_gcrs_m_s", "obs_vy_gcrs_m_s", "obs_vz_gcrs_m_s",
    "tle_x_gcrs_m", "tle_y_gcrs_m", "tle_z_gcrs_m", "tle_vx_gcrs_m_s", "tle_vy_gcrs_m_s", "tle_vz_gcrs_m_s",
    "obs_x_ecef_m", "obs_y_ecef_m", "obs_z_ecef_m", "tle_x_ecef_m", "tle_y_ecef_m", "tle_z_ecef_m",
    "dr_x_m", "dr_y_m", "dr_z_m", "position_residual_3d_m",
    "dv_x_m_s", "dv_y_m_s", "dv_z_m_s", "velocity_residual_3d_m_s",
    "dr_radial_m", "dr_along_track_m", "dr_cross_track_m",
    "dv_radial_m_s", "dv_along_track_m_s", "dv_cross_track_m_s",
    "sigma_radial_obs_m", "sigma_along_track_obs_m", "sigma_cross_track_obs_m",
    "cov_radial_along_obs_m2", "cov_radial_cross_obs_m2", "cov_along_cross_obs_m2",
    "sigma_dv_radial_obs_m_s", "sigma_dv_along_track_obs_m_s", "sigma_dv_cross_track_obs_m_s",
    "sigma_radial_ref_m", "sigma_along_track_ref_m", "sigma_cross_track_ref_m",
    "geocentric_angular_residual_deg", "geocentric_range_residual_m", "equivalent_along_track_time_offset_s",
    "position_mahalanobis_obs_cov", "state_mahalanobis_obs_cov",
    "position_mahalanobis_total_cov", "state_mahalanobis_total_cov",
    "n_tle_ensemble", "stage5_state_source", "stage5_valid", "stage5_rejection_reason", "precision_validation",
    "valid",
    "quality_flags", "stage5_quality_flags",
    *covariance_columns("tle_dispersion"),
]

TLE_TRAJECTORY_COLUMNS = [
    "time_utc", "norad_id", "tle_epoch_utc", "tle_age_days",
    "x_gcrs_m", "y_gcrs_m", "z_gcrs_m", "vx_gcrs_m_s", "vy_gcrs_m_s", "vz_gcrs_m_s",
    "x_ecef_m", "y_ecef_m", "z_ecef_m",
    "r_hat_x", "r_hat_y", "r_hat_z", "t_hat_x", "t_hat_y", "t_hat_z", "n_hat_x", "n_hat_y", "n_hat_z",
]

FLAG_COLUMNS = ("quality_flags", "stage5_quality_flags")


def _frame(rows: Sequence[Dict[str, object]], columns):
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
        c.to_csv(p, index=False, float_format="%.15g")
        paths["csv"] = str(p)
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        if not csv:
            p = out_dir / f"{stem}.csv"
            df.to_csv(p, index=False, float_format="%.15g")
            paths["csv"] = str(p)
        return paths
    p = out_dir / f"{stem}.parquet"
    df.to_parquet(p, index=False)
    paths["parquet"] = str(p)
    return paths


def write_tables(residuals: Sequence[Dict[str, object]], tle_rows: Sequence[Dict[str, object]],
                 out_dir: Path) -> Dict[str, Optional[str]]:
    """Write ``residuals.{csv,parquet}`` and ``tle_trajectory.parquet`` (CSV only without pyarrow)."""
    r = _write(_frame(residuals, RESIDUAL_COLUMNS), out_dir, "residuals")
    t = _write(_frame(tle_rows, TLE_TRAJECTORY_COLUMNS), out_dir, "tle_trajectory", csv=False)
    return {"residuals_csv": r["csv"], "residuals_parquet": r["parquet"],
            "tle_trajectory_parquet": t["parquet"], "tle_trajectory_csv": t["csv"]}


def read_table(path):
    """Read a Stage-6 table; flag columns become lists, ``time_utc`` UTC datetimes."""
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
