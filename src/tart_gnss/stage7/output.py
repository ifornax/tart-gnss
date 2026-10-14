"""Stage-7 tables (CSV/Parquet), ``summary.json`` and ``stage7_manifest.json``.

``diagnosis``          one row per Stage-6 epoch (evidence tags, level, confidence)
``station_influence``  long format, one row per (epoch, station): leave-one-out solution and
                       shift, direction residuals, line residual, metadata, flags
``pairwise``           one row per (epoch, station pair): two-line diagnostic solution
``correlations.csv``   every correlation test with its interpretation tag

Every table carries ``norad_id`` (and the long tables ``station_id``), so
runs for several satellites can later be concatenated without schema
changes. List columns are ``;``-joined strings in CSV and lists in Parquet.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional, Sequence

from tart_gnss.stage4.output import _clean, write_manifest  # noqa: F401  (same JSON conventions)

STAGE7_SCHEMA_VERSION = 1

DIAGNOSIS_COLUMNS = [
    "time_utc", "norad_id", "diagnosed", "position_residual_3d_m", "dr_radial_m", "dr_along_track_m",
    "dr_cross_track_m", "position_mahalanobis_obs_cov", "position_mahalanobis_total_cov", "anomaly_statistic",
    "anomaly_threshold", "statistical_anomaly", "stage4_d2_tle", "strong", "persistent", "anomaly_run_length",
    "anomaly_duration_s", "stage5_nis", "stage5_flags", "stage4_condition_number", "stage4_min_crossing_angle_deg",
    "stage4_n_stations", "stage4_residual_3d_m", "max_station_line_residual_m", "max_station_direction_residual_deg",
    "max_leave_one_out_shift_m", "max_leave_one_out_mahalanobis", "most_influential_station", "station_pick",
    "station_pick_method", "max_pair_separation_z", "predicted_condition_number", "predicted_min_crossing_angle_deg",
    "mean_station_time_offset_s", "max_abs_station_time_offset_s", "tle_age_days", "precision_validation",
    "diagnostic_tags",
    "evidence_level", "diagnostic_confidence", "quality_flags",
]
INFLUENCE_COLUMNS = [
    "time_utc", "norad_id", "station_id", "station_time_offset_s", "n_stations", "n_stations_remaining",
    "loo_x_ecef_m", "loo_y_ecef_m", "loo_z_ecef_m", "loo_shift_3d_m", "loo_shift_radial_m",
    "loo_shift_along_track_m", "loo_shift_cross_track_m", "loo_mahalanobis", "loo_dof", "loo_condition_number",
    "loo_min_crossing_angle_deg", "degraded_geometry", "loo_d2_tle", "direction_residual_angle_deg",
    "direction_residual_az_deg", "direction_residual_el_deg", "direction_residual_along_motion_deg",
    "direction_residual_across_motion_deg", "line_residual_m", "normalised_line_residual", "snr", "reduced_chi2",
    "flagged_fraction", "vis_image_separation_deg", "zero_gain_count", "gain_solution_index",
    "calibration_warning", "timing_warning", "stage3_flags", "stage4_flags",
]
PAIR_COLUMNS = ["time_utc", "norad_id", "pair_station_1", "pair_station_2", "pair_x_ecef_m", "pair_y_ecef_m",
                "pair_z_ecef_m", "pair_line_separation_m", "pair_separation_z", "pair_geometry_angle_deg",
                "n_stations"]
CORRELATION_COLUMNS = ["scope", "station_id", "family", "x_variable", "y_variable", "method", "coefficient",
                       "partial_coefficient", "p_value_naive", "p_value", "p_value_holm", "n_samples", "n_effective", "x_unique",
                       "interpretation_tag"]
LIST_COLUMNS = ("stage5_flags", "diagnostic_tags", "quality_flags", "stage3_flags", "stage4_flags")


def _frame(rows, columns):
    import pandas as pd

    df = pd.DataFrame(list(rows))
    for c in columns:
        if c not in df.columns:
            df[c] = None
    return df[columns]


def _write(df, out_dir: Path, stem: str) -> Dict[str, Optional[str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    c = df.copy()
    for col in LIST_COLUMNS:
        if col in c.columns:
            c[col] = c[col].map(lambda f: ";".join(map(str, f)) if isinstance(f, list) else (f or ""))
    p = out_dir / f"{stem}.csv"
    c.to_csv(p, index=False, float_format="%.12g")
    paths: Dict[str, Optional[str]] = {"csv": str(p), "parquet": None}
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return paths
    p = out_dir / f"{stem}.parquet"
    df.to_parquet(p, index=False)
    paths["parquet"] = str(p)
    return paths


def write_tables(diagnosis: Sequence[dict], influence: Sequence[dict], pairs: Sequence[dict], correlations,
                 out_dir: Path) -> Dict[str, Optional[str]]:
    """Write the diagnosis, station-influence, pairwise and correlation tables."""
    d = _write(_frame(diagnosis, DIAGNOSIS_COLUMNS), out_dir, "diagnosis")
    i = _write(_frame(influence, INFLUENCE_COLUMNS), out_dir, "station_influence")
    p = _write(_frame(pairs, PAIR_COLUMNS), out_dir, "pairwise")
    cp = out_dir / "correlations.csv"
    _frame(correlations.to_dict("records") if len(correlations) else [], CORRELATION_COLUMNS).to_csv(
        cp, index=False, float_format="%.6g")
    return {"diagnosis_csv": d["csv"], "diagnosis_parquet": d["parquet"], "station_influence_csv": i["csv"],
            "station_influence_parquet": i["parquet"], "pairwise_csv": p["csv"], "pairwise_parquet": p["parquet"],
            "correlations_csv": str(cp)}


def read_table(path):
    """Read a Stage-7 table; list columns become lists."""
    import pandas as pd

    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
        for col in LIST_COLUMNS:
            if col in df.columns:
                df[col] = df[col].map(lambda f: list(f) if f is not None else [])
    else:
        df = pd.read_csv(path)
        for col in LIST_COLUMNS:
            if col in df.columns:
                df[col] = df[col].map(lambda f: [x for x in str(f).split(";") if x] if isinstance(f, str) else [])
    if "time_utc" in df.columns:
        df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    return df


def write_json(doc, path: Path) -> str:
    """Indented JSON with NaN/Inf as null."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(doc), indent=2) + "\n", encoding="utf-8")
    return str(path)
