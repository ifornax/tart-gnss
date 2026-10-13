"""Stage-3 tables (CSV/Parquet) and the ``stage3_manifest.json``.

Table schema (version :data:`STAGE3_SCHEMA_VERSION`): one row per station and
non-empty time bin. Angles in degrees, ``l``/``m`` dimensionless direction
cosines (East, North), ``time_mjd`` in days (UTC MJD). Columns are listed in
:data:`COLUMNS`. ``az_tle_deg``/``el_tle_deg``/``tle_separation_deg`` are TLE
*reference predictions* for convenience, never fitted quantities.
``quality_flags`` is a ``;``-separated string in CSV and a list in Parquet.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

STAGE3_SCHEMA_VERSION = 1

COLUMNS: List[str] = [
    "time_utc", "time_mjd", "station_id", "tart_id", "norad_id", "primary_method",
    "n_times", "n_samples",
    "l", "m", "sigma_l", "sigma_m", "cov_lm",
    "az_obs_deg", "el_obs_deg", "sigma_az_deg", "sigma_el_deg", "cov_az_el",
    "complex_amp_real", "complex_amp_imag", "apparent_amplitude",
    "chi2", "reduced_chi2", "rms_complex_residual", "median_abs_residual", "max_abs_residual",
    "sigma_vis", "residual_tail_ratio", "snr",
    "n_baselines_total", "n_baselines_used", "flagged_fraction",
    "az_image_deg", "el_image_deg", "image_peak_value", "image_snr", "vis_image_separation_deg",
    "az_boot_p16", "az_boot_p50", "az_boot_p84", "el_boot_p16", "el_boot_p50", "el_boot_p84",
    "fit_success", "optimizer_status", "usable", "quality_flags",
    "stage2_three_station_valid",
    "az_tle_deg", "el_tle_deg", "range_tle_km", "tle_separation_deg",
]

#: Columns that hold TLE reference predictions (not measurements).
REFERENCE_COLUMNS = ("az_tle_deg", "el_tle_deg", "range_tle_km", "tle_separation_deg")


def _frame(rows: Sequence[Dict[str, object]]):
    import pandas as pd

    df = pd.DataFrame(list(rows))
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = None
    return df[COLUMNS]


def write_tables(rows: Sequence[Dict[str, object]], out_dir: Path) -> Dict[str, Optional[str]]:
    """Write ``directions.csv`` and (if pyarrow is installed) ``directions.parquet``.

    Returns the written paths (``parquet`` is ``None`` when unavailable).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _frame(rows)
    csv = df.copy()
    csv["quality_flags"] = csv["quality_flags"].map(lambda f: ";".join(f) if isinstance(f, list) else (f or ""))
    csv_path = out_dir / "directions.csv"
    csv.to_csv(csv_path, index=False, float_format="%.10g")
    parquet_path: Optional[str] = None
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pass
    else:
        p = out_dir / "directions.parquet"
        df.to_parquet(p, index=False)
        parquet_path = str(p)
    return {"csv": str(csv_path), "parquet": parquet_path}


def read_directions(path):
    """Read a Stage-3 table (CSV or Parquet); ``quality_flags`` becomes a list."""
    import pandas as pd

    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
        df["quality_flags"] = df["quality_flags"].map(lambda f: list(f) if f is not None else [])
    else:
        df = pd.read_csv(path, keep_default_na=True)
        df["quality_flags"] = df["quality_flags"].map(
            lambda f: [x for x in str(f).split(";") if x] if isinstance(f, str) else [])
    df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    return df


def write_manifest(doc: Dict[str, object], path: Path) -> str:
    """Write the Stage-3 manifest as indented JSON (NaN written as null)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(doc), indent=2) + "\n", encoding="utf-8")
    return str(path)


def _clean(obj):
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, float) and obj != obj:
        return None
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        return _clean(obj.item())
    return obj
