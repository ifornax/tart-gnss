"""Read a Stage-5 run: the smoothed trajectory and its full 6x6 covariances.

The manifest must be a Stage-5 manifest (``stage == 5``, ``schema_version ==
1``, ``state_frame == "GCRS"``) and is hashed, as is the trajectory table.
The table is looked for next to the manifest first (Parquet, then CSV), then
at the recorded output paths. Stage 6 only reads these files; it never
writes to or modifies the Stage-5 run.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from tart_gnss.stage4.loader import sha256_file
from tart_gnss.stage5.output import covariance_columns, state_covariance
from tart_gnss.stage6.errors import InvalidStage5ManifestError, MissingTrajectoryError

STATE_COLUMNS = ["x_gcrs_smoothed_m", "y_gcrs_smoothed_m", "z_gcrs_smoothed_m",
                 "vx_gcrs_smoothed_m_s", "vy_gcrs_smoothed_m_s", "vz_gcrs_smoothed_m_s"]
FILTERED_COLUMNS = ["x_gcrs_filtered_m", "y_gcrs_filtered_m", "z_gcrs_filtered_m",
                    "vx_gcrs_filtered_m_s", "vy_gcrs_filtered_m_s", "vz_gcrs_filtered_m_s"]
REQUIRED_COLUMNS = (
    "time_utc", "time_mjd", "norad_id", *STATE_COLUMNS,
    "x_ecef_smoothed_m", "y_ecef_smoothed_m", "z_ecef_smoothed_m", "valid", "quality_flags",
    *covariance_columns("smoothed"),
)


@dataclass
class Stage5Trajectory:
    """A loaded Stage-5 run."""

    manifest_path: Path
    manifest: Dict[str, object]
    sha256: str
    trajectory_path: Path
    trajectory_sha256: str
    table: object
    norad_id: Optional[int]
    times_s: np.ndarray
    innovations_path: Optional[Path] = None

    @property
    def midpoint_s(self) -> Optional[float]:
        """Midpoint of the epochs with a finite smoothed state (Unix s)."""
        st = self.table[STATE_COLUMNS].to_numpy(dtype=float)
        ok = np.all(np.isfinite(st), axis=1)
        t = self.times_s[ok] if ok.any() else self.times_s
        return float(0.5 * (t.min() + t.max())) if len(t) else None


def load_stage5_manifest(path) -> Tuple[Dict[str, object], str]:
    """Read and validate ``stage5_manifest.json``; return ``(manifest, sha256)``."""
    path = Path(path)
    if not path.is_file():
        raise InvalidStage5ManifestError(f"Stage-5 manifest not found: {path}")
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise InvalidStage5ManifestError(f"cannot read Stage-5 manifest {path}: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("stage") != 5:
        raise InvalidStage5ManifestError(f"{path} is not a Stage-5 manifest (stage != 5)")
    if doc.get("schema_version") != 1:
        raise InvalidStage5ManifestError(f"unsupported Stage-5 schema_version {doc.get('schema_version')!r}")
    if doc.get("state_frame", "GCRS") != "GCRS":
        raise InvalidStage5ManifestError(f"{path}: unsupported state_frame {doc.get('state_frame')!r}")
    return doc, sha256_file(path)


def _resolve(manifest: Dict[str, object], manifest_path: Path, stem: str) -> Optional[Path]:
    outputs = manifest.get("outputs") or {}
    base = manifest_path.parent
    cands: List[Path] = [base / f"{stem}.parquet", base / f"{stem}.csv"]
    for ext in ("parquet", "csv"):
        rec = outputs.get(f"{stem}_{ext}")
        if rec:
            cands += [Path(rec), base / Path(rec).name]
    for c in cands:
        if c.is_file():
            if c.suffix == ".parquet":
                try:
                    import pyarrow  # noqa: F401
                except ImportError:
                    continue
            return c
    return None


def load_stage5_trajectory(manifest_path) -> Stage5Trajectory:
    """Load the Stage-5 smoothed trajectory and full covariance."""
    from tart_gnss.stage4.timing import time_seconds
    from tart_gnss.stage5.output import read_table

    manifest_path = Path(manifest_path)
    manifest, sha = load_stage5_manifest(manifest_path)
    tpath = _resolve(manifest, manifest_path, "trajectory")
    if tpath is None:
        raise MissingTrajectoryError(f"no Stage-5 trajectory table found for {manifest_path}")
    try:
        df = read_table(tpath)
    except Exception as exc:  # noqa: BLE001 - any reader failure is a data error
        raise MissingTrajectoryError(f"cannot read Stage-5 trajectory {tpath}: {exc}") from exc
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise MissingTrajectoryError(f"{tpath}: missing required columns {missing[:6]}"
                                     f"{' ...' if len(missing) > 6 else ''}")
    df["valid"] = df["valid"].map(lambda v: bool(v) if v == v and v is not None else False)
    norad = manifest.get("norad_id")
    if norad is None and df["norad_id"].notna().any():
        norad = df["norad_id"].dropna().iloc[0]
    return Stage5Trajectory(
        manifest_path=manifest_path, manifest=manifest, sha256=sha, trajectory_path=tpath,
        trajectory_sha256=sha256_file(tpath), table=df, norad_id=int(norad) if norad is not None else None,
        times_s=time_seconds(df["time_utc"]) if len(df) else np.zeros(0),
        innovations_path=_resolve(manifest, manifest_path, "innovations"),
    )


def row_state(row, kind: str = "smoothed") -> Tuple[np.ndarray, np.ndarray]:
    """``(state (6,), covariance (6, 6))`` of one trajectory row (``kind`` smoothed or filtered)."""
    cols = STATE_COLUMNS if kind == "smoothed" else FILTERED_COLUMNS
    x = np.array([row[c] for c in cols], dtype=float)
    try:
        p = state_covariance(row, kind)
    except (KeyError, TypeError, ValueError):
        p = np.full((6, 6), np.nan)
    return x, p
