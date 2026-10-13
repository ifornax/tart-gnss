"""Read Stage-4 outputs and turn them into GCRS position measurements.

The Stage-4 manifest is validated (``stage == 4``, ``schema_version == 1``)
and hashed. The position table is looked for next to the manifest first
(Parquet, then CSV), then at the path recorded in ``outputs``, so a moved or
copied result directory still loads. Only Stage-4 data enter: the
``satellite_reference`` block of the Stage-4 manifest (TLE information kept
for Stage 6) is never read.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from tart_gnss.stage4.loader import sha256_file
from tart_gnss.stage5 import quality as q
from tart_gnss.stage5.errors import InvalidStage4ManifestError, MissingPositionDataError
from tart_gnss.stage5.filter import TrajectoryMeasurement
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices, rotate_measurement

REQUIRED_COLUMNS = (
    "time_utc", "time_mjd", "norad_id", "x_ecef_m", "y_ecef_m", "z_ecef_m",
    "cov_xx", "cov_xy", "cov_xz", "cov_yy", "cov_yz", "cov_zz", "sigma_x_m", "sigma_y_m", "sigma_z_m",
    "n_stations", "condition_number", "min_crossing_angle_deg", "max_crossing_angle_deg",
    "chi2", "reduced_chi2", "fit_success", "quality_flags",
)
DUPLICATE_TIMESTAMP = "DUPLICATE_TIMESTAMP"
DUPLICATE_TOL_S = 1e-6


@dataclass
class Stage4TrajectoryInput:
    """A loaded Stage-4 run."""

    manifest_path: Path
    manifest: Dict[str, object]
    sha256: str
    positions_path: Path
    positions_sha256: str
    positions: object
    norad_id: Optional[int]
    contributions_path: Optional[Path] = None


@dataclass
class MeasurementSet:
    """Time-ordered measurements, their ECEF->GCRS matrices and the rejection bookkeeping."""

    measurements: List[TrajectoryMeasurement]
    matrices: np.ndarray
    rows: object
    rejected: Dict[str, int]
    n_input: int
    sorted_input: bool
    covariance_notes: List[str]


def load_stage4_manifest(path) -> Tuple[Dict[str, object], str]:
    """Read and validate ``stage4_manifest.json``; return ``(manifest, sha256)``."""
    path = Path(path)
    if not path.is_file():
        raise InvalidStage4ManifestError(f"Stage-4 manifest not found: {path}")
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise InvalidStage4ManifestError(f"cannot read Stage-4 manifest {path}: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("stage") != 4:
        raise InvalidStage4ManifestError(f"{path} is not a Stage-4 manifest (stage != 4)")
    if doc.get("schema_version") != 1:
        raise InvalidStage4ManifestError(f"unsupported Stage-4 schema_version {doc.get('schema_version')!r}")
    frame = doc.get("coordinate_frame")
    if frame is not None and frame != "ITRS/ECEF":
        raise InvalidStage4ManifestError(f"{path}: unsupported coordinate_frame {frame!r} (expected ITRS/ECEF)")
    return doc, sha256_file(path)


def _resolve(manifest: Dict[str, object], manifest_path: Path, stem: str) -> Optional[Path]:
    outputs = manifest.get("outputs") or {}
    base = manifest_path.parent
    cands: List[Path] = []
    for ext in ("parquet", "csv"):
        cands.append(base / f"{stem}.{ext}")
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


def resolve_positions_path(manifest: Dict[str, object], manifest_path: Path) -> Path:
    """Locate ``positions.{parquet,csv}`` (next to the manifest first)."""
    p = _resolve(manifest, manifest_path, "positions")
    if p is None:
        raise MissingPositionDataError(f"no Stage-4 position table found for {manifest_path}")
    return p


def read_positions(path: Path):
    """Read a Stage-4 position table and check the required columns."""
    from tart_gnss.stage4.output import read_table

    try:
        df = read_table(path)
    except Exception as exc:  # noqa: BLE001 - any reader failure is a data error
        raise MissingPositionDataError(f"cannot read Stage-4 positions {path}: {exc}") from exc
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise MissingPositionDataError(f"{path}: missing required columns {missing}")
    df["fit_success"] = df["fit_success"].map(lambda v: bool(v) if v == v and v is not None else False)
    return df


def load_stage4_positions(manifest_path) -> Stage4TrajectoryInput:
    """Load validated Stage-4 positions, covariances and provenance."""
    manifest_path = Path(manifest_path)
    manifest, sha = load_stage4_manifest(manifest_path)
    ppath = resolve_positions_path(manifest, manifest_path)
    df = read_positions(ppath)
    norad = manifest.get("norad_id")
    if norad is None and len(df):
        norad = df["norad_id"].dropna().iloc[0] if df["norad_id"].notna().any() else None
    return Stage4TrajectoryInput(
        manifest_path=manifest_path, manifest=manifest, sha256=sha, positions_path=ppath,
        positions_sha256=sha256_file(ppath), positions=df, norad_id=int(norad) if norad is not None else None,
        contributions_path=_resolve(manifest, manifest_path, "station_contributions"),
    )


def covariance_matrix(row) -> np.ndarray:
    """Full 3x3 ECEF covariance from the six Stage-4 components (never reduced to the diagonal)."""
    xx, xy, xz = float(row["cov_xx"]), float(row["cov_xy"]), float(row["cov_xz"])
    yy, yz, zz = float(row["cov_yy"]), float(row["cov_yz"]), float(row["cov_zz"])
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def build_measurements(df, policy: str) -> MeasurementSet:
    """Sort, apply the quality policy, validate covariances and rotate into GCRS."""
    from tart_gnss.stage4.timing import time_seconds

    n_input = len(df)
    t_raw = time_seconds(df["time_utc"]) if n_input else np.zeros(0)
    sorted_input = bool(np.all(np.diff(t_raw) >= 0)) if n_input else True
    order = np.argsort(t_raw, kind="stable")
    df = df.iloc[order].reset_index(drop=True)
    t = t_raw[order]
    mats = itrs_to_gcrs_matrices(t) if n_input else np.zeros((0, 3, 3))
    out: List[TrajectoryMeasurement] = []
    rejected: Dict[str, int] = {}
    notes: List[str] = []
    prev_t = -np.inf
    for k, row in df.iterrows():
        flags = list(row["quality_flags"])
        y = np.array([row["x_ecef_m"], row["y_ecef_m"], row["z_ecef_m"]], dtype=float)
        r = covariance_matrix(row)
        reason = q.stage4_rejection_reason(flags, bool(row["fit_success"]), policy)
        if not reason and not np.all(np.isfinite(y)):
            reason = q.NO_MEASUREMENT
        cov_reason = ""
        if np.all(np.isfinite(y)):
            r, cov_reason = q.validate_covariance(r)
            if cov_reason:
                notes.append(f"{row['time_utc']}: {cov_reason}")
                reason = reason or q.INVALID_STAGE4_COVARIANCE
        if not reason and t[k] - prev_t <= DUPLICATE_TOL_S:
            reason = DUPLICATE_TIMESTAMP
        if not reason:
            prev_t = t[k]
        if cov_reason or not np.all(np.isfinite(y)):
            yg, rg = np.full(3, np.nan), np.full((3, 3), np.nan)
        else:
            yg, rg = rotate_measurement(mats[k], y, r)
        if reason:
            rejected[reason] = rejected.get(reason, 0) + 1
        m = TrajectoryMeasurement(time_s=float(t[k]), position_gcrs_m=yg, covariance_gcrs_m2=rg,
                                  used=not reason, reason=reason, index=int(k), stage4_flags=flags,
                                  covariance_problem=cov_reason)
        out.append(m)
    return MeasurementSet(measurements=out, matrices=mats, rows=df, rejected=rejected, n_input=n_input,
                          sorted_input=sorted_input, covariance_notes=notes)
