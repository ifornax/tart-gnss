"""Read Stage-3 outputs and the frozen station coordinates for Stage 4.

Station coordinates come from the run's own provenance: the ``stations``
section of ``stage3_manifest.json``, which Stage 3 copied from the Stage-2
manifest. When that Stage-2 manifest is still reachable and its SHA-256 matches
the hash recorded by Stage 3, the coordinates are cross-checked against it
(mismatch is an error). The live station registry is never consulted, so a
later registry update cannot change a past result.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage4.errors import InvalidStage3ManifestError, MissingDirectionDataError, TriangulationError
from tart_gnss.stage4.geometry import (
    azel_covariance_to_ecef,
    azel_to_ecef_unit_vector,
    geodetic_to_ecef,
    lm_covariance_to_ecef,
    lm_to_ecef_unit_vector,
)
from tart_gnss.stage4.quality import stage3_rejection_reason

REQUIRED_COLUMNS = (
    "time_utc", "station_id", "norad_id", "az_obs_deg", "el_obs_deg", "sigma_az_deg", "sigma_el_deg",
    "cov_az_el", "l", "m", "sigma_l", "sigma_m", "cov_lm", "fit_success", "quality_flags",
)
#: Coordinates agreeing to this tolerance are considered identical (deg, m).
COORD_TOL_DEG = 1e-9
COORD_TOL_M = 1e-6


@dataclass
class StationCoordinates:
    """Frozen geodetic and ECEF coordinates of one station."""

    station_id: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    ecef_m: np.ndarray
    source: str
    cross_check: str
    calibration_warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        return {
            "latitude_deg": self.latitude_deg, "longitude_deg": self.longitude_deg, "altitude_m": self.altitude_m,
            "ecef_m": [float(v) for v in self.ecef_m], "datum": "WGS84", "source": self.source,
            "cross_check": self.cross_check, "calibration_warnings": list(self.calibration_warnings),
        }


@dataclass
class Stage3Run:
    """A loaded Stage-3 run."""

    manifest_path: Path
    manifest: Dict[str, object]
    sha256: str
    directions_path: Path
    directions_sha256: str
    directions: object
    stations: Dict[str, StationCoordinates]
    norad_id: Optional[int]


def sha256_file(path: Path) -> str:
    """Hex SHA-256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_stage3_manifest(path) -> Tuple[Dict[str, object], str]:
    """Read and validate ``stage3_manifest.json``; return ``(manifest, sha256)``."""
    path = Path(path)
    if not path.is_file():
        raise InvalidStage3ManifestError(f"Stage-3 manifest not found: {path}")
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise InvalidStage3ManifestError(f"cannot read Stage-3 manifest {path}: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("stage") != 3:
        raise InvalidStage3ManifestError(f"{path} is not a Stage-3 manifest (stage != 3)")
    if doc.get("schema_version") != 1:
        raise InvalidStage3ManifestError(f"unsupported Stage-3 schema_version {doc.get('schema_version')!r}")
    if not isinstance(doc.get("stations"), dict) or not doc["stations"]:
        raise InvalidStage3ManifestError(f"{path}: no stations recorded")
    return doc, sha256_file(path)


def resolve_directions_path(manifest: Dict[str, object], manifest_path: Path) -> Path:
    """Locate the Stage-3 direction table (Parquet preferred, then CSV).

    Tried in order: the file next to the manifest, then the path recorded in
    ``outputs`` (as given, then relative to the manifest directory).
    """
    outputs = manifest.get("outputs") or {}
    base = manifest_path.parent
    cands: List[Path] = []
    for key, name in (("directions_parquet", "directions.parquet"), ("directions_csv", "directions.csv")):
        cands.append(base / name)
        rec = outputs.get(key)
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
    raise MissingDirectionDataError(f"no Stage-3 direction table found for {manifest_path}")


def read_directions(path: Path):
    """Read a Stage-3 table and check the columns Stage 4 needs."""
    from tart_gnss.stage3.output import read_directions as read3

    try:
        df = read3(path)
    except Exception as exc:  # noqa: BLE001 - any reader failure is a data error
        raise MissingDirectionDataError(f"cannot read Stage-3 directions {path}: {exc}") from exc
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise MissingDirectionDataError(f"{path}: missing required columns {missing}")
    df["station_id"] = df["station_id"].astype(str)
    df["fit_success"] = df["fit_success"].map(lambda v: bool(v) if v == v and v is not None else False)
    return df


def _stage2_cross_check(manifest: Dict[str, object], manifest_path: Path) -> Tuple[Optional[Dict], str]:
    rec = manifest.get("input_manifest")
    if not rec:
        return None, "no Stage-2 manifest recorded"
    cands = [Path(rec), manifest_path.parent / Path(rec)]
    found = next((c for c in cands if c.is_file()), None)
    if found is None:
        return None, f"Stage-2 manifest {rec} not reachable; coordinates taken from the Stage-3 manifest only"
    want = manifest.get("input_manifest_sha256")
    got = sha256_file(found)
    if want and got != want:
        return None, f"Stage-2 manifest {found} changed since Stage 3 (SHA-256 mismatch); not used"
    try:
        return json.loads(found.read_text(encoding="utf-8")), f"verified against Stage-2 manifest {found} (sha256 {got})"
    except ValueError:
        return None, f"Stage-2 manifest {found} unreadable; not used"


def station_coordinates(manifest: Dict[str, object], manifest_path: Path) -> Dict[str, StationCoordinates]:
    """Frozen station coordinates from the Stage-3 manifest (cross-checked against Stage 2).

    Raises
    ------
    InvalidStage3ManifestError
        Missing coordinates, or disagreement with the hash-verified Stage-2 manifest.
    """
    stage2, check = _stage2_cross_check(manifest, manifest_path)
    out = {}
    for sid, st in manifest["stations"].items():
        try:
            lat, lon = float(st["latitude_deg"]), float(st["longitude_deg"])
            alt = float(st.get("altitude_m") or 0.0)
        except (KeyError, TypeError, ValueError) as exc:
            raise InvalidStage3ManifestError(f"station {sid}: missing or invalid coordinates") from exc
        if not (np.isfinite(lat) and np.isfinite(lon) and np.isfinite(alt)) or abs(lat) > 90:
            raise InvalidStage3ManifestError(f"station {sid}: invalid coordinates")
        station_check = check
        if stage2 is not None:
            s2 = (stage2.get("stations") or {}).get(sid)
            if s2 is None:
                station_check = f"station absent from Stage-2 manifest ({check})"
            else:
                d = (abs(float(s2["latitude_deg"]) - lat), abs(float(s2["longitude_deg"]) - lon),
                     abs(float(s2.get("altitude_m") or 0.0) - alt))
                if d[0] > COORD_TOL_DEG or d[1] > COORD_TOL_DEG or d[2] > COORD_TOL_M:
                    raise InvalidStage3ManifestError(
                        f"station {sid}: Stage-3 coordinates disagree with the Stage-2 manifest")
        warn = list(((st.get("calibration") or {}).get("warnings")) or [])
        out[sid] = StationCoordinates(
            station_id=sid, latitude_deg=lat, longitude_deg=lon, altitude_m=alt,
            ecef_m=geodetic_to_ecef(lat, lon, alt),
            source=f"{manifest_path}: stations.{sid} (frozen Stage-1/2 run provenance)",
            cross_check=station_check, calibration_warnings=warn,
        )
    return out


def load_stage3_run(manifest_path, stations: Optional[Sequence[str]] = None) -> Stage3Run:
    """Load the Stage-3 manifest, direction table and frozen station coordinates.

    ``stations`` restricts the run to those station ids (unknown ids raise
    :class:`InvalidStage3ManifestError`).
    """
    manifest_path = Path(manifest_path)
    manifest, sha = load_stage3_manifest(manifest_path)
    coords = station_coordinates(manifest, manifest_path)
    if stations:
        unknown = [s for s in stations if s not in coords]
        if unknown:
            raise InvalidStage3ManifestError(f"unknown station(s) {unknown}; available: {sorted(coords)}")
        coords = {s: coords[s] for s in stations}
    dpath = resolve_directions_path(manifest, manifest_path)
    df = read_directions(dpath)
    unknown_rows = sorted(set(df["station_id"]) - set(manifest["stations"]))
    if unknown_rows:
        raise MissingDirectionDataError(f"direction rows for stations not in the manifest: {unknown_rows}")
    df = df[df["station_id"].isin(list(coords))].reset_index(drop=True)
    norad = manifest.get("norad_id")
    return Stage3Run(manifest_path=manifest_path, manifest=manifest, sha256=sha, directions_path=dpath,
                     directions_sha256=sha256_file(dpath), directions=df, stations=coords,
                     norad_id=int(norad) if norad is not None else None)


def sight_vector(row, coords: StationCoordinates) -> np.ndarray:
    """ECEF unit sight vector of one Stage-3 row (``(l, m)`` preferred, else Az/El)."""
    l, m = float(row["l"]), float(row["m"])
    if np.isfinite(l) and np.isfinite(m):
        return lm_to_ecef_unit_vector(l, m, coords.latitude_deg, coords.longitude_deg)
    return azel_to_ecef_unit_vector(float(row["az_obs_deg"]), float(row["el_obs_deg"]), coords.latitude_deg,
                                    coords.longitude_deg)


def sight_line(row, coords: StationCoordinates) -> Tuple[np.ndarray, np.ndarray, str]:
    """ECEF unit sight vector and 3x3 tangent covariance (rad^2) for one Stage-3 row.

    Uses ``(l, m)`` and its covariance when finite (well defined at the
    zenith); otherwise Az/El and ``cov_az_el``. Returns ``(s, cov3, source)``.

    Raises
    ------
    TriangulationError
        No finite direction.
    CovarianceError
        No finite, positive-definite covariance.
    """
    lat, lon = coords.latitude_deg, coords.longitude_deg
    l, m = float(row["l"]), float(row["m"])
    if np.isfinite(l) and np.isfinite(m):
        s = lm_to_ecef_unit_vector(l, m, lat, lon)
        sl, sm, c = float(row["sigma_l"]), float(row["sigma_m"]), float(row["cov_lm"])
        if np.isfinite(sl) and np.isfinite(sm) and np.isfinite(c):
            return s, lm_covariance_to_ecef(l, m, [[sl * sl, c], [c, sm * sm]], lat, lon), "lm"
    az, el = float(row["az_obs_deg"]), float(row["el_obs_deg"])
    if not (np.isfinite(az) and np.isfinite(el)):
        raise TriangulationError(f"{coords.station_id}: no finite direction")
    s = azel_to_ecef_unit_vector(az, el, lat, lon)
    sa, se, c = float(row["sigma_az_deg"]), float(row["sigma_el_deg"]), float(row["cov_az_el"])
    return s, azel_covariance_to_ecef(az, el, [[sa * sa, c], [c, se * se]], lat, lon), "azel"


def apply_quality_policy(df, policy: str):
    """Split Stage-3 rows into kept rows and a ``reason -> count`` map of rejections.

    Rows without ``fit_success`` or without a finite direction are rejected
    (``FIT_FAILED`` / ``NO_DIRECTION``) whatever their flags say.
    """
    reasons = []
    for _, r in df.iterrows():
        why = stage3_rejection_reason(r["quality_flags"], policy)
        if not why and not r["fit_success"]:
            why = "FIT_FAILED"
        if not why and not ((np.isfinite(r["l"]) and np.isfinite(r["m"]))
                            or (np.isfinite(r["az_obs_deg"]) and np.isfinite(r["el_obs_deg"]))):
            why = "NO_DIRECTION"
        reasons.append(why)
    reasons = np.array(reasons, dtype=object)
    kept = df[reasons == ""]
    counts: Dict[str, int] = {}
    for why in reasons[reasons != ""]:
        counts[why] = counts.get(why, 0) + 1
    return kept, counts
