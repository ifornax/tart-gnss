"""Load a Stage-6 run and its upstream provenance chain (Stages 5, 4, 3, 2).

Every link follows a path recorded in the downstream manifest, resolved as
recorded, then relative to that run's ``provenance.working_directory``, then
next to the downstream manifest. Where the downstream manifest recorded the
upstream SHA-256, it is verified. Each link is reported as ``verified``,
``unverified`` (no hash recorded), ``mismatch`` or ``missing``. Unrelated
files are never discovered heuristically.

Policy: under ``strict`` a hash mismatch is fatal (:class:`ProvenanceChainError`);
under ``permissive`` it is recorded and the data are used. Stages 5 and 4 are
required; a missing Stage 3 disables the re-triangulation analyses and a
missing Stage 2 is recorded as ``CALIBRATION_METADATA_MISSING``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from tart_gnss.stage4.loader import StationCoordinates, sha256_file
from tart_gnss.stage7.errors import InvalidStage6ManifestError, ProvenanceChainError


@dataclass
class StageLink:
    """One provenance link."""

    stage: int
    recorded: Optional[str]
    path: Optional[Path]
    sha256: Optional[str]
    expected_sha256: Optional[str]
    status: str

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        return {"stage": self.stage, "recorded": self.recorded, "path": str(self.path) if self.path else None,
                "sha256": self.sha256, "expected_sha256": self.expected_sha256, "status": self.status}


@dataclass
class DiagnosticInputs:
    """All loaded stage outputs and their provenance."""

    stage6_manifest_path: Path
    manifests: Dict[int, Dict[str, object]]
    links: List[StageLink]
    norad_id: Optional[int]
    residuals: object
    tle_trajectory: object
    summary6: Dict[str, object]
    trajectory: object
    innovations: Optional[object]
    positions: object
    contributions: Optional[object]
    directions: Optional[object]
    stations: Dict[str, StationCoordinates]
    stage2: Optional[Dict[str, object]]
    table_hashes: Dict[str, str] = field(default_factory=dict)

    @property
    def calibration_metadata_available(self) -> bool:
        """A Stage-2 manifest with per-file calibration was loaded."""
        return self.stage2 is not None

    def link(self, stage: int) -> Optional[StageLink]:
        """The link to ``stage`` (None if not attempted)."""
        return next((lk for lk in self.links if lk.stage == stage), None)


def _read_json(path: Path) -> Dict[str, object]:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ProvenanceChainError(f"cannot read {path}: {exc}") from exc


def resolve(recorded: Optional[str], downstream: Path, working_directory: Optional[str]) -> Optional[Path]:
    """Resolve a recorded path (as recorded, relative to the run's working directory, next to the manifest)."""
    if not recorded:
        return None
    cands = [Path(recorded)]
    if working_directory and not Path(recorded).is_absolute():
        cands.append(Path(working_directory) / recorded)
    cands.append(Path(downstream).parent / Path(recorded).name)
    return next((c for c in cands if c.is_file()), None)


def _link(stage: int, recorded, downstream_path: Path, downstream: Dict[str, object], expected: Optional[str],
          policy: str) -> StageLink:
    wd = (downstream.get("provenance") or {}).get("working_directory")
    p = resolve(recorded, downstream_path, wd)
    if p is None:
        return StageLink(stage, recorded, None, None, expected, "missing")
    got = sha256_file(p)
    status = "unverified" if not expected else ("verified" if got == expected else "mismatch")
    if status == "mismatch" and policy == "strict":
        raise ProvenanceChainError(f"Stage-{stage} manifest {p} does not match the SHA-256 recorded downstream")
    return StageLink(stage, recorded, p, got, expected, status)


def _table(manifest: Dict[str, object], manifest_path: Path, stem: str, required: bool = True):
    from tart_gnss.stage6.output import read_table

    outputs = manifest.get("outputs") or {}
    base = Path(manifest_path).parent
    cands = [base / f"{stem}.parquet", base / f"{stem}.csv"]
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
            df = read_table(c)
            for col in ("quality_flags", "stage3_quality_flags", "stage4_station_flags", "stage4_quality_flags",
                        "stage5_quality_flags"):
                if col in df.columns:
                    df[col] = df[col].map(_as_list)
            for col in ("reference_time_utc", "station_time_utc"):
                if col in df.columns:
                    import pandas as pd

                    df[col] = pd.to_datetime(df[col], utc=True)
            return df, c
    if required:
        raise ProvenanceChainError(f"no {stem} table found for {manifest_path}")
    return None, None


def _as_list(f):
    if isinstance(f, list):
        return f
    if f is None or (isinstance(f, float) and f != f):
        return []
    if isinstance(f, str):
        return [x for x in f.split(";") if x]
    return list(f)


def station_coordinates(manifest4: Dict[str, object]) -> Dict[str, StationCoordinates]:
    """Frozen station coordinates recorded by Stage 4."""
    out = {}
    for sid, st in (manifest4.get("stations") or {}).items():
        try:
            out[sid] = StationCoordinates(
                station_id=sid, latitude_deg=float(st["latitude_deg"]), longitude_deg=float(st["longitude_deg"]),
                altitude_m=float(st.get("altitude_m") or 0.0), ecef_m=np.asarray(st["ecef_m"], float),
                source="stage4_manifest.stations", cross_check=str(st.get("cross_check", "")),
                calibration_warnings=list(st.get("calibration_warnings") or []))
        except (KeyError, TypeError, ValueError) as exc:
            raise ProvenanceChainError(f"Stage-4 manifest: invalid coordinates for station {sid}") from exc
    return out


def load_diagnostic_chain(stage6_manifest, policy: str = "strict") -> DiagnosticInputs:
    """Load and verify linked outputs from Stages 2-6."""
    p6 = Path(stage6_manifest)
    if not p6.is_file():
        raise InvalidStage6ManifestError(f"Stage-6 manifest not found: {p6}")
    try:
        m6 = json.loads(p6.read_text(encoding="utf-8"))
    except ValueError as exc:
        raise InvalidStage6ManifestError(f"cannot read {p6}: {exc}") from exc
    if not isinstance(m6, dict) or m6.get("stage") != 6 or m6.get("schema_version") != 1:
        raise InvalidStage6ManifestError(f"{p6} is not a schema-1 Stage-6 manifest")
    links = [StageLink(6, str(p6), p6, sha256_file(p6), None, "input")]
    residuals, rpath = _table(m6, p6, "residuals")
    tle_traj, tpath = _table(m6, p6, "tle_trajectory")
    summary_path = p6.parent / "summary.json"
    summary6 = _read_json(summary_path) if summary_path.is_file() else {}

    l5 = _link(5, m6.get("input_stage5_manifest"), p6, m6, m6.get("input_stage5_sha256"), policy)
    links.append(l5)
    if l5.path is None:
        raise ProvenanceChainError(f"Stage-5 manifest {l5.recorded} not reachable")
    m5 = _read_json(l5.path)
    trajectory, _ = _table(m5, l5.path, "trajectory")
    innovations, _ = _table(m5, l5.path, "innovations", required=False)

    l4 = _link(4, m5.get("input_stage4_manifest"), l5.path, m5, m5.get("input_sha256"), policy)
    links.append(l4)
    if l4.path is None:
        raise ProvenanceChainError(f"Stage-4 manifest {l4.recorded} not reachable")
    m4 = _read_json(l4.path)
    positions, _ = _table(m4, l4.path, "positions")
    contributions, _ = _table(m4, l4.path, "station_contributions", required=False)

    l3 = _link(3, m4.get("input_stage3_manifest"), l4.path, m4, m4.get("input_sha256"), policy)
    links.append(l3)
    m3, directions = None, None
    if l3.path is not None:
        m3 = _read_json(l3.path)
        directions, _ = _table(m3, l3.path, "directions", required=False)
    stage2 = None
    m2 = None
    if m3 is not None:
        l2 = _link(2, m3.get("input_manifest"), l3.path, m3, m3.get("input_manifest_sha256"), policy)
        links.append(l2)
        if l2.path is not None:
            m2 = _read_json(l2.path)
            stage2 = m2
    manifests = {6: m6, 5: m5, 4: m4}
    if m3 is not None:
        manifests[3] = m3
    if m2 is not None:
        manifests[2] = m2
    norad = m6.get("norad_id")
    return DiagnosticInputs(
        stage6_manifest_path=p6, manifests=manifests, links=links, norad_id=int(norad) if norad is not None else None,
        residuals=residuals, tle_trajectory=tle_traj, summary6=summary6, trajectory=trajectory,
        innovations=innovations, positions=positions, contributions=contributions, directions=directions,
        stations=station_coordinates(m4), stage2=stage2,
        table_hashes={"residuals": sha256_file(rpath), "tle_trajectory": sha256_file(tpath)},
    )


def _parse_utc(text: str) -> float:
    return datetime.fromisoformat(str(text).replace("Z", "+00:00")).timestamp()


def stage2_calibration_table(stage2: Optional[Dict[str, object]]) -> Dict[str, List[Dict[str, object]]]:
    """Per-station native-file calibration records sorted by time.

    Each record: ``start_s``, ``end_s``, ``zero_gain_count``, ``fingerprint``,
    ``solution_index`` (0 for the first gain solution seen, +1 at each change).
    """
    out: Dict[str, List[Dict[str, object]]] = {}
    if not stage2:
        return out
    for sid, st in (stage2.get("stations") or {}).items():
        files = ((st.get("observation") or {}).get("files")) or []
        recs = []
        for f in files:
            try:
                recs.append({"start_s": _parse_utc(f["first_timestamp_utc"]), "end_s": _parse_utc(f["last_timestamp_utc"]),
                             "zero_gain_count": len(f.get("zero_gain_antennas") or []),
                             "fingerprint": f.get("calibration_fingerprint")})
            except (KeyError, ValueError):
                continue
        recs.sort(key=lambda r: r["start_s"])
        idx, last = -1, object()
        for r in recs:
            if r["fingerprint"] != last:
                idx += 1
                last = r["fingerprint"]
            r["solution_index"] = idx
        out[sid] = recs
    return out


def calibration_at(table: Dict[str, List[Dict[str, object]]], station_id: str, t_s: float) -> Optional[Dict[str, object]]:
    """The native-file calibration record covering ``t_s`` (nearest file within 5 s of its span)."""
    recs = table.get(station_id) or []
    best, best_d = None, np.inf
    for r in recs:
        d = 0.0 if r["start_s"] <= t_s <= r["end_s"] else min(abs(t_s - r["start_s"]), abs(t_s - r["end_s"]))
        if d < best_d:
            best, best_d = r, d
    return best if best_d <= 5.0 else None
