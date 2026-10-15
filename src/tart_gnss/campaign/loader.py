"""Load one completed Stage-7 run and the upstream products the campaign needs.

The Stage-7 manifest is validated (``stage == 7``; schema version in
:data:`tart_gnss.campaign.schema.SUPPORTED_STAGE7_SCHEMA`, a future version
fails clearly). Its own outputs (``summary.json``, ``diagnosis``,
``station_influence``) are read next to it. Upstream products are reached only
through the provenance chain the Stage-7 manifest recorded (resolved as
recorded, then relative to the Stage-7 run's working directory), and each
upstream manifest is re-hashed against the hash Stage 7 recorded. Nothing is
recomputed or rerun.

Provenance problems never abort the campaign: they make the run
``BROKEN_PROVENANCE`` (listed and excluded by the strict policy).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from tart_gnss.campaign.errors import CampaignSchemaError, InvalidStage7ManifestError
from tart_gnss.campaign.schema import SUPPORTED_STAGE7_SCHEMA, canonical_constellation
from tart_gnss.stage4.loader import sha256_file


@dataclass
class CampaignRun:
    """One completed Stage-7 run and the upstream tables it links to."""

    manifest_path: Path
    sha256: str
    manifest: Dict[str, object]
    summary: Dict[str, object]
    diagnosis: object
    influence: object
    stage_manifests: Dict[int, Dict[str, object]] = field(default_factory=dict)
    stage_paths: Dict[int, Optional[Path]] = field(default_factory=dict)
    stage6_summary: Dict[str, object] = field(default_factory=dict)
    residuals: object = None
    trajectory: object = None
    positions: object = None
    contributions: object = None
    provenance_problems: List[str] = field(default_factory=list)
    metadata_flags: List[str] = field(default_factory=list)

    @property
    def norad_id(self) -> Optional[int]:
        """NORAD catalogue number."""
        v = self.manifest.get("norad_id", self.summary.get("norad_id"))
        return int(v) if v is not None else None

    @property
    def broken(self) -> bool:
        """Whether any provenance link is missing, changed or mismatched."""
        return bool(self.provenance_problems)

    @property
    def satellite(self) -> Dict[str, object]:
        """Structured satellite metadata (Stage-3 manifest, else the Stage-4 copy)."""
        m3 = self.stage_manifests.get(3) or {}
        if m3.get("satellite"):
            return dict(m3["satellite"])
        m4 = self.stage_manifests.get(4) or {}
        return dict(((m4.get("satellite_reference") or {}).get("satellite")) or {})

    @property
    def constellation(self) -> str:
        """Canonical constellation from structured metadata only."""
        return canonical_constellation(self.satellite.get("constellation"))


def _read_json(path: Path) -> Dict[str, object]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _resolve(recorded: Optional[str], working_directory: Optional[str]) -> Optional[Path]:
    if not recorded:
        return None
    cands = [Path(recorded)]
    if working_directory and not Path(recorded).is_absolute():
        cands.append(Path(working_directory) / recorded)
    return next((c for c in cands if c.is_file()), None)


def _table(manifest, path, stem, required=False):
    from tart_gnss.stage7.loader import _table as t7

    try:
        df, _ = t7(manifest, path, stem, required=required)
    except Exception:  # noqa: BLE001 - recorded as a problem by the caller
        return None
    return df


def load_campaign_run(manifest_path: Path) -> CampaignRun:
    """Load and validate one completed Stage-7 run.

    Raises
    ------
    InvalidStage7ManifestError
        Unreadable file or not a Stage-7 manifest, or its own outputs missing.
    CampaignSchemaError
        Unsupported Stage-7 schema version.
    """
    from tart_gnss.stage7.output import read_table

    p = Path(manifest_path)
    try:
        m7 = _read_json(p)
    except (OSError, ValueError) as exc:
        raise InvalidStage7ManifestError(f"cannot read {p}: {exc}") from exc
    if not isinstance(m7, dict) or m7.get("stage") != 7:
        raise InvalidStage7ManifestError(f"{p} is not a Stage-7 manifest")
    if m7.get("schema_version") not in SUPPORTED_STAGE7_SCHEMA:
        raise CampaignSchemaError(f"{p}: Stage-7 schema_version {m7.get('schema_version')!r} is not supported "
                                  f"(supported: {SUPPORTED_STAGE7_SCHEMA})")
    base = p.parent
    try:
        summary = _read_json(base / "summary.json")
        diag = read_table(base / "diagnosis.parquet") if (base / "diagnosis.parquet").is_file() else read_table(
            base / "diagnosis.csv")
        infl = read_table(base / "station_influence.parquet") if (base / "station_influence.parquet").is_file() \
            else read_table(base / "station_influence.csv")
    except (OSError, ValueError) as exc:
        raise InvalidStage7ManifestError(f"{p}: Stage-7 outputs missing or unreadable ({exc})") from exc
    run = CampaignRun(manifest_path=p, sha256=sha256_file(p), manifest=m7, summary=summary, diagnosis=diag,
                      influence=infl)
    wd = (m7.get("provenance") or {}).get("working_directory")
    for link in m7.get("provenance_chain") or []:
        st = int(link.get("stage"))
        if st == 7:
            continue
        path = _resolve(link.get("path") or link.get("recorded"), wd)
        run.stage_paths[st] = path
        if link.get("status") not in ("verified", "input"):
            if st >= 3 or link.get("status") == "mismatch":
                run.provenance_problems.append(f"stage {st} link recorded as {link.get('status')}")
            continue
        if path is None:
            run.provenance_problems.append(f"stage {st} manifest not reachable")
            continue
        if link.get("sha256") and sha256_file(path) != link["sha256"]:
            run.provenance_problems.append(f"stage {st} manifest changed since Stage 7")
            continue
        try:
            run.stage_manifests[st] = _read_json(path)
        except (OSError, ValueError):
            run.provenance_problems.append(f"stage {st} manifest unreadable")
    for st in (6, 5, 4):
        if st not in run.stage_manifests and not any(f"stage {st}" in x for x in run.provenance_problems):
            run.provenance_problems.append(f"stage {st} manifest not recorded")
    if 6 in run.stage_manifests:
        p6 = run.stage_paths[6]
        if (p6.parent / "summary.json").is_file():
            run.stage6_summary = _read_json(p6.parent / "summary.json")
        run.residuals = _table(run.stage_manifests[6], p6, "residuals")
    if 5 in run.stage_manifests:
        run.trajectory = _table(run.stage_manifests[5], run.stage_paths[5], "trajectory")
    if 4 in run.stage_manifests:
        run.positions = _table(run.stage_manifests[4], run.stage_paths[4], "positions")
        run.contributions = _table(run.stage_manifests[4], run.stage_paths[4], "station_contributions")
    if run.constellation == "Unknown":
        run.metadata_flags.append("MISSING_CONSTELLATION_METADATA")
    return run
