"""Selection of the reference TLE (and an optional ensemble) with provenance.

Primary reference: exactly one of

``--tle-file PATH``
    A TLE text file (2- or 3-line sets). Every set of the target NORAD ID is
    read (unlike the Stage-1 catalogue loader, duplicates are kept: a history
    file may hold several epochs). If several match, the one whose epoch is
    nearest the midpoint of the Stage-5 trajectory is selected; the rule and
    the candidates are recorded. A ``.meta.json`` sidecar, when present, must
    carry the file's SHA-256.
``--tle-source stage1``
    The exact element set the pipeline used in Stage 1, recovered through the
    run provenance: Stage-5 manifest -> ``input_stage4_manifest`` ->
    ``satellite_reference.satellite`` (copied unchanged from Stages 1-3).
``--tle-source download``
    An explicit CelesTrak download, archived with its provenance sidecar
    before use. CelesTrak serves only *current* elements, so this is never
    appropriate for historical observations; the age checks flag it.

Nothing is downloaded unless ``download`` is requested explicitly.

Ensemble (optional): every distinct set (by line content) of the same NORAD
ID found in ``--tle-ensemble-files`` / ``--tle-ensemble-dir`` (``*.tle``,
``*.txt``) whose epoch is within ``--tle-ensemble-window-days`` of the
trajectory midpoint. The primary set is included when it falls in the window.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from tart_gnss.stage4.loader import sha256_file
from tart_gnss.stage6.errors import TLENoradMismatchError, TLEReferenceError
from tart_gnss.tle import TLEError, format_utc, parse_tle_text, tle_epoch_utc

SOURCES = ("stage1", "download")


@dataclass
class TLEReference:
    """One element set with its provenance."""

    norad_id: int
    name: str
    line1: str
    line2: str
    epoch_utc: datetime
    source_kind: str
    path: Optional[str] = None
    sha256: Optional[str] = None
    source: str = ""
    selection_rule: str = ""
    provenance: Dict[str, object] = field(default_factory=dict)

    @property
    def epoch_unix_s(self) -> float:
        """TLE epoch as Unix seconds."""
        return self.epoch_utc.timestamp()

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable description."""
        return {
            "norad_id": self.norad_id, "name": self.name, "line1": self.line1, "line2": self.line2,
            "epoch_utc": format_utc(self.epoch_utc, 6), "source_kind": self.source_kind, "file": self.path,
            "sha256": self.sha256, "source": self.source, "selection_rule": self.selection_rule,
            "provenance": self.provenance,
        }


def _parse(text: str, where: str) -> List[Tuple[str, str, str]]:
    try:
        return parse_tle_text(text)
    except TLEError as exc:
        raise TLEReferenceError(f"{where}: {exc}") from exc


def _sidecar_source(path: Path, digest: str) -> str:
    meta = path.with_name(path.name + ".meta.json")
    if not meta.is_file():
        return f"local file {path}"
    doc = json.loads(meta.read_text(encoding="utf-8"))
    if doc.get("sha256") and doc["sha256"] != digest:
        raise TLEReferenceError(f"SHA-256 of {path} does not match its metadata {meta}; the file was modified")
    return str(doc.get("source", f"local file {path}"))


def _records(path: Path) -> Tuple[List[Tuple[str, str, str]], str, str]:
    if not path.is_file():
        raise TLEReferenceError(f"TLE file not found: {path}")
    digest = sha256_file(path)
    triples = _parse(path.read_text(encoding="utf-8"), str(path))
    if not triples:
        raise TLEReferenceError(f"no TLEs in {path}")
    return triples, digest, _sidecar_source(path, digest)


def load_tle_reference(norad_id: int, tle_file, reference_time_s: Optional[float] = None) -> TLEReference:
    """Load the selected NORAD element set from a TLE file with provenance.

    Raises
    ------
    TLENoradMismatchError
        The file holds a single element set, of another object.
    TLEReferenceError
        File missing/modified/malformed, or no set for ``norad_id``.
    """
    path = Path(tle_file)
    triples, digest, source = _records(path)
    cands = [(n, l1, l2) for n, l1, l2 in triples if int(l1[2:7]) == int(norad_id)]
    if not cands:
        if len(triples) == 1:
            raise TLENoradMismatchError(f"{path} holds NORAD {int(triples[0][1][2:7])}, not {norad_id}")
        raise TLEReferenceError(f"TLE_NOT_FOUND: no element set for NORAD {norad_id} in {path} "
                                f"({len(triples)} sets)")
    epochs = [tle_epoch_utc(l1, l2) for _, l1, l2 in cands]
    if len(cands) == 1 or reference_time_s is None:
        k = 0
        rule = "only element set for this NORAD ID in the file" if len(cands) == 1 else "first in file"
    else:
        k = min(range(len(cands)), key=lambda i: abs(epochs[i].timestamp() - reference_time_s))
        rule = f"epoch nearest the Stage-5 trajectory midpoint among {len(cands)} sets for this NORAD ID"
    name, l1, l2 = cands[k]
    return TLEReference(
        norad_id=int(norad_id), name=name, line1=l1, line2=l2, epoch_utc=epochs[k], source_kind="tle-file",
        path=str(path), sha256=digest, source=source, selection_rule=rule,
        provenance={"candidates": [format_utc(e, 6) for e in epochs]},
    )


def _resolve(recorded: str, manifest_path: Path, working_directory: Optional[str]) -> Optional[Path]:
    cands = [Path(recorded)]
    if working_directory and not Path(recorded).is_absolute():
        cands.append(Path(working_directory) / recorded)
    cands.append(manifest_path.parent / Path(recorded).name)
    return next((c for c in cands if c.is_file()), None)


def stage1_tle_reference(stage5_manifest_path: Path, stage5_manifest: Dict[str, object],
                         norad_id: int) -> TLEReference:
    """The Stage-1 element set recovered from the run provenance (Stage 5 -> Stage 4 manifest)."""
    rec = stage5_manifest.get("input_stage4_manifest")
    if not rec:
        raise TLEReferenceError("Stage-5 manifest records no input_stage4_manifest; use --tle-file")
    wd = (stage5_manifest.get("provenance") or {}).get("working_directory")
    p4 = _resolve(str(rec), Path(stage5_manifest_path), wd)
    if p4 is None:
        raise TLEReferenceError(f"Stage-4 manifest {rec} not reachable; use --tle-file")
    try:
        m4 = json.loads(p4.read_text(encoding="utf-8"))
    except ValueError as exc:
        raise TLEReferenceError(f"cannot read Stage-4 manifest {p4}: {exc}") from exc
    sat = ((m4.get("satellite_reference") or {}).get("satellite")) or {}
    l1, l2 = sat.get("tle_line1"), sat.get("tle_line2")
    if not l1 or not l2:
        raise TLEReferenceError(f"{p4}: no Stage-1 TLE recorded in satellite_reference")
    _parse(f"{l1}\n{l2}\n", f"Stage-1 TLE in {p4}")
    if int(l1[2:7]) != int(norad_id):
        raise TLENoradMismatchError(f"Stage-1 TLE in {p4} is NORAD {int(l1[2:7])}, not {norad_id}")
    return TLEReference(
        norad_id=int(norad_id), name=str(sat.get("satellite_name") or f"NORAD {norad_id}"), line1=l1, line2=l2,
        epoch_utc=tle_epoch_utc(l1, l2), source_kind="stage1", path=None, sha256=None,
        source=str(sat.get("tle_source") or ""),
        selection_rule="exact Stage-1 element set from the run provenance",
        provenance={"chain": [str(stage5_manifest_path), str(p4), m4.get("input_stage3_manifest")],
                    "stage4_manifest_sha256": sha256_file(p4)},
    )


def download_tle_reference(norad_id: int, out_dir: Path, reference_time_s: Optional[float], fetch=None
                           ) -> TLEReference:
    """Explicit CelesTrak download (archived with a provenance sidecar), then :func:`load_tle_reference`."""
    from tart_gnss.tle import download_gnss_catalogue

    kw = {"fetch": fetch} if fetch is not None else {}
    try:
        path = download_gnss_catalogue(out_dir, **kw)
    except (TLEError, OSError) as exc:
        raise TLEReferenceError(f"TLE download failed: {exc}") from exc
    ref = load_tle_reference(norad_id, path, reference_time_s)
    ref.source_kind = "download"
    ref.provenance["downloaded"] = True
    ref.provenance["retrieved_utc"] = format_utc(datetime.now(timezone.utc))
    return ref


def load_ensemble(norad_id: int, files: Sequence = (), dirs: Sequence = (), centre_time_s: Optional[float] = None,
                  window_days: float = 7.0) -> Tuple[List[TLEReference], Dict[str, object]]:
    """Distinct element sets of ``norad_id`` within ``window_days`` of ``centre_time_s``."""
    paths: List[Path] = [Path(f) for f in files]
    for d in dirs:
        d = Path(d)
        if not d.is_dir():
            raise TLEReferenceError(f"TLE ensemble directory not found: {d}")
        paths += sorted(p for p in d.iterdir() if p.suffix in (".tle", ".txt") and p.is_file())
    seen, out = set(), []
    n_found = n_outside = 0
    for p in paths:
        triples, digest, source = _records(p)
        for name, l1, l2 in triples:
            if int(l1[2:7]) != int(norad_id):
                continue
            n_found += 1
            key = (l1.strip(), l2.strip())
            if key in seen:
                continue
            seen.add(key)
            ep = tle_epoch_utc(l1, l2)
            if centre_time_s is not None and abs(ep.timestamp() - centre_time_s) > window_days * 86400.0:
                n_outside += 1
                continue
            out.append(TLEReference(norad_id=int(norad_id), name=name, line1=l1, line2=l2, epoch_utc=ep,
                                    source_kind="ensemble", path=str(p), sha256=digest, source=source,
                                    selection_rule=f"within {window_days:g} d of the trajectory midpoint"))
    out.sort(key=lambda r: r.epoch_utc)
    info = {"files": [str(p) for p in paths], "n_sets_for_norad": n_found, "n_distinct_in_window": len(out),
            "n_outside_window": n_outside, "window_days": window_days,
            "deduplication": "by exact line content"}
    return out, info
