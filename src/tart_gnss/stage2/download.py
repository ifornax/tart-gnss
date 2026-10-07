"""Downloading native files, with integrity checks and caching.

Each downloaded file gets a provenance sidecar ``<file>.prov.json`` recording
the remote identifier, URL, size, ETag, retrieval time and SHA-256. An existing
local file is reused only if its sidecar matches the archive record and its
current SHA-256 matches the sidecar. Otherwise the run stops. Files are never
silently overwritten.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from tart_gnss.stage2.archive import ArchiveBackend, ObservationRecord
from tart_gnss.stage2.errors import ArchiveUnavailableError, DownloadError
from tart_gnss.tle import format_utc, sha256_file

logger = logging.getLogger(__name__)

_MD5_ETAG = re.compile(r"^[0-9a-f]{32}$")


@dataclass
class RawFile:
    """A native file on local disk, with its provenance.

    Attributes
    ----------
    sha256, size_bytes
        Integrity of the local file.
    md5_matches_etag
        True/False if the S3 ETag is a plain MD5 and was compared, ``None``
        if no comparison was possible (e.g. multipart ETag).
    reused
        True if an existing, verified local copy was used.
    """

    original_filename: str
    local_path: str
    remote_identifier: str
    archive_url: str
    retrieval_utc: str
    sha256: str
    size_bytes: int
    etag: Optional[str]
    md5_matches_etag: Optional[bool]
    reused: bool = False

    def to_dict(self) -> dict:
        """JSON-serialisable representation."""
        return asdict(self)


def _md5(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 (S3 ETag comparison, not security)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sidecar(path: Path) -> Path:
    return path.with_name(path.name + ".prov.json")


def _verify_against_record(path: Path, record: ObservationRecord) -> Optional[bool]:
    size = path.stat().st_size
    if record.size_bytes is not None and size != record.size_bytes:
        raise DownloadError(
            f"{path}: size {size} B differs from archive size {record.size_bytes} B"
        )
    if record.etag and _MD5_ETAG.match(record.etag):
        if _md5(path) != record.etag:
            raise DownloadError(f"{path}: MD5 does not match the archive ETag {record.etag}")
        return True
    return None


def download_observation(
    record: ObservationRecord, raw_dir: Path, archive: ArchiveBackend
) -> RawFile:
    """Retrieve one archive file into ``raw_dir``, or reuse a verified local copy.

    Parameters
    ----------
    record
        Archive record from :func:`tart_gnss.stage2.archive.find_observations`.
    raw_dir
        Destination directory. The original filename is kept.
    archive
        Archive backend used for the download.

    Returns
    -------
    RawFile

    Raises
    ------
    DownloadError
        If the transfer fails, the file fails integrity checks, or an existing
        local file cannot be shown to be the same archive object.
    """
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    local = raw_dir / record.filename
    sidecar = _sidecar(local)

    if local.exists():
        if not sidecar.is_file():
            raise DownloadError(
                f"{local} exists without provenance; refusing to overwrite or reuse it"
            )
        prov = json.loads(sidecar.read_text(encoding="utf-8"))
        if prov.get("remote_identifier") != record.remote_identifier:
            raise DownloadError(f"{local} was retrieved from a different archive object")
        if sha256_file(local) != prov.get("sha256"):
            raise DownloadError(f"{local} no longer matches its recorded SHA-256")
        md5_ok = _verify_against_record(local, record)
        logger.info("Reusing verified local file %s (sha256 %s)", local, prov["sha256"])
        prov["reused"] = True
        prov["md5_matches_etag"] = md5_ok
        return RawFile(**{k: prov[k] for k in RawFile.__dataclass_fields__})

    partial = local.with_name(local.name + ".part")
    if partial.exists():
        partial.unlink()  # leftover from an interrupted transfer of this same object
    retrieval = datetime.now(timezone.utc)
    try:
        archive.download(record.remote_identifier, partial)
    except ArchiveUnavailableError as exc:
        if partial.exists():
            partial.unlink()
        raise DownloadError(str(exc)) from exc
    if not partial.is_file():
        raise DownloadError(f"Download of {record.remote_identifier} produced no file")
    try:
        md5_ok = _verify_against_record(partial, record)
    except DownloadError:
        partial.unlink()
        raise
    os.replace(partial, local)
    raw = RawFile(
        original_filename=record.filename,
        local_path=str(local),
        remote_identifier=record.remote_identifier,
        archive_url=record.archive_url,
        retrieval_utc=format_utc(retrieval, 3),
        sha256=sha256_file(local),
        size_bytes=local.stat().st_size,
        etag=record.etag,
        md5_matches_etag=md5_ok,
    )
    sidecar.write_text(json.dumps(raw.to_dict(), indent=2), encoding="utf-8")
    logger.info("Downloaded %s -> %s (sha256 %s)", record.remote_identifier, local, raw.sha256)
    return raw
