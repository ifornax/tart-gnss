"""Persistent, versioned station/day inventory snapshots.

Layout: ``<cache_dir>/<host>/<bucket>/<prefix>/listing_<completed>.json``,
where ``<prefix>`` is the archive day prefix (``<archive_name>/vis/<Y>/<M>/<D>``),
so the cache key distinguishes host, bucket, station archive and day. Every
successful live listing writes a *new* snapshot (earlier ones are kept as
lineage); writes are atomic (temporary file, then ``os.replace``). Failed
listings are never cached, so a failure can never be read back as an empty
day. Snapshots are verified (schema version and checksum) on reading.
"""

from __future__ import annotations

import json
import logging
import os
import re
import tempfile
from datetime import date
from pathlib import Path
from typing import List, Optional, Tuple, Union

from tart_gnss.archive_scan.errors import InventoryCacheError
from tart_gnss.archive_scan.inventory import (
    COMPLETE_STATUSES,
    DayListing,
    InventoryObject,
    listing_checksum,
)
from tart_gnss.stage2.selection import parse_utc
from tart_gnss.tle import format_utc

logger = logging.getLogger(__name__)

#: Schema version of a cached listing; other versions are rejected.
CACHE_SCHEMA_VERSION = 1

_SAFE = re.compile(r"[^A-Za-z0-9._-]")


def _safe(part: str) -> str:
    return _SAFE.sub("_", part) or "_"


def _stamp(listing: DayListing) -> str:
    return listing.listing_completed_utc.strftime("%Y%m%dT%H%M%S%fZ")


def listing_to_dict(listing: DayListing) -> dict:
    """Serialise a listing into the cache document."""
    return {
        "schema_version": CACHE_SCHEMA_VERSION,
        "station_id": listing.station_id,
        "archive_name": listing.archive_name,
        "day": listing.day.isoformat(),
        "bucket": listing.bucket,
        "host": listing.host,
        "prefix": listing.prefix,
        "listing_started_utc": format_utc(listing.listing_started_utc, 6),
        "listing_completed_utc": format_utc(listing.listing_completed_utc, 6),
        "listing_status": listing.status,
        "n_duplicate_keys": listing.n_duplicate_keys,
        "ignored_keys": list(listing.ignored_keys),
        "objects": [o.to_dict() for o in listing.objects],
        "source": "live",
        "software_version": listing.software_version,
        "checksum_of_canonical_listing": listing.checksum,
    }


def listing_from_dict(doc: dict, path: Optional[Path] = None) -> DayListing:
    """Rebuild a listing from a cache document, verifying schema and checksum.

    Raises
    ------
    InventoryCacheError
        On a different schema version, a non-successful status, missing
        fields or a checksum mismatch.
    """
    if doc.get("schema_version") != CACHE_SCHEMA_VERSION:
        raise InventoryCacheError(f"{path}: cache schema {doc.get('schema_version')!r}, "
                                  f"expected {CACHE_SCHEMA_VERSION}")
    try:
        status = doc["listing_status"]
        if status not in COMPLETE_STATUSES:
            raise InventoryCacheError(f"{path}: cached status {status!r} is not a successful listing")
        objects = [
            InventoryObject(key=o["key"], size_bytes=o["size_bytes"], etag=o["etag"],
                            last_modified_utc=parse_utc(o["last_modified_utc"]) if o["last_modified_utc"] else None,
                            filename_utc=parse_utc(o["filename_utc"]))
            for o in doc["objects"]
        ]
        ignored = list(doc["ignored_keys"])
        checksum = listing_checksum(status, objects, ignored)
        if checksum != doc["checksum_of_canonical_listing"]:
            raise InventoryCacheError(f"{path}: checksum mismatch (cache file corrupt or edited)")
        if (status == "FOUND") != bool(objects):
            raise InventoryCacheError(f"{path}: status {status} inconsistent with {len(objects)} objects")
        return DayListing(
            station_id=doc["station_id"], archive_name=doc["archive_name"], day=date.fromisoformat(doc["day"]),
            prefix=doc["prefix"], host=doc["host"], bucket=doc["bucket"], status=status, objects=objects,
            ignored_keys=ignored, n_duplicate_keys=int(doc.get("n_duplicate_keys", 0)), source="cache",
            listing_started_utc=parse_utc(doc["listing_started_utc"]),
            listing_completed_utc=parse_utc(doc["listing_completed_utc"]), checksum=checksum,
            cache_path=str(path) if path else None, software_version=doc.get("software_version"),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise InventoryCacheError(f"{path}: malformed cache document ({exc})") from exc


class InventoryCache:
    """Read and write station/day listing snapshots under ``root``."""

    def __init__(self, root: Union[str, Path]):
        self.root = Path(root)

    def directory(self, host: str, bucket: str, prefix: str) -> Path:
        """Snapshot directory for one host/bucket/day prefix."""
        parts = [_safe(p) for p in prefix.strip("/").split("/")]
        return self.root.joinpath(_safe(host), _safe(bucket), *parts)

    def write(self, listing: DayListing) -> Path:
        """Atomically write a new snapshot of a successful listing; returns its path.

        Raises
        ------
        InventoryCacheError
            If the listing did not succeed (failures are never cached).
        """
        if not listing.complete:
            raise InventoryCacheError(f"refusing to cache a {listing.status} listing of {listing.prefix}")
        d = self.directory(listing.host, listing.bucket, listing.prefix)
        d.mkdir(parents=True, exist_ok=True)
        path = d / f"listing_{_stamp(listing)}.json"
        text = json.dumps(listing_to_dict(listing), indent=1, sort_keys=True) + "\n"
        fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp_", suffix=".json")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(text)
            os.replace(tmp, path)
        except BaseException:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise
        return path

    def snapshots(self, host: str, bucket: str, prefix: str) -> List[Path]:
        """All snapshot files for a prefix, oldest first."""
        d = self.directory(host, bucket, prefix)
        return sorted(d.glob("listing_*.json")) if d.is_dir() else []

    def latest(self, host: str, bucket: str, prefix: str) -> Tuple[Optional[DayListing], int]:
        """The newest valid snapshot of a prefix, and the number of corrupt ones skipped.

        A corrupt or incompatible newer snapshot is reported (warning) and the
        next older valid one is used.
        """
        corrupt = 0
        for path in reversed(self.snapshots(host, bucket, prefix)):
            try:
                doc = json.loads(path.read_text(encoding="utf-8"))
                listing = listing_from_dict(doc, path)
            except (OSError, json.JSONDecodeError, InventoryCacheError) as exc:
                logger.warning("Ignoring invalid inventory cache file %s: %s", path, exc)
                corrupt += 1
                continue
            if listing.prefix != prefix or listing.host != host or listing.bucket != bucket:
                logger.warning("Ignoring cache file %s: it describes %s/%s/%s", path, listing.host,
                               listing.bucket, listing.prefix)
                corrupt += 1
                continue
            return listing, corrupt
        return None, corrupt
