"""TART network registry: retrieval, snapshotting, normalisation and diffing.

The TART network monitor publishes every registered telescope at
:data:`DEFAULT_REGISTRY_URL`. This module is the ONLY place that talks to that
endpoint or knows its JSON layout. It turns a download into:

1. a raw, never-overwritten snapshot (``registry_<UTC>.json``) with a metadata
   sidecar recording the URL, retrieval time and SHA-256;
2. validated :class:`RegistryStation` records plus diagnostics for malformed
   entries (coordinates are validated by
   :func:`tart_gnss.stations.validate_coordinates`; nothing is guessed);
3. a generated, normalised station catalogue that
   :func:`tart_gnss.stations.load_stations` reads like any other.

The live registry is for *discovery*. Stage 1 reads a frozen catalogue by
default, so a later registry update never changes an old analysis.

Registry presence says nothing about archived data: ``archive_name`` is set
only when :func:`verify_archives` found the prefix in the visibility archive
(or such a verification is carried over from a previous catalogue).
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Union

import yaml

from tart_gnss.stations import Station, StationConfigError, validate_coordinates

logger = logging.getLogger(__name__)

#: Official TART network registry endpoint (JSON). The interactive map at
#: https://map.elec.ac.nz is built on it; its HTML is never scraped.
DEFAULT_REGISTRY_URL = "https://map.elec.ac.nz/api/v1/telescopes"

#: Base of the per-telescope TART REST API (convention used by tart2ms).
TART_API_BASE = "https://api.elec.ac.nz/tart"

#: Version of the generated catalogue layout written by :func:`build_catalogue_document`.
CATALOGUE_FORMAT_VERSION = 1

#: Default flag level (metres) above which a station displacement between two
#: catalogues is logged as a warning. It is a reporting threshold chosen to
#: exceed the metre-level rounding seen in the registry, NOT a validity bound:
#: no change is ever rejected on this basis.
DEFAULT_DISPLACEMENT_WARNING_M = 100.0


class RegistryError(Exception):
    """Base error for TART registry operations."""


class RegistryFetchError(RegistryError):
    """Raised when the remote TART registry cannot be retrieved."""


class RegistrySchemaError(RegistryError):
    """Raised when the registry payload has an unsupported schema."""


class RegistryStationError(RegistryError):
    """Raised when an individual registry station is invalid."""


def utc_stamp(dt: datetime) -> str:
    """Format an aware datetime as compact UTC, e.g. ``20261009T154800Z``."""
    return dt.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def utc_iso(dt: datetime) -> str:
    """Format an aware datetime as ISO-8601 UTC to the second, ``Z`` suffix."""
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(data: bytes) -> str:
    """Return the SHA-256 hex digest of ``data`` (deterministic for equal bytes)."""
    return hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------
# Retrieval
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RegistryDownload:
    """One retrieval of the registry, exactly as received.

    Attributes
    ----------
    url
        Endpoint the payload was fetched from.
    retrieved_utc
        Aware UTC datetime at which the response was received.
    raw
        Response body bytes, unmodified (hashed and archived as-is).
    payload
        ``raw`` parsed as JSON.
    """

    url: str
    retrieved_utc: datetime
    raw: bytes
    payload: Any

    @property
    def sha256(self) -> str:
        """SHA-256 of :attr:`raw`."""
        return sha256_bytes(self.raw)


def download_registry(url: str = DEFAULT_REGISTRY_URL, timeout_s: float = 30.0) -> RegistryDownload:
    """Download the TART registry, keeping the raw bytes for archiving.

    Parameters
    ----------
    url
        TART network registry endpoint.
    timeout_s
        HTTP connect/read timeout in seconds.

    Returns
    -------
    RegistryDownload
        Raw bytes, parsed JSON and retrieval time (UTC).

    Raises
    ------
    RegistryFetchError
        If the endpoint cannot be reached, returns an HTTP error status, or
        returns a body that is not valid JSON.
    """
    import requests

    logger.info("Fetching TART registry from %s", url)
    try:
        response = requests.get(url, timeout=timeout_s, headers={"Accept": "application/json"})
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RegistryFetchError(f"Cannot retrieve TART registry {url}: {exc}") from exc
    retrieved = datetime.now(timezone.utc)
    raw = response.content
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise RegistryFetchError(f"TART registry {url} did not return valid JSON: {exc}") from exc
    return RegistryDownload(url=url, retrieved_utc=retrieved, raw=raw, payload=payload)


def fetch_registry(url: str = DEFAULT_REGISTRY_URL, timeout_s: float = 30.0) -> Any:
    """Download the current TART telescope registry and return the parsed JSON.

    Parameters
    ----------
    url
        TART network registry endpoint.
    timeout_s
        HTTP request timeout in seconds.

    Returns
    -------
    dict
        Parsed registry JSON.

    Raises
    ------
    RegistryFetchError
        If the endpoint cannot be reached or returns invalid JSON.
    """
    return download_registry(url, timeout_s).payload


# --------------------------------------------------------------------------
# Snapshots
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SnapshotFiles:
    """Paths written by :func:`write_snapshot`."""

    raw: Path
    meta: Path
    latest: Path


def write_snapshot(
    download: RegistryDownload, snapshot_dir: Union[str, Path], station_count: int
) -> SnapshotFiles:
    """Archive a raw registry payload with a provenance sidecar.

    Writes ``registry_<UTC>.json`` (the exact bytes received),
    ``registry_<UTC>.meta.json`` (``source_url``, ``retrieved_utc``,
    ``sha256``, ``station_count``) and refreshes ``latest.json`` (a copy of
    the newest raw payload).

    Parameters
    ----------
    download
        Output of :func:`download_registry`.
    snapshot_dir
        Directory for snapshots (created if needed).
    station_count
        Number of entries that validated into stations.

    Returns
    -------
    SnapshotFiles

    Raises
    ------
    RegistryError
        If a snapshot with the same timestamp already exists (historical
        snapshots are never overwritten).
    """
    directory = Path(snapshot_dir)
    directory.mkdir(parents=True, exist_ok=True)
    stem = f"registry_{utc_stamp(download.retrieved_utc)}"
    raw_path = directory / f"{stem}.json"
    meta_path = directory / f"{stem}.meta.json"
    if raw_path.exists() or meta_path.exists():
        raise RegistryError(f"Snapshot {raw_path} already exists; refusing to overwrite it")
    raw_path.write_bytes(download.raw)
    meta = {
        "source_url": download.url,
        "retrieved_utc": utc_iso(download.retrieved_utc),
        "sha256": download.sha256,
        "station_count": station_count,
        "snapshot_file": raw_path.name,
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    latest = directory / "latest.json"
    latest.write_bytes(download.raw)
    logger.info("Wrote registry snapshot %s (sha256 %s)", raw_path, download.sha256)
    return SnapshotFiles(raw=raw_path, meta=meta_path, latest=latest)


# --------------------------------------------------------------------------
# Parsing and normalisation
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RegistryStation:
    """A validated, normalised registry entry.

    Attributes
    ----------
    tart_id
        Canonical TART telescope identifier (registry ``nodeName``).
    name
        Display name (registry ``telescopeName``, else the TART id).
    latitude_deg, longitude_deg
        WGS84 geodetic latitude/longitude, degrees (east positive).
    altitude_m
        Height in metres as published (vertical datum undocumented by TART).
    api_url
        TART API URL published by the registry, or ``None`` (the registry
        does not currently publish one; see :func:`default_api_url`).
    archive_name
        Archive prefix published by the registry, or ``None``.
    registry_source
        Registry URL the entry came from.
    registry_retrieved_utc
        UTC ISO-8601 retrieval time of the registry snapshot.
    online
        Registry ``online`` flag at retrieval (informational only; not a
        verified statement about the telescope API).
    """

    tart_id: str
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    api_url: Optional[str]
    archive_name: Optional[str]
    registry_source: str
    registry_retrieved_utc: str
    online: Optional[bool] = None


@dataclass(frozen=True)
class RegistryDiagnostic:
    """Why a registry entry was not turned into a station.

    Attributes
    ----------
    index
        Position of the entry in the registry ``telescopes`` list.
    tart_id
        The entry's TART id, if it had one.
    reason
        Human-readable explanation.
    """

    index: int
    tart_id: Optional[str]
    reason: str


@dataclass
class ParsedRegistry:
    """Result of :func:`parse_registry`.

    Attributes
    ----------
    entry_count
        Number of entries in the registry.
    stations
        Valid stations keyed by TART id, in registry order.
    diagnostics
        One record per skipped entry.
    """

    entry_count: int
    stations: Dict[str, RegistryStation] = field(default_factory=dict)
    diagnostics: List[RegistryDiagnostic] = field(default_factory=list)


def _get_path(entry: Mapping[str, Any], *keys: str) -> Any:
    node: Any = entry
    for key in keys:
        if not isinstance(node, Mapping):
            return None
        node = node.get(key)
    return node


def _optional_str(value: Any) -> Optional[str]:
    return value.strip() if isinstance(value, str) and value.strip() else None


def parse_registry(
    payload: Any,
    retrieved_utc: str,
    source: str = DEFAULT_REGISTRY_URL,
    strict: bool = False,
) -> ParsedRegistry:
    """Normalise the TART registry into validated station records.

    The supported layout is ``{"telescopes": [{"nodeName": ..., "telescopeName":
    ..., "telescopeInfo": {"info": {"location": {"lat", "lon", "alt"}}}}, ...]}``.
    Optional ``apiUrl``/``archiveName`` entry fields are used if the registry
    ever publishes them (it does not as of 2026-10). This is the only
    function that reads registry JSON fields.

    Parameters
    ----------
    payload
        Parsed registry JSON.
    retrieved_utc
        UTC ISO-8601 retrieval time, recorded on every station.
    source
        Registry URL, recorded on every station.
    strict
        If true, any malformed entry raises instead of being diagnosed.

    Returns
    -------
    ParsedRegistry
        Valid stations plus one diagnostic per skipped entry. Entries with a
        missing, non-finite or out-of-range coordinate are skipped (never set
        to zero). All entries sharing a duplicated TART id are skipped, since
        the registry is then ambiguous about that station.

    Raises
    ------
    RegistrySchemaError
        If the top-level structure is not a mapping with a ``telescopes`` list.
    RegistryStationError
        In ``strict`` mode, for the first malformed or duplicated entry.
    """
    if not isinstance(payload, Mapping) or not isinstance(payload.get("telescopes"), list):
        found = sorted(payload) if isinstance(payload, Mapping) else type(payload).__name__
        raise RegistrySchemaError(
            f"Unsupported TART registry schema: expected a mapping with a 'telescopes' list, got {found}"
        )
    entries = payload["telescopes"]
    logger.info("Registry contains %d telescope entries", len(entries))
    result = ParsedRegistry(entry_count=len(entries))

    def reject(index: int, tart_id: Optional[str], reason: str) -> None:
        if strict:
            raise RegistryStationError(f"Registry entry {index} ({tart_id or 'no id'}): {reason}")
        logger.warning("Skipping station %s: %s", tart_id or f"#{index}", reason)
        result.diagnostics.append(RegistryDiagnostic(index, tart_id, reason))

    candidates: List[RegistryStation] = []
    indices: List[int] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            reject(index, None, "entry is not a JSON object")
            continue
        tart_id = _optional_str(entry.get("nodeName"))
        if tart_id is None:
            reject(index, None, "missing 'nodeName' (TART id)")
            continue
        location = _get_path(entry, "telescopeInfo", "info", "location")
        if not isinstance(location, Mapping):
            reject(index, tart_id, "missing valid coordinates: no telescopeInfo.info.location")
            continue
        lat, lon, alt = location.get("lat"), location.get("lon"), location.get("alt")
        try:
            validate_coordinates(tart_id, lat, lon, alt)
        except StationConfigError as exc:
            reject(index, tart_id, f"missing valid coordinates: {exc}")
            continue
        online = entry.get("online")
        candidates.append(
            RegistryStation(
                tart_id=tart_id,
                name=_optional_str(entry.get("telescopeName")) or tart_id,
                latitude_deg=float(lat),
                longitude_deg=float(lon),
                altitude_m=float(alt),
                api_url=_optional_str(entry.get("apiUrl")),
                archive_name=_optional_str(entry.get("archiveName")),
                registry_source=source,
                registry_retrieved_utc=retrieved_utc,
                online=online if isinstance(online, bool) else None,
            )
        )
        indices.append(index)

    counts: Dict[str, int] = {}
    for station in candidates:
        counts[station.tart_id] = counts.get(station.tart_id, 0) + 1
    for index, station in zip(indices, candidates):
        if counts[station.tart_id] > 1:
            reject(index, station.tart_id,
                   f"duplicate TART id ({counts[station.tart_id]} entries); ambiguous, all skipped")
        else:
            result.stations[station.tart_id] = station
    return result


def default_api_url(tart_id: str) -> str:
    """Return the conventional TART REST API URL for a telescope.

    This is *convention-derived* metadata (``https://api.elec.ac.nz/tart/<id>``,
    the pattern used by tart2ms), not something the registry publishes, and
    says nothing about whether the API is online.
    """
    return f"{TART_API_BASE}/{tart_id}"


# --------------------------------------------------------------------------
# Archive verification
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ArchiveStatus:
    """Archive availability for one TART id.

    ``archive_name`` is set only when ``archive_verified`` is true.
    """

    archive_name: Optional[str]
    archive_verified: bool
    archive_verified_utc: Optional[str]


def s3_prefix_exists(prefix: str, timeout_s: float = 30.0) -> bool:
    """Return whether any object exists under ``prefix`` in the TART archive.

    Uses an anonymous S3 ``ListObjectsV2`` request limited to one key.

    Raises
    ------
    RegistryFetchError
        If the archive cannot be queried.
    """
    import requests

    from tart_gnss.stage2.archive import ARCHIVE_BUCKET, ARCHIVE_HOST

    url = f"https://{ARCHIVE_HOST}/{ARCHIVE_BUCKET}"
    try:
        response = requests.get(url, params={"list-type": "2", "prefix": prefix, "max-keys": "1"},
                                timeout=timeout_s)
        response.raise_for_status()
        root = ET.fromstring(response.content)
    except (requests.RequestException, ET.ParseError) as exc:
        raise RegistryFetchError(f"Cannot list {url}?prefix={prefix}: {exc}") from exc
    count = next((el.text for el in root.iter() if el.tag.endswith("KeyCount")), None)
    if count is not None:
        return int(count) > 0
    return any(el.tag.endswith("Contents") for el in root.iter())


def verify_archives(
    tart_ids: Sequence[str],
    prefix_exists: Callable[[str], bool] = s3_prefix_exists,
    now: Optional[datetime] = None,
) -> Dict[str, ArchiveStatus]:
    """Check which TART ids have a ``<tart_id>/vis/`` prefix in the archive.

    Parameters
    ----------
    tart_ids
        Telescopes to check.
    prefix_exists
        Callable answering whether an archive prefix holds any object
        (injectable for tests).
    now
        Verification time (default: current UTC time).

    Returns
    -------
    dict
        TART id -> :class:`ArchiveStatus`. A failed query is logged and the
        station is recorded as unverified (never assumed available).
    """
    stamp = utc_iso(now or datetime.now(timezone.utc))
    out: Dict[str, ArchiveStatus] = {}
    for tart_id in tart_ids:
        try:
            found = prefix_exists(f"{tart_id}/vis/")
        except RegistryFetchError as exc:
            logger.warning("Archive check for %s failed: %s", tart_id, exc)
            out[tart_id] = ArchiveStatus(None, False, None)
            continue
        logger.info("Archive prefix %s/vis/: %s", tart_id, "found" if found else "not found")
        out[tart_id] = ArchiveStatus(tart_id if found else None, found, stamp)
    return out


def archive_metadata_from_catalogue(stations: Mapping[str, Station]) -> Dict[str, ArchiveStatus]:
    """Carry forward archive verifications recorded in an earlier catalogue.

    Only entries with ``archive_verified`` true and a verification time are
    kept, so the provenance of each carried-over value stays explicit.
    """
    out: Dict[str, ArchiveStatus] = {}
    for st in stations.values():
        if st.archive_verified and st.archive_name and st.archive_verified_utc:
            out[st.tart_id or st.id] = ArchiveStatus(st.archive_name, True, st.archive_verified_utc)
    return out


# --------------------------------------------------------------------------
# Conversion to Station objects and the generated catalogue
# --------------------------------------------------------------------------


def _station_source(rs: RegistryStation, snapshot_sha256: Optional[str]) -> str:
    text = f"TART network registry {rs.registry_source} (node {rs.tart_id}), retrieved {rs.registry_retrieved_utc}"
    if snapshot_sha256:
        text += f", snapshot sha256 {snapshot_sha256}"
    return text


def registry_to_stations(
    registry_stations: Mapping[str, RegistryStation],
    archive_metadata: Optional[Mapping[str, ArchiveStatus]] = None,
    snapshot_sha256: Optional[str] = None,
) -> Dict[str, Station]:
    """Convert normalised registry entries into validated :class:`Station` objects.

    Parameters
    ----------
    registry_stations
        ``ParsedRegistry.stations``.
    archive_metadata
        TART id -> :class:`ArchiveStatus`; stations without an entry get
        ``archive_name=None, archive_verified=False``.
    snapshot_sha256
        SHA-256 of the raw registry snapshot, recorded in each ``source``.

    Returns
    -------
    dict
        TART id -> :class:`Station` whose ``id`` *is* the canonical TART id.
        ``api_url`` is the registry value, else :func:`default_api_url`.

    Raises
    ------
    StationConfigError
        If a station fails validation (cannot happen for parser output).
    """
    archive_metadata = archive_metadata or {}
    out: Dict[str, Station] = {}
    for tart_id, rs in registry_stations.items():
        status = archive_metadata.get(tart_id)
        if status is None and rs.archive_name:
            status = ArchiveStatus(rs.archive_name, False, None)
        out[tart_id] = Station(
            id=tart_id,
            name=rs.name,
            latitude_deg=rs.latitude_deg,
            longitude_deg=rs.longitude_deg,
            altitude_m=rs.altitude_m,
            source=_station_source(rs, snapshot_sha256),
            tart_id=tart_id,
            api_url=rs.api_url or default_api_url(tart_id),
            archive_name=status.archive_name if status else None,
            archive_verified=status.archive_verified if status else False,
            archive_verified_utc=status.archive_verified_utc if status else None,
        )
    return out


def build_catalogue_document(
    parsed: ParsedRegistry,
    stations: Mapping[str, Station],
    *,
    source_url: str,
    retrieved_utc: str,
    sha256: str,
    snapshot_file: Optional[str],
) -> Dict[str, Any]:
    """Assemble the generated station catalogue (``stations_registry.yaml``).

    The document is readable by :func:`tart_gnss.stations.load_stations`.
    Its ``metadata`` block records where and when the registry was
    retrieved, the SHA-256 of the raw snapshot, and the skipped entries.
    Unknown values are written as ``null``, never guessed. Stations are
    sorted by TART id so successive catalogues diff cleanly.
    """
    api_from_registry = {k for k, v in parsed.stations.items() if v.api_url}
    return {
        "metadata": {
            "format_version": CATALOGUE_FORMAT_VERSION,
            "generated_by": "tart-stations-update",
            "source": source_url,
            "retrieved_utc": retrieved_utc,
            "sha256": sha256,
            "snapshot_file": snapshot_file,
            "registry_entry_count": parsed.entry_count,
            "station_count": len(stations),
            "api_url_note": (
                "api_url is convention-derived (https://api.elec.ac.nz/tart/<tart_id>) unless the "
                "registry publishes one; it is not a statement that the API is online"
            ),
            "api_url_from_registry": sorted(api_from_registry),
            "archive_note": (
                "archive_name is set only for a verified archive prefix; registry presence does "
                "not imply archived data"
            ),
            "skipped_entries": [asdict(d) for d in parsed.diagnostics],
        },
        "stations": {
            sid: {
                "name": st.name,
                "tart_id": st.tart_id,
                "latitude_deg": st.latitude_deg,
                "longitude_deg": st.longitude_deg,
                "altitude_m": st.altitude_m,
                "api_url": st.api_url,
                "archive_name": st.archive_name,
                "archive_verified": st.archive_verified,
                "archive_verified_utc": st.archive_verified_utc,
                "source": st.source,
            }
            for sid, st in sorted(stations.items())
        },
    }


CATALOGUE_HEADER = """\
# GENERATED FILE - do not edit by hand. Regenerate with `tart-stations-update`.
#
# Normalised TART station catalogue built from a frozen snapshot of the TART
# network registry (see `metadata`). Stations are keyed by canonical TART id.
# Coordinates are copied verbatim from the registry and validated; entries
# with missing/invalid coordinates are listed under metadata.skipped_entries.
# Hand-maintained convenience names live in station_aliases.yaml.

"""


def write_catalogue(path: Union[str, Path], document: Mapping[str, Any]) -> Path:
    """Write a generated catalogue as YAML (atomically, via a temporary file)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(CATALOGUE_HEADER + yaml.safe_dump(dict(document), sort_keys=False, allow_unicode=True),
                   encoding="utf-8")
    tmp.replace(path)
    logger.info("Wrote normalised station catalogue %s", path)
    return path


# --------------------------------------------------------------------------
# Diffing
# --------------------------------------------------------------------------

#: Station attributes compared as "metadata" by :func:`diff_catalogues`.
DIFF_METADATA_FIELDS = ("api_url", "archive_name", "archive_verified")


def displacement_m(a: Station, b: Station) -> float:
    """Straight-line distance (m) between two stations' WGS84 positions.

    Both geodetic positions are converted to Earth-centred Cartesian
    coordinates with Skyfield's WGS84 model, so latitude, longitude and
    height changes all contribute.
    """
    from skyfield.api import wgs84

    pa = wgs84.latlon(a.latitude_deg, a.longitude_deg, elevation_m=a.altitude_m).itrs_xyz.m
    pb = wgs84.latlon(b.latitude_deg, b.longitude_deg, elevation_m=b.altitude_m).itrs_xyz.m
    return float(math.sqrt(sum((x - y) ** 2 for x, y in zip(pa, pb))))


@dataclass
class RegistryDiff:
    """Differences between two station catalogues, keyed by TART id.

    Attributes
    ----------
    added, removed
        TART ids present only in the new / old catalogue.
    coordinate_changes
        One dict per moved station with old/new values, deltas
        (``delta_latitude_deg``, ``delta_longitude_deg``, ``delta_altitude_m``)
        and ``displacement_m``.
    name_changes
        ``{"tart_id", "old", "new"}`` records.
    metadata_changes
        ``{"tart_id", "field", "old", "new"}`` records for
        :data:`DIFF_METADATA_FIELDS`.
    """

    added: List[str] = field(default_factory=list)
    removed: List[str] = field(default_factory=list)
    coordinate_changes: List[Dict[str, Any]] = field(default_factory=list)
    name_changes: List[Dict[str, Any]] = field(default_factory=list)
    metadata_changes: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        """True if the catalogues are equivalent."""
        return not (self.added or self.removed or self.coordinate_changes
                    or self.name_changes or self.metadata_changes)

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serialisable form."""
        return asdict(self)

    def format(self) -> str:
        """Human-readable multi-line summary."""
        if self.is_empty:
            return "No changes."
        lines: List[str] = []
        lines += [f"+ new station {t}" for t in self.added]
        lines += [f"- removed station {t}" for t in self.removed]
        for c in self.coordinate_changes:
            lines.append(
                f"~ {c['tart_id']}: coordinates moved {c['displacement_m']:.2f} m "
                f"(dlat {c['delta_latitude_deg']:+.7f} deg, dlon {c['delta_longitude_deg']:+.7f} deg, "
                f"dalt {c['delta_altitude_m']:+.2f} m)"
            )
        lines += [f"~ {c['tart_id']}: name {c['old']!r} -> {c['new']!r}" for c in self.name_changes]
        lines += [f"~ {c['tart_id']}: {c['field']} {c['old']!r} -> {c['new']!r}" for c in self.metadata_changes]
        return "\n".join(lines)


def _by_tart_id(stations: Mapping[str, Station]) -> Dict[str, Station]:
    return {(st.tart_id or st.id): st for st in stations.values()}


def diff_catalogues(
    old: Mapping[str, Station],
    new: Mapping[str, Station],
    displacement_warning_m: float = DEFAULT_DISPLACEMENT_WARNING_M,
) -> RegistryDiff:
    """Compare two station catalogues by TART id.

    Stations are matched on ``tart_id`` (falling back to ``id``), so a legacy
    alias-keyed catalogue can be compared with a generated one.

    Parameters
    ----------
    old, new
        Catalogues as returned by :func:`tart_gnss.stations.load_stations`.
    displacement_warning_m
        Log a warning for displacements above this many metres (a reporting
        flag level only; nothing is rejected).

    Returns
    -------
    RegistryDiff
    """
    a, b = _by_tart_id(old), _by_tart_id(new)
    diff = RegistryDiff(added=sorted(set(b) - set(a)), removed=sorted(set(a) - set(b)))
    for tart_id in sorted(set(a) & set(b)):
        o, n = a[tart_id], b[tart_id]
        if (o.latitude_deg, o.longitude_deg, o.altitude_m) != (n.latitude_deg, n.longitude_deg, n.altitude_m):
            moved = displacement_m(o, n)
            diff.coordinate_changes.append({
                "tart_id": tart_id,
                "old": {"latitude_deg": o.latitude_deg, "longitude_deg": o.longitude_deg, "altitude_m": o.altitude_m},
                "new": {"latitude_deg": n.latitude_deg, "longitude_deg": n.longitude_deg, "altitude_m": n.altitude_m},
                "delta_latitude_deg": n.latitude_deg - o.latitude_deg,
                "delta_longitude_deg": n.longitude_deg - o.longitude_deg,
                "delta_altitude_m": n.altitude_m - o.altitude_m,
                "displacement_m": moved,
            })
            if moved > displacement_warning_m:
                logger.warning("Station %s moved %.1f m between catalogues (flag level %.0f m)",
                               tart_id, moved, displacement_warning_m)
        if o.name != n.name:
            diff.name_changes.append({"tart_id": tart_id, "old": o.name, "new": n.name})
        for attr in DIFF_METADATA_FIELDS:
            ov, nv = getattr(o, attr), getattr(n, attr)
            if ov != nv:
                diff.metadata_changes.append({"tart_id": tart_id, "field": attr, "old": ov, "new": nv})
    return diff
