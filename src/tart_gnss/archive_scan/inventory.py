"""Station/day inventory of the TART visibility archive.

Each required ``(station, UTC day)`` prefix is listed at most once per scan,
through the Stage-2 :class:`~tart_gnss.stage2.archive.ArchiveBackend` (the
public S3 archive by default), and only object *metadata* is read: no native
file is ever downloaded. Prefixes come from
:func:`tart_gnss.stage2.archive.day_prefixes` and filenames are parsed with
:func:`tart_gnss.stage2.archive.parse_filename_utc`, so the layout and the
filename convention exist in one place only.

Listing statuses:

* ``FOUND``: the listing succeeded and holds at least one ``vis_*`` file;
* ``EMPTY``: the listing succeeded and holds no ``vis_*`` file;
* ``ERROR``: the listing failed (an S3 error is never empty data);
* ``UNKNOWN``: the prefix was not listed (``offline`` without a cache entry).
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta, timezone
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from tart_gnss.stage2.archive import ArchiveBackend, ArchiveObject, day_prefixes, parse_filename_utc
from tart_gnss.stage2.errors import ArchiveUnavailableError
from tart_gnss.stations import Station
from tart_gnss.tle import format_utc

logger = logging.getLogger(__name__)

STATUS_FOUND, STATUS_EMPTY, STATUS_ERROR, STATUS_UNKNOWN = "FOUND", "EMPTY", "ERROR", "UNKNOWN"
#: Statuses of a listing that actually succeeded.
COMPLETE_STATUSES = (STATUS_FOUND, STATUS_EMPTY)

#: Default bounded retry for transient listing errors: attempts and first backoff (s).
DEFAULT_LIST_ATTEMPTS = 3
DEFAULT_BACKOFF_S = 2.0


@dataclass(frozen=True)
class InventoryObject:
    """One ``vis_*`` archive object (metadata only)."""

    key: str
    size_bytes: Optional[int]
    etag: Optional[str]
    last_modified_utc: Optional[datetime]
    filename_utc: datetime

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable form (UTC strings with microseconds)."""
        return {
            "key": self.key, "size_bytes": self.size_bytes, "etag": self.etag,
            "last_modified_utc": format_utc(self.last_modified_utc, 6) if self.last_modified_utc else None,
            "filename_utc": format_utc(self.filename_utc, 6),
        }


@dataclass
class DayListing:
    """The listing of one ``(station, UTC day)`` archive prefix.

    Attributes
    ----------
    status
        ``FOUND``, ``EMPTY``, ``ERROR`` or ``UNKNOWN`` (module docstring).
    objects
        ``vis_*`` objects, sorted by ``(filename_utc, key)`` and unique by key.
    ignored_keys
        Objects whose last path component is not a ``vis_*`` filename.
    n_duplicate_keys
        Repeated keys dropped while canonicalising the listing.
    source
        ``live`` (listed in this scan), ``cache`` (an earlier snapshot) or
        ``none`` (not listed).
    checksum
        SHA-256 of the canonical object list (:func:`listing_checksum`).
    """

    station_id: str
    archive_name: Optional[str]
    day: date
    prefix: Optional[str]
    host: str
    bucket: str
    status: str
    objects: List[InventoryObject] = field(default_factory=list)
    ignored_keys: List[str] = field(default_factory=list)
    n_duplicate_keys: int = 0
    error: Optional[str] = None
    source: str = "none"
    listing_started_utc: Optional[datetime] = None
    listing_completed_utc: Optional[datetime] = None
    checksum: Optional[str] = None
    cache_path: Optional[str] = None
    software_version: Optional[str] = None
    attempts: int = 0

    @property
    def complete(self) -> bool:
        """True if the listing succeeded (``FOUND`` or ``EMPTY``)."""
        return self.status in COMPLETE_STATUSES

    @property
    def n_files(self) -> int:
        """Number of ``vis_*`` files."""
        return len(self.objects)

    @property
    def total_bytes(self) -> int:
        """Aggregate object size (objects without a size count as 0)."""
        return int(sum(o.size_bytes or 0 for o in self.objects))

    @property
    def earliest_utc(self) -> Optional[datetime]:
        """Earliest filename time, if any."""
        return self.objects[0].filename_utc if self.objects else None

    @property
    def latest_utc(self) -> Optional[datetime]:
        """Latest filename time, if any."""
        return self.objects[-1].filename_utc if self.objects else None

    @property
    def covers_full_day(self) -> Optional[bool]:
        """True if the listing started after the day had settled (see :func:`day_settled`)."""
        if self.listing_started_utc is None:
            return None
        return day_settled(self.day, self.listing_started_utc)


#: A day is "settled" this long after its UTC end; earlier listings may miss late uploads.
DAY_SETTLE_MARGIN = timedelta(hours=1)


def day_settled(day: date, at_utc: datetime) -> bool:
    """True if a listing made at ``at_utc`` can hold the whole of UTC ``day``.

    Files are uploaded minutes after they are recorded, so a day counts as
    settled :data:`DAY_SETTLE_MARGIN` after its UTC end.
    """
    end = datetime(day.year, day.month, day.day, tzinfo=timezone.utc) + timedelta(days=1)
    return at_utc >= end + DAY_SETTLE_MARGIN


def required_days(start_utc: datetime, end_utc: datetime, file_span_s: float) -> List[date]:
    """UTC days whose prefixes must be listed for the half-open query ``[start, end)``.

    A file named ``T`` is assumed to hold samples in ``[T - span, T]``, so it
    can overlap the query only if ``start < T < end + span``. The days of
    ``start - span`` to ``end + span`` are therefore listed: the day after
    ``end`` when ``end + span`` crosses midnight, and the day before ``start``
    when ``start - span`` does (an explicit edge allowance in case a file is
    stored under an earlier day directory than its filename time).
    """
    if end_utc <= start_utc:
        raise ValueError("end must be after start")
    span = timedelta(seconds=file_span_s)
    lo = (start_utc - span).astimezone(timezone.utc).date()
    hi = (end_utc + span).astimezone(timezone.utc).date()
    out = []
    day = lo
    while day <= hi:
        out.append(day)
        day += timedelta(days=1)
    return out


def day_prefix(archive_name: str, day: date) -> str:
    """The archive prefix of one UTC day, exactly as Stage 2 builds it."""
    t = datetime(day.year, day.month, day.day, tzinfo=timezone.utc)
    return day_prefixes(archive_name, t, t)[0]


def canonicalise(raw: Iterable[ArchiveObject]) -> Tuple[List[InventoryObject], List[str], int]:
    """Parse, filter and deduplicate raw archive objects.

    Keys whose last component is not a ``vis_*`` filename are ignored (and
    returned); repeated keys are dropped (the first occurrence is kept).
    The result is sorted by ``(filename_utc, key)``, so it does not depend on
    the listing order.

    Returns
    -------
    (objects, ignored_keys, n_duplicates)
    """
    seen: Dict[str, InventoryObject] = {}
    ignored: List[str] = []
    duplicates = 0
    for obj in raw:
        if obj.key in seen:
            duplicates += 1
            continue
        t = parse_filename_utc(obj.key.rsplit("/", 1)[-1])
        if t is None:
            if obj.key not in ignored:
                ignored.append(obj.key)
            continue
        seen[obj.key] = InventoryObject(
            key=obj.key, size_bytes=obj.size_bytes, etag=obj.etag,
            last_modified_utc=obj.last_modified_utc, filename_utc=t,
        )
    objects = sorted(seen.values(), key=lambda o: (o.filename_utc, o.key))
    return objects, sorted(ignored), duplicates


def listing_checksum(status: str, objects: Sequence[InventoryObject], ignored_keys: Sequence[str]) -> str:
    """SHA-256 of the canonical JSON of a listing's status, objects and ignored keys."""
    doc = {"status": status, "objects": [o.to_dict() for o in objects], "ignored_keys": list(ignored_keys)}
    text = json.dumps(doc, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def list_day_live(
    backend: ArchiveBackend,
    station: Station,
    day: date,
    host: str,
    bucket: str,
    software_version: str,
    attempts: int = DEFAULT_LIST_ATTEMPTS,
    backoff_s: float = DEFAULT_BACKOFF_S,
    sleep: Callable[[float], None] = time.sleep,
    now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
) -> DayListing:
    """List one station/day prefix live, with bounded retry and exponential backoff.

    Exactly one successful ``backend.list_objects`` call is made; failed calls
    are retried up to ``attempts`` times in total. A failure after the last
    attempt gives an ``ERROR`` listing (never an empty one).
    """
    if not station.archive_name:
        return DayListing(station_id=station.id, archive_name=None, day=day, prefix=None, host=host,
                          bucket=bucket, status=STATUS_ERROR, source="none",
                          error=f"station '{station.id}' has no 'archive_name' in the station catalogue")
    prefix = day_prefix(station.archive_name, day)
    started = now()
    error: Optional[str] = None
    for attempt in range(1, attempts + 1):
        try:
            logger.debug("Listing s3://%s/%s (attempt %d)", bucket, prefix, attempt)
            objects, ignored, dups = canonicalise(backend.list_objects(prefix))
        except ArchiveUnavailableError as exc:
            error = str(exc)
            logger.warning("Listing %s failed (attempt %d/%d): %s", prefix, attempt, attempts, exc)
            if attempt < attempts:
                sleep(backoff_s * 2 ** (attempt - 1))
            continue
        status = STATUS_FOUND if objects else STATUS_EMPTY
        return DayListing(
            station_id=station.id, archive_name=station.archive_name, day=day, prefix=prefix,
            host=host, bucket=bucket, status=status, objects=objects, ignored_keys=ignored,
            n_duplicate_keys=dups, source="live", listing_started_utc=started,
            listing_completed_utc=now(), checksum=listing_checksum(status, objects, ignored),
            software_version=software_version, attempts=attempt,
        )
    return DayListing(station_id=station.id, archive_name=station.archive_name, day=day, prefix=prefix,
                      host=host, bucket=bucket, status=STATUS_ERROR, error=error, source="live",
                      listing_started_utc=started, listing_completed_utc=now(),
                      software_version=software_version, attempts=attempts)


@dataclass
class InventoryStats:
    """Counts reported in the manifest and the log."""

    prefixes_required: int = 0
    live_listings: int = 0
    live_list_calls: int = 0
    cache_hits: int = 0
    cache_stale: int = 0
    cache_corrupt: int = 0
    failures: int = 0
    unknown: int = 0

    def to_dict(self) -> Dict[str, int]:
        """Plain dictionary."""
        return dict(self.__dict__)


class _CountingBackend:
    """Wrap a backend to count ``list_objects`` calls (each call is one remote listing)."""

    def __init__(self, backend: ArchiveBackend, stats: InventoryStats):
        self._backend = backend
        self._stats = stats

    def list_objects(self, prefix: str):
        self._stats.live_list_calls += 1
        return self._backend.list_objects(prefix)


#: Log a progress line every this many station/day prefixes.
PROGRESS_EVERY = 100


def build_inventory(
    stations: Sequence[Station],
    days: Sequence[date],
    backend_factory: Callable[[], ArchiveBackend],
    cache,
    policy: str,
    max_age_hours: float,
    host: str,
    bucket: str,
    software_version: str,
    now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
    sleep: Callable[[float], None] = time.sleep,
    attempts: int = DEFAULT_LIST_ATTEMPTS,
) -> Tuple[Dict[Tuple[str, date], DayListing], InventoryStats]:
    """List (or load from cache) every required station/day prefix, once each.

    Parameters
    ----------
    stations, days
        Stations (with ``archive_name``) and the UTC days to cover.
    backend_factory
        Called at most once, and only if a live listing is needed, so the
        ``offline`` policy and full cache reuse never construct an S3 client.
    cache
        :class:`~tart_gnss.archive_scan.cache.InventoryCache`.
    policy
        ``refresh`` (always list live), ``reuse`` (use a valid cached snapshot
        younger than ``max_age_hours`` that was made after the day settled;
        list the rest live) or ``offline`` (cache only; missing days are
        ``UNKNOWN``).

    Returns
    -------
    (dict, InventoryStats)
        ``(station_id, day)`` -> :class:`DayListing`, and counts.
    """
    if policy not in ("refresh", "reuse", "offline"):
        raise ValueError(f"unknown cache policy '{policy}'")
    stats = InventoryStats(prefixes_required=len(stations) * len(days))
    backend: Optional[_CountingBackend] = None
    out: Dict[Tuple[str, date], DayListing] = {}
    done = 0
    for st in stations:
        for day in days:
            done += 1
            if done % PROGRESS_EVERY == 0:
                logger.info("Inventory progress: %d/%d prefixes (%d live, %d cache hits, %d failures)", done,
                            stats.prefixes_required, stats.live_listings, stats.cache_hits, stats.failures)
            key = (st.id, day)
            if not st.archive_name:
                out[key] = list_day_live(None, st, day, host, bucket, software_version)  # type: ignore[arg-type]
                stats.failures += 1
                continue
            prefix = day_prefix(st.archive_name, day)
            cached = None
            if policy in ("reuse", "offline"):
                cached, n_corrupt = cache.latest(host, bucket, prefix)
                stats.cache_corrupt += n_corrupt
                if cached is not None and policy == "reuse":
                    age_h = (now() - cached.listing_completed_utc).total_seconds() / 3600.0
                    if age_h > max_age_hours or not cached.covers_full_day:
                        logger.debug("Cache for %s is stale (age %.1f h, full day %s)", prefix, age_h,
                                     cached.covers_full_day)
                        stats.cache_stale += 1
                        cached = None
            if cached is not None:
                cached.station_id = st.id
                out[key] = cached
                stats.cache_hits += 1
                continue
            if policy == "offline":
                out[key] = DayListing(station_id=st.id, archive_name=st.archive_name, day=day, prefix=prefix,
                                      host=host, bucket=bucket, status=STATUS_UNKNOWN, source="none",
                                      error="offline: no valid cached listing")
                stats.unknown += 1
                continue
            if backend is None:
                backend = _CountingBackend(backend_factory(), stats)
            listing = list_day_live(backend, st, day, host, bucket, software_version,
                                    attempts=attempts, sleep=sleep, now=now)
            stats.live_listings += 1
            if listing.complete:
                listing.cache_path = str(cache.write(listing))
            else:
                stats.failures += 1
            out[key] = listing
            logger.debug("%s %s: %s (%d files)", st.id, day, listing.status, listing.n_files)
    logger.info("Inventory: %d prefixes, %d listed live (%d list calls), %d cache hits, %d stale, "
                "%d corrupt cache files, %d failures, %d unknown", stats.prefixes_required, stats.live_listings,
                stats.live_list_calls, stats.cache_hits, stats.cache_stale, stats.cache_corrupt,
                stats.failures, stats.unknown)
    return out, stats
