"""Finding observations in the TART visibility archive.

The TART project stores native HDF5 visibility files in a public S3 bucket.
This module reuses the bucket location from the maintained client
:mod:`tart_tools.archive_handler` (``MINIO_API_HOST``, ``BUCKET_NAME``) and the
same ``minio`` client library, so it does not hard-code a second copy of the
archive location when ``tart_tools`` is installed. The archive layout is
``<archive_name>/vis/<year>/<month>/<day>/vis_<YYYY-MM-DD_HH_MM_SS.ffffff>.hdf``.

``tart_tools.handle_archive_request`` itself is not called, because:

1. it selects objects by S3 ``last_modified``, which is the *upload* time and
   lags the observation by minutes;
2. it renames files to ``obs_00000.hdf``, which loses the original identifiers
   needed for provenance.

Instead, files are selected by the observation time encoded in the filename.
In every file inspected, the filename time equals the timestamp of the last
sample, and each file holds about 60 one-second samples. A file named ``T``
therefore holds data from roughly ``T - 61 s`` to ``T``. The search keeps
files with ``start <= T <= end + margin``. The actual timestamps are always
checked afterwards by reading the HDF5 file (:mod:`tart_gnss.stage2.native`).

The rest of the pipeline sees only :class:`ObservationRecord` objects, never
raw HTTP or S3 responses.
"""

from __future__ import annotations

import logging
import re
from dataclasses import asdict, dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable, List, Optional, Protocol

from tart_gnss.stage2.errors import ArchiveUnavailableError
from tart_gnss.stations import Station
from tart_gnss.tle import format_utc

logger = logging.getLogger(__name__)

try:  # reuse the maintained TART archive location when available
    from tart_tools.archive_handler import BUCKET_NAME as _TT_BUCKET
    from tart_tools.archive_handler import MINIO_API_HOST as _TT_HOST
except ImportError:  # tart_tools not installed: same values, copied from tart_tools 1.4.6
    _TT_HOST, _TT_BUCKET = "s3.us-west-2.amazonaws.com", "tart"

ARCHIVE_HOST: str = _TT_HOST
ARCHIVE_BUCKET: str = _TT_BUCKET

#: Default margin (s) after the query end for selecting files by filename time.
#: Each file holds ~61 s of data ending at its filename time, so 120 s keeps
#: every file that can overlap the query interval.
DEFAULT_FILE_SPAN_MARGIN_S = 120.0

_FILENAME_RE = re.compile(
    r"vis_(\d{4})-(\d{2})-(\d{2})_(\d{2})_(\d{2})_(\d{2})(?:\.(\d{1,6}))?\.(hdf|h5|hdf5)$"
)


@dataclass(frozen=True)
class ArchiveObject:
    """An object listed in the archive (backend-neutral)."""

    key: str
    size_bytes: Optional[int]
    etag: Optional[str]
    last_modified_utc: Optional[datetime]


class ArchiveBackend(Protocol):
    """The interface an archive backend must provide (real S3 or a test mock)."""

    def list_objects(self, prefix: str) -> Iterable[ArchiveObject]:
        """List every object under ``prefix``, recursively."""

    def download(self, key: str, destination: Path) -> None:
        """Download object ``key`` to ``destination``."""

    def url(self, key: str) -> str:
        """Return a stable reference (URL) for ``key``."""


class TartS3Archive:
    """Anonymous, read-only access to the public TART S3 visibility archive.

    Parameters
    ----------
    host, bucket
        Archive endpoint and bucket. The defaults come from ``tart_tools``.
    client
        Optional pre-built ``minio.Minio`` client (for testing).

    Raises
    ------
    ArchiveUnavailableError
        If the ``minio`` package is missing or the archive cannot be reached.
    """

    def __init__(self, host: str = ARCHIVE_HOST, bucket: str = ARCHIVE_BUCKET, client=None):
        self.host = host
        self.bucket = bucket
        if client is None:
            try:
                from minio import Minio
            except ImportError as exc:
                raise ArchiveUnavailableError(
                    "Archive access needs the 'minio' package: pip install 'tart-gnss[stage2]'"
                ) from exc
            client = Minio(endpoint=host, secure=True)
        self._client = client

    def _errors(self):
        errors = [OSError]
        try:
            from minio.error import MinioException
            errors.append(MinioException)
        except ImportError:
            pass
        try:
            from urllib3.exceptions import HTTPError
            errors.append(HTTPError)
        except ImportError:
            pass
        return tuple(errors)

    def list_objects(self, prefix: str) -> List[ArchiveObject]:
        """List objects under ``prefix``. Wraps backend errors in ArchiveUnavailableError."""
        try:
            objects = self._client.list_objects(self.bucket, prefix=prefix, recursive=True)
            return [
                ArchiveObject(
                    key=o.object_name,
                    size_bytes=o.size,
                    etag=(o.etag or "").strip('"') or None,
                    last_modified_utc=o.last_modified.astimezone(timezone.utc) if o.last_modified else None,
                )
                for o in objects
            ]
        except self._errors() as exc:
            raise ArchiveUnavailableError(
                f"Cannot list s3://{self.bucket}/{prefix} at {self.host}: {exc}"
            ) from exc

    def download(self, key: str, destination: Path) -> None:
        """Download ``key`` to ``destination``. Wraps backend errors in ArchiveUnavailableError."""
        try:
            self._client.fget_object(self.bucket, key, str(destination))
        except self._errors() as exc:
            raise ArchiveUnavailableError(f"Cannot download s3://{self.bucket}/{key}: {exc}") from exc

    def url(self, key: str) -> str:
        """HTTPS URL of an object."""
        return f"https://{self.host}/{self.bucket}/{key}"


@dataclass
class ObservationRecord:
    """A native TART observation file found in an archive.

    Attributes
    ----------
    filename_utc
        Observation time encoded in the filename (the last-sample time, UTC).
    start_utc, end_utc
        First and last sample times. These stay ``None`` until the file has
        been read; they are never guessed from the filename.
    """

    station_id: str
    tart_id: Optional[str]
    archive_name: str
    filename: str
    remote_identifier: str
    archive_url: str
    filename_utc: datetime
    format: str
    size_bytes: Optional[int]
    etag: Optional[str]
    last_modified_utc: Optional[datetime]
    start_utc: Optional[datetime] = None
    end_utc: Optional[datetime] = None

    def to_dict(self) -> dict:
        """JSON-serialisable representation."""
        d = asdict(self)
        for key in ("filename_utc", "last_modified_utc", "start_utc", "end_utc"):
            d[key] = format_utc(d[key], 6) if d[key] else None
        return d


def parse_filename_utc(filename: str) -> Optional[datetime]:
    """Return the UTC observation time encoded in a TART ``vis_*.hdf`` filename.

    Returns ``None`` if the name does not follow the archive convention.
    """
    m = _FILENAME_RE.search(filename)
    if not m:
        return None
    y, mo, d, h, mi, s, frac, _ = m.groups()
    micro = int((frac or "0").ljust(6, "0"))
    return datetime(int(y), int(mo), int(d), int(h), int(mi), int(s), micro, tzinfo=timezone.utc)


def day_prefixes(archive_name: str, start_utc: datetime, end_utc: datetime) -> List[str]:
    """Archive day prefixes ``<name>/vis/<Y>/<M>/<D>/`` covering ``[start, end]``.

    Months and days are not zero-padded, matching the archive.
    """
    out = []
    day: date = start_utc.astimezone(timezone.utc).date()
    last: date = end_utc.astimezone(timezone.utc).date()
    while day <= last:
        out.append(f"{archive_name}/vis/{day.year}/{day.month}/{day.day}/")
        day += timedelta(days=1)
    return out


def find_observations(
    station: Station,
    start_utc: datetime,
    end_utc: datetime,
    archive: Optional[ArchiveBackend] = None,
    file_span_margin_s: float = DEFAULT_FILE_SPAN_MARGIN_S,
) -> List[ObservationRecord]:
    """Find TART observations that may intersect a UTC interval.

    Parameters
    ----------
    station
        Validated station. ``station.archive_name`` must be set.
    start_utc, end_utc
        UTC query boundaries (aware datetimes).
    archive
        Archive backend. Defaults to the public TART S3 archive.
    file_span_margin_s
        Files whose filename time is up to this many seconds after
        ``end_utc`` are also returned, because a file ends at its filename time.

    Returns
    -------
    list of ObservationRecord
        Candidate files in time order. Sample times are confirmed later by
        reading the files. An empty list means nothing was found.

    Raises
    ------
    ArchiveUnavailableError
        If the station has no archive name, or the archive cannot be queried.
    """
    if not station.archive_name:
        raise ArchiveUnavailableError(
            f"Station '{station.id}' has no 'archive_name' in the station configuration"
        )
    if end_utc <= start_utc:
        raise ValueError("end_utc must be after start_utc")
    archive = archive or TartS3Archive()
    upper = end_utc + timedelta(seconds=file_span_margin_s)
    records: List[ObservationRecord] = []
    for prefix in day_prefixes(station.archive_name, start_utc, upper):
        logger.debug("Listing %s", prefix)
        for obj in archive.list_objects(prefix):
            filename = obj.key.rsplit("/", 1)[-1]
            t = parse_filename_utc(filename)
            if t is None:
                logger.debug("Ignoring non-visibility object %s", obj.key)
                continue
            if start_utc <= t <= upper:
                records.append(
                    ObservationRecord(
                        station_id=station.id, tart_id=station.tart_id,
                        archive_name=station.archive_name, filename=filename,
                        remote_identifier=obj.key, archive_url=archive.url(obj.key),
                        filename_utc=t, format="hdf5", size_bytes=obj.size_bytes,
                        etag=obj.etag, last_modified_utc=obj.last_modified_utc,
                    )
                )
    records.sort(key=lambda r: r.filename_utc)
    logger.info(
        "%s (%s): %d archive files between %s and %s (+%.0f s margin)",
        station.id, station.archive_name, len(records),
        format_utc(start_utc), format_utc(end_utc), file_span_margin_s,
    )
    return records
