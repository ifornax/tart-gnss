"""TLE catalogue handling: parsing, constellation classification, provenance.

All network access of the package is isolated in :func:`download_gnss_catalogue`.

A downloaded catalogue is written once to ``<tle_dir>/gnss_<retrieval>.tle``
together with a sidecar ``<file>.meta.json`` recording the source URLs,
retrieval time, SHA-256 hash and the constellation of every NORAD ID. Files
are never overwritten, so a historical catalogue can always be re-used with
``--tle-file``.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from sgp4.api import Satrec

logger = logging.getLogger(__name__)

#: Constellations known to the classifier (canonical spelling).
CONSTELLATIONS: Tuple[str, ...] = ("GPS", "Galileo", "BeiDou", "GLONASS")

#: CelesTrak GP groups used for each constellation.
CELESTRAK_GROUPS: Dict[str, str] = {
    "GPS": "gps-ops",
    "Galileo": "galileo",
    "BeiDou": "beidou",
    "GLONASS": "glo-ops",
}

#: CelesTrak GP query URL template (TLE format).
CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=tle"

UNKNOWN_CONSTELLATION = "UNKNOWN"

# Name heuristics, used only when no sidecar metadata is available.
# The GLONASS pattern matches CelesTrak's "COSMOS <n> (<GLONASS slot no.>)"
# naming; plain "COSMOS" names are NOT assumed to be GLONASS.
_NAME_PATTERNS: Tuple[Tuple[str, "re.Pattern[str]"], ...] = (
    ("GPS", re.compile(r"\bGPS\b|NAVSTAR", re.IGNORECASE)),
    ("Galileo", re.compile(r"GALILEO|^GSAT0\d{3}", re.IGNORECASE)),
    ("BeiDou", re.compile(r"BEIDOU|COMPASS", re.IGNORECASE)),
    ("GLONASS", re.compile(r"GLONASS|^COSMOS\s+\d+\s+\(\d+[A-Z]?\)\s*$", re.IGNORECASE)),
)


class TLEError(ValueError):
    """Raised for malformed TLE input or inconsistent TLE provenance."""


@dataclass(frozen=True)
class TLERecord:
    """A single satellite's two-line element set and its provenance.

    Attributes
    ----------
    name
        Satellite name (title line, stripped).
    norad_id
        NORAD catalogue number (from TLE line 1).
    line1, line2
        The TLE lines exactly as read.
    constellation
        GNSS constellation (``GPS``, ``Galileo``, ``BeiDou``, ``GLONASS``) or
        ``UNKNOWN``.
    epoch_utc
        TLE epoch as a timezone-aware UTC datetime.
    source
        Free-text TLE source (URL or file path).
    """

    name: str
    norad_id: int
    line1: str
    line2: str
    constellation: str
    epoch_utc: datetime
    source: str

    def to_dict(self) -> Dict[str, object]:
        """Return a JSON-serialisable representation."""
        return {
            "satellite_name": self.name,
            "norad_id": self.norad_id,
            "constellation": self.constellation,
            "epoch_utc": format_utc(self.epoch_utc),
            "line1": self.line1,
            "line2": self.line2,
            "source": self.source,
        }


@dataclass
class TLECatalogue:
    """A set of TLE records with file-level provenance.

    Attributes
    ----------
    records
        Satellites in file order (duplicate NORAD IDs removed).
    path
        Path of the TLE file actually used.
    sha256
        SHA-256 hex digest of that file.
    source
        Origin of the data (URL(s) or ``"local file"``).
    retrieved_utc
        Retrieval time of a downloaded catalogue (``None`` if unknown).
    metadata_path
        Sidecar metadata file, if one was found/written.
    """

    records: List[TLERecord]
    path: Path
    sha256: str
    source: str
    retrieved_utc: Optional[datetime] = None
    metadata_path: Optional[Path] = None
    extra: Dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        """Return the file-level provenance as a JSON-serialisable dict."""
        return {
            "tle_file": str(self.path),
            "tle_sha256": self.sha256,
            "tle_source": self.source,
            "tle_retrieved_utc": format_utc(self.retrieved_utc) if self.retrieved_utc else None,
            "tle_metadata_file": str(self.metadata_path) if self.metadata_path else None,
            "n_satellites": len(self.records),
        }


def format_utc(dt: datetime, decimals: int = 3) -> str:
    """Format an aware datetime as ISO-8601 UTC with a ``Z`` suffix.

    Parameters
    ----------
    dt
        Timezone-aware datetime (naive datetimes are rejected).
    decimals
        Number of fractional-second digits (0-6).
    """
    if dt.tzinfo is None:
        raise ValueError("Refusing to format a timezone-naive datetime.")
    dt = dt.astimezone(timezone.utc)
    text = dt.strftime("%Y-%m-%dT%H:%M:%S")
    if decimals > 0:
        text += "." + f"{dt.microsecond:06d}"[:decimals]
    return text + "Z"


def sha256_file(path: Union[str, Path]) -> str:
    """Return the SHA-256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tle_checksum(line: str) -> int:
    """Compute the modulo-10 checksum of a TLE line (digits + 1 per minus sign)."""
    total = 0
    for char in line[:68]:
        if char.isdigit():
            total += int(char)
        elif char == "-":
            total += 1
    return total % 10


def tle_epoch_utc(line1: str, line2: str) -> datetime:
    """Return the epoch of a TLE as an aware UTC datetime (via ``sgp4``)."""
    sat = Satrec.twoline2rv(line1, line2)
    jd = sat.jdsatepoch + sat.jdsatepochF
    return datetime(2000, 1, 1, 12, tzinfo=timezone.utc) + timedelta(days=jd - 2451545.0)


def classify_constellation(name: str) -> str:
    """Classify a satellite into a GNSS constellation from its name.

    This is a heuristic used only when no explicit constellation metadata is
    available (e.g. a user-supplied TLE file without sidecar). Returns
    ``"UNKNOWN"`` if no pattern matches.
    """
    for constellation, pattern in _NAME_PATTERNS:
        if pattern.search(name.strip()):
            return constellation
    return UNKNOWN_CONSTELLATION


def normalise_constellation(name: str) -> str:
    """Map a user-supplied constellation name to its canonical spelling."""
    for canonical in CONSTELLATIONS + (UNKNOWN_CONSTELLATION,):
        if name.strip().lower() == canonical.lower():
            return canonical
    raise TLEError(
        f"Unknown constellation '{name}'. Choose from: {', '.join(CONSTELLATIONS)}"
    )


def parse_tle_text(text: str) -> List[Tuple[str, str, str]]:
    """Parse TLE text into ``(name, line1, line2)`` tuples.

    Accepts three-line (title + two lines) and bare two-line sets. Lines are
    checked for structure, checksum and NORAD-ID consistency.

    Raises
    ------
    TLEError
        On any malformed element set.
    """
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    out: List[Tuple[str, str, str]] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            name, l1, l2 = None, lines[i], lines[i + 1]
            i += 2
        elif i + 2 < len(lines) and lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
            name, l1, l2 = lines[i].strip(), lines[i + 1], lines[i + 2]
            i += 3
        else:
            raise TLEError(f"Malformed TLE near line {i + 1}: {lines[i]!r}")
        for label, line in (("1", l1), ("2", l2)):
            if len(line) < 69:
                raise TLEError(f"TLE line {label} too short ({len(line)} chars): {line!r}")
            if not line[68].isdigit() or tle_checksum(line) != int(line[68]):
                raise TLEError(f"TLE line {label} checksum mismatch: {line!r}")
        if l1[2:7].strip() != l2[2:7].strip():
            raise TLEError(f"NORAD IDs differ between lines: {l1!r} / {l2!r}")
        if name is None:
            name = f"NORAD {int(l1[2:7])}"
        out.append((name, l1, l2))
    return out


def _metadata_path(tle_path: Path) -> Path:
    return tle_path.with_name(tle_path.name + ".meta.json")


def load_tle_file(path: Union[str, Path]) -> TLECatalogue:
    """Load a TLE catalogue exactly as stored on disk.

    If a sidecar ``<file>.meta.json`` exists (written by
    :func:`download_gnss_catalogue`), its source, retrieval time and
    per-satellite constellation labels are used; its recorded SHA-256 must
    match the file. Otherwise constellations are inferred from names with
    :func:`classify_constellation`.

    Parameters
    ----------
    path
        TLE text file.

    Returns
    -------
    TLECatalogue
    """
    path = Path(path)
    if not path.is_file():
        raise TLEError(f"TLE file not found: {path}")
    digest = sha256_file(path)
    triples = parse_tle_text(path.read_text(encoding="utf-8", errors="strict"))
    if not triples:
        raise TLEError(f"No TLEs found in {path}")

    meta_path = _metadata_path(path)
    meta: Dict[str, object] = {}
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("sha256") and meta["sha256"] != digest:
            raise TLEError(
                f"SHA-256 of {path} does not match its metadata {meta_path}; "
                "the catalogue has been modified since retrieval."
            )
        logger.info("Using TLE provenance metadata %s", meta_path)
    else:
        meta_path = None
        logger.info(
            "No sidecar metadata for %s; constellations inferred from satellite names",
            path,
        )

    constellation_map: Mapping[str, str] = meta.get("constellations", {})  # type: ignore[assignment]
    source = str(meta.get("source", f"local file {path}"))
    retrieved = meta.get("retrieved_utc")
    retrieved_dt = (
        datetime.fromisoformat(str(retrieved).replace("Z", "+00:00")) if retrieved else None
    )

    records: List[TLERecord] = []
    seen: Dict[int, str] = {}
    for name, l1, l2 in triples:
        norad = int(l1[2:7])
        if norad in seen:
            logger.warning(
                "Duplicate NORAD ID %d (%s / %s) in %s; keeping the first entry",
                norad, seen[norad], name, path,
            )
            continue
        seen[norad] = name
        constellation = constellation_map.get(str(norad)) or classify_constellation(name)
        records.append(
            TLERecord(
                name=name,
                norad_id=norad,
                line1=l1,
                line2=l2,
                constellation=constellation,
                epoch_utc=tle_epoch_utc(l1, l2),
                source=source,
            )
        )
    return TLECatalogue(
        records=records,
        path=path,
        sha256=digest,
        source=source,
        retrieved_utc=retrieved_dt,
        metadata_path=meta_path,
    )


def filter_by_constellation(
    records: Iterable[TLERecord], constellations: Optional[Sequence[str]]
) -> List[TLERecord]:
    """Keep only records whose constellation is in ``constellations``.

    ``None`` keeps every record. Names are matched case-insensitively.
    """
    records = list(records)
    if not constellations:
        return records
    wanted = {normalise_constellation(c) for c in constellations}
    return [r for r in records if r.constellation in wanted]


def _http_get(url: str, timeout: float, attempts: int = 4, backoff_s: float = 10.0) -> str:
    """HTTP GET with retries on transient failures (exponential back-off)."""
    request = urllib.request.Request(url, headers={"User-Agent": "tart-gnss (TART research)"})
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == attempts:
                raise
            wait = backoff_s * 2 ** (attempt - 1)
            logger.warning("GET %s failed (%s); retry %d/%d in %.0f s",
                           url, exc, attempt, attempts - 1, wait)
            time.sleep(wait)
    raise AssertionError("unreachable")


def download_gnss_catalogue(
    out_dir: Union[str, Path],
    constellations: Sequence[str] = CONSTELLATIONS,
    timeout: float = 60.0,
    fetch=_http_get,
) -> Path:
    """Download current GNSS TLEs from CelesTrak and archive them locally.

    Parameters
    ----------
    out_dir
        Directory receiving ``gnss_<YYYY-MM-DDTHHMMSSZ>.tle`` and its
        ``.meta.json`` sidecar.
    constellations
        Constellations to download (see :data:`CELESTRAK_GROUPS`).
    timeout
        Per-request timeout in seconds.
    fetch
        Callable ``(url, timeout) -> str``; injectable for testing.

    Returns
    -------
    pathlib.Path
        Path to the new TLE file.

    Notes
    -----
    CelesTrak serves only the *current* element sets. For reconstructing
    past observations, archived catalogues (or historical TLEs from e.g.
    Space-Track) should be supplied with ``--tle-file``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    retrieved = datetime.now(timezone.utc).replace(microsecond=0)
    stamp = retrieved.strftime("%Y-%m-%dT%H%M%SZ")
    tle_path = out_dir / f"gnss_{stamp}.tle"
    if tle_path.exists():
        raise TLEError(f"Refusing to overwrite existing TLE catalogue {tle_path}")

    chunks: List[str] = []
    constellation_map: Dict[str, str] = {}
    urls: Dict[str, str] = {}
    for constellation in constellations:
        canonical = normalise_constellation(constellation)
        url = CELESTRAK_URL.format(group=CELESTRAK_GROUPS[canonical])
        logger.info("Downloading %s TLEs from %s", canonical, url)
        text = fetch(url, timeout)
        triples = parse_tle_text(text)
        if not triples:
            raise TLEError(f"No TLEs returned by {url}: {text[:200]!r}")
        for name, l1, l2 in triples:
            constellation_map.setdefault(str(int(l1[2:7])), canonical)
            chunks.append(f"{name}\n{l1}\n{l2}\n")
        urls[canonical] = url
        logger.info("  %d %s element sets", len(triples), canonical)

    tle_path.write_text("".join(chunks), encoding="utf-8")
    meta = {
        "source": "CelesTrak GP data: " + ", ".join(urls.values()),
        "source_urls": urls,
        "retrieved_utc": format_utc(retrieved, decimals=0),
        "sha256": sha256_file(tle_path),
        "constellations": constellation_map,
        "notes": "Current element sets at retrieval time. Do not edit the .tle file; "
        "its SHA-256 is verified on load.",
    }
    _metadata_path(tle_path).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    logger.info("Saved TLE catalogue to %s", tle_path)
    return tle_path
