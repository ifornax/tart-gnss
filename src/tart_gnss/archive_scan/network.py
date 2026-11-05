"""Network-wide archive discovery: eligible stations, physical sites, baselines, sweep and triads.

Used by ``tart-archive-scan --all-stations``. The station set is derived
from the registry at run time; nothing about the network is hard-coded.

Four concepts are kept apart (see ``docs/archive_scan.md``): a *registered*
station (a registry entry), a *filename-derived archive candidate* (S3 keys
suggest an interval), *native-validated* simultaneous data (Stage 2) and a
*detected* satellite (Stage 3). Everything here is the second kind only.

Algorithm
---------
1. Eligibility: a registry entry needs valid coordinates and an
   ``archive_name``; duplicated archive names are listed once (the first
   station keeps it). Excluded entries are reported with a reason.
2. Physical sites: stations whose geodesic separation is at most the site
   tolerance (default 1 km) are joined by single linkage into one *site
   group*; co-located receivers never count as independent sites.
3. Sweep: one event sweep over all stations' candidate segments
   (:func:`network_active_segments`) gives maximal intervals with a constant
   active station set; those with at least ``min_sites`` distinct site
   groups are the *network active segments*.
4. Station sets ("triads" for size 3): only sets of stations from distinct
   site groups that are co-active in some network active segment are
   considered. Each set's contiguous candidate intervals are the
   intersection of its members' own segments (:func:`intersect_all`), so a
   triad stays one segment while a fourth station joins or leaves, and a
   set is never credited with coverage that only a co-located sibling has.
5. Baselines are WGS84 ellipsoidal geodesics (Vincenty inverse); the ECEF
   chord is reported separately. A baseline cut is a coarse site-diversity
   filter, not a triangulation-quality metric.
"""

from __future__ import annotations

import itertools
import logging
import math
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import yaml

from tart_gnss.archive_scan.intervals import Interval, intersect_all, k_station_segments
from tart_gnss.stage2.synchronise import total_seconds
from tart_gnss.stage3.geometry import geodetic_to_ecef
from tart_gnss.stations import Station, StationConfigError, parse_station_config

logger = logging.getLogger(__name__)

#: WGS84 ellipsoid.
WGS84_A_M = 6378137.0
WGS84_F = 1.0 / 298.257223563

#: Default separation (km) below which stations are one physical site.
DEFAULT_SITE_TOLERANCE_KM = 1.0

BASELINE_METHOD = "WGS84 ellipsoidal geodesic (Vincenty inverse), km; ECEF chord reported separately"

#: Network ranking rule, recorded in the manifest. Revise the version when the rule changes.
NETWORK_RANKING_VERSION = "1"
NETWORK_RANKING_RULE = [
    "1. inventory_complete true first",
    "2. possible_shared_duration_s descending (union of the set's contiguous segments >= min duration)",
    "3. longest_contiguous_duration_s descending",
    "4. min_pairwise_baseline_km descending (secondary; NOT a triangulation-quality proxy)",
    "5. station_ids lexicographic (stable tie-break)",
]

# Exclusion reasons.
NO_ARCHIVE_NAME = "NO_ARCHIVE_NAME"
INVALID_COORDINATES = "INVALID_COORDINATES"
DUPLICATE_ARCHIVE_NAME = "DUPLICATE_ARCHIVE_NAME"
NOT_IN_INCLUDE_LIST = "NOT_IN_INCLUDE_LIST"
EXCLUDED_BY_USER = "EXCLUDED_BY_USER"


class GeodesicError(ArithmeticError):
    """The Vincenty iteration did not converge (nearly antipodal points)."""


# --------------------------------------------------------------------------- geodesy

def geodesic_distance_km(lat1: float, lon1: float, lat2: float, lon2: float,
                         a: float = WGS84_A_M, f: float = WGS84_F, tol: float = 1e-12, max_iter: int = 200) -> float:
    """Ellipsoidal surface distance between two geodetic points (Vincenty 1975 inverse), km.

    Parameters
    ----------
    lat1, lon1, lat2, lon2
        Geodetic latitude/longitude in degrees.
    a, f
        Ellipsoid semi-major axis (m) and flattening (default WGS84).

    Raises
    ------
    GeodesicError
        If the iteration does not converge (only for nearly antipodal
        points); never returns a silently wrong value.
    """
    if lat1 == lat2 and lon1 == lon2:
        return 0.0
    b = a * (1.0 - f)
    L = math.radians(lon2 - lon1)
    U1 = math.atan((1.0 - f) * math.tan(math.radians(lat1)))
    U2 = math.atan((1.0 - f) * math.tan(math.radians(lat2)))
    sinU1, cosU1, sinU2, cosU2 = math.sin(U1), math.cos(U1), math.sin(U2), math.cos(U2)
    lam = L
    for _ in range(max_iter):
        sin_lam, cos_lam = math.sin(lam), math.cos(lam)
        sin_sigma = math.hypot(cosU2 * sin_lam, cosU1 * sinU2 - sinU1 * cosU2 * cos_lam)
        if sin_sigma == 0.0:
            return 0.0
        cos_sigma = sinU1 * sinU2 + cosU1 * cosU2 * cos_lam
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cosU1 * cosU2 * sin_lam / sin_sigma
        cos2_alpha = 1.0 - sin_alpha ** 2
        cos_2sm = cos_sigma - 2.0 * sinU1 * sinU2 / cos2_alpha if cos2_alpha != 0.0 else 0.0
        C = f / 16.0 * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha))
        lam_prev = lam
        lam = L + (1.0 - C) * f * sin_alpha * (
            sigma + C * sin_sigma * (cos_2sm + C * cos_sigma * (-1.0 + 2.0 * cos_2sm ** 2)))
        if abs(lam - lam_prev) < tol:
            break
    else:
        raise GeodesicError(f"Vincenty inverse did not converge for ({lat1}, {lon1}) -> ({lat2}, {lon2})")
    u2 = cos2_alpha * (a * a - b * b) / (b * b)
    A = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)))
    B = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)))
    d_sigma = B * sin_sigma * (cos_2sm + B / 4.0 * (
        cos_sigma * (-1.0 + 2.0 * cos_2sm ** 2)
        - B / 6.0 * cos_2sm * (-3.0 + 4.0 * sin_sigma ** 2) * (-3.0 + 4.0 * cos_2sm ** 2)))
    return b * A * (sigma - d_sigma) / 1000.0


def chord_distance_km(s1: Station, s2: Station) -> float:
    """Straight-line (3-D ECEF) distance between two stations, km (includes heights)."""
    p1 = geodetic_to_ecef(s1.latitude_deg, s1.longitude_deg, s1.altitude_m)
    p2 = geodetic_to_ecef(s2.latitude_deg, s2.longitude_deg, s2.altitude_m)
    return float(math.dist(p1, p2)) / 1000.0


class Baselines:
    """Pairwise geodesic and chord distances, each pair computed once."""

    def __init__(self, stations: Mapping[str, Station]):
        self._stations = stations
        self._geo: Dict[Tuple[str, str], float] = {}
        self._chord: Dict[Tuple[str, str], float] = {}

    @staticmethod
    def _key(a: str, b: str) -> Tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    def geodesic_km(self, a: str, b: str) -> float:
        """Geodesic distance between stations ``a`` and ``b`` (km)."""
        k = self._key(a, b)
        if k not in self._geo:
            s1, s2 = self._stations[k[0]], self._stations[k[1]]
            self._geo[k] = geodesic_distance_km(s1.latitude_deg, s1.longitude_deg, s2.latitude_deg, s2.longitude_deg)
        return self._geo[k]

    def chord_km(self, a: str, b: str) -> float:
        """ECEF chord distance between stations ``a`` and ``b`` (km)."""
        k = self._key(a, b)
        if k not in self._chord:
            self._chord[k] = chord_distance_km(self._stations[k[0]], self._stations[k[1]])
        return self._chord[k]

    def pairs(self, ids: Sequence[str]) -> List[Tuple[str, str, float, float]]:
        """``(a, b, geodesic_km, chord_km)`` for every pair of ``ids``."""
        return [(a, b, self.geodesic_km(a, b), self.chord_km(a, b)) for a, b in itertools.combinations(ids, 2)]


# --------------------------------------------------------------------------- eligibility

@dataclass
class Exclusion:
    """A registry entry that cannot take part in network discovery."""

    station_id: str
    tart_id: Optional[str]
    reason: str
    detail: str = ""
    archive_name: Optional[str] = None
    archive_verified: Optional[bool] = None


def load_registry_tolerant(path: Path) -> Tuple[Dict[str, Station], List[Exclusion]]:
    """Load every valid station of a catalogue; invalid entries become ``INVALID_COORDINATES`` exclusions.

    Entries already skipped by ``tart-stations-update`` (``metadata.skipped_entries``
    of a generated registry) are reported as ``INVALID_COORDINATES`` too.
    Coordinates are never synthesised.

    Raises
    ------
    StationConfigError
        If the file is missing or has no ``stations`` mapping.
    """
    path = Path(path)
    if not path.is_file():
        raise StationConfigError(f"Station configuration file not found: {path}")
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise StationConfigError(f"Cannot parse {path}: {exc}") from exc
    if not isinstance(data, Mapping) or not isinstance(data.get("stations"), Mapping):
        raise StationConfigError(f"{path} must contain a 'stations' mapping")
    stations: Dict[str, Station] = {}
    excluded: List[Exclusion] = []
    for sid, entry in data["stations"].items():
        try:
            stations.update(parse_station_config({"stations": {sid: entry}}))
        except StationConfigError as exc:
            entry = entry if isinstance(entry, Mapping) else {}
            excluded.append(Exclusion(str(sid), entry.get("tart_id"), INVALID_COORDINATES, str(exc),
                                      entry.get("archive_name"), entry.get("archive_verified")))
    meta = data.get("metadata") if isinstance(data.get("metadata"), Mapping) else {}
    for skipped in meta.get("skipped_entries") or []:
        tid = str(skipped.get("tart_id") or f"index {skipped.get('index')}")
        excluded.append(Exclusion(tid, skipped.get("tart_id"), INVALID_COORDINATES,
                                  "skipped by tart-stations-update: " + str(skipped.get("reason", ""))))
    return stations, excluded


def select_eligible(catalogue: Mapping[str, Station], include: Optional[Sequence[Station]] = None,
                    exclude: Optional[Sequence[Station]] = None) -> Tuple[List[Station], List[Exclusion]]:
    """Stations usable for network discovery, in catalogue-id order, and the excluded ones with reasons.

    ``archive_verified`` is recorded but never treated as evidence of data.
    """
    include_ids = None if include is None else {s.id for s in include}
    exclude_ids = {s.id for s in exclude or []}
    eligible: List[Station] = []
    excluded: List[Exclusion] = []
    archive_owner: Dict[str, str] = {}
    for sid in sorted(catalogue):
        st = catalogue[sid]

        def drop(reason: str, detail: str = "") -> None:
            excluded.append(Exclusion(st.id, st.tart_id, reason, detail, st.archive_name, st.archive_verified))

        if include_ids is not None and sid not in include_ids:
            drop(NOT_IN_INCLUDE_LIST)
        elif sid in exclude_ids:
            drop(EXCLUDED_BY_USER)
        elif not st.archive_name:
            drop(NO_ARCHIVE_NAME, "archive_name is null in the registry")
        elif st.archive_name in archive_owner:
            drop(DUPLICATE_ARCHIVE_NAME, f"archive '{st.archive_name}' is already listed for {archive_owner[st.archive_name]}")
        else:
            archive_owner[st.archive_name] = sid
            eligible.append(st)
    return eligible, excluded


# --------------------------------------------------------------------------- physical sites

def site_groups(stations: Sequence[Station], baselines: Baselines, tolerance_km: float) -> Dict[str, str]:
    """Station id -> site-group id, by single linkage at ``tolerance_km`` (geodesic).

    Group ids are ``SITE_<lowest station id in the group>``, so they are
    deterministic and readable. A tolerance of 0 groups only identical
    coordinates.
    """
    ids = sorted(s.id for s in stations)
    parent = {i: i for i in ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(ids, 2):
        if baselines.geodesic_km(a, b) <= tolerance_km:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    return {i: f"SITE_{find(i)}" for i in ids}


# --------------------------------------------------------------------------- sweep

@dataclass
class ActiveSegment:
    """A maximal interval with a constant set of active stations and >= ``min_sites`` site groups."""

    segment_id: int
    start_utc: datetime
    end_utc: datetime
    station_ids: Tuple[str, ...]
    site_group_ids: Tuple[str, ...]
    inventory_complete: bool = True

    @property
    def duration_s(self) -> float:
        """Possible (candidate) duration, seconds."""
        return (self.end_utc - self.start_utc).total_seconds()

    @property
    def n_sites(self) -> int:
        """Number of distinct site groups active."""
        return len(self.site_group_ids)


def network_active_segments(station_segments: Mapping[str, Sequence[Interval]], site_of: Mapping[str, str],
                            min_sites: int) -> List[ActiveSegment]:
    """One event sweep over every station's candidate segments.

    Uses :func:`~tart_gnss.archive_scan.intervals.k_station_segments` with
    ``k = 1`` (maximal intervals of constant active station set, half-open,
    ends before starts at equal instants, zero-length pieces dropped) and
    keeps the intervals whose active stations span at least ``min_sites``
    distinct site groups. Complexity O(E log E) for E boundary events.
    """
    out: List[ActiveSegment] = []
    for (a, b), ids in k_station_segments(station_segments, 1):
        sites = tuple(sorted({site_of[s] for s in ids}))
        if len(sites) >= min_sites:
            out.append(ActiveSegment(len(out) + 1, a, b, ids, sites))
    return out


# --------------------------------------------------------------------------- station sets

@dataclass
class TriadSegment:
    """One contiguous candidate interval during which every member of a station set is active."""

    start_utc: datetime
    end_utc: datetime
    file_counts: Dict[str, int] = field(default_factory=dict)
    co_active: Tuple[str, ...] = ()
    inventory_complete: bool = True
    segment_index: int = 0

    @property
    def duration_s(self) -> float:
        """Possible (candidate) duration, seconds."""
        return (self.end_utc - self.start_utc).total_seconds()


@dataclass
class Triad:
    """A set of stations at distinct physical sites with shared candidate archive coverage."""

    station_ids: Tuple[str, ...]
    site_group_ids: Tuple[str, ...]
    segments: List[TriadSegment]
    pair_baselines: List[Tuple[str, str, float, float]]
    n_short_segments: int = 0
    rank: int = 0
    passes_baseline: bool = True

    @property
    def triad_id(self) -> str:
        """Stable identifier: the sorted station ids joined by ``+``."""
        return "+".join(self.station_ids)

    @property
    def possible_shared_s(self) -> float:
        """Union of the qualifying segments (they are disjoint by construction), seconds."""
        return total_seconds([(s.start_utc, s.end_utc) for s in self.segments])

    @property
    def longest_s(self) -> float:
        """Longest qualifying contiguous segment, seconds."""
        return max((s.duration_s for s in self.segments), default=0.0)

    @property
    def inventory_complete(self) -> bool:
        """All qualifying segments have complete station/day listings."""
        return all(s.inventory_complete for s in self.segments)

    def baseline_stats(self) -> Tuple[float, float, float]:
        """(min, median, max) pairwise geodesic baseline, km."""
        vals = sorted(p[2] for p in self.pair_baselines)
        n = len(vals)
        median = vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])
        return vals[0], median, vals[-1]

    @property
    def longest_segment(self) -> Optional[TriadSegment]:
        """The longest qualifying segment (earliest on ties)."""
        return max(self.segments, key=lambda s: (s.duration_s, -s.start_utc.timestamp()), default=None)


def candidate_station_sets(active: Sequence[ActiveSegment], site_of: Mapping[str, str], size: int) -> List[Tuple[str, ...]]:
    """Station sets of ``size`` members at distinct sites that are co-active in some network active segment."""
    sets: Set[Tuple[str, ...]] = set()
    for seg in active:
        for combo in itertools.combinations(seg.station_ids, size):
            if len({site_of[s] for s in combo}) == size:
                sets.add(combo)
    return sorted(sets)


def build_triads(sets: Iterable[Tuple[str, ...]], station_segments: Mapping[str, Sequence[Interval]],
                 site_of: Mapping[str, str], baselines: Baselines, min_duration_s: float,
                 file_counter, completeness, co_active_of) -> List[Triad]:
    """Contiguous shared segments of each station set (intersection of its own members' segments).

    Parameters
    ----------
    sets
        Output of :func:`candidate_station_sets`.
    min_duration_s
        Minimum contiguous candidate duration of a segment; sets with no
        qualifying segment are dropped (their short segments are counted).
    file_counter
        ``(station_id, start, end) -> int``, files in the Stage-2 lookup envelope.
    completeness
        ``(station_ids, start, end) -> bool``, all required listings complete.
    co_active_of
        ``(start, end) -> tuple`` of every station active at some time in the interval.
    """
    triads: List[Triad] = []
    for ids in sets:
        shared = intersect_all([station_segments[s] for s in ids])
        segs, n_short = [], 0
        for a, b in shared:
            if (b - a).total_seconds() < min_duration_s:
                n_short += 1
                continue
            segs.append(TriadSegment(a, b, {s: file_counter(s, a, b) for s in ids},
                                     tuple(s for s in co_active_of(a, b) if s not in ids),
                                     completeness(ids, a, b), len(segs) + 1))
        if segs:
            triads.append(Triad(tuple(ids), tuple(site_of[s] for s in ids), segs, baselines.pairs(ids), n_short))
    return triads


def rank_triads(triads: Sequence[Triad]) -> List[Triad]:
    """Order by :data:`NETWORK_RANKING_RULE` and assign ``rank`` (1-based)."""
    def key(t: Triad):
        return (not t.inventory_complete, -t.possible_shared_s, -t.longest_s, -t.baseline_stats()[0], t.station_ids)

    ranked = sorted(triads, key=key)
    for i, t in enumerate(ranked, start=1):
        t.rank = i
    return ranked
