"""Matching every Stage-1 visibility window with the shared archive candidate segments.

The Stage-1 JSON is read with :func:`tart_gnss.stage2.selection.load_stage1_results`
and hashed once. Windows are grouped by NORAD ID and ordered by their
``start_utc`` string exactly as :func:`tart_gnss.stage2.selection.select_window`
does, so ``window_index`` agrees with Stage 2. Signal relevance comes from
:func:`tart_gnss.stage2.selection.signal_relevance`.

Matching bisects each window into the sorted, disjoint segment list
(``O(W log S + K)`` for ``W`` windows, ``S`` segments and ``K`` matches); it
never scans windows against individual files. A long window that meets
several separated segments gives one candidate per segment: they are never
joined into a fabricated continuous overlap.
"""

from __future__ import annotations

import bisect
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from tart_gnss.archive_scan.intervals import CommonSegment, count_files
from tart_gnss.stage2.selection import parse_utc, signal_relevance

logger = logging.getLogger(__name__)

#: Signal relevance values that place the satellite inside the TART band.
IN_BAND = ("GPS_L1", "GALILEO_E1", "BEIDOU_B1C")


@dataclass
class VisibilityWindow:
    """One Stage-1 window, reduced to what matching and ranking need (copied verbatim)."""

    rank: int
    window_index: int
    n_windows_for_satellite: int
    satellite_name: str
    norad_id: int
    constellation: str
    start_utc: datetime
    end_utc: datetime
    signal_relevance: str
    signal_reason: str
    tle_epoch_utc: Optional[datetime]
    tle_epoch_text: Optional[str]
    min_mid_elevation_deg: Optional[float]
    clipped: bool

    @property
    def duration_s(self) -> float:
        """Predicted window duration, seconds."""
        return (self.end_utc - self.start_utc).total_seconds()


def _min_mid_elevation(window: Mapping[str, Any]) -> Optional[float]:
    try:
        values = [float(s["mid"]["elevation_deg"]) for s in window["stations"].values()]
    except (KeyError, TypeError, ValueError, AttributeError):
        return None
    return min(values) if values else None


def load_windows(doc: Mapping[str, Any]) -> List[VisibilityWindow]:
    """All Stage-1 windows of a loaded results document, in Stage-1 rank order.

    The TLE epoch is the per-satellite element-set epoch recorded with each
    window (``window["tle"]["epoch"]``), never the catalogue retrieval time;
    it is ``None`` when missing or unparsable.
    """
    by_norad: Dict[int, List[Mapping[str, Any]]] = defaultdict(list)
    for w in doc["windows"]:
        by_norad[int(w["norad_id"])].append(w)
    out: List[VisibilityWindow] = []
    for norad, wins in by_norad.items():
        wins = sorted(wins, key=lambda w: w["start_utc"])
        for index, w in enumerate(wins, start=1):
            name = str(w["satellite"]).strip()
            relevance, reason = signal_relevance(str(w["constellation"]), name)
            epoch_text = (w.get("tle") or {}).get("epoch")
            try:
                epoch = parse_utc(epoch_text) if epoch_text else None
            except ValueError:
                epoch = None
            out.append(VisibilityWindow(
                rank=int(w["rank"]), window_index=index, n_windows_for_satellite=len(wins),
                satellite_name=name, norad_id=norad, constellation=str(w["constellation"]),
                start_utc=parse_utc(w["start_utc"]), end_utc=parse_utc(w["end_utc"]),
                signal_relevance=relevance, signal_reason=reason, tle_epoch_utc=epoch,
                tle_epoch_text=epoch_text, min_mid_elevation_deg=_min_mid_elevation(w),
                clipped=bool(w.get("start_clipped_by_search_interval") or w.get("end_clipped_by_search_interval")),
            ))
    out.sort(key=lambda v: (v.rank, v.norad_id, v.start_utc))
    return out


def filter_windows(
    windows: Sequence[VisibilityWindow],
    relevance: Optional[Sequence[str]],
    constellations: Optional[Sequence[str]],
) -> Tuple[List[VisibilityWindow], Dict[str, int]]:
    """Apply the signal-relevance and constellation filters; returns kept windows and exclusion counts."""
    rel = set(relevance) if relevance else None
    cons = {c.lower() for c in constellations} if constellations else None
    kept: List[VisibilityWindow] = []
    excluded: Dict[str, int] = defaultdict(int)
    for w in windows:
        if rel is not None and w.signal_relevance not in rel:
            excluded[f"signal_relevance:{w.signal_relevance}"] += 1
        elif cons is not None and w.constellation.lower() not in cons:
            excluded[f"constellation:{w.constellation}"] += 1
        else:
            kept.append(w)
    return kept, dict(sorted(excluded.items()))


@dataclass
class Candidate:
    """The intersection of one Stage-1 window with one shared archive candidate segment."""

    window: VisibilityWindow
    segment: CommonSegment
    overlap_start_utc: datetime
    overlap_end_utc: datetime
    file_counts: Dict[str, int]
    tle_age_days: Optional[float]
    flags: List[str] = field(default_factory=list)
    candidate_id: str = ""
    order: int = 0
    ranking_components: Dict[str, object] = field(default_factory=dict)

    @property
    def duration_s(self) -> float:
        """Possible overlap duration (candidate, not verified), seconds."""
        return (self.overlap_end_utc - self.overlap_start_utc).total_seconds()

    @property
    def mid_utc(self) -> datetime:
        """Overlap midpoint."""
        return self.overlap_start_utc + (self.overlap_end_utc - self.overlap_start_utc) / 2


def tle_age_days(epoch: Optional[datetime], at_utc: datetime) -> Optional[float]:
    """Signed TLE age ``at - epoch`` in days (the Stage-1 ``age_*_days`` convention)."""
    if epoch is None:
        return None
    return (at_utc - epoch).total_seconds() / 86400.0


def match_windows(
    windows: Sequence[VisibilityWindow],
    segments: Sequence[CommonSegment],
    file_times: Mapping[str, Sequence[datetime]],
    span_s: float,
    min_duration_s: float,
) -> Tuple[List[Candidate], int]:
    """Intersect every window with the sorted, disjoint shared segments.

    Parameters
    ----------
    windows
        Filtered Stage-1 windows (any order).
    segments
        Shared candidate segments, sorted and disjoint.
    file_times
        Station id -> sorted filename times (for per-overlap file counts).
    span_s
        File-span assumption, used for the file-count lookup envelope.
    min_duration_s
        Minimum possible overlap duration to keep (seconds, inclusive).

    Returns
    -------
    (candidates, n_short)
        Candidates (unranked) and the number of positive-length intersections
        shorter than ``min_duration_s``.
    """
    ends = [s.end_utc for s in segments]
    out: List[Candidate] = []
    n_short = 0
    for w in windows:
        i = bisect.bisect_right(ends, w.start_utc)
        while i < len(segments) and segments[i].start_utc < w.end_utc:
            seg = segments[i]
            a, b = max(seg.start_utc, w.start_utc), min(seg.end_utc, w.end_utc)
            i += 1
            if b <= a:
                continue
            if (b - a).total_seconds() < min_duration_s:
                n_short += 1
                continue
            counts = {sid: count_files(file_times[sid], a, b, span_s) for sid in seg.station_ids}
            age = tle_age_days(w.tle_epoch_utc, a + (b - a) / 2)
            flags: List[str] = []
            if age is None:
                flags.append("TLE_AGE_MISSING")
            if w.signal_relevance not in IN_BAND:
                flags.append(f"SIGNAL_{w.signal_relevance}")
            if w.clipped:
                flags.append("STAGE1_WINDOW_CLIPPED")
            if not seg.inventory_complete:
                flags.append("INVENTORY_INCOMPLETE")
            if not seg.all_stations:
                flags.append("PARTIAL_NETWORK_DIAGNOSTIC")
            out.append(Candidate(window=w, segment=seg, overlap_start_utc=a, overlap_end_utc=b,
                                 file_counts=counts, tle_age_days=age, flags=flags))
    return out, n_short
