"""Filename-based candidate coverage and multi-station intersection.

A file named ``T`` gives the *candidate* interval ``[T - span, T]``, where
``span`` is the explicit ``--file-span-assumption-s``. The span is a search
envelope, not a verified exposure time: filename times alone never prove
that any intermediate second was sampled. Candidate intervals of a station
are unioned **without filling gaps** (two intervals merge only if they
overlap or touch), and all coverage is clipped to the half-open scan
interval ``[start, end)``. Zero-duration intersections are discarded.

Station x station intersections reuse
:func:`tart_gnss.stage2.synchronise.intersect_segments`, which needs sorted,
disjoint lists; :func:`merge_intervals` guarantees that.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from tart_gnss.stage2.synchronise import intersect_segments, total_seconds

Interval = Tuple[datetime, datetime]


def merge_intervals(intervals: Iterable[Interval]) -> List[Interval]:
    """Sort and union intervals; only overlapping or touching ones merge (gaps stay gaps).

    Zero- and negative-length intervals are dropped.
    """
    out: List[Interval] = []
    for a, b in sorted((a, b) for a, b in intervals if b > a):
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def clip(segments: Sequence[Interval], window: Interval) -> List[Interval]:
    """Clip disjoint sorted segments to ``window``, dropping empty pieces."""
    lo, hi = window
    return [(max(a, lo), min(b, hi)) for a, b in segments if min(b, hi) > max(a, lo)]


def station_candidate_segments(file_times: Sequence[datetime], span_s: float, window: Interval) -> List[Interval]:
    """Candidate coverage of one station: ``∪ [T - span, T]``, gaps kept, clipped to ``window``.

    Parameters
    ----------
    file_times
        Filename (last-sample) times of the station's files, any order.
    span_s
        Assumed maximum prior duration of a native file (seconds, > 0).
    window
        Scan interval ``[start, end)``.
    """
    if span_s <= 0:
        raise ValueError("file span assumption must be positive")
    span = timedelta(seconds=span_s)
    return clip(merge_intervals((t - span, t) for t in file_times), window)


def intersect_all(segment_lists: Sequence[Sequence[Interval]]) -> List[Interval]:
    """Intersection of several stations' disjoint sorted segment lists (all stations present)."""
    if not segment_lists:
        return []
    common: List[Interval] = list(segment_lists[0])
    for segs in segment_lists[1:]:
        common = intersect_segments(common, list(segs))
    return common


def k_station_segments(station_segments: Mapping[str, Sequence[Interval]], k: int) -> List[Tuple[Interval, Tuple[str, ...]]]:
    """Maximal intervals in which at least ``k`` stations have candidate coverage.

    A sweep over all segment boundaries gives elementary intervals with a
    constant set of covering stations; adjacent pieces with the same set are
    joined. Each returned interval therefore has one fixed station set
    (sorted ids). With ``k`` equal to the number of stations the result
    equals :func:`intersect_all`.
    """
    events: List[Tuple[datetime, int, str]] = []
    for sid, segs in station_segments.items():
        for a, b in segs:
            events.append((a, 1, sid))
            events.append((b, -1, sid))
    events.sort(key=lambda e: (e[0], e[1]))  # ends (-1) before starts (+1) at equal times
    active: Dict[str, int] = {}
    out: List[Tuple[Interval, Tuple[str, ...]]] = []
    i = 0
    prev_t: Optional[datetime] = None
    while i < len(events):
        t = events[i][0]
        if prev_t is not None and t > prev_t:
            ids = tuple(sorted(s for s, n in active.items() if n > 0))
            if len(ids) >= k:
                if out and out[-1][0][1] == prev_t and out[-1][1] == ids:
                    out[-1] = ((out[-1][0][0], t), ids)
                else:
                    out.append(((prev_t, t), ids))
        while i < len(events) and events[i][0] == t:
            _, delta, sid = events[i]
            active[sid] = active.get(sid, 0) + delta
            i += 1
        prev_t = t
    return out


def count_files(sorted_times: Sequence[datetime], start: datetime, end: datetime, span_s: float) -> int:
    """Files whose filename time lies in ``[start, end + span]`` (the Stage-2 selection rule)."""
    upper = end + timedelta(seconds=span_s)
    return bisect.bisect_right(sorted_times, upper) - bisect.bisect_left(sorted_times, start)


def gap_groups(segments: Sequence[Interval], tolerance_s: float) -> List[Dict[str, object]]:
    """Diagnostic grouping of segments separated by gaps of at most ``tolerance_s``.

    Gaps are never counted as covered: each group reports its span *and* the
    covered seconds of its member segments separately.
    """
    groups: List[List[Interval]] = []
    tol = timedelta(seconds=tolerance_s)
    for seg in segments:
        if groups and seg[0] - groups[-1][-1][1] <= tol:
            groups[-1].append(seg)
        else:
            groups.append([seg])
    return [{"start": g[0][0], "end": g[-1][1], "n_segments": len(g), "covered_s": total_seconds(g),
             "span_s": (g[-1][1] - g[0][0]).total_seconds()} for g in groups]


@dataclass
class CommonSegment:
    """One shared archive-candidate interval with a fixed covering station set."""

    segment_id: int
    start_utc: datetime
    end_utc: datetime
    station_ids: Tuple[str, ...]
    file_counts: Dict[str, int] = field(default_factory=dict)
    inventory_complete: bool = True
    all_stations: bool = True

    @property
    def duration_s(self) -> float:
        """Possible (candidate) duration, seconds; not a verified integration time."""
        return (self.end_utc - self.start_utc).total_seconds()
