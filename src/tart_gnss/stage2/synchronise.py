"""Temporal overlap, coverage of the predicted window, and cross-station timing.

Definitions (all times UTC):

* station data interval: ``[first_timestamp_i, last_timestamp_i]``;
* actual common interval: ``[max_i first_i, min_i last_i]``, valid only if
  its end is after its start;
* coverage fraction of a predicted window ``W``:
  ``T(data segments ∩ W) / T(W)``. Data segments are runs of samples without
  gaps (see :func:`tart_gnss.stage2.native.summarise_coverage`), so gaps do not
  count as covered;
* three-station coverage: ``T(∩_i segments_i ∩ W) / T(W)``.

Exact dump synchronisation is not assumed. For each station pair, the
separation between each sample and the other station's nearest sample is
characterised instead.
"""

from __future__ import annotations

from datetime import datetime
from itertools import combinations
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage2.errors import NoCommonTimeError
from tart_gnss.stage2.native import StationCoverage
from tart_gnss.tle import format_utc

Interval = Tuple[datetime, datetime]


def common_interval(intervals: Mapping[str, Interval]) -> Interval:
    """Intersection ``[max_i start_i, min_i end_i]`` of the station data intervals.

    Raises
    ------
    NoCommonTimeError
        If no interval is given, or the intersection has zero or negative length.
    """
    if not intervals:
        raise NoCommonTimeError("No station intervals to intersect")
    start = max(a for a, _ in intervals.values())
    end = min(b for _, b in intervals.values())
    if end <= start:
        detail = ", ".join(f"{k}: {format_utc(a, 0)}-{format_utc(b, 0)}" for k, (a, b) in intervals.items())
        raise NoCommonTimeError(f"Station data intervals have no common overlap ({detail})")
    return start, end


def intersect_segments(a: Sequence[Interval], b: Sequence[Interval]) -> List[Interval]:
    """Pairwise intersection of two lists of disjoint, time-ordered intervals."""
    out: List[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        start = max(a[i][0], b[j][0])
        end = min(a[i][1], b[j][1])
        if end > start:
            out.append((start, end))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def total_seconds(segments: Sequence[Interval]) -> float:
    """Total length of disjoint intervals, in seconds."""
    return float(sum((b - a).total_seconds() for a, b in segments))


def coverage_fraction(segments: Sequence[Interval], window: Interval) -> float:
    """Fraction of ``window`` covered by ``segments``, in [0, 1]."""
    length = (window[1] - window[0]).total_seconds()
    if length <= 0:
        raise ValueError("window must have positive duration")
    return total_seconds(intersect_segments(list(segments), [window])) / length


def nearest_separation(a: Sequence[datetime], b: Sequence[datetime]) -> np.ndarray:
    """For each time in ``a``, the signed offset (s) to the nearest time in ``b``.

    The value is ``t_b - t_a``. Both inputs must be sorted.
    """
    ta = np.array([t.timestamp() for t in a])
    tb = np.array([t.timestamp() for t in b])
    if ta.size == 0 or tb.size == 0:
        return np.array([])
    idx = np.clip(np.searchsorted(tb, ta), 1, tb.size - 1) if tb.size > 1 else np.zeros(ta.size, int)
    if tb.size > 1:
        left, right = tb[idx - 1], tb[idx]
        nearest = np.where(np.abs(ta - left) <= np.abs(right - ta), left, right)
    else:
        nearest = np.full(ta.size, tb[0])
    return nearest - ta


def synchronise(
    coverages: Mapping[str, StationCoverage],
    predicted_window: Interval,
    retrieved_window: Optional[Interval] = None,
) -> Dict[str, object]:
    """Compute the actual three-station overlap, coverage and timing statistics.

    Parameters
    ----------
    coverages
        Station id -> coverage, for the stations with valid data.
    predicted_window
        The Stage-1 predicted common window.
    retrieved_window
        The part of the predicted window that retrieval was asked to cover
        (shorter than ``predicted_window`` when ``--max-duration`` is used).
        Coverage fractions relative to it are reported as well.

    Returns
    -------
    dict
        Manifest-ready synchronisation section. ``actual_common_*`` are
        ``None`` and ``no_common_time_reason`` is set when there is no overlap.
    """
    ids = list(coverages)
    t_pred = (predicted_window[1] - predicted_window[0]).total_seconds()
    retrieved = retrieved_window or predicted_window
    result: Dict[str, object] = {
        "predicted_window_duration_s": t_pred,
        "retrieved_window_utc": [format_utc(retrieved[0], 3), format_utc(retrieved[1], 3)],
        "retrieved_window_duration_s": (retrieved[1] - retrieved[0]).total_seconds(),
        "station_retrieved_window_coverage_fraction": {
            sid: coverage_fraction(coverages[sid].segments, retrieved) for sid in ids
        },
        "stations": ids,
        "station_coverage_fraction": {
            sid: coverage_fraction(coverages[sid].segments, predicted_window) for sid in ids
        },
        "station_median_cadence_s": {sid: coverages[sid].median_time_spacing_seconds for sid in ids},
    }
    try:
        start, end = common_interval(
            {sid: (c.first_timestamp_utc, c.last_timestamp_utc) for sid, c in coverages.items()}
        )
    except NoCommonTimeError as exc:
        result.update({
            "actual_common_start_utc": None, "actual_common_end_utc": None,
            "actual_common_duration_seconds": 0.0, "three_station_coverage_fraction": 0.0,
            "three_station_retrieved_window_coverage_fraction": 0.0, "common_segments": [], "station_pair_time_offset_statistics": {},
            "no_common_time_reason": str(exc),
        })
        return result

    common_segments: Optional[List[Interval]] = None
    for sid in ids:
        segs = list(coverages[sid].segments)
        common_segments = segs if common_segments is None else intersect_segments(common_segments, segs)
    common_segments = common_segments or []
    pair_stats: Dict[str, Dict[str, object]] = {}
    all_abs: List[np.ndarray] = []
    for a, b in combinations(ids, 2):
        ta = [t for t in coverages[a].timestamps if start <= t <= end]
        dt = nearest_separation(ta, coverages[b].timestamps)
        if dt.size:
            all_abs.append(np.abs(dt))
            pair_stats[f"{a}__{b}"] = {
                "n_samples": int(dt.size),
                "median_offset_s": float(np.median(dt)),
                "median_abs_separation_s": float(np.median(np.abs(dt))),
                "max_abs_separation_s": float(np.max(np.abs(dt))),
            }
    joined = np.concatenate(all_abs) if all_abs else np.array([])
    result.update({
        "actual_common_start_utc": format_utc(start, 6),
        "actual_common_end_utc": format_utc(end, 6),
        "actual_common_duration_seconds": (end - start).total_seconds(),
        "common_segments": [[format_utc(a, 6), format_utc(b, 6)] for a, b in common_segments],
        "three_station_coverage_fraction": coverage_fraction(common_segments, predicted_window)
        if common_segments else 0.0,
        "three_station_retrieved_window_coverage_fraction": coverage_fraction(common_segments, retrieved)
        if common_segments else 0.0,
        "median_nearest_timestamp_separation_s": float(np.median(joined)) if joined.size else None,
        "max_nearest_timestamp_separation_s": float(np.max(joined)) if joined.size else None,
        "station_pair_time_offset_statistics": pair_stats,
        "no_common_time_reason": None,
    })
    return result
