"""Transparent, deterministic ranking of archive candidates.

The sort key, in order (earlier components dominate):

1. complete inventory for all required stations first;
2. all selected stations present (a partial-network diagnostic ranks last);
3. in-band signal relevance (``GPS_L1``, ``GALILEO_E1``, ``BEIDOU_B1C``) first;
   ``UNKNOWN`` is never preferred;
4. possible overlap duration **capped** at ``duration_cap_s`` (default 600 s),
   descending, so a very long (e.g. geostationary) window cannot win on
   length alone;
5. absolute TLE age at the overlap midpoint, ascending (missing last);
6. Stage-1 geometry score: the minimum over stations of the Stage-1
   mid-window elevation (``windows[].stations[].mid.elevation_deg``),
   descending (missing last);
7. Stage-1 rank, NORAD ID and overlap start, ascending (stable tie-breakers).
"""

from __future__ import annotations

import math
from typing import Dict, List, Sequence

from tart_gnss.archive_scan.matching import IN_BAND, Candidate

DEFAULT_DURATION_CAP_S = 600.0

RANKING_RULE = [
    "inventory_complete (True first)",
    "all_selected_stations (True first)",
    "in_band_signal (GPS_L1/GALILEO_E1/BEIDOU_B1C first)",
    "min(possible_overlap_duration_s, duration_cap_s) descending",
    "abs(tle_age_days) ascending (missing last)",
    "stage1_min_mid_elevation_deg descending (missing last)",
    "rank_stage1, norad_id, overlap_start_utc ascending",
]


def ranking_components(c: Candidate, duration_cap_s: float) -> Dict[str, object]:
    """The published values of every ranking component for one candidate."""
    return {
        "inventory_complete": c.segment.inventory_complete,
        "all_selected_stations": c.segment.all_stations,
        "in_band_signal": c.window.signal_relevance in IN_BAND,
        "capped_duration_s": min(c.duration_s, duration_cap_s),
        "abs_tle_age_days": abs(c.tle_age_days) if c.tle_age_days is not None else None,
        "stage1_min_mid_elevation_deg": c.window.min_mid_elevation_deg,
    }


def _key(c: Candidate):
    r = c.ranking_components
    age = r["abs_tle_age_days"]
    elev = r["stage1_min_mid_elevation_deg"]
    return (
        not r["inventory_complete"],
        not r["all_selected_stations"],
        not r["in_band_signal"],
        -r["capped_duration_s"],
        math.inf if age is None else age,
        math.inf if elev is None else -elev,
        c.window.rank, c.window.norad_id, c.overlap_start_utc,
    )


def rank_candidates(candidates: Sequence[Candidate], duration_cap_s: float = DEFAULT_DURATION_CAP_S) -> List[Candidate]:
    """Sort candidates by :data:`RANKING_RULE`; set ``order`` (1 = best) and ``candidate_id``.

    ``candidate_id`` is ``R<stage1 rank>-S<segment id>``, unique because a
    window meets each disjoint segment at most once.
    """
    if duration_cap_s <= 0:
        raise ValueError("duration cap must be positive")
    for c in candidates:
        c.ranking_components = ranking_components(c, duration_cap_s)
        c.candidate_id = f"R{c.window.rank}-S{c.segment.segment_id}"
    ranked = sorted(candidates, key=_key)
    for i, c in enumerate(ranked, start=1):
        c.order = i
    return ranked
