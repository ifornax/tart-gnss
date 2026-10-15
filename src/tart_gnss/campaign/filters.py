"""Run filters, the campaign quality policy and exclusion reasons.

A run is either *included* (enters the campaign tables) or *excluded* with an
explicit reason; runs are never dropped silently. Reasons:

``SCHEMA_MISMATCH``          unsupported Stage-7 schema version
``MISSING_STAGE7_OUTPUT``    unreadable manifest or missing Stage-7 tables
``BROKEN_PROVENANCE``        an upstream link missing, changed or mismatched (``strict`` only)
``INSUFFICIENT_EPOCHS``      fewer than ``--min-epochs-per-run`` diagnosed epochs
``NO_PRECISION_VALIDATION``  ``--precision-validation-only`` and the run is not a precision validation
``FILTERED_<NAME>``          removed by a user filter (NORAD ID, constellation, station, time range,
                             elevation, TLE age, evidence level)
``DUPLICATE_RUN``            the same Stage-7 manifest content (SHA-256) was already included

Included runs then fall into one *cohort*: ``precision_validation`` (Stage-6
precision validation true) or ``descriptive_only``. Scientific summaries use
the precision cohort by default (``--precision-only``); descriptive runs are
still counted and reported.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Sequence

import numpy as np

POLICIES = ("strict", "permissive")


@dataclass
class RunFilters:
    """User filters (all optional) and thresholds."""

    norad_ids: Sequence[int] = field(default_factory=list)
    constellations: Sequence[str] = field(default_factory=list)
    stations: Sequence[str] = field(default_factory=list)
    start: Optional[str] = None
    end: Optional[str] = None
    min_elevation_deg: Optional[float] = None
    max_tle_age_days: Optional[float] = None
    max_evidence_level: Optional[int] = None
    precision_validation_only: bool = False
    min_epochs_per_run: int = 3

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable form."""
        return {k: (list(v) if isinstance(v, (list, tuple)) else v) for k, v in asdict(self).items()}


def _time(text: Optional[str]) -> Optional[float]:
    if not text:
        return None
    return datetime.fromisoformat(str(text).replace("Z", "+00:00").replace(" ", "T")).timestamp()


def exclusion_reason(row: Dict[str, object], filters: RunFilters, policy: str) -> str:
    """Why an aggregated run row is excluded ("" if it is included)."""
    if policy not in POLICIES:
        raise ValueError(f"unknown quality policy {policy!r}")
    flags = set(row.get("quality_flags") or [])
    if policy == "strict" and "BROKEN_PROVENANCE" in flags:
        return "BROKEN_PROVENANCE"
    if (row.get("n_valid_epochs") or 0) < filters.min_epochs_per_run:
        return "INSUFFICIENT_EPOCHS"
    if filters.precision_validation_only and not row.get("precision_validation"):
        return "NO_PRECISION_VALIDATION"
    if filters.norad_ids and row.get("norad_id") not in set(int(x) for x in filters.norad_ids):
        return "FILTERED_NORAD_ID"
    if filters.constellations and str(row.get("constellation")).lower() not in {c.lower() for c in filters.constellations}:
        return "FILTERED_CONSTELLATION"
    if filters.stations and not set(filters.stations) & set(row.get("stations") or []):
        return "FILTERED_STATION"
    t0, t1 = _time(filters.start), _time(filters.end)
    ts = _time(row.get("start_utc"))
    if ts is not None and ((t0 is not None and ts < t0) or (t1 is not None and ts > t1)):
        return "FILTERED_TIME_RANGE"
    if filters.min_elevation_deg is not None:
        el = row.get("median_elevation_deg")
        if el is None or not np.isfinite(el) or el < filters.min_elevation_deg:
            return "FILTERED_ELEVATION"
    if filters.max_tle_age_days is not None:
        age = row.get("max_abs_tle_age_days")
        if age is None or not np.isfinite(age) or age > filters.max_tle_age_days:
            return "FILTERED_TLE_AGE"
    if filters.max_evidence_level is not None:
        lv = row.get("overall_evidence_level")
        if lv is not None and lv > filters.max_evidence_level:
            return "FILTERED_EVIDENCE_LEVEL"
    return ""


def split_runs(run_rows: List[Dict[str, object]], filters: RunFilters, policy: str):
    """``(included_rows, excluded_rows)``; excluded rows carry ``reason``."""
    inc, exc = [], []
    seen = set()
    for r in run_rows:
        why = exclusion_reason(r, filters, policy)
        if not why and r.get("manifest_sha256") in seen:
            why = "DUPLICATE_RUN"
        if why:
            exc.append({**r, "reason": why})
        else:
            seen.add(r.get("manifest_sha256"))
            inc.append(r)
    return inc, exc
