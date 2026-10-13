"""Time matching of independent station direction estimates, and timing diagnostics.

Station timestamps are independent (each TART station dumps on its own clock
and Stage 3 bins each station separately), so rows are never joined by row
number and nothing is interpolated.

Matching
    One *reference station* supplies the reference epochs: by default the
    station with the most valid estimates (ties broken by station id). For
    every reference epoch ``t`` each other station contributes its nearest
    valid estimate ``t_i`` if ``|t_i - t| <= tolerance``. An estimate is used
    in at most one epoch: if two reference epochs pick the same estimate, the
    closer epoch keeps it.

Diagnostics
    The angular rate of the target at each station is estimated from that
    station's own Stage-3 track (a straight-line fit of the measured unit
    vectors over a sliding window; never from a TLE). The angular motion during a station's time
    offset is ``omega_i |t_i - t|``. ``TIME_MISMATCH`` is raised when this
    motion exceeds ``time_mismatch_sigma`` times the station's 1-sigma angular
    uncertainty, i.e. when the timing error is comparable with the
    measurement error. The light-travel time ``rho_i / c`` (about 70 ms for
    GNSS) is recorded but not corrected: Stage 4 does no retarded-time
    modelling.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from tart_gnss.stage4.errors import TimeMatchError


@dataclass
class MatchedDirection:
    """One station's contribution to a matched epoch."""

    station_id: str
    time_utc: object
    time_offset_s: float
    station_position_ecef_m: np.ndarray
    sight_vector_ecef: np.ndarray
    angular_covariance: np.ndarray
    quality_flags: List[str]
    row: Dict[str, object] = field(default_factory=dict)
    angular_rate_deg_s: float = float("nan")


@dataclass
class MatchedEpoch:
    """A reference epoch and the row index of each matched station's estimate."""

    reference_time_s: float
    reference_station: str
    rows: Dict[str, int]
    offsets_s: Dict[str, float]

    @property
    def station_ids(self) -> List[str]:
        """Matched station ids (sorted)."""
        return sorted(self.rows)

    @property
    def max_offset_s(self) -> float:
        """Spread of the matched station timestamps (max minus min, seconds)."""
        t = [self.reference_time_s + o for o in self.offsets_s.values()]
        return float(max(t) - min(t)) if t else float("nan")


def time_seconds(times) -> np.ndarray:
    """Seconds since the Unix epoch for an array of tz-aware/naive UTC datetimes."""
    import pandas as pd

    t = pd.to_datetime(pd.Series(times), utc=True)
    return (t - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().to_numpy(dtype=float)


def choose_reference_station(direction_tables: Dict[str, object]) -> str:
    """The station with the most rows (ties: alphabetical first)."""
    if not direction_tables:
        raise TimeMatchError("no station has valid direction estimates")
    return sorted(direction_tables, key=lambda k: (-len(direction_tables[k]), k))[0]


def match_station_epochs(
    direction_tables: Dict[str, object],
    tolerance_s: float,
    reference_station: Optional[str] = None,
    reference_time_s: Optional[float] = None,
) -> List[MatchedEpoch]:
    """Match independent station direction estimates in time.

    Parameters
    ----------
    direction_tables
        ``station_id -> DataFrame`` of the station's valid estimates, with a
        ``time_utc`` column. Row labels are returned in :class:`MatchedEpoch`.
    tolerance_s
        Maximum ``|t_i - t|`` for a station to contribute to epoch ``t``.
    reference_station
        Station supplying the reference epochs (default: most estimates).
    reference_time_s
        If given (Unix seconds), only the reference epoch nearest this time
        is returned.

    Returns
    -------
    list of MatchedEpoch
        In time order. An epoch may contain a single station (the reference)
        when nothing else matches; callers decide what is sufficient.
    """
    if not tolerance_s >= 0:
        raise TimeMatchError("time-match tolerance must be >= 0")
    tables = {k: v for k, v in direction_tables.items() if len(v)}
    if reference_station is None:
        reference_station = choose_reference_station(tables)
    if reference_station not in tables:
        raise TimeMatchError(f"reference station {reference_station!r} has no valid estimates")
    ref = tables[reference_station]
    ref_t = time_seconds(ref["time_utc"])
    ref_idx = list(ref.index)
    order = np.argsort(ref_t, kind="stable")
    if reference_time_s is not None:
        order = order[[int(np.argmin(np.abs(ref_t[order] - reference_time_s)))]]
    epochs = [MatchedEpoch(float(ref_t[k]), reference_station, {reference_station: ref_idx[k]},
                           {reference_station: 0.0}) for k in order]
    for sid, df in sorted(tables.items()):
        if sid == reference_station:
            continue
        t = time_seconds(df["time_utc"])
        idx = list(df.index)
        srt = np.argsort(t, kind="stable")
        ts = t[srt]
        best: Dict[int, tuple] = {}
        for e_i, ep in enumerate(epochs):
            j = int(np.searchsorted(ts, ep.reference_time_s))
            cands = [c for c in (j - 1, j) if 0 <= c < len(ts)]
            if not cands:
                continue
            c = min(cands, key=lambda c: abs(ts[c] - ep.reference_time_s))
            off = float(ts[c] - ep.reference_time_s)
            if abs(off) > tolerance_s + 1e-9:
                continue
            k = int(srt[c])
            if k not in best or abs(off) < abs(best[k][1]):
                best[k] = (e_i, off)
        for k, (e_i, off) in best.items():
            epochs[e_i].rows[sid] = idx[k]
            epochs[e_i].offsets_s[sid] = off
    return epochs


def angular_rates_deg_s(times_s: np.ndarray, unit_vectors: np.ndarray, window_s: float = 120.0,
                        min_samples: int = 3) -> np.ndarray:
    """Angular rate (deg/s) of a measured track at each sample.

    The unit vectors within ``+-window_s/2`` of each sample (widened to the
    ``min_samples`` nearest if needed) are fitted with a straight line in time;
    the rate is the norm of the slope component perpendicular to the
    direction. Fitting over a window rather than differencing neighbours keeps
    the angular measurement noise from dominating the rate. NaN when fewer
    than ``min_samples`` finite samples exist.
    """
    t = np.asarray(times_s, dtype=float)
    s = np.asarray(unit_vectors, dtype=float)
    out = np.full(len(t), np.nan)
    ok = np.isfinite(t) & np.all(np.isfinite(s), axis=1)
    if ok.sum() < min_samples:
        return out
    tk, sk = t[ok], s[ok]
    for i in np.flatnonzero(ok):
        dt = tk - t[i]
        sel = np.abs(dt) <= 0.5 * window_s
        if sel.sum() < min_samples:
            sel = np.zeros(len(tk), dtype=bool)
            sel[np.argsort(np.abs(dt))[:min_samples]] = True
        x = dt[sel]
        if np.ptp(x) <= 0:
            continue
        xc = x - x.mean()
        slope = (xc[:, None] * sk[sel]).sum(axis=0) / np.sum(xc * xc)
        slope -= np.dot(slope, s[i]) * s[i]
        out[i] = float(np.degrees(np.linalg.norm(slope)))
    return out


def motion_during_offset_deg(angular_rate_deg_s: float, offset_s: float) -> float:
    """Angular motion (deg) of the target during a station's time offset."""
    return float(abs(angular_rate_deg_s) * abs(offset_s))
