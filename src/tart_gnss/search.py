"""Orchestration of the Stage-1 common-visibility search.

Combines propagation, visibility and window extraction for every satellite
and characterises each resulting window (geometry at start/mid/end, ITRS
satellite position, line-of-sight vectors, TLE age). No file I/O here.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, List, Optional, Sequence

import numpy as np

from tart_gnss import propagation
from tart_gnss.stations import Station
from tart_gnss.tle import TLERecord, format_utc
from tart_gnss.visibility import common_visibility, make_common_visibility_predicate
from tart_gnss.windows import DEFAULT_TOLERANCE_S, Window, extract_windows

logger = logging.getLogger(__name__)

SECONDS_PER_DAY = 86400.0


@dataclass
class SearchConfig:
    """Parameters of a common-visibility search.

    Attributes
    ----------
    start, end
        Aware UTC datetimes bounding the search.
    min_elevation_deg
        Elevation threshold, degrees.
    step_s
        Coarse sampling interval, seconds.
    tolerance_s
        Boundary refinement tolerance, seconds.
    min_duration_s
        Discard windows shorter than this (seconds).
    """

    start: datetime
    end: datetime
    min_elevation_deg: float = 20.0
    step_s: float = 10.0
    tolerance_s: float = DEFAULT_TOLERANCE_S
    min_duration_s: float = 0.0

    def __post_init__(self) -> None:
        for label, value in (("start", self.start), ("end", self.end)):
            if value.tzinfo is None:
                raise ValueError(f"{label} must be timezone-aware UTC")
        self.start = self.start.astimezone(timezone.utc)
        self.end = self.end.astimezone(timezone.utc)
        if self.end <= self.start:
            raise ValueError("end must be later than start")
        if self.step_s <= 0:
            raise ValueError("step must be positive")
        if self.tolerance_s <= 0:
            raise ValueError("tolerance must be positive")
        if not -90.0 <= self.min_elevation_deg <= 90.0:
            raise ValueError("min elevation must be within [-90, 90] degrees")
        if self.min_duration_s < 0:
            raise ValueError("min duration must be non-negative")

    @property
    def span_s(self) -> float:
        """Search interval length in seconds."""
        return (self.end - self.start).total_seconds()

    def offsets(self) -> np.ndarray:
        """Coarse sample offsets (s), always including both interval ends."""
        grid = np.arange(0.0, self.span_s, self.step_s)
        return np.append(grid, self.span_s)

    def to_utc(self, offset_s: float) -> datetime:
        """Convert an elapsed-SI-second offset to a UTC datetime.

        Uses Skyfield so that the conversion is consistent with the
        propagation epochs (offsets are added on the TT scale).
        """
        t = propagation.times_from_offsets(self.start, float(offset_s))
        return t.utc_datetime()


@dataclass
class StationGeometry:
    """Predicted geometry of a satellite from one station at one epoch."""

    azimuth_deg: float
    elevation_deg: float
    range_km: float
    range_rate_km_s: float

    def to_dict(self) -> Dict[str, float]:
        """Return a JSON-serialisable dict."""
        return {
            "azimuth_deg": self.azimuth_deg,
            "elevation_deg": self.elevation_deg,
            "range_km": self.range_km,
            "range_rate_km_s": self.range_rate_km_s,
        }


@dataclass
class CommonWindow:
    """One satellite x common-visibility-window result.

    ``geometry[station_id][epoch]`` holds :class:`StationGeometry` for
    epochs ``"start"``, ``"mid"`` and ``"end"``. ``satellite_itrs`` holds the
    satellite ITRS position (km) and velocity (km/s) at the same epochs.
    ``line_of_sight_mid`` maps station id to the ITRS line-of-sight vector at
    the window midpoint. TLE ages are signed (observation - TLE epoch), days.
    """

    tle: TLERecord
    window: Window
    start_utc: datetime
    end_utc: datetime
    mid_utc: datetime
    station_ids: List[str]
    geometry: Dict[str, Dict[str, StationGeometry]]
    satellite_itrs: Dict[str, Dict[str, List[float]]]
    line_of_sight_mid: Dict[str, Dict[str, object]]
    tle_age_start_days: float
    tle_age_mid_days: float
    tle_age_end_days: float
    rank: Optional[int] = None

    @property
    def duration_s(self) -> float:
        """Window duration in seconds."""
        return self.window.duration_s


@dataclass
class SearchResult:
    """All windows found by :func:`find_common_windows` plus diagnostics."""

    windows: List[CommonWindow]
    n_satellites: int
    failed_satellites: List[Dict[str, object]] = field(default_factory=list)


def tle_age_days(observation: datetime, epoch: datetime) -> float:
    """Signed TLE age (observation - epoch) in days; negative = before epoch."""
    return (observation - epoch).total_seconds() / SECONDS_PER_DAY


def _characterise(
    record: TLERecord,
    satellite,
    stations: Sequence[Station],
    positions,
    station_itrs: Dict[str, np.ndarray],
    window: Window,
    config: SearchConfig,
) -> CommonWindow:
    epochs = {"start": window.start_s, "mid": window.mid_s, "end": window.end_s}
    ts = propagation.get_timescale()
    times = propagation.times_from_offsets(config.start, list(epochs.values()), ts)

    geometry: Dict[str, Dict[str, StationGeometry]] = {}
    for station, pos in zip(stations, positions):
        az, el, rng, rr = propagation.topocentric(satellite, pos, times)
        geometry[station.id] = {
            key: StationGeometry(float(az[k]), float(el[k]), float(rng[k]), float(rr[k]))
            for k, key in enumerate(epochs)
        }

    xyz, vel = propagation.satellite_itrs(satellite, times)
    satellite_itrs = {
        key: {"position_km": xyz[:, k].tolist(), "velocity_km_s": vel[:, k].tolist()}
        for k, key in enumerate(epochs)
    }
    los = {
        station.id: propagation.line_of_sight(xyz[:, 1], station_itrs[station.id])
        for station in stations
    }
    start_utc = config.to_utc(window.start_s)
    mid_utc = config.to_utc(window.mid_s)
    end_utc = config.to_utc(window.end_s)
    return CommonWindow(
        tle=record,
        window=window,
        start_utc=start_utc,
        end_utc=end_utc,
        mid_utc=mid_utc,
        station_ids=[s.id for s in stations],
        geometry=geometry,
        satellite_itrs=satellite_itrs,
        line_of_sight_mid=los,
        tle_age_start_days=tle_age_days(start_utc, record.epoch_utc),
        tle_age_mid_days=tle_age_days(mid_utc, record.epoch_utc),
        tle_age_end_days=tle_age_days(end_utc, record.epoch_utc),
    )


def find_common_windows(
    records: Sequence[TLERecord], stations: Sequence[Station], config: SearchConfig
) -> SearchResult:
    """Find every continuous window of simultaneous visibility.

    Parameters
    ----------
    records
        Satellites to search.
    stations
        Any number (>= 1) of stations; all must see the satellite.
    config
        Search parameters.

    Returns
    -------
    SearchResult
        Unranked windows (one per satellite x window) and diagnostics.
    """
    if not stations:
        raise ValueError("At least one station is required")
    ts = propagation.get_timescale()
    offsets = config.offsets()
    times = propagation.times_from_offsets(config.start, offsets, ts)
    positions = [propagation.station_position(s) for s in stations]
    station_itrs = {s.id: propagation.station_itrs_km(s) for s in stations}

    logger.info(
        "Searching %d GNSS satellites between %s and %s (%d samples, step %.3g s)",
        len(records), format_utc(config.start), format_utc(config.end),
        offsets.size, config.step_s,
    )
    results: List[CommonWindow] = []
    failed: List[Dict[str, object]] = []
    for record in records:
        satellite = propagation.make_satellite(record, ts)
        el = propagation.elevations(satellite, positions, times)
        if np.isnan(el).any():
            n_bad = int(np.isnan(el).any(axis=0).sum())
            logger.warning(
                "SGP4 propagation failed for %s (NORAD %d) at %d epochs; "
                "those epochs are treated as not visible",
                record.name, record.norad_id, n_bad,
            )
            failed.append({"norad_id": record.norad_id, "satellite_name": record.name,
                           "n_failed_epochs": n_bad})
        mask = common_visibility(el, config.min_elevation_deg)
        if not mask.any():
            continue
        predicate = make_common_visibility_predicate(
            satellite, positions, config.start, config.min_elevation_deg
        )
        for window in extract_windows(offsets, mask, predicate, config.tolerance_s):
            if window.duration_s < config.min_duration_s:
                continue
            results.append(
                _characterise(record, satellite, stations, positions, station_itrs, window, config)
            )
        logger.debug("%s: %d common samples", record.name, int(mask.sum()))
    logger.info("Found %d common visibility windows", len(results))
    return SearchResult(windows=results, n_satellites=len(records), failed_satellites=failed)


def assess_tle_ages(
    records: Sequence[TLERecord], config: SearchConfig, warn_days: float
) -> Dict[str, object]:
    """Diagnose separation between TLE epochs and the observation interval.

    Logs a warning when any satellite is propagated backwards (observation
    before TLE epoch -- e.g. newer TLEs used for an older observation) or
    when ``|t_obs - t_epoch|`` exceeds ``warn_days`` anywhere in the interval.

    ``warn_days`` is a user-configurable flag level, *not* an accuracy bound:
    SGP4 error growth depends on the orbit and on the TLE quality.

    Returns
    -------
    dict
        Summary statistics stored in the run metadata.
    """
    if not records:
        return {}
    start_age = np.array([tle_age_days(config.start, r.epoch_utc) for r in records])
    end_age = np.array([tle_age_days(config.end, r.epoch_utc) for r in records])
    max_abs = np.maximum(np.abs(start_age), np.abs(end_age))
    backwards = [r for r, a in zip(records, start_age) if a < 0]
    stale = [r for r, a in zip(records, max_abs) if a > warn_days]
    epochs = [r.epoch_utc for r in records]
    summary = {
        "tle_epoch_min_utc": format_utc(min(epochs)),
        "tle_epoch_max_utc": format_utc(max(epochs)),
        "max_abs_tle_age_days": float(max_abs.max()),
        "median_abs_tle_age_days": float(np.median(max_abs)),
        "n_backward_propagated": len(backwards),
        "n_exceeding_warning_threshold": len(stale),
        "warning_threshold_days": warn_days,
    }
    logger.info(
        "TLE epochs span %s .. %s; max |t_obs - t_epoch| = %.2f d",
        summary["tle_epoch_min_utc"], summary["tle_epoch_max_utc"], summary["max_abs_tle_age_days"],
    )
    if backwards:
        logger.warning(
            "%d satellites have TLE epochs AFTER the observation start: they are "
            "propagated backwards in time. For historical observations use a TLE "
            "catalogue contemporaneous with the observation (--tle-file).",
            len(backwards),
        )
    if stale:
        logger.warning(
            "%d satellites have |t_obs - t_TLE epoch| > %.1f days (max %.1f d). SGP4 "
            "accuracy degrades with propagation span; consider TLEs closer to the "
            "observation.",
            len(stale), warn_days, summary["max_abs_tle_age_days"],
        )
    return summary
