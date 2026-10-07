"""Loading Stage-1 results and selecting a single common-visibility window.

Stage-1 geometry is never recomputed. Everything predicted in Stage 1
(times, az/el, ITRS positions, TLE provenance) is carried forward verbatim.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from tart_gnss.stage2.errors import Stage1SelectionError
from tart_gnss.tle import sha256_file

logger = logging.getLogger(__name__)

#: Centre frequency of the TART receivers (GPS L1 / Galileo E1 / BeiDou B1C), Hz.
TART_CENTRE_FREQUENCY_HZ = 1575.42e6

#: Allowed ``signal_relevance`` values.
SIGNAL_RELEVANCE_VALUES = ("GPS_L1", "GALILEO_E1", "BEIDOU_B1C", "OUT_OF_BAND", "UNKNOWN")


def parse_utc(text: str) -> datetime:
    """Parse an ISO-8601 UTC string (with ``Z`` or an offset) into an aware datetime.

    Raises
    ------
    ValueError
        If the string has no timezone or cannot be parsed.
    """
    dt = datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    if dt.tzinfo is None:
        raise ValueError(f"Timestamp '{text}' has no timezone")
    return dt.astimezone(timezone.utc)


def signal_relevance(constellation: str, satellite_name: str) -> Tuple[str, str]:
    """Classify whether a satellite transmits inside the TART band at 1575.42 MHz.

    Only well-established, constellation-level facts are used:

    * GPS satellites all broadcast L1 C/A at 1575.42 MHz -> ``GPS_L1``.
    * Galileo satellites broadcast E1 at 1575.42 MHz -> ``GALILEO_E1``.
    * BeiDou-3 MEO/IGSO satellites broadcast B1C at 1575.42 MHz -> ``BEIDOU_B1C``.
    * BeiDou-2 satellites broadcast B1I at 1561.098 MHz -> ``OUT_OF_BAND``.
    * Anything else (BeiDou-3 GEO, experimental BeiDou-3S, GLONASS, unknown)
      -> ``UNKNOWN``, because in-band transmission is not established here.

    The health or operational status of an individual satellite is not checked.

    Parameters
    ----------
    constellation
        Stage-1 constellation label.
    satellite_name
        Satellite name as written in the TLE.

    Returns
    -------
    (str, str)
        The relevance code and a short justification.
    """
    name = satellite_name.upper()
    if constellation == "GPS":
        return "GPS_L1", "GPS L1 C/A at 1575.42 MHz"
    if constellation == "Galileo":
        return "GALILEO_E1", "Galileo E1 at 1575.42 MHz"
    if constellation == "BeiDou":
        if name.startswith("BEIDOU-2"):
            return "OUT_OF_BAND", "BeiDou-2 B1I at 1561.098 MHz is outside the TART band"
        if re.match(r"BEIDOU-3 (M|IGSO)", name):
            return "BEIDOU_B1C", "BeiDou-3 MEO/IGSO B1C at 1575.42 MHz"
        return "UNKNOWN", "BeiDou satellite without an established B1C (1575.42 MHz) signal"
    if constellation == "GLONASS":
        return "UNKNOWN", "GLONASS L1 FDMA (1598-1606 MHz) is outside the TART band; in-band signals not established"
    return "UNKNOWN", "constellation unknown"


@dataclass
class StationInfo:
    """A station exactly as recorded in the Stage-1 results."""

    id: str
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    tart_id: Optional[str]
    raw: Dict[str, Any] = field(repr=False, default_factory=dict)


@dataclass
class Stage1Window:
    """One Stage-1 satellite x common-visibility window, carried forward unchanged.

    Attributes
    ----------
    rank
        Stage-1 rank (1 = longest window).
    window_index
        1-based index of this window among the same satellite's windows, in time order.
    n_windows_for_satellite
        Number of Stage-1 windows found for this satellite.
    start_utc, end_utc
        Predicted common-visibility interval (aware UTC datetimes).
    stations
        Stage-1 stations, in Stage-1 order.
    window
        The complete Stage-1 JSON record for the window.
    results_path, results_sha256
        The Stage-1 results file and its SHA-256.
    run
        The Stage-1 ``run`` section (threshold, step, ...).
    """

    rank: int
    window_index: int
    n_windows_for_satellite: int
    satellite_name: str
    norad_id: int
    constellation: str
    start_utc: datetime
    end_utc: datetime
    stations: List[StationInfo]
    window: Dict[str, Any]
    results_path: Path
    results_sha256: str
    run: Dict[str, Any]

    @property
    def duration_s(self) -> float:
        """Predicted window duration, seconds."""
        return (self.end_utc - self.start_utc).total_seconds()

    @property
    def mid_utc(self) -> datetime:
        """Predicted window midpoint."""
        return self.start_utc + (self.end_utc - self.start_utc) / 2

    @property
    def clipped(self) -> bool:
        """True if Stage 1 flagged either boundary as clipped by its search interval."""
        return bool(self.window.get("start_clipped_by_search_interval")
                    or self.window.get("end_clipped_by_search_interval"))


def load_stage1_results(path: Union[str, Path]) -> Tuple[Dict[str, Any], Path]:
    """Load the detailed Stage-1 JSON results.

    A Stage-1 CSV path is also accepted: the JSON written alongside it (same
    stem) is used, because only the JSON holds the full provenance.

    Returns
    -------
    (dict, pathlib.Path)
        The parsed document and the JSON path actually read.

    Raises
    ------
    Stage1SelectionError
        If the file is missing, is not valid JSON, or lacks the required sections.
    """
    path = Path(path)
    if path.suffix.lower() == ".csv":
        json_path = path.with_suffix(".json")
        if not json_path.is_file():
            raise Stage1SelectionError(
                f"{path} is a Stage-1 CSV but the detailed JSON {json_path} was not found; "
                "pass the JSON file"
            )
        logger.info("Using detailed Stage-1 JSON %s for CSV %s", json_path, path)
        path = json_path
    if not path.is_file():
        raise Stage1SelectionError(f"Stage-1 results file not found: {path}")
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise Stage1SelectionError(f"Cannot parse Stage-1 results {path}: {exc}") from exc
    for key in ("run", "stations", "windows"):
        if key not in doc:
            raise Stage1SelectionError(f"{path} is not a Stage-1 results file (missing '{key}')")
    return doc, path


def _windows_by_satellite(doc: Dict[str, Any], norad_id: int) -> List[Dict[str, Any]]:
    windows = [w for w in doc["windows"] if int(w["norad_id"]) == int(norad_id)]
    return sorted(windows, key=lambda w: w["start_utc"])


def select_window(
    doc: Dict[str, Any],
    results_path: Union[str, Path],
    rank: Optional[int] = None,
    norad_id: Optional[int] = None,
    window_index: Optional[int] = None,
) -> Stage1Window:
    """Select exactly one Stage-1 window.

    Selection rules:

    * ``rank`` alone selects that rank.
    * ``norad_id`` alone is valid only if the satellite has exactly one window.
    * ``norad_id`` with ``window_index`` (1-based, time order) selects that window.
    * ``rank`` with ``norad_id`` must refer to the same window.

    Raises
    ------
    Stage1SelectionError
        If the selection is missing, out of range, inconsistent or ambiguous.
    """
    results_path = Path(results_path)
    windows = doc["windows"]
    if not windows:
        raise Stage1SelectionError("The Stage-1 results contain no common visibility windows")
    if rank is None and norad_id is None:
        raise Stage1SelectionError("Select a window with --rank or --norad-id")
    if window_index is not None and norad_id is None:
        raise Stage1SelectionError("--window-index requires --norad-id")

    chosen: Optional[Dict[str, Any]] = None
    if rank is not None:
        matches = [w for w in windows if int(w["rank"]) == int(rank)]
        if not matches:
            raise Stage1SelectionError(
                f"Rank {rank} not found; valid ranks are 1..{len(windows)}"
            )
        chosen = matches[0]
        if norad_id is not None and int(chosen["norad_id"]) != int(norad_id):
            raise Stage1SelectionError(
                f"Rank {rank} is NORAD {chosen['norad_id']}, not {norad_id}"
            )
    if norad_id is not None:
        sat_windows = _windows_by_satellite(doc, norad_id)
        if not sat_windows:
            raise Stage1SelectionError(f"NORAD {norad_id} has no Stage-1 common windows")
        if window_index is not None:
            if not 1 <= window_index <= len(sat_windows):
                raise Stage1SelectionError(
                    f"NORAD {norad_id} has {len(sat_windows)} window(s); "
                    f"--window-index {window_index} is out of range"
                )
            candidate = sat_windows[window_index - 1]
            if chosen is not None and candidate is not chosen:
                raise Stage1SelectionError("--rank and --window-index select different windows")
            chosen = candidate
        elif chosen is None:
            if len(sat_windows) > 1:
                listing = "; ".join(
                    f"window {i} (rank {w['rank']}): {w['start_utc']} - {w['end_utc']}"
                    for i, w in enumerate(sat_windows, start=1)
                )
                raise Stage1SelectionError(
                    f"NORAD {norad_id} has {len(sat_windows)} windows; add --window-index. {listing}"
                )
            chosen = sat_windows[0]

    assert chosen is not None
    sat_windows = _windows_by_satellite(doc, int(chosen["norad_id"]))
    index = next(i for i, w in enumerate(sat_windows, start=1) if w is chosen)
    stations = [
        StationInfo(
            id=s["id"], name=s.get("name", s["id"]), latitude_deg=float(s["latitude_deg"]),
            longitude_deg=float(s["longitude_deg"]), altitude_m=float(s["altitude_m"]),
            tart_id=s.get("tart_id"), raw=s,
        )
        for s in doc["stations"]
    ]
    return Stage1Window(
        rank=int(chosen["rank"]),
        window_index=index,
        n_windows_for_satellite=len(sat_windows),
        satellite_name=str(chosen["satellite"]).strip(),
        norad_id=int(chosen["norad_id"]),
        constellation=str(chosen["constellation"]),
        start_utc=parse_utc(chosen["start_utc"]),
        end_utc=parse_utc(chosen["end_utc"]),
        stations=stations,
        window=chosen,
        results_path=results_path,
        results_sha256=sha256_file(results_path),
        run=doc["run"],
    )


def query_interval(
    start_utc: datetime, end_utc: datetime, padding_s: float, max_duration_s: Optional[float] = None
) -> Tuple[datetime, datetime]:
    """Padded archive query interval ``[t_start - P, t_end + P]``.

    Parameters
    ----------
    start_utc, end_utc
        Predicted window (aware UTC).
    padding_s
        Padding ``P`` in seconds (must be >= 0).
    max_duration_s
        If given, the padded interval is shortened to at most this many seconds,
        centred on the predicted window midpoint. Use this for a first, short
        retrieval.

    Returns
    -------
    (datetime, datetime)
        Query start and end, UTC.
    """
    if padding_s < 0:
        raise ValueError("padding must be non-negative")
    q0 = start_utc - timedelta(seconds=padding_s)
    q1 = end_utc + timedelta(seconds=padding_s)
    if max_duration_s is not None:
        if max_duration_s <= 0:
            raise ValueError("max duration must be positive")
        if (q1 - q0).total_seconds() > max_duration_s:
            mid = start_utc + (end_utc - start_utc) / 2
            half = timedelta(seconds=max_duration_s / 2)
            q0, q1 = mid - half, mid + half
    return q0, q1
