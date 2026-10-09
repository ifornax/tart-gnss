"""Ranking and serialisation of common-visibility results (CSV, JSON, summary)."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Union

from tart_gnss.propagation import COORDINATE_FRAME, COORDINATE_FRAME_DESCRIPTION
from tart_gnss.search import CommonWindow, SearchConfig
from tart_gnss.stations import Station
from tart_gnss.tle import TLECatalogue, TLERecord, format_utc

logger = logging.getLogger(__name__)

EPOCHS = ("start", "mid", "end")

#: Version of the Stage-1 results/metadata JSON layout. Version 1 (files
#: without this field) identified stations by catalogue names such as
#: ``mauritius``; version 2 identifies them by canonical TART id (``mu-udm``)
#: whenever the default registry catalogue is used, and records the station
#: registry snapshot in the run metadata.
RESULTS_SCHEMA_VERSION = 2


def rank_windows(windows: Iterable[CommonWindow]) -> List[CommonWindow]:
    """Sort windows by decreasing duration and assign 1-based ranks.

    Ties are broken by earlier start time, then by NORAD ID, so the ranking
    is deterministic.
    """
    ordered = sorted(
        windows, key=lambda w: (-round(w.duration_s, 6), w.start_utc, w.tle.norad_id)
    )
    for i, w in enumerate(ordered, start=1):
        w.rank = i
    return ordered


def csv_columns(n_stations: int) -> List[str]:
    """Return the CSV header for ``n_stations`` stations."""
    cols = [
        "rank", "satellite_name", "norad_id", "constellation",
        "common_start_utc", "common_end_utc", "duration_seconds", "duration_minutes",
        "start_clipped", "end_clipped", "min_elevation_deg",
    ]
    for i in range(1, n_stations + 1):
        cols.append(f"station_{i}")
        for epoch in EPOCHS:
            cols += [f"station_{i}_az_{epoch}_deg", f"station_{i}_el_{epoch}_deg"]
        cols.append(f"station_{i}_range_mid_km")
    cols += [
        "mid_utc", "coordinate_frame",
        "sat_x_mid_km", "sat_y_mid_km", "sat_z_mid_km",
        "tle_epoch", "tle_age_start_days", "tle_age_mid_days", "tle_source",
    ]
    return cols


def window_row(w: CommonWindow, min_elevation_deg: float) -> Dict[str, object]:
    """Flatten one window into a CSV row (angles in degrees, distances in km)."""
    row: Dict[str, object] = {
        "rank": w.rank,
        "satellite_name": w.tle.name,
        "norad_id": w.tle.norad_id,
        "constellation": w.tle.constellation,
        "common_start_utc": format_utc(w.start_utc),
        "common_end_utc": format_utc(w.end_utc),
        "duration_seconds": f"{w.duration_s:.3f}",
        "duration_minutes": f"{w.duration_s / 60.0:.3f}",
        "start_clipped": w.window.start_clipped,
        "end_clipped": w.window.end_clipped,
        "min_elevation_deg": min_elevation_deg,
    }
    for i, sid in enumerate(w.station_ids, start=1):
        row[f"station_{i}"] = sid
        for epoch in EPOCHS:
            g = w.geometry[sid][epoch]
            row[f"station_{i}_az_{epoch}_deg"] = f"{g.azimuth_deg:.4f}"
            row[f"station_{i}_el_{epoch}_deg"] = f"{g.elevation_deg:.4f}"
        row[f"station_{i}_range_mid_km"] = f"{w.geometry[sid]['mid'].range_km:.3f}"
    x, y, z = w.satellite_itrs["mid"]["position_km"]
    row.update(
        {
            "mid_utc": format_utc(w.mid_utc),
            "coordinate_frame": COORDINATE_FRAME,
            "sat_x_mid_km": f"{x:.4f}",
            "sat_y_mid_km": f"{y:.4f}",
            "sat_z_mid_km": f"{z:.4f}",
            "tle_epoch": format_utc(w.tle.epoch_utc, decimals=6),
            "tle_age_start_days": f"{w.tle_age_start_days:.5f}",
            "tle_age_mid_days": f"{w.tle_age_mid_days:.5f}",
            "tle_source": w.tle.source,
        }
    )
    return row


def write_csv(
    path: Union[str, Path], windows: Sequence[CommonWindow], n_stations: int, min_elevation_deg: float
) -> Path:
    """Write ranked windows to CSV (one row per satellite x window).

    Full TLE lines are not written here; they are in the JSON output.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_columns(n_stations))
        writer.writeheader()
        for w in windows:
            writer.writerow(window_row(w, min_elevation_deg))
    logger.info("Wrote %s", path)
    return path


def window_to_dict(w: CommonWindow) -> Dict[str, object]:
    """Full JSON representation of one window, including TLE provenance."""
    return {
        "rank": w.rank,
        "satellite": w.tle.name,
        "norad_id": w.tle.norad_id,
        "constellation": w.tle.constellation,
        "start_utc": format_utc(w.start_utc),
        "mid_utc": format_utc(w.mid_utc),
        "end_utc": format_utc(w.end_utc),
        "duration_seconds": w.duration_s,
        "start_clipped_by_search_interval": w.window.start_clipped,
        "end_clipped_by_search_interval": w.window.end_clipped,
        "stations": {
            sid: {epoch: w.geometry[sid][epoch].to_dict() for epoch in EPOCHS}
            for sid in w.station_ids
        },
        "satellite_position": {
            "coordinate_frame": COORDINATE_FRAME,
            "units": {"position": "km", "velocity": "km/s"},
            **w.satellite_itrs,
        },
        "line_of_sight_mid": {"coordinate_frame": COORDINATE_FRAME, **w.line_of_sight_mid},
        "tle": {
            "epoch": format_utc(w.tle.epoch_utc, decimals=6),
            "line1": w.tle.line1,
            "line2": w.tle.line2,
            "source": w.tle.source,
            "age_start_days": w.tle_age_start_days,
            "age_mid_days": w.tle_age_mid_days,
            "age_end_days": w.tle_age_end_days,
        },
    }


def build_results_json(
    windows: Sequence[CommonWindow],
    stations: Sequence[Station],
    station_itrs_km: Dict[str, Sequence[float]],
    config: SearchConfig,
    catalogue: TLECatalogue,
    searched: Sequence[TLERecord],
) -> Dict[str, object]:
    """Assemble the detailed machine-readable results document."""
    return {
        "schema_version": RESULTS_SCHEMA_VERSION,
        "run": {
            "start_utc": format_utc(config.start),
            "end_utc": format_utc(config.end),
            "min_elevation_deg": config.min_elevation_deg,
            "step_seconds": config.step_s,
            "boundary_tolerance_seconds": config.tolerance_s,
            "min_duration_seconds": config.min_duration_s,
            "visibility_definition": "all stations simultaneously have elevation >= min_elevation_deg",
        },
        "coordinate_frame": {"name": COORDINATE_FRAME, "description": COORDINATE_FRAME_DESCRIPTION},
        "stations": [
            {**s.to_dict(), "itrs_position_km": list(map(float, station_itrs_km[s.id]))}
            for s in stations
        ],
        "tle_catalogue": catalogue.to_dict(),
        "satellites_searched": [r.to_dict() for r in searched],
        "windows": [window_to_dict(w) for w in windows],
    }


def write_json(path: Union[str, Path], document: Dict[str, object]) -> Path:
    """Write a JSON document with stable formatting."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    logger.info("Wrote %s", path)
    return path


def format_summary(
    windows: Sequence[CommonWindow],
    stations: Sequence[Station],
    config: SearchConfig,
    max_rows: int = 20,
) -> str:
    """Human-readable terminal summary of the ranked windows."""
    lines = ["COMMON GNSS VISIBILITY", "", "Stations:"]
    lines += [f"  {s.name} [{s.id}]" for s in stations]
    lines += [
        "",
        "Interval:",
        f"  {format_utc(config.start, 0)}",
        f"  {format_utc(config.end, 0)}",
        "",
        f"Minimum elevation: {config.min_elevation_deg:g} deg",
        "",
    ]
    if not windows:
        lines.append("No common GNSS visibility windows found.")
        return "\n".join(lines)
    header = f"{'Rank':<5}{'Satellite':<28}{'NORAD':<7}{'Const.':<9}{'Start UTC':<22}{'End UTC':<22}{'Duration':>10}"
    lines += [header, "-" * len(header)]
    for w in windows[:max_rows]:
        flag = "*" if (w.window.start_clipped or w.window.end_clipped) else " "
        lines.append(
            f"{w.rank:<5}{w.tle.name[:27]:<28}{w.tle.norad_id:<7}{w.tle.constellation:<9}"
            f"{format_utc(w.start_utc, 0):<22}{format_utc(w.end_utc, 0):<22}"
            f"{w.duration_s / 60.0:>8.1f} min{flag}"
        )
    if len(windows) > max_rows:
        lines.append(f"... {len(windows) - max_rows} more windows in the CSV/JSON output")
    if any(w.window.start_clipped or w.window.end_clipped for w in windows):
        lines.append("* window truncated by the search interval")
    lines.append(f"Total: {len(windows)} common visibility windows")
    return "\n".join(lines)
