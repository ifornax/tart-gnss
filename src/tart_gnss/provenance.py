"""Reproducibility metadata: software versions, hashes and run parameters."""

from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone
from importlib import metadata
from typing import Dict, Optional, Sequence

from tart_gnss import __version__
from tart_gnss.propagation import COORDINATE_FRAME, COORDINATE_FRAME_DESCRIPTION
from tart_gnss.stations import Station
from tart_gnss.tle import TLECatalogue, format_utc, sha256_file  # noqa: F401  (re-export)


def _version(dist: str) -> Optional[str]:
    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return None


def software_versions() -> Dict[str, Optional[str]]:
    """Return versions of this package, Python and key dependencies."""
    return {
        "tart_gnss": __version__,
        "python": sys.version.split()[0],
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "skyfield": _version("skyfield"),
        "sgp4": _version("sgp4"),
        "numpy": _version("numpy"),
        "pyyaml": _version("PyYAML"),
    }


def utc_now() -> datetime:
    """Current time as an aware UTC datetime."""
    return datetime.now(timezone.utc)


def build_run_metadata(
    *,
    run_timestamp: datetime,
    start: datetime,
    end: datetime,
    stations: Sequence[Station],
    station_file: str,
    min_elevation_deg: float,
    step_s: float,
    tolerance_s: float,
    min_duration_s: float,
    constellations: Optional[Sequence[str]],
    catalogue: TLECatalogue,
    n_satellites_searched: int,
    n_windows: int,
    tle_age_summary: Dict[str, object],
    command_line: Sequence[str],
    output_files: Dict[str, str],
) -> Dict[str, object]:
    """Assemble the run metadata dictionary written to ``run_metadata_*.json``.

    All times are UTC ISO-8601 strings with a ``Z`` suffix.
    """
    return {
        "software": software_versions(),
        "run_timestamp_utc": format_utc(run_timestamp),
        "command_line": list(command_line),
        "start_utc": format_utc(start),
        "end_utc": format_utc(end),
        "station_ids": [s.id for s in stations],
        "stations": [s.to_dict() for s in stations],
        "station_file": station_file,
        "station_file_sha256": sha256_file(station_file),
        "min_elevation_deg": min_elevation_deg,
        "elevation_definition": "geometric (unrefracted) topocentric elevation above the WGS84 local horizon",
        "sampling_interval_s": step_s,
        "boundary_tolerance_s": tolerance_s,
        "min_duration_s": min_duration_s,
        "constellations": list(constellations) if constellations else "all",
        "timescale": "Skyfield built-in leap-second and Delta-T tables (load.timescale(builtin=True))",
        "coordinate_frame": COORDINATE_FRAME,
        "coordinate_frame_description": COORDINATE_FRAME_DESCRIPTION,
        **catalogue.to_dict(),
        "n_satellites_searched": n_satellites_searched,
        "n_common_windows": n_windows,
        "tle_age_summary": tle_age_summary,
        "output_files": output_files,
    }
