"""The Stage-2 manifest, the single interface from Stage 2 to Stage 3.

``manifest.json`` records the selected Stage-1 window (unchanged), the
retrieval, the native files with their hashes, the validated observation
metadata, the Measurement Sets, the synchronisation diagnostics and the
software versions. ``three_station_valid`` is true only if all three stations
have real, validated visibility data with a non-zero common interval (and,
when conversion was requested, a validated MS each).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

from tart_gnss.provenance import software_versions

logger = logging.getLogger(__name__)

MANIFEST_VERSION = "1.0"

REQUIRED_TOP_LEVEL = ("stage", "manifest_version", "three_station_valid", "validity_reasons",
                      "satellite", "prediction", "stations", "synchronisation", "software", "provenance")
REQUIRED_SATELLITE = ("satellite_name", "norad_id", "constellation", "signal_relevance",
                      "tle_epoch", "tle_source")
REQUIRED_PREDICTION = ("predicted_common_start_utc", "predicted_common_end_utc",
                       "predicted_duration_seconds", "query_start_utc", "query_end_utc",
                       "padding_seconds", "stage1_result_file", "stage1_result_sha256")
REQUIRED_STATION = ("station_id", "tart_id", "station_name", "latitude_deg", "longitude_deg",
                    "altitude_m", "archive_query_status", "archive_query_start",
                    "archive_query_end", "raw_files", "observation", "measurement_set")
REQUIRED_SYNC = ("actual_common_start_utc", "actual_common_end_utc",
                 "actual_common_duration_seconds", "three_station_coverage_fraction",
                 "station_pair_time_offset_statistics", "three_station_valid")

#: Number of stations a localisation dataset needs.
N_STATIONS = 3


def validate_manifest(doc: Dict[str, object]) -> List[str]:
    """Return a list of problems with a manifest (an empty list means it is valid).

    Checks the required fields of every section, that there are exactly
    three stations, and that the top-level and synchronisation
    ``three_station_valid`` flags agree.
    """
    problems = [f"missing '{k}'" for k in REQUIRED_TOP_LEVEL if k not in doc]
    if problems:
        return problems
    for section, keys in (("satellite", REQUIRED_SATELLITE), ("prediction", REQUIRED_PREDICTION),
                          ("synchronisation", REQUIRED_SYNC)):
        problems += [f"missing '{section}.{k}'" for k in keys if k not in doc[section]]  # type: ignore[operator]
    stations = doc["stations"]
    if not isinstance(stations, dict) or len(stations) != N_STATIONS:
        problems.append(f"expected {N_STATIONS} stations, found {len(stations) if isinstance(stations, dict) else 'none'}")
    else:
        for sid, entry in stations.items():
            problems += [f"missing 'stations.{sid}.{k}'" for k in REQUIRED_STATION if k not in entry]
    if doc["three_station_valid"] != doc["synchronisation"].get("three_station_valid"):  # type: ignore[union-attr]
        problems.append("three_station_valid differs between top level and synchronisation")
    if doc["three_station_valid"] and doc["validity_reasons"]:
        problems.append("three_station_valid is true but validity_reasons is not empty")
    return problems


def build_manifest(
    *,
    three_station_valid: bool,
    validity_reasons: List[str],
    satellite: Dict[str, object],
    prediction: Dict[str, object],
    stations: Dict[str, Dict[str, object]],
    synchronisation: Dict[str, object],
    provenance: Dict[str, object],
    stage1_window: Dict[str, object],
) -> Dict[str, object]:
    """Assemble the manifest dictionary (no I/O)."""
    synchronisation = dict(synchronisation)
    synchronisation["three_station_valid"] = three_station_valid
    return {
        "stage": 2,
        "manifest_version": MANIFEST_VERSION,
        "three_station_valid": three_station_valid,
        "validity_reasons": validity_reasons,
        "satellite": satellite,
        "prediction": prediction,
        "stations": stations,
        "synchronisation": synchronisation,
        "stage1_window": stage1_window,
        "software": software_versions_stage2(),
        "provenance": provenance,
    }


def software_versions_stage2() -> Dict[str, Optional[str]]:
    """Stage-1 software versions plus the Stage-2 dependencies."""
    from importlib import metadata

    versions = software_versions()
    for dist in ("h5py", "minio", "tart2ms", "tart-tools", "python-casacore"):
        try:
            versions[dist] = metadata.version(dist)
        except metadata.PackageNotFoundError:
            versions[dist] = None
    return versions


def write_manifest(path, doc: Dict[str, object], run_tag: str) -> Path:
    """Write ``manifest.json``. An existing manifest is renamed, never overwritten.

    Parameters
    ----------
    path
        Destination (normally ``<run dir>/manifest.json``).
    doc
        Manifest from :func:`build_manifest`.
    run_tag
        Compact UTC timestamp used to rename a previous manifest
        (``manifest_<previous tag>.json``).

    Raises
    ------
    ValueError
        If the manifest fails :func:`validate_manifest`.
    """
    problems = validate_manifest(doc)
    if problems:
        raise ValueError("Refusing to write an invalid manifest: " + "; ".join(problems))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        try:
            old_tag = json.loads(path.read_text())["provenance"]["run_tag"]
        except (ValueError, KeyError, TypeError):
            old_tag = "previous"
        archived = path.with_name(f"manifest_{old_tag}.json")
        n = 1
        while archived.exists():
            archived = path.with_name(f"manifest_{old_tag}_{n}.json")
            n += 1
        path.rename(archived)
        logger.info("Previous manifest kept as %s", archived)
    path.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    logger.info("Wrote manifest %s", path)
    return path
