"""TART station configuration: loading, validation and identifier resolution.

Station coordinates are never hard-coded in the software. They are read from a
YAML or JSON configuration file whose entries record the provenance of every
coordinate (see ``config/stations.yaml``).
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict, dataclass
from importlib import resources
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Union

import yaml

logger = logging.getLogger(__name__)

#: Plausible geodetic-height range for a ground station (metres). Values
#: outside this range almost certainly indicate a units or data-entry error.
MIN_ALTITUDE_M = -500.0
MAX_ALTITUDE_M = 9000.0


class StationConfigError(ValueError):
    """Raised when the station configuration is missing, malformed or invalid."""


@dataclass(frozen=True)
class Station:
    """A ground station (TART telescope) at a fixed geodetic position.

    Parameters
    ----------
    id
        Short configuration identifier used on the command line
        (e.g. ``"mauritius"``).
    name
        Human-readable station name.
    latitude_deg
        WGS84 geodetic latitude in degrees, north positive, in [-90, 90].
    longitude_deg
        WGS84 longitude in degrees, east positive, in [-180, 180].
    altitude_m
        Height above the WGS84 ellipsoid in metres (as published by the
        source; the vertical datum of TART-published heights is undocumented).
    source
        Free-text provenance of the coordinates.
    tart_id
        TART telescope identifier (e.g. ``"mu-udm"``), if known.
    api_url
        TART REST API base URL, if known. (The configuration key
        ``tart_api`` is accepted as a legacy alias.)
    archive_name
        Prefix of the telescope in the TART visibility archive (used by
        Stage 2). ``None`` means the archive location is unknown.

    Raises
    ------
    StationConfigError
        If any coordinate is missing, non-finite or out of range.
    """

    id: str
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    source: str
    tart_id: Optional[str] = None
    api_url: Optional[str] = None
    archive_name: Optional[str] = None

    def __post_init__(self) -> None:
        validate_coordinates(
            self.id, self.latitude_deg, self.longitude_deg, self.altitude_m
        )
        if not self.source or not str(self.source).strip():
            raise StationConfigError(
                f"Station '{self.id}': a 'source' recording the provenance of "
                "the coordinates is required."
            )

    def to_dict(self) -> Dict[str, object]:
        """Return the station as a plain JSON-serialisable dictionary."""
        return asdict(self)


def _as_float(station_id: str, field: str, value: object) -> float:
    if value is None:
        raise StationConfigError(
            f"Station '{station_id}': '{field}' is missing (null). Coordinates "
            "must come from an authoritative source and are never guessed."
        )
    if isinstance(value, bool):
        raise StationConfigError(f"Station '{station_id}': '{field}' must be a number.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise StationConfigError(
            f"Station '{station_id}': '{field}'={value!r} is not a number."
        ) from exc
    if not math.isfinite(number):
        raise StationConfigError(f"Station '{station_id}': '{field}' must be finite.")
    return number


def validate_coordinates(
    station_id: str, latitude_deg: object, longitude_deg: object, altitude_m: object
) -> None:
    """Validate geodetic station coordinates.

    Parameters
    ----------
    station_id
        Identifier used in error messages.
    latitude_deg, longitude_deg
        Geodetic latitude in [-90, 90] deg and longitude in [-180, 180] deg.
    altitude_m
        Height in metres, within [MIN_ALTITUDE_M, MAX_ALTITUDE_M].

    Raises
    ------
    StationConfigError
        If a value is missing, non-numeric, non-finite or out of range.
    """
    lat = _as_float(station_id, "latitude_deg", latitude_deg)
    lon = _as_float(station_id, "longitude_deg", longitude_deg)
    alt = _as_float(station_id, "altitude_m", altitude_m)
    if not -90.0 <= lat <= 90.0:
        raise StationConfigError(
            f"Station '{station_id}': latitude {lat} deg is outside [-90, 90]."
        )
    if not -180.0 <= lon <= 180.0:
        raise StationConfigError(
            f"Station '{station_id}': longitude {lon} deg is outside [-180, 180]."
        )
    if not MIN_ALTITUDE_M <= alt <= MAX_ALTITUDE_M:
        raise StationConfigError(
            f"Station '{station_id}': altitude {alt} m is outside "
            f"[{MIN_ALTITUDE_M}, {MAX_ALTITUDE_M}]."
        )


def default_station_file() -> Path:
    """Return the path of the station configuration bundled with the package."""
    return Path(str(resources.files("tart_gnss").joinpath("config/stations.yaml")))


def parse_station_config(data: Mapping[str, object]) -> Dict[str, Station]:
    """Build validated :class:`Station` objects from a parsed configuration.

    Parameters
    ----------
    data
        Mapping with a top-level ``stations`` key whose value maps station
        identifiers to dictionaries with ``name``, ``latitude_deg``,
        ``longitude_deg``, ``altitude_m`` and ``source`` (plus optional
        ``tart_id``, ``api_url`` and ``archive_name``).

    Returns
    -------
    dict
        Station identifier -> :class:`Station`.
    """
    if not isinstance(data, Mapping) or "stations" not in data:
        raise StationConfigError("Station configuration must contain a 'stations' mapping.")
    entries = data["stations"]
    if not isinstance(entries, Mapping) or not entries:
        raise StationConfigError("'stations' must be a non-empty mapping.")
    stations: Dict[str, Station] = {}
    for station_id, entry in entries.items():
        station_id = str(station_id)
        if not isinstance(entry, Mapping):
            raise StationConfigError(f"Station '{station_id}' must be a mapping.")
        stations[station_id] = Station(
            id=station_id,
            name=str(entry.get("name") or station_id),
            latitude_deg=_as_float(station_id, "latitude_deg", entry.get("latitude_deg")),
            longitude_deg=_as_float(station_id, "longitude_deg", entry.get("longitude_deg")),
            altitude_m=_as_float(station_id, "altitude_m", entry.get("altitude_m")),
            source=str(entry.get("source") or ""),
            tart_id=entry.get("tart_id"),
            api_url=entry.get("api_url") or entry.get("tart_api"),
            archive_name=entry.get("archive_name"),
        )
    return stations


def load_stations(path: Optional[Union[str, Path]] = None) -> Dict[str, Station]:
    """Load and validate a YAML or JSON station configuration file.

    Parameters
    ----------
    path
        Configuration file (``.yaml``/``.yml``/``.json``). ``None`` selects
        the configuration bundled with the package.

    Returns
    -------
    dict
        Station identifier -> validated :class:`Station`.
    """
    path = Path(path) if path is not None else default_station_file()
    if not path.is_file():
        raise StationConfigError(f"Station configuration file not found: {path}")
    text = path.read_text(encoding="utf-8")
    try:
        if path.suffix.lower() == ".json":
            data = json.loads(text)
        else:
            data = yaml.safe_load(text)
    except (json.JSONDecodeError, yaml.YAMLError) as exc:
        raise StationConfigError(f"Cannot parse station configuration {path}: {exc}") from exc
    stations = parse_station_config(data)
    logger.debug("Loaded %d stations from %s", len(stations), path)
    return stations


def resolve_stations(
    station_ids: Iterable[str], catalogue: Mapping[str, Station]
) -> List[Station]:
    """Resolve station identifiers against a loaded configuration.

    Identifiers are matched case-insensitively; hyphens and underscores are
    treated as equivalent. The TART identifier (``tart_id``) is also accepted.

    Parameters
    ----------
    station_ids
        Requested identifiers, in the order they should be used.
    catalogue
        Output of :func:`load_stations`.

    Returns
    -------
    list of Station
        Stations in the requested order.

    Raises
    ------
    StationConfigError
        For unknown or duplicated stations.
    """

    def norm(s: str) -> str:
        return s.strip().lower().replace("-", "_")

    lookup: Dict[str, Station] = {}
    for station in catalogue.values():
        lookup[norm(station.id)] = station
        if station.tart_id:
            lookup.setdefault(norm(station.tart_id), station)

    resolved: List[Station] = []
    for requested in station_ids:
        station = lookup.get(norm(requested))
        if station is None:
            known = ", ".join(sorted(catalogue))
            raise StationConfigError(
                f"Unknown station '{requested}'. Known stations: {known}"
            )
        if station in resolved:
            raise StationConfigError(f"Station '{requested}' was requested more than once.")
        resolved.append(station)
    return resolved
