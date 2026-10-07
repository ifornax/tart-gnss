"""SGP4 orbit propagation and station-relative geometry (via Skyfield).

No orbital mechanics or frame transformations are implemented here by hand:
SGP4 propagation, TEME->ITRS conversion and topocentric angles are all
delegated to Skyfield/sgp4.

Conventions
-----------
* Times are Skyfield ``Time`` objects built from UTC; offsets are elapsed SI
  seconds from the run start (added on the TT scale, so leap seconds are
  handled correctly).
* Azimuth is measured from geographic north through east, in degrees [0, 360).
* Elevation is the geometric (unrefracted) altitude above the WGS84 local
  horizontal, in degrees.
* Earth-fixed vectors are given in Skyfield's ITRS frame, in km. Skyfield's
  ITRS ignores polar motion unless polar-motion tables are loaded (they are
  not, by default), an approximation at the ~10 m level.
"""

from __future__ import annotations

from datetime import datetime
from functools import lru_cache
from typing import Dict, Sequence, Tuple

import numpy as np
from skyfield.api import EarthSatellite, Timescale, load, wgs84
from skyfield.framelib import itrs
from skyfield.toposlib import GeographicPosition
from skyfield.timelib import Time

from tart_gnss.stations import Station
from tart_gnss.tle import TLERecord

#: Identifier of the Earth-fixed frame used for all Cartesian outputs.
COORDINATE_FRAME = "ITRS"
COORDINATE_FRAME_DESCRIPTION = (
    "ITRS (International Terrestrial Reference System) Earth-centred, "
    "Earth-fixed Cartesian coordinates as implemented by Skyfield "
    "(skyfield.framelib.itrs; TEME->ITRS via SGP4/Skyfield; polar motion "
    "neglected). Station positions from WGS84 geodetic coordinates."
)


@lru_cache(maxsize=1)
def get_timescale() -> Timescale:
    """Return a Skyfield timescale using Skyfield's built-in UT1/leap-second tables.

    The built-in tables avoid network downloads, which keeps runs reproducible
    for a given Skyfield version.
    """
    return load.timescale(builtin=True)


def make_satellite(record: TLERecord, ts: Timescale = None) -> EarthSatellite:
    """Instantiate a Skyfield SGP4 satellite from a TLE record."""
    ts = ts or get_timescale()
    return EarthSatellite(record.line1, record.line2, record.name, ts)


def station_position(station: Station) -> GeographicPosition:
    """Return the WGS84 geographic position of a station for use with Skyfield."""
    return wgs84.latlon(
        station.latitude_deg, station.longitude_deg, elevation_m=station.altitude_m
    )


def station_itrs_km(station: Station) -> np.ndarray:
    """Return the station's ITRS Cartesian position in km, shape (3,)."""
    return np.asarray(station_position(station).itrs_xyz.km, dtype=float)


def times_from_offsets(start: datetime, offsets_s, ts: Timescale = None) -> Time:
    """Build Skyfield times at ``start`` + ``offsets_s`` elapsed SI seconds.

    Parameters
    ----------
    start
        Timezone-aware UTC datetime of the run start.
    offsets_s
        Scalar or array of elapsed seconds since ``start``.
    """
    ts = ts or get_timescale()
    t0 = ts.from_datetime(start)
    return ts.tt_jd(t0.tt + np.asarray(offsets_s, dtype=float) / 86400.0)


def topocentric(
    satellite: EarthSatellite, station: GeographicPosition, times: Time
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute satellite azimuth, elevation, range and range rate for a station.

    Parameters
    ----------
    satellite
        Skyfield ``EarthSatellite``.
    station
        Station geographic position (:func:`station_position`).
    times
        Skyfield time (scalar or array).

    Returns
    -------
    azimuth_deg : numpy.ndarray
        Topocentric azimuth, degrees east of north.
    elevation_deg : numpy.ndarray
        Topocentric geometric elevation, degrees.
    range_km : numpy.ndarray
        Station-to-satellite distance, km.
    range_rate_km_s : numpy.ndarray
        Time derivative of the range, km/s (positive = receding).

    Notes
    -----
    Geometric positions (no light-time, aberration or refraction), which
    are appropriate for line-of-sight geometry at the ~arcsecond level.
    """
    topo = (satellite - station).at(times)
    el, az, distance, _, _, range_rate = topo.frame_latlon_and_rates(station)
    return (
        np.asarray(az.degrees, dtype=float),
        np.asarray(el.degrees, dtype=float),
        np.asarray(distance.km, dtype=float),
        np.asarray(range_rate.km_per_s, dtype=float),
    )


def elevations(
    satellite: EarthSatellite, stations: Sequence[GeographicPosition], times: Time
) -> np.ndarray:
    """Elevation (degrees) of a satellite at every station and time.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n_stations, n_times)`` (or ``(n_stations,)`` for a
        scalar time).
    """
    rows = []
    for station in stations:
        el, _, _ = (satellite - station).at(times).altaz()
        rows.append(np.asarray(el.degrees, dtype=float))
    return np.vstack(rows) if np.ndim(rows[0]) else np.asarray(rows)


def satellite_itrs(satellite: EarthSatellite, times: Time) -> Tuple[np.ndarray, np.ndarray]:
    """Satellite position (km) and velocity (km/s) in the ITRS frame.

    Returns
    -------
    position_km, velocity_km_s : numpy.ndarray
        Arrays of shape ``(3,)`` or ``(3, n_times)``. The velocity is the
        Earth-fixed (rotating-frame) velocity.
    """
    position, velocity = satellite.at(times).frame_xyz_and_velocity(itrs)
    return np.asarray(position.km, dtype=float), np.asarray(velocity.km_per_s, dtype=float)


def line_of_sight(sat_itrs_km: np.ndarray, station_itrs: np.ndarray) -> Dict[str, object]:
    """Geometric line of sight d = r_s - r_i from a station to the satellite.

    Parameters
    ----------
    sat_itrs_km
        Satellite ITRS position, km, shape (3,).
    station_itrs
        Station ITRS position, km, shape (3,).

    Returns
    -------
    dict
        ``vector_km`` (d), ``unit_vector`` (d/|d|) and ``range_km`` (|d|),
        all in the ITRS frame.
    """
    d = np.asarray(sat_itrs_km, dtype=float) - np.asarray(station_itrs, dtype=float)
    rng = float(np.linalg.norm(d))
    return {"vector_km": d.tolist(), "unit_vector": (d / rng).tolist(), "range_km": rng}
