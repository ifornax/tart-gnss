"""Physical sanity checks of propagation and frame consistency."""

from datetime import datetime, timezone

import numpy as np
import pytest

from tart_gnss import propagation
from tart_gnss.stations import Station
from tart_gnss.tle import load_tle_file

START = datetime(2026, 10, 1, tzinfo=timezone.utc)
STATION = Station("test", "test", -20.0, 57.0, 100.0, "synthetic test value")


@pytest.fixture
def records(sample_tle_path):
    return {r.constellation: r for r in load_tle_file(sample_tle_path).records}


@pytest.mark.parametrize(
    "constellation,radius_km",
    [("GPS", 26560.0), ("Galileo", 29600.0), ("BeiDou", 42164.0), ("GLONASS", 25510.0)],
)
def test_orbit_radius(records, constellation, radius_km):
    sat = propagation.make_satellite(records[constellation])
    t = propagation.times_from_offsets(START, np.arange(0, 86400, 3600.0))
    xyz, _ = propagation.satellite_itrs(sat, t)
    r = np.linalg.norm(xyz, axis=0)
    assert np.all(np.abs(r - radius_km) < 0.02 * radius_km)


def test_altaz_consistent_with_itrs_line_of_sight(records):
    """Az/el from Skyfield must match an ENU projection of the ITRS line of sight."""
    sat = propagation.make_satellite(records["GPS"])
    t = propagation.times_from_offsets(START, np.arange(0, 86400, 600.0))
    pos = propagation.station_position(STATION)
    az, el, rng, _ = propagation.topocentric(sat, pos, t)
    xyz, _ = propagation.satellite_itrs(sat, t)
    d = xyz - propagation.station_itrs_km(STATION)[:, None]
    lat, lon = np.radians(STATION.latitude_deg), np.radians(STATION.longitude_deg)
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    e, n, u = east @ d, north @ d, up @ d
    el2 = np.degrees(np.arctan2(u, np.hypot(e, n)))
    az2 = np.degrees(np.arctan2(e, n)) % 360.0
    assert np.allclose(np.linalg.norm(d, axis=0), rng, rtol=1e-6)
    assert np.max(np.abs(el2 - el)) < 1e-3
    daz = (az2 - az + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(daz[el > 0])) < 1e-3


def test_line_of_sight_unit_vector():
    los = propagation.line_of_sight(np.array([3.0, 4.0, 0.0]), np.zeros(3))
    assert los["range_km"] == pytest.approx(5.0)
    assert los["unit_vector"] == pytest.approx([0.6, 0.8, 0.0])


def test_geostationary_elevation_nearly_constant(records):
    sat = propagation.make_satellite(records["BeiDou"])
    pos = propagation.station_position(STATION)
    t = propagation.times_from_offsets(START, np.arange(0, 86400, 1800.0))
    el = propagation.elevations(sat, [pos], t)
    assert el.shape == (1, t.shape[0])
    assert np.ptp(el) < 10.0
