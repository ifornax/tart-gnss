import json

import pytest

from tart_gnss.stations import (
    Station,
    StationConfigError,
    default_station_file,
    legacy_station_file,
    load_stations,
    parse_station_config,
    resolve_stations,
)


def make(**overrides):
    entry = dict(name="x", latitude_deg=0.0, longitude_deg=0.0, altitude_m=0.0, source="test")
    entry.update(overrides)
    return {"stations": {"s": entry}}


@pytest.mark.parametrize("lat", [100.0, -90.5, float("nan")])
def test_invalid_latitude_rejected(lat):
    with pytest.raises(StationConfigError):
        parse_station_config(make(latitude_deg=lat))


@pytest.mark.parametrize("lon", [300.0, -180.1, float("inf")])
def test_invalid_longitude_rejected(lon):
    with pytest.raises(StationConfigError):
        parse_station_config(make(longitude_deg=lon))


@pytest.mark.parametrize("alt", [20000.0, -1000.0, "high"])
def test_invalid_altitude_rejected(alt):
    with pytest.raises(StationConfigError):
        parse_station_config(make(altitude_m=alt))


@pytest.mark.parametrize("field", ["latitude_deg", "longitude_deg", "altitude_m"])
def test_null_coordinates_rejected(field):
    with pytest.raises(StationConfigError, match="missing"):
        parse_station_config(make(**{field: None}))


def test_missing_source_rejected():
    with pytest.raises(StationConfigError, match="source"):
        parse_station_config(make(source=""))


def test_direct_construction_validates():
    with pytest.raises(StationConfigError):
        Station("s", "x", 100.0, 0.0, 0.0, "test")


def test_load_yaml_and_resolve(station_yaml):
    catalogue = load_stations(station_yaml)
    assert set(catalogue) == {"alpha", "beta", "gamma"}
    resolved = resolve_stations(["gamma", "ALPHA", "beta"], catalogue)
    assert [s.id for s in resolved] == ["gamma", "alpha", "beta"]
    # TART id with hyphen also resolves
    assert resolve_stations(["zz-gamma"], catalogue)[0].id == "gamma"


def test_load_json(tmp_path):
    path = tmp_path / "stations.json"
    path.write_text(json.dumps(make(latitude_deg=10.0)))
    assert load_stations(path)["s"].latitude_deg == 10.0


def test_unknown_station_clear_error(station_yaml):
    with pytest.raises(StationConfigError, match="Unknown station 'atlantis'"):
        resolve_stations(["alpha", "atlantis"], load_stations(station_yaml))


def test_duplicate_station_rejected(station_yaml):
    with pytest.raises(StationConfigError, match="more than once"):
        resolve_stations(["alpha", "alpha"], load_stations(station_yaml))


def test_missing_file():
    with pytest.raises(StationConfigError, match="not found"):
        load_stations("/nonexistent/stations.yaml")


def test_bundled_configuration_is_valid_and_sourced():
    catalogue = load_stations(default_station_file())
    for sid in ("mu-udm", "mg-tana", "za-rhodes"):
        station = catalogue[sid]
        assert station.source and "elec.ac.nz" in station.source
        assert station.tart_id == sid


def test_legacy_configuration_is_valid_and_sourced():
    catalogue = load_stations(legacy_station_file())
    for sid in ("mauritius", "madagascar", "south_africa"):
        station = catalogue[sid]
        assert station.source and "elec.ac.nz" in station.source
        assert station.tart_id


def test_legacy_and_registry_coordinates_agree():
    legacy = {s.tart_id: s for s in load_stations(legacy_station_file()).values()}
    current = load_stations(default_station_file())
    for tart_id, old in legacy.items():
        new = current[tart_id]
        assert (old.latitude_deg, old.longitude_deg, old.altitude_m) == (
            new.latitude_deg, new.longitude_deg, new.altitude_m)
