"""TART registry tests. All offline: payloads are local fixtures, HTTP is mocked.

Fixture coordinates are synthetic test values in the shape of the real
registry (https://map.elec.ac.nz/api/v1/telescopes), not authoritative data.
"""

import copy
import json
from dataclasses import replace
from datetime import datetime, timezone

import pytest

from tart_gnss import registry, registry_cli
from tart_gnss.registry import (
    RegistryDownload,
    RegistryFetchError,
    RegistrySchemaError,
    RegistryStationError,
    diff_catalogues,
    parse_registry,
    registry_to_stations,
)
from tart_gnss.stations import (
    StationConfigError,
    default_alias_file,
    default_station_file,
    load_aliases,
    load_catalogue_metadata,
    load_stations,
    resolve_stations,
)

RETRIEVED = "2026-10-09T12:00:00Z"


def entry(node, lat=-20.0, lon=57.0, alt=20.0, name=None):
    location = {}
    for key, value in (("lat", lat), ("lon", lon), ("alt", alt)):
        if value is not ...:
            location[key] = value
    return {"nodeName": node, "telescopeName": name or f"Site {node}", "online": True,
            "telescopeInfo": {"info": {"location": location, "name": node}}}


def payload(*entries):
    return {"telescopes": list(entries), "totalCount": len(entries), "onlineCount": 0,
            "timestamp": "2026-10-09T12:00:00Z"}


VALID = payload(entry("mu-udm", -20.25, 57.75, 20.0), entry("mg-tana", -18.9, 47.55, 1280.0),
                entry("za-rhodes", -33.3, 26.5, 642.26))


def download(doc, when=datetime(2026, 10, 9, 12, 0, tzinfo=timezone.utc)):
    raw = json.dumps(doc).encode()
    return RegistryDownload(url=registry.DEFAULT_REGISTRY_URL, retrieved_utc=when, raw=raw, payload=doc)


# --- parsing -------------------------------------------------------------


def test_valid_registry_normalises_all_entries():
    parsed = parse_registry(VALID, RETRIEVED)
    assert parsed.entry_count == 3 and not parsed.diagnostics
    assert list(parsed.stations) == ["mu-udm", "mg-tana", "za-rhodes"]
    rs = parsed.stations["za-rhodes"]
    assert (rs.latitude_deg, rs.longitude_deg, rs.altitude_m) == (-33.3, 26.5, 642.26)
    assert rs.name == "Site za-rhodes" and rs.registry_retrieved_utc == RETRIEVED


@pytest.mark.parametrize("field", ["lat", "lon", "alt"])
def test_null_coordinate_is_never_zero(field):
    doc = payload(entry("aa-null", **{field: None}), entry("bb-ok"))
    parsed = parse_registry(doc, RETRIEVED)
    assert "aa-null" not in parsed.stations and "bb-ok" in parsed.stations
    assert parsed.diagnostics[0].tart_id == "aa-null" and "missing" in parsed.diagnostics[0].reason


def test_absent_altitude_is_skipped():
    # bw-biust in the real registry publishes lat/lon but no alt
    parsed = parse_registry(payload(entry("aa-noalt", alt=...)), RETRIEVED)
    assert not parsed.stations and parsed.diagnostics[0].tart_id == "aa-noalt"


def test_invalid_longitude_rejected_by_existing_validation():
    parsed = parse_registry(payload(entry("aa-lon", lon=200.0)), RETRIEVED)
    assert not parsed.stations
    assert "outside [-180, 180]" in parsed.diagnostics[0].reason


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), "north"])
def test_non_finite_or_non_numeric_rejected(bad):
    assert not parse_registry(payload(entry("aa-bad", lat=bad)), RETRIEVED).stations


def test_missing_location_and_id_diagnosed():
    no_loc = {"nodeName": "aa-noloc", "telescopeInfo": None}
    parsed = parse_registry(payload(no_loc, {"telescopeName": "x"}, "junk", entry("bb-ok")), RETRIEVED)
    assert list(parsed.stations) == ["bb-ok"]
    assert len(parsed.diagnostics) == 3


def test_strict_mode_raises():
    with pytest.raises(RegistryStationError, match="aa-null"):
        parse_registry(payload(entry("aa-null", lat=None)), RETRIEVED, strict=True)


def test_duplicate_tart_id_diagnosed_and_excluded():
    doc = payload(entry("aa-dup", -10.0), entry("aa-dup", -11.0), entry("bb-ok"))
    parsed = parse_registry(doc, RETRIEVED)
    assert list(parsed.stations) == ["bb-ok"]
    assert [d.tart_id for d in parsed.diagnostics] == ["aa-dup", "aa-dup"]
    assert all("duplicate" in d.reason for d in parsed.diagnostics)
    with pytest.raises(RegistryStationError, match="duplicate"):
        parse_registry(doc, RETRIEVED, strict=True)


@pytest.mark.parametrize("doc", [[], {"nodes": []}, {"telescopes": {}}, "text", None])
def test_schema_change_raises(doc):
    with pytest.raises(RegistrySchemaError):
        parse_registry(doc, RETRIEVED)


# --- conversion ----------------------------------------------------------


def test_registry_to_stations_keeps_identity_and_provenance():
    parsed = parse_registry(VALID, RETRIEVED)
    stations = registry_to_stations(parsed.stations, snapshot_sha256="abc")
    st = stations["mu-udm"]
    assert st.id == st.tart_id == "mu-udm"
    assert st.api_url == registry.default_api_url("mu-udm") == "https://api.elec.ac.nz/tart/mu-udm"
    assert st.archive_name is None and st.archive_verified is False
    assert "map.elec.ac.nz" in st.source and RETRIEVED in st.source and "abc" in st.source


def test_archive_status_only_when_verified():
    parsed = parse_registry(VALID, RETRIEVED)
    seen = []

    def exists(prefix):
        seen.append(prefix)
        return prefix.startswith("za-rhodes/")

    status = registry.verify_archives(list(parsed.stations), exists,
                                      now=datetime(2026, 10, 9, tzinfo=timezone.utc))
    assert seen == ["mu-udm/vis/", "mg-tana/vis/", "za-rhodes/vis/"]
    stations = registry_to_stations(parsed.stations, status)
    assert stations["za-rhodes"].archive_name == "za-rhodes" and stations["za-rhodes"].archive_verified
    assert stations["za-rhodes"].archive_verified_utc == "2026-10-09T00:00:00Z"
    assert stations["mu-udm"].archive_name is None and stations["mu-udm"].archive_verified is False


def test_archive_check_failure_is_unverified():
    def broken(prefix):
        raise RegistryFetchError("offline")

    status = registry.verify_archives(["mu-udm"], broken)
    assert status["mu-udm"] == registry.ArchiveStatus(None, False, None)


def test_catalogue_round_trip(tmp_path):
    parsed = parse_registry(payload(*VALID["telescopes"], entry("bw-bad", alt=None)), RETRIEVED)
    stations = registry_to_stations(parsed.stations)
    doc = registry.build_catalogue_document(parsed, stations, source_url=registry.DEFAULT_REGISTRY_URL,
                                            retrieved_utc=RETRIEVED, sha256="f" * 64, snapshot_file="s.json")
    path = registry.write_catalogue(tmp_path / "cat.yaml", doc)
    assert load_stations(path) == stations
    meta = load_catalogue_metadata(path)
    assert meta["sha256"] == "f" * 64 and meta["retrieved_utc"] == RETRIEVED
    assert meta["skipped_entries"][0]["tart_id"] == "bw-bad"
    assert path.read_text().startswith("# GENERATED FILE")


# --- aliases and resolution ---------------------------------------------


@pytest.fixture
def catalogue():
    return registry_to_stations(parse_registry(VALID, RETRIEVED).stations)


ALIASES = {"mauritius": "mu-udm", "madagascar": "mg-tana", "south_africa": "za-rhodes"}


def test_alias_and_tart_id_resolve_to_same_station(catalogue):
    assert resolve_stations(["mauritius"], catalogue, ALIASES) == resolve_stations(["mu-udm"], catalogue, ALIASES)
    same = [resolve_stations([x], catalogue, ALIASES)[0] for x in ("za-rhodes", "za_rhodes", "ZA-Rhodes", "south_africa")]
    assert all(s.id == "za-rhodes" for s in same)


def test_alias_and_id_duplicate_rejected(catalogue):
    with pytest.raises(StationConfigError, match="more than once"):
        resolve_stations(["mauritius", "mu-udm"], catalogue, ALIASES)


def test_unknown_alias_fails_clearly(catalogue):
    with pytest.raises(StationConfigError, match="Unknown station 'atlantis'"):
        resolve_stations(["atlantis"], catalogue, ALIASES)
    with pytest.raises(StationConfigError, match="points to 'xx-gone'"):
        resolve_stations(["gone"], catalogue, {"gone": "xx-gone"})


def test_bundled_aliases_resolve_in_bundled_catalogue():
    catalogue = load_stations(default_station_file())
    aliases = load_aliases(default_alias_file())
    resolved = resolve_stations(["mauritius", "madagascar", "south_africa"], catalogue, aliases)
    assert [s.id for s in resolved] == ["mu-udm", "mg-tana", "za-rhodes"]
    assert all(s.archive_verified and s.archive_name == s.id for s in resolved)


def test_legacy_catalogue_still_resolves_by_name_and_tart_id():
    from tart_gnss.stations import legacy_station_file

    legacy = load_stations(legacy_station_file())
    aliases = load_aliases()
    assert resolve_stations(["mauritius"], legacy, aliases)[0].id == "mauritius"
    assert resolve_stations(["mu-udm"], legacy, aliases)[0].id == "mauritius"


# --- hashing and snapshots ----------------------------------------------


def test_snapshot_hash_is_deterministic(tmp_path):
    d = download(VALID)
    assert d.sha256 == registry.sha256_bytes(json.dumps(VALID).encode()) == download(copy.deepcopy(VALID)).sha256
    assert download(payload(entry("other"))).sha256 != d.sha256
    files = registry.write_snapshot(d, tmp_path, 3)
    assert files.raw.name == "registry_20261009T120000Z.json"
    assert registry.sha256_bytes(files.raw.read_bytes()) == d.sha256
    meta = json.loads(files.meta.read_text())
    assert meta == {"source_url": registry.DEFAULT_REGISTRY_URL, "retrieved_utc": "2026-10-09T12:00:00Z",
                    "sha256": d.sha256, "station_count": 3, "snapshot_file": files.raw.name}
    assert files.latest.read_bytes() == d.raw
    with pytest.raises(registry.RegistryError, match="refusing to overwrite"):
        registry.write_snapshot(d, tmp_path, 3)


# --- diff ----------------------------------------------------------------


def test_registry_diff(catalogue):
    new = dict(catalogue)
    del new["mg-tana"]
    rhodes = new["za-rhodes"]
    new["za-rhodes"] = replace(rhodes, latitude_deg=rhodes.latitude_deg + 0.001, name="Renamed")
    new["zz-new"] = replace(rhodes, id="zz-new", tart_id="zz-new")
    new["mu-udm"] = replace(new["mu-udm"], archive_name="mu-udm", archive_verified=True)
    diff = diff_catalogues(catalogue, new)
    assert diff.added == ["zz-new"] and diff.removed == ["mg-tana"]
    (move,) = diff.coordinate_changes
    assert move["tart_id"] == "za-rhodes"
    assert move["delta_latitude_deg"] == pytest.approx(0.001) and move["delta_longitude_deg"] == 0.0
    assert move["displacement_m"] == pytest.approx(110.9, abs=0.5)  # 0.001 deg of latitude near 33 S
    assert diff.name_changes == [{"tart_id": "za-rhodes", "old": rhodes.name, "new": "Renamed"}]
    assert {c["field"] for c in diff.metadata_changes} == {"archive_name", "archive_verified"}
    assert "moved 110" in diff.format()
    assert diff_catalogues(catalogue, catalogue).is_empty


# --- fetch / refresh -----------------------------------------------------


class FakeResponse:
    def __init__(self, content, status=200):
        self.content, self.status_code = content, status

    def raise_for_status(self):
        import requests

        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} error")


def test_fetch_success(monkeypatch):
    import requests

    monkeypatch.setattr(requests, "get", lambda *a, **k: FakeResponse(json.dumps(VALID).encode()))
    assert registry.fetch_registry() == VALID


@pytest.mark.parametrize("failure", ["timeout", "http", "json"])
def test_fetch_failure_leaves_catalogue_untouched(monkeypatch, tmp_path, failure, capsys):
    import requests

    def fake_get(*args, **kwargs):
        if failure == "timeout":
            raise requests.Timeout("timed out")
        if failure == "http":
            return FakeResponse(b"bad gateway", 502)
        return FakeResponse(b"<html>not json</html>")

    monkeypatch.setattr(requests, "get", fake_get)
    with pytest.raises(RegistryFetchError):
        registry.download_registry()

    catalogue = tmp_path / "stations_registry.yaml"
    catalogue.write_text("original")
    snapdir = tmp_path / "snaps"
    code = registry_cli.main_update(["--output", str(catalogue), "--snapshot-dir", str(snapdir)])
    assert code == 1
    assert catalogue.read_text() == "original" and not snapdir.exists()
    assert "remains available" in capsys.readouterr().err


def test_update_cli_writes_snapshot_catalogue_and_diff(monkeypatch, tmp_path):
    first = download(VALID)
    second = download(payload(entry("mu-udm", -20.26, 57.75, 20.0), entry("zz-new")),
                      when=datetime(2026, 10, 10, tzinfo=timezone.utc))
    downloads = iter([first, second])
    monkeypatch.setattr(registry, "download_registry", lambda *a, **k: next(downloads))
    out, snaps = tmp_path / "cat.yaml", tmp_path / "snaps"
    argv = ["--output", str(out), "--snapshot-dir", str(snaps), "--diff", str(out)]

    assert registry_cli.main_update(argv[:4]) == 0
    assert set(load_stations(out)) == {"mu-udm", "mg-tana", "za-rhodes"}
    assert registry_cli.main_update(argv) == 0
    assert set(load_stations(out)) == {"mu-udm", "zz-new"}
    assert sorted(p.name for p in snaps.iterdir()) == [
        "latest.json",
        "registry_20261009T120000Z.diff.json",  # first run: diffed against the packaged catalogue
        "registry_20261009T120000Z.json", "registry_20261009T120000Z.meta.json",
        "registry_20261010T000000Z.diff.json", "registry_20261010T000000Z.json",
        "registry_20261010T000000Z.meta.json",
    ]
    change = json.loads((snaps / "registry_20261010T000000Z.diff.json").read_text())
    assert change["added"] == ["zz-new"] and sorted(change["removed"]) == ["mg-tana", "za-rhodes"]
    assert change["coordinate_changes"][0]["tart_id"] == "mu-udm"
    assert load_catalogue_metadata(out)["changes_from_previous"]["added"] == ["zz-new"]


def test_update_cli_dry_run_writes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr(registry, "download_registry", lambda *a, **k: download(VALID))
    out = tmp_path / "cat.yaml"
    assert registry_cli.main_update(["--output", str(out), "--snapshot-dir", str(tmp_path / "s"), "--dry-run"]) == 0
    assert not out.exists() and not (tmp_path / "s").exists()


def test_list_cli(capsys):
    assert registry_cli.main_list([]) == 0
    out = capsys.readouterr().out
    assert "mu-udm" in out and "mauritius" in out
    assert registry_cli.main_list(["--archive-only", "--format", "json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert rows and all(r["archive_name"] for r in rows)
    assert next(r for r in rows if r["id"] == "za-rhodes")["aliases"] == ["south_africa"]


def test_failed_archive_check_keeps_previous_verification(monkeypatch, tmp_path):
    parsed = parse_registry(VALID, RETRIEVED)
    old = registry_to_stations(parsed.stations, {"mu-udm": registry.ArchiveStatus("mu-udm", True, "2026-10-01T00:00:00Z")})
    prev = tmp_path / "prev.yaml"
    registry.write_catalogue(prev, registry.build_catalogue_document(
        parsed, old, source_url="u", retrieved_utc=RETRIEVED, sha256="0" * 64, snapshot_file=None))

    def flaky(prefix):
        raise RegistryFetchError("timeout")

    monkeypatch.setattr(registry, "download_registry", lambda *a, **k: download(VALID))
    monkeypatch.setattr(registry, "s3_prefix_exists", flaky)
    monkeypatch.setattr(registry.verify_archives, "__defaults__", (flaky, None))
    result = registry_cli.refresh_station_catalogue(output=tmp_path / "new.yaml", snapshot_dir=tmp_path / "s",
                                                    previous=prev, verify_archives=True, dry_run=True)
    assert result.stations["mu-udm"].archive_verified_utc == "2026-10-01T00:00:00Z"
    assert result.stations["mg-tana"].archive_verified is False
