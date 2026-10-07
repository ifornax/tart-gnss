import json
from datetime import datetime, timezone

import pytest

from tart_gnss import tle

GPS_L1 = "1 39533U 14008A   26278.86890240  .00000025  00000+0  00000+0 0  9993"
GPS_L2 = "2 39533  53.6663 144.8014 0083568 231.5298 127.8020  2.00566942 91909"


def test_checksum():
    assert tle.tle_checksum(GPS_L1) == int(GPS_L1[68])
    assert tle.tle_checksum(GPS_L2) == int(GPS_L2[68])


def test_parse_three_and_two_line():
    text = f"GPS BIIF-5  (PRN 30)\n{GPS_L1}\n{GPS_L2}\n{GPS_L1}\n{GPS_L2}\n"
    parsed = tle.parse_tle_text(text)
    assert parsed[0][0] == "GPS BIIF-5  (PRN 30)"
    assert parsed[1][0] == "NORAD 39533"


def test_bad_checksum_rejected():
    bad = GPS_L1[:68] + str((int(GPS_L1[68]) + 1) % 10)
    with pytest.raises(tle.TLEError, match="checksum"):
        tle.parse_tle_text(f"X\n{bad}\n{GPS_L2}\n")


def test_epoch():
    epoch = tle.tle_epoch_utc(GPS_L1, GPS_L2)
    # day-of-year 278.86890240 of 2026 = 2026-10-05 20:51:13.167 UTC
    assert epoch.tzinfo is not None
    assert abs((epoch - datetime(2026, 10, 5, 20, 51, 13, 167000, tzinfo=timezone.utc)).total_seconds()) < 1e-3


@pytest.mark.parametrize(
    "name,expected",
    [
        ("GPS BIIF-5  (PRN 30)", "GPS"),
        ("NAVSTAR 81 (USA 319)", "GPS"),
        ("GSAT0219 (GALILEO 23)", "Galileo"),
        ("BEIDOU-2 G6 (C02)", "BeiDou"),
        ("COSMOS 2433 (720)", "GLONASS"),
        ("COSMOS 2251", "UNKNOWN"),
        ("ISS (ZARYA)", "UNKNOWN"),
    ],
)
def test_classification(name, expected):
    assert tle.classify_constellation(name) == expected


def test_load_file_and_filter(sample_tle_path):
    cat = tle.load_tle_file(sample_tle_path)
    assert len(cat.records) == 4
    assert {r.constellation for r in cat.records} == {"GPS", "Galileo", "BeiDou", "GLONASS"}
    assert len(cat.sha256) == 64
    only = tle.filter_by_constellation(cat.records, ["gps", "Galileo"])
    assert {r.constellation for r in only} == {"GPS", "Galileo"}
    with pytest.raises(tle.TLEError):
        tle.filter_by_constellation(cat.records, ["Iridium"])


def test_download_archives_with_provenance(tmp_path, sample_tle_path):
    text = sample_tle_path.read_text()
    calls = []

    def fake_fetch(url, timeout):
        calls.append(url)
        return text if "gps-ops" in url else ""

    path = tle.download_gnss_catalogue(tmp_path, ["GPS"], fetch=fake_fetch)
    meta = json.loads(path.with_name(path.name + ".meta.json").read_text())
    assert calls and "gps-ops" in calls[0]
    assert meta["sha256"] == tle.sha256_file(path)
    assert meta["retrieved_utc"].endswith("Z")
    # constellation recorded from the download group, not from the name
    cat = tle.load_tle_file(path)
    assert all(r.constellation == "GPS" for r in cat.records)
    assert "celestrak" in cat.source.lower()


def test_modified_catalogue_detected(tmp_path, sample_tle_path):
    path = tle.download_gnss_catalogue(
        tmp_path, ["GPS"], fetch=lambda url, timeout: sample_tle_path.read_text()
    )
    path.write_text(path.read_text() + "\n")
    with pytest.raises(tle.TLEError, match="SHA-256"):
        tle.load_tle_file(path)
