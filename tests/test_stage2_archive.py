from datetime import timedelta

import pytest

from stage2_helpers import MockArchive, archive_key, make_tart_hdf5, utc
from tart_gnss.stage2.archive import day_prefixes, find_observations, parse_filename_utc
from tart_gnss.stage2.errors import ArchiveUnavailableError
from tart_gnss.stations import Station

ALPHA = Station("alpha", "a", -20.0, 57.0, 0.0, "test", tart_id="aa-alpha", archive_name="aa-alpha")


def test_parse_filename():
    assert parse_filename_utc("vis_2026-10-01_00_00_36.068519.hdf") == utc(2026, 10, 1, 0, 0, 36, 68519)
    assert parse_filename_utc("vis_2026-10-01_00_00_36.hdf") == utc(2026, 10, 1, 0, 0, 36)
    assert parse_filename_utc("README.txt") is None


def test_day_prefixes_cross_midnight_unpadded():
    assert day_prefixes("za-rhodes", utc(2026, 6, 30, 23, 50), utc(2026, 7, 1, 0, 10)) == [
        "za-rhodes/vis/2026/6/30/", "za-rhodes/vis/2026/7/1/"]


def test_find_observations_uses_filename_time_and_margin(tmp_path):
    arch = MockArchive()
    start = utc(2026, 6, 25, 10, 0)
    for minutes in (-3, 0.5, 5, 10.5, 12, 30):
        t = start + timedelta(minutes=minutes)
        arch.add_file(archive_key("aa-alpha", t), make_tart_hdf5(tmp_path / f"{minutes}.hdf", "alpha", t, n=2))
    recs = find_observations(ALPHA, start, start + timedelta(minutes=10), arch, file_span_margin_s=120)
    times = [r.filename_utc for r in recs]
    # 10:00:30 .. 10:12:00 kept (12:00 lies within the 120 s margin); 09:57 and 10:30 rejected
    assert times == [start + timedelta(minutes=m) for m in (0.5, 5, 10.5, 12)]
    assert all(r.start_utc is None for r in recs), "sample times are not guessed before reading files"
    assert recs[0].remote_identifier.startswith("aa-alpha/vis/2026/6/25/")


def test_no_archive_name():
    station = Station("x", "x", 0, 0, 0, "test")
    with pytest.raises(ArchiveUnavailableError, match="archive_name"):
        find_observations(station, utc(2026, 1, 1), utc(2026, 1, 2), MockArchive())


def test_archive_failure_is_explicit():
    arch = MockArchive()
    arch.fail_list.add("aa-alpha")
    with pytest.raises(ArchiveUnavailableError):
        find_observations(ALPHA, utc(2026, 1, 1), utc(2026, 1, 1, 1), arch)
