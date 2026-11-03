"""Station/day listing: prefixes, edge days, one call per prefix, parsing, statuses, retry."""

from datetime import date, timedelta

import pytest

from tart_gnss.archive_scan.cache import InventoryCache
from tart_gnss.archive_scan.inventory import (
    STATUS_EMPTY,
    STATUS_ERROR,
    STATUS_FOUND,
    STATUS_UNKNOWN,
    build_inventory,
    canonicalise,
    day_prefix,
    list_day_live,
    required_days,
)
from tart_gnss.stage2.archive import ArchiveObject, day_prefixes
from tart_gnss.stations import load_stations

from .helpers import NOW, TRIAD, MockArchive, run_of_files, scan, t, vis_key, window, write_stage1

STATIONS = [load_stations()[s] for s in TRIAD]


def _inventory(tmp_path, archive, days, policy="refresh", **kw):
    return build_inventory(STATIONS, days, lambda: archive, InventoryCache(tmp_path / "cache"), policy, 24.0,
                           "host", "tart", "test", now=lambda: NOW, sleep=lambda s: None, **kw)


def test_day_prefix_matches_stage2_and_is_not_zero_padded():
    assert day_prefix("mg-tana", date(2026, 1, 5)) == "mg-tana/vis/2026/1/5/"
    assert day_prefix("mu-udm", date(2026, 10, 31)) == "mu-udm/vis/2026/10/31/"
    assert day_prefix("za-rhodes", date(2026, 11, 3)) == day_prefixes("za-rhodes", t("2026-11-03T05:00:00"),
                                                                     t("2026-11-03T06:00:00"))[0]


def test_required_days_month_and_year_transitions():
    days = required_days(t("2025-12-31T12:00:00"), t("2026-01-01T12:00:00"), 120)
    assert days == [date(2025, 12, 31), date(2026, 1, 1)]
    days = required_days(t("2026-09-30T12:00:00"), t("2026-10-01T12:00:00"), 120)
    assert [day_prefix("x", d) for d in days] == ["x/vis/2026/9/30/", "x/vis/2026/10/1/"]


def test_required_days_edge_allowance():
    # start at midnight: start - span crosses into the previous day; end + span crosses into the next one.
    days = required_days(t("2026-10-09T00:00:00"), t("2026-11-03T00:00:00"), 120)
    assert days[0] == date(2026, 10, 8) and days[-1] == date(2026, 11, 3) and len(days) == 27
    # well inside a day: only that day.
    assert required_days(t("2026-10-09T06:00:00"), t("2026-10-09T07:00:00"), 120) == [date(2026, 10, 9)]


def test_exactly_one_listing_per_station_day_over_25_days(tmp_path):
    start, end = t("2026-10-09T00:00:00"), t("2026-11-03T00:00:00")
    files = {s: run_of_files(t("2026-10-15T10:01:00"), 20) for s in TRIAD}
    wins = [window(i + 1, 1000 + i, f"GPS TEST (PRN {i})", "GPS", start + timedelta(hours=7 * i),
                   start + timedelta(hours=7 * i + 5), epoch=start) for i in range(80)]
    s1 = write_stage1(tmp_path / "s1.json", wins, start, end)
    archive = MockArchive(files)
    code, m, _ = scan(tmp_path, archive, "--visibility-results", str(s1))
    days = required_days(start, end, 120)
    assert len(archive.calls) == len(set(archive.calls)) == 3 * len(days) == 81
    assert m["listing"]["live_list_calls"] == 81 and m["listing"]["prefixes_required"] == 81
    assert code == 0


def test_canonicalise_parses_filters_and_deduplicates():
    a = t("2026-10-09T12:00:00.5")
    good = ArchiveObject(vis_key("mg-tana", a), 10, "e", None)
    later = ArchiveObject(vis_key("mg-tana", a + timedelta(minutes=1)), 10, "e", None)
    junk = [ArchiveObject("mg-tana/vis/2026/10/9/readme.txt", 1, None, None),
            ArchiveObject("mg-tana/vis/2026/10/9/vis_2026-10-09_12_00.hdf", 1, None, None)]
    objs, ignored, dups = canonicalise([later, good, junk[0], good, junk[1]])
    assert [o.key for o in objs] == [good.key, later.key]
    assert objs[0].filename_utc == a and dups == 1 and len(ignored) == 2


def test_found_empty_error_unknown_are_distinct(tmp_path):
    d = date(2026, 10, 9)
    archive = MockArchive({"mg-tana": run_of_files(t("2026-10-09T01:00:00"), 3)},
                          errors={day_prefix("mu-udm", d): -1})
    listings, stats = _inventory(tmp_path, archive, [d])
    assert listings[("mg-tana", d)].status == STATUS_FOUND
    assert listings[("za-rhodes", d)].status == STATUS_EMPTY
    err = listings[("mu-udm", d)]
    assert err.status == STATUS_ERROR and err.n_files == 0 and "mock failure" in err.error and not err.complete
    assert stats.failures == 1
    # offline: the failed day was never cached, so it is UNKNOWN (not EMPTY); successes come from cache.
    listings2, stats2 = _inventory(tmp_path, MockArchive(), [d], policy="offline")
    assert listings2[("mu-udm", d)].status == STATUS_UNKNOWN
    assert listings2[("za-rhodes", d)].status == STATUS_EMPTY and listings2[("za-rhodes", d)].source == "cache"
    assert stats2.cache_hits == 2 and stats2.unknown == 1


def test_bounded_retry_then_success():
    st = STATIONS[0]
    d = date(2026, 10, 9)
    archive = MockArchive({"mg-tana": run_of_files(t("2026-10-09T01:00:00"), 2)}, errors={day_prefix("mg-tana", d): 2})
    sleeps = []
    L = list_day_live(archive, st, d, "h", "b", "v", attempts=3, backoff_s=1.0, sleep=sleeps.append)
    assert L.status == STATUS_FOUND and L.attempts == 3 and sleeps == [1.0, 2.0]
    archive = MockArchive(errors={day_prefix("mg-tana", d): -1})
    L = list_day_live(archive, st, d, "h", "b", "v", attempts=2, sleep=lambda s: None)
    assert L.status == STATUS_ERROR and len(archive.calls) == 2


def test_station_without_archive_name_is_an_error(tmp_path):
    from dataclasses import replace

    st = replace(STATIONS[0], archive_name=None)
    listings, stats = build_inventory([st], [date(2026, 10, 9)], lambda: pytest.fail("no listing expected"),
                                      InventoryCache(tmp_path), "refresh", 24, "h", "b", "v")
    L = listings[(st.id, date(2026, 10, 9))]
    assert L.status == STATUS_ERROR and "archive_name" in L.error and stats.failures == 1
