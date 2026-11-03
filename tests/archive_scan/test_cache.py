"""Inventory cache: refresh lineage, reuse TTL and partial days, offline, checksum validation."""

import json
from datetime import date, timedelta

import pytest

from tart_gnss.archive_scan.cache import InventoryCache, listing_from_dict
from tart_gnss.archive_scan.errors import InventoryCacheError
from tart_gnss.archive_scan.inventory import STATUS_ERROR, STATUS_UNKNOWN, build_inventory, day_prefix
from tart_gnss.stations import load_stations

from .helpers import NOW, MockArchive, run_of_files, t

ST = [load_stations()["mg-tana"]]
DAY = date(2026, 10, 9)
PREFIX = day_prefix("mg-tana", DAY)


def _run(cache, archive, policy, now=NOW, max_age=24.0):
    def factory():
        if archive is None:
            pytest.fail("no S3 client may be constructed")
        return archive
    return build_inventory(ST, [DAY], factory, cache, policy, max_age, "host", "tart", "test",
                           now=lambda: now, sleep=lambda s: None)


def _archive():
    return MockArchive({"mg-tana": run_of_files(t("2026-10-09T03:00:00"), 5)})


def test_refresh_writes_new_snapshot_and_keeps_history(tmp_path):
    cache = InventoryCache(tmp_path)
    _run(cache, _archive(), "refresh")
    _run(cache, _archive(), "refresh", now=NOW + timedelta(minutes=5))
    snaps = cache.snapshots("host", "tart", PREFIX)
    assert len(snaps) == 2
    doc = json.loads(snaps[-1].read_text())
    for key in ("schema_version", "station_id", "archive_name", "bucket", "host", "prefix", "listing_started_utc",
                "listing_completed_utc", "listing_status", "objects", "source", "software_version",
                "checksum_of_canonical_listing"):
        assert key in doc
    assert set(doc["objects"][0]) == {"key", "size_bytes", "etag", "last_modified_utc", "filename_utc"}


def test_reuse_within_ttl_and_refresh_when_stale(tmp_path):
    cache = InventoryCache(tmp_path)
    _run(cache, _archive(), "refresh")
    listings, stats = _run(cache, None, "reuse", now=NOW + timedelta(hours=2))
    assert stats.cache_hits == 1 and stats.live_list_calls == 0 and listings[("mg-tana", DAY)].n_files == 5
    archive = _archive()
    listings, stats = _run(cache, archive, "reuse", now=NOW + timedelta(hours=30))
    assert stats.cache_stale == 1 and len(archive.calls) == 1 and listings[("mg-tana", DAY)].source == "live"


def test_reuse_relists_a_day_cached_before_it_had_settled(tmp_path):
    cache = InventoryCache(tmp_path)
    during = t("2026-10-09T18:00:00")
    _run(cache, _archive(), "refresh", now=during)
    archive = _archive()
    listings, stats = _run(cache, archive, "reuse", now=during + timedelta(hours=8))
    assert stats.cache_stale == 1 and len(archive.calls) == 1
    assert listings[("mg-tana", DAY)].covers_full_day is True


def test_offline_uses_cache_only_and_missing_is_unknown(tmp_path):
    cache = InventoryCache(tmp_path)
    listings, stats = _run(cache, None, "offline")
    assert listings[("mg-tana", DAY)].status == STATUS_UNKNOWN and stats.unknown == 1
    _run(cache, _archive(), "refresh")
    listings, stats = _run(cache, None, "offline", now=NOW + timedelta(days=400))  # TTL ignored offline
    assert listings[("mg-tana", DAY)].n_files == 5 and stats.cache_hits == 1


def test_checksum_mismatch_is_rejected_and_falls_back(tmp_path):
    cache = InventoryCache(tmp_path)
    _run(cache, _archive(), "refresh")
    _run(cache, MockArchive({"mg-tana": run_of_files(t("2026-10-09T03:00:00"), 7)}), "refresh",
         now=NOW + timedelta(minutes=1))
    newest = cache.snapshots("host", "tart", PREFIX)[-1]
    doc = json.loads(newest.read_text())
    doc["objects"].pop()
    newest.write_text(json.dumps(doc))
    with pytest.raises(InventoryCacheError, match="checksum"):
        listing_from_dict(doc, newest)
    listings, stats = _run(cache, None, "offline")
    assert stats.cache_corrupt == 1 and listings[("mg-tana", DAY)].n_files == 5  # older valid snapshot


def test_incompatible_schema_and_garbage_are_rejected(tmp_path):
    cache = InventoryCache(tmp_path)
    _run(cache, _archive(), "refresh")
    path = cache.snapshots("host", "tart", PREFIX)[-1]
    doc = json.loads(path.read_text())
    doc["schema_version"] = 99
    path.write_text(json.dumps(doc))
    (path.parent / "listing_99999999T999999999999Z.json").write_text("{not json")
    listings, stats = _run(cache, None, "offline")
    assert listings[("mg-tana", DAY)].status == STATUS_UNKNOWN and stats.cache_corrupt == 2


def test_failed_listing_is_never_cached(tmp_path):
    cache = InventoryCache(tmp_path)
    listings, _ = _run(cache, MockArchive(errors={PREFIX: -1}), "refresh")
    assert listings[("mg-tana", DAY)].status == STATUS_ERROR
    assert cache.snapshots("host", "tart", PREFIX) == []
    with pytest.raises(InventoryCacheError):
        cache.write(listings[("mg-tana", DAY)])


def test_cache_key_separates_host_bucket_and_station(tmp_path):
    cache = InventoryCache(tmp_path)
    assert cache.directory("h1", "tart", PREFIX) != cache.directory("h2", "tart", PREFIX)
    assert cache.directory("h1", "tart", PREFIX) != cache.directory("h1", "other", PREFIX)
    assert cache.directory("h1", "tart", PREFIX) != cache.directory("h1", "tart", day_prefix("mu-udm", DAY))
