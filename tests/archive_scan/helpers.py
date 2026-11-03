"""Mocked archive backend and synthetic Stage-1 results for the archive-scan tests (no network)."""

from __future__ import annotations

import csv
import json
import random
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from tart_gnss.archive_scan import cli
from tart_gnss.stage2.archive import ArchiveObject
from tart_gnss.stage2.errors import ArchiveUnavailableError
from tart_gnss.stations import load_stations
from tart_gnss.tle import format_utc

UTC = timezone.utc
TRIAD = ("mg-tana", "za-rhodes", "mu-udm")
#: Fixed "now" well after the synthetic days, so every listing covers settled days.
NOW = datetime(2026, 11, 10, 12, tzinfo=UTC)


def t(text: str) -> datetime:
    """Parse ``YYYY-MM-DDTHH:MM:SS[.f]`` as UTC."""
    if "." in text:  # Python 3.9 fromisoformat needs 3 or 6 fractional digits
        head, frac = text.split(".")
        text = f"{head}.{frac.ljust(6, '0')}"
    return datetime.fromisoformat(text).replace(tzinfo=UTC)


def vis_key(archive_name: str, when: datetime) -> str:
    """Archive key of a file whose filename time is ``when`` (the real layout)."""
    return (f"{archive_name}/vis/{when.year}/{when.month}/{when.day}/"
            f"vis_{when:%Y-%m-%d_%H_%M_%S}.{when.microsecond:06d}.hdf")


def run_of_files(first_end: datetime, n: int, step_s: float = 60.0) -> List[datetime]:
    """``n`` consecutive one-minute files, the first ending at ``first_end``."""
    return [first_end + timedelta(seconds=i * step_s) for i in range(n)]


class MockArchive:
    """In-memory ArchiveBackend: records every listing call and refuses downloads."""

    def __init__(self, files: Optional[Dict[str, Iterable[datetime]]] = None,
                 errors: Optional[Dict[str, int]] = None, extra: Optional[Dict[str, List[ArchiveObject]]] = None,
                 shuffle_seed: Optional[int] = None):
        self.objects: List[ArchiveObject] = []
        for name, times in (files or {}).items():
            for when in times:
                self.objects.append(ArchiveObject(key=vis_key(name, when), size_bytes=1_000_000,
                                                  etag="abc", last_modified_utc=when + timedelta(minutes=3)))
        for objs in (extra or {}).values():
            self.objects.extend(objs)
        self.errors = dict(errors or {})  # prefix -> number of failures before success (-1: always)
        self.calls: List[str] = []
        self.seed = shuffle_seed

    def list_objects(self, prefix: str):
        self.calls.append(prefix)
        n = self.errors.get(prefix, 0)
        if n != 0:
            if n > 0:
                self.errors[prefix] = n - 1
            raise ArchiveUnavailableError(f"mock failure listing {prefix}")
        out = [o for o in self.objects if o.key.startswith(prefix)]
        if self.seed is not None:
            random.Random(self.seed).shuffle(out)
        return out

    def download(self, key, destination):
        raise AssertionError("the archive scanner must never download")

    def url(self, key):
        return f"mock://{key}"


def registry_station(tart_id: str) -> dict:
    """A Stage-1 station entry with the packaged registry coordinates."""
    st = load_stations()[tart_id]
    return {"id": st.id, "name": st.name, "latitude_deg": st.latitude_deg, "longitude_deg": st.longitude_deg,
            "altitude_m": st.altitude_m, "source": st.source, "tart_id": st.tart_id, "api_url": st.api_url,
            "archive_name": st.archive_name}


def window(rank: int, norad: int, name: str, constellation: str, start: datetime, end: datetime,
           epoch: Optional[datetime] = None, elevation: float = 45.0, stations: Sequence[str] = TRIAD) -> dict:
    """A Stage-1 window record with the fields the scanner reads."""
    w = {
        "rank": rank, "satellite": name, "norad_id": norad, "constellation": constellation,
        "start_utc": format_utc(start, 3), "mid_utc": format_utc(start + (end - start) / 2, 3),
        "end_utc": format_utc(end, 3), "duration_seconds": (end - start).total_seconds(),
        "start_clipped_by_search_interval": False, "end_clipped_by_search_interval": False,
        "stations": {sid: {"mid": {"azimuth_deg": 0.0, "elevation_deg": elevation}} for sid in stations},
    }
    if epoch is not None:
        w["tle"] = {"epoch": format_utc(epoch, 6), "line1": "", "line2": "", "source": "test"}
    return w


def write_stage1(path: Path, windows: List[dict], start: datetime, end: datetime,
                 stations: Optional[List[dict]] = None) -> Path:
    """Write a synthetic Stage-1 JSON."""
    doc = {
        "run": {"start_utc": format_utc(start, 3), "end_utc": format_utc(end, 3), "min_elevation_deg": 20.0},
        "stations": stations or [registry_station(s) for s in TRIAD],
        "tle_catalogue": {"tle_retrieved_utc": "2026-11-03T13:15:15.000Z"},
        "windows": windows,
    }
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


def scan(tmp_path: Path, archive: MockArchive, *argv: str, now: datetime = NOW, sleep=None, cache_dir=None):
    """Run the scanner in-process; returns (exit code, manifest dict, output paths)."""
    out = tmp_path / "out" / "archive_candidates.csv"
    args = ["--output", str(out), "--cache-dir", str(cache_dir or tmp_path / "cache"), *argv]
    ns = cli.build_parser().parse_args(args)
    code = cli.run(ns, ["tart-archive-scan", *args], archive=archive, now=lambda: now,
                   sleep=sleep or (lambda s: None))
    from tart_gnss.archive_scan.output import output_paths

    paths = output_paths(out)
    manifest = json.loads(paths["manifest"].read_text()) if paths["manifest"].exists() else None
    return code, manifest, paths


def read_csv(path: Path) -> List[Dict[str, str]]:
    """Rows of a CSV as dictionaries."""
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))
