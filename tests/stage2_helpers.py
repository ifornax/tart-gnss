"""Shared fixtures for Stage-2 tests: synthetic TART HDF5 files and a mock archive.

Synthetic files follow the layout of real TART archive files (config JSON,
timestamp, vis, baselines, antenna_positions, gains, phases, phase_elaz).
Their contents are random numbers, NOT real observations.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np

from tart_gnss.stage2.archive import ArchiveObject

STATIONS_YAML = """
stations:
  alpha:
    name: "Test station alpha"
    tart_id: "aa-alpha"
    archive_name: "aa-alpha"
    latitude_deg: -20.0
    longitude_deg: 57.0
    altitude_m: 0.0
    source: "synthetic test value"
  beta:
    name: "Test station beta"
    tart_id: "bb-beta"
    archive_name: "bb-beta"
    latitude_deg: -19.0
    longitude_deg: 47.0
    altitude_m: 1000.0
    source: "synthetic test value"
  gamma:
    name: "Test station gamma"
    tart_id: "cc-gamma"
    archive_name: "cc-gamma"
    latitude_deg: -33.0
    longitude_deg: 26.0
    altitude_m: 500.0
    source: "synthetic test value"
"""

STATION_COORDS = {"alpha": (-20.0, 57.0, 0.0), "beta": (-19.0, 47.0, 1000.0), "gamma": (-33.0, 26.0, 500.0)}


def utc(*args) -> datetime:
    return datetime(*args, tzinfo=timezone.utc)


def iso(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


def make_tart_hdf5(path: Path, station: str, t_end: datetime, n: int = 60, cadence: float = 1.0,
                   frequency: float = 1575.42e6, n_ant: int = 24, lat=None, zero_gains=(),
                   timestamps: List[datetime] = None) -> Path:
    """Write a synthetic native TART file whose last sample is ``t_end``."""
    import h5py

    lat0, lon0, alt0 = STATION_COORDS[station]
    if timestamps is None:
        timestamps = [t_end - timedelta(seconds=cadence * (n - 1 - k)) for k in range(n)]
    n_bl = n_ant * (n_ant - 1) // 2
    rng = np.random.default_rng(0)
    gains = np.ones(n_ant, dtype=np.float32)
    gains[list(zero_gains)] = 0.0
    # same keys as the config of a real za-rhodes archive file
    config = {"name": station, "frequency": frequency, "L0_frequency": 1571328000.0,
              "baseband_frequency": 4092000.0, "sampling_frequency": 16368000.0, "bandwidth": 2.5e6,
              "num_antenna": n_ant, "lat": lat0 if lat is None else lat, "lon": lon0, "alt": alt0,
              "orientation": 0.0, "axes": ["East", "North", "Up"]}
    with h5py.File(path, "w") as f:
        f.create_dataset("config", data=np.array([json.dumps(config)], dtype=object), dtype=h5py.string_dtype())
        f.create_dataset("timestamp", data=np.array([t.isoformat() for t in timestamps], dtype=object),
                         dtype=h5py.string_dtype())
        vis = (rng.normal(size=(len(timestamps), n_bl)) + 1j * rng.normal(size=(len(timestamps), n_bl))).astype(np.complex64)
        f.create_dataset("vis", data=vis)
        f.create_dataset("baselines", data=np.array([(i, j) for i in range(n_ant) for j in range(i + 1, n_ant)]))
        f.create_dataset("antenna_positions", data=rng.uniform(-1, 1, size=(n_ant, 3)).astype(np.float32))
        f.create_dataset("gains", data=gains)
        f.create_dataset("phases", data=np.zeros(n_ant, dtype=np.float32))
        f.create_dataset("phase_elaz", data=np.array([90.0, 0.0]))
    return path


def archive_key(archive_name: str, t_end: datetime) -> str:
    name = "vis_" + t_end.strftime("%Y-%m-%d_%H_%M_%S.%f") + ".hdf"
    return f"{archive_name}/vis/{t_end.year}/{t_end.month}/{t_end.day}/{name}"


class MockArchive:
    """In-memory archive backend serving bytes by key."""

    def __init__(self) -> None:
        self.objects: Dict[str, bytes] = {}
        self.downloads: List[str] = []
        self.fail_list: set = set()

    def add_file(self, key: str, path: Path) -> None:
        self.objects[key] = Path(path).read_bytes()

    def list_objects(self, prefix: str):
        from tart_gnss.stage2.errors import ArchiveUnavailableError

        if any(prefix.startswith(p) for p in self.fail_list):
            raise ArchiveUnavailableError(f"mock failure for {prefix}")
        return [ArchiveObject(k, len(v), hashlib.md5(v).hexdigest(), None)  # noqa: S324
                for k, v in sorted(self.objects.items()) if k.startswith(prefix)]

    def download(self, key: str, destination: Path) -> None:
        self.downloads.append(key)
        Path(destination).write_bytes(self.objects[key])

    def url(self, key: str) -> str:
        return f"mock://{key}"


def populate(archive: MockArchive, tmp: Path, station: str, start: datetime, n_files: int, offset_s: float = 0.0):
    """Add ``n_files`` consecutive 60-s synthetic files for a station."""
    tmp.mkdir(parents=True, exist_ok=True)
    for k in range(n_files):
        t_end = start + timedelta(seconds=offset_s + 60.0 * (k + 1) - 1.0)
        p = make_tart_hdf5(tmp / f"{station}_{k}.hdf", station, t_end)
        archive.add_file(archive_key(f"{dict(alpha='aa-alpha', beta='bb-beta', gamma='cc-gamma')[station]}", t_end), p)


def stage1_doc(windows=None) -> dict:
    """A minimal Stage-1 results document using the synthetic stations."""
    stations = [
        {"id": sid, "name": sid, "latitude_deg": c[0], "longitude_deg": c[1], "altitude_m": c[2],
         "source": "synthetic test value", "tart_id": None}
        for sid, c in STATION_COORDS.items()
    ]
    tle = {"epoch": "2026-06-24T00:00:00.000000Z", "line1": "1 ...", "line2": "2 ...",
           "source": "test", "age_start_days": 1.0, "age_mid_days": 1.0, "age_end_days": 1.0}
    geom = {sid: {e: {"azimuth_deg": 0.0, "elevation_deg": 45.0, "range_km": 2e4, "range_rate_km_s": 0.0}
                  for e in ("start", "mid", "end")} for sid in STATION_COORDS}
    if windows is None:
        windows = [
            (1, "GPS BIIF-9  (PRN 26)", 40534, "GPS", utc(2026, 6, 25, 10, 0), utc(2026, 6, 25, 10, 20)),
            (2, "GSAT0211 (GALILEO 14)", 41549, "Galileo", utc(2026, 6, 25, 12, 0), utc(2026, 6, 25, 12, 15)),
            (3, "GSAT0211 (GALILEO 14)", 41549, "Galileo", utc(2026, 6, 25, 20, 0), utc(2026, 6, 25, 20, 10)),
            (4, "BEIDOU-2 G6 (C02)", 38953, "BeiDou", utc(2026, 6, 25, 0, 0), utc(2026, 6, 25, 0, 5)),
        ]
    out = []
    for rank, name, norad, const, a, b in windows:
        out.append({"rank": rank, "satellite": name, "norad_id": norad, "constellation": const,
                    "start_utc": iso(a), "mid_utc": iso(a + (b - a) / 2), "end_utc": iso(b),
                    "duration_seconds": (b - a).total_seconds(),
                    "start_clipped_by_search_interval": False, "end_clipped_by_search_interval": False,
                    "stations": geom, "satellite_position": {"coordinate_frame": "ITRS"},
                    "line_of_sight_mid": {"coordinate_frame": "ITRS"}, "tle": tle})
    return {"run": {"start_utc": "2026-06-25T00:00:00.000Z", "end_utc": "2026-06-26T00:00:00.000Z",
                    "min_elevation_deg": 20.0, "step_seconds": 10.0},
            "stations": stations, "windows": out}
