"""Synthetic TART-like observations for tests and the end-to-end demonstration.

Visibilities are generated with the same forward model as the estimator
(:func:`tart_gnss.stage3.visibility_model.predict_visibilities`), and UVW is
written in the ``tart2ms`` convention (``UVW = r_ANTENNA1 - r_ANTENNA2``,
zenith phase centre), so a synthetic Measurement Set exercises exactly the
code path used for real data. Nothing here is real TART data.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage3.geometry import azel_to_unit, enu_to_ecef, ecef_to_enu
from tart_gnss.stage3.loader import VisibilityDataset, dataset_from_arrays
from tart_gnss.stage3.visibility_model import wavenumber

GPS_L1_HZ = 1575.42e6


def tart_like_layout(n_ant: int = 24, radius_m: float = 2.0, seed: int = 1) -> np.ndarray:
    """Random, roughly planar antenna layout (ENU metres) inside a disc of ``radius_m``.

    A synthetic stand-in for a TART array (24 antennas, baselines up to ~4 m);
    antenna heights scatter by ~1 cm.
    """
    rng = np.random.default_rng(seed)
    r = radius_m * np.sqrt(rng.uniform(0, 1, n_ant))
    th = rng.uniform(0, 2 * np.pi, n_ant)
    return np.stack([r * np.sin(th), r * np.cos(th), 0.01 * rng.standard_normal(n_ant)], axis=1)


@dataclass
class SyntheticObservation:
    """MS-like arrays of a synthetic observation."""

    time_mjd_seconds: np.ndarray
    antenna1: np.ndarray
    antenna2: np.ndarray
    uvw_m: np.ndarray
    data: np.ndarray
    flag: np.ndarray
    weight: np.ndarray
    frequency_hz: np.ndarray
    antenna_positions_ecef_m: np.ndarray
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    true_az_deg: np.ndarray
    true_el_deg: np.ndarray
    amplitude: complex
    noise_sigma: float


def simulate_observation(
    antenna_enu_m: np.ndarray,
    latitude_deg: float,
    longitude_deg: float,
    altitude_m: float,
    times_mjd_seconds: Sequence[float],
    az_deg,
    el_deg,
    amplitude: complex = 1.0 + 0.0j,
    noise_sigma: float = 0.0,
    frequency_hz: Sequence[float] = (GPS_L1_HZ,),
    flag_fraction: float = 0.0,
    flagged_antennas: Sequence[int] = (),
    extra_sources: Sequence[Tuple[float, float, complex]] = (),
    rng: Optional[np.random.Generator] = None,
) -> SyntheticObservation:
    """Simulate a point source moving along ``(az_deg[t], el_deg[t])``.

    ``az_deg``/``el_deg`` are scalars or one value per time. ``noise_sigma`` is
    the per-component (real, imaginary) Gaussian noise. ``flag_fraction`` of
    the cross-correlation samples are flagged at random; baselines touching
    ``flagged_antennas`` are flagged entirely. ``extra_sources`` adds fixed
    ``(az, el, amplitude)`` interferers.
    """
    rng = rng or np.random.default_rng(0)
    times = np.asarray(times_mjd_seconds, dtype=float)
    az = np.broadcast_to(np.asarray(az_deg, dtype=float), times.shape)
    el = np.broadcast_to(np.asarray(el_deg, dtype=float), times.shape)
    freq = np.asarray(frequency_hz, dtype=float)
    n_ant = len(antenna_enu_m)
    p, q = np.triu_indices(n_ant, k=1)
    nbl = p.size
    a1, a2 = np.tile(p, times.size), np.tile(q, times.size)
    t = np.repeat(times, nbl)
    ecef = enu_to_ecef(antenna_enu_m, latitude_deg, longitude_deg, altitude_m)
    enu = ecef_to_enu(ecef, latitude_deg, longitude_deg)  # re-centred, as the loader does
    b = enu[a2] - enu[a1]
    k = wavenumber(freq)
    s = azel_to_unit(np.repeat(az, nbl), np.repeat(el, nbl))
    proj = np.sum(b * s, axis=1)
    data = amplitude * np.exp(-1j * proj[:, None] * k[None, :])
    for e_az, e_el, e_amp in extra_sources:
        data = data + e_amp * np.exp(-1j * (b @ azel_to_unit(e_az, e_el))[:, None] * k[None, :])
    if noise_sigma > 0:
        data = data + noise_sigma * (rng.standard_normal(data.shape) + 1j * rng.standard_normal(data.shape))
    flag = rng.uniform(size=data.shape) < flag_fraction
    if len(flagged_antennas):
        bad = np.isin(a1, flagged_antennas) | np.isin(a2, flagged_antennas)
        flag[bad] = True
    return SyntheticObservation(
        time_mjd_seconds=t, antenna1=a1, antenna2=a2, uvw_m=-b, data=data[:, :, None],
        flag=flag[:, :, None], weight=np.ones((t.size, 1)), frequency_hz=freq,
        antenna_positions_ecef_m=ecef, latitude_deg=latitude_deg, longitude_deg=longitude_deg,
        altitude_m=altitude_m, true_az_deg=np.asarray(az), true_el_deg=np.asarray(el),
        amplitude=complex(amplitude), noise_sigma=noise_sigma,
    )


def to_dataset(obs: SyntheticObservation, station_id: str = "synthetic") -> VisibilityDataset:
    """Build a :class:`VisibilityDataset` directly (no Measurement Set)."""
    return dataset_from_arrays(
        station_id=station_id, time_mjd_seconds=obs.time_mjd_seconds, antenna1=obs.antenna1,
        antenna2=obs.antenna2, uvw_m=obs.uvw_m, data=obs.data, flag=obs.flag, weight=obs.weight,
        frequency_hz=obs.frequency_hz, antenna_positions_ecef_m=obs.antenna_positions_ecef_m,
        latitude_deg=obs.latitude_deg, longitude_deg=obs.longitude_deg,
    )


def write_ms(obs: SyntheticObservation, path) -> Path:
    """Write a minimal CASA Measurement Set (requires python-casacore)."""
    from casacore.tables import default_ms, default_ms_subtable, makearrcoldesc, maketabdesc, table

    path = Path(path)
    n_chan = obs.frequency_hz.size
    desc = maketabdesc([makearrcoldesc("DATA", 0j, ndim=2, valuetype="complex", shape=[n_chan, 1])])
    with default_ms(str(path), desc) as t:
        t.addrows(obs.time_mjd_seconds.size)
        t.putcol("TIME", obs.time_mjd_seconds)
        t.putcol("TIME_CENTROID", obs.time_mjd_seconds)
        t.putcol("ANTENNA1", obs.antenna1)
        t.putcol("ANTENNA2", obs.antenna2)
        t.putcol("UVW", obs.uvw_m)
        t.putcol("DATA", obs.data)
        t.putcol("FLAG", obs.flag)
        t.putcol("WEIGHT", obs.weight)
        t.putcol("SIGMA", np.ones_like(obs.weight))
        t.putcol("INTERVAL", np.ones(obs.time_mjd_seconds.size))
        t.putcol("EXPOSURE", np.ones(obs.time_mjd_seconds.size))
    with table(str(path / "ANTENNA"), readonly=False, ack=False) as ant:
        n_ant = obs.antenna_positions_ecef_m.shape[0]
        ant.addrows(n_ant)
        ant.putcol("POSITION", obs.antenna_positions_ecef_m)
        ant.putcol("NAME", [f"ANTENNA-{i}" for i in range(n_ant)])
        ant.putcol("DISH_DIAMETER", np.full(n_ant, 0.1))
    with default_ms_subtable("SPECTRAL_WINDOW", str(path / "SPECTRAL_WINDOW")) as spw:
        spw.addrows(1)
        spw.putcell("NUM_CHAN", 0, n_chan)
        spw.putcell("CHAN_FREQ", 0, obs.frequency_hz)
        spw.putcell("CHAN_WIDTH", 0, np.full(n_chan, 2.5e6))
        spw.putcell("REF_FREQUENCY", 0, float(obs.frequency_hz[0]))
        spw.putcell("TOTAL_BANDWIDTH", 0, 2.5e6 * n_chan)
    with table(str(path / "OBSERVATION"), readonly=False, ack=False) as o:
        o.addrows(1)
        o.putcell("TELESCOPE_NAME", 0, "SYNTHETIC-TART")
    return path


def write_manifest(
    path,
    stations: Dict[str, Dict[str, object]],
    norad_id: int,
    tle_lines: Optional[Tuple[str, str]] = None,
    three_station_valid: bool = True,
    query_start_utc: str = "2026-06-25T09:14:02.460Z",
) -> Path:
    """Write a minimal Stage-2-style manifest for synthetic Measurement Sets.

    ``stations`` maps station_id to a dict with ``latitude_deg``,
    ``longitude_deg``, ``altitude_m``, ``ms_path`` (relative to the manifest
    directory), ``n_antennas``, ``first_time_utc``, ``last_time_utc`` and
    optional ``zero_gain_antennas``/``calibration_fingerprints``/``tart_id``.
    """
    path = Path(path)
    doc: Dict[str, object] = {
        "stage": 2,
        "manifest_version": "1.0",
        "three_station_valid": three_station_valid,
        "validity_reasons": [] if three_station_valid else ["synthetic: marked invalid"],
        "satellite": {
            "satellite_name": "SYNTHETIC", "norad_id": norad_id, "constellation": "GPS",
            "signal_relevance": "GPS_L1",
            "tle_line1": tle_lines[0] if tle_lines else None,
            "tle_line2": tle_lines[1] if tle_lines else None,
            "tle_epoch": None, "tle_source": "synthetic test input",
        },
        "prediction": {"query_start_utc": query_start_utc},
        "stations": {},
        "provenance": {"run_directory": "synthetic", "note": "SYNTHETIC DATA - not a TART observation"},
    }
    for sid, st in stations.items():
        doc["stations"][sid] = {
            "station_id": sid, "tart_id": st.get("tart_id", f"zz-{sid}"),
            "station_name": f"Synthetic station {sid}",
            "latitude_deg": st["latitude_deg"], "longitude_deg": st["longitude_deg"],
            "altitude_m": st["altitude_m"],
            "observation": {
                "coverage": {"zero_gain_antennas": list(st.get("zero_gain_antennas", [])),
                             "calibration_fingerprints": list(st.get("calibration_fingerprints", ["synthetic"]))},
                "calibration_status": "synthetic (no calibration)",
            },
            "measurement_set_status": st.get("measurement_set_status", "VALID"),
            "measurement_set": {
                "path": st["ms_path"], "valid": True, "station_id": sid,
                "number_of_antennas": st["n_antennas"], "frequency_hz": [GPS_L1_HZ],
                "first_time_utc": st["first_time_utc"], "last_time_utc": st["last_time_utc"],
                "calibration_applied_by_tart2ms": True,
            },
            "errors": [],
        }
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path


def utc_to_mjd_seconds(year: int, month: int, day: int, hour: int, minute: int, second: float) -> float:
    """UTC calendar time to MS ``TIME`` (MJD seconds)."""
    from datetime import datetime, timezone

    from tart_gnss.stage3.loader import MJD_UNIX_OFFSET_S

    whole = int(second)
    dt = datetime(year, month, day, hour, minute, whole, tzinfo=timezone.utc)
    return dt.timestamp() + (second - whole) + MJD_UNIX_OFFSET_S


def baseline_pairs(n_ant: int) -> List[Tuple[int, int]]:
    """All ``(p, q)`` antenna pairs with ``p < q``."""
    p, q = np.triu_indices(n_ant, k=1)
    return list(zip(p.tolist(), q.tolist()))
