#!/usr/bin/env python3
"""Inspect tab-sim Measurement Sets and validate saved-orbit baseline phases.

Inputs:
    A tab-sim MS; for the orbit subcommand, its used_orbits.json.

Outputs:
    MS component statistics or three near-field visibility phase scores.

Logic:
    Read the first integration, and optionally propagate the exact saved
    OMM via Skyfield at sub-integration times. Calculate antenna positions
    in GCRS with Astropy and tab-sim's corrected path-length phasors.
    This is a forward-model check, not orbit-blind localisation.
"""
from __future__ import annotations

import argparse
import json
import logging
from datetime import timezone

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from casacore.tables import table
from scipy.constants import c
from skyfield.api import EarthSatellite, load

LOGGER = logging.getLogger(__name__)


def read_first(ms_path: str) -> dict:
    """Read the first MS integration with geometry and frequency metadata.

    Args:
        ms_path: Measurement Set directory.

    Returns:
        A dictionary of visibility components and antenna metadata.
    """
    with table(ms_path, readonly=True, ack=False) as ms:
        times = ms.getcol("TIME")
        rows = np.flatnonzero(times == times.min())
        components = {}
        for name in ("DATA", "RFI_MODEL_DATA", "RFI_DATA", "NOISE_DATA"):
            if name in ms.colnames():
                components[name] = ms.getcol(name)[rows, 0, 0]
        result = dict(
            components=components, a1=ms.getcol("ANTENNA1")[rows],
            a2=ms.getcol("ANTENNA2")[rows], uvw=ms.getcol("UVW")[rows],
            flag=ms.getcol("FLAG")[rows, 0, 0],
            time=float(times.min()),
            interval=float(ms.getcol("INTERVAL")[rows[0]]),
            n_rows=int(ms.nrows()),
            n_times=int(len(np.unique(times))),
        )
    with table(f"{ms_path}/ANTENNA", readonly=True, ack=False) as ant:
        result["ants"] = np.asarray(ant.getcol("POSITION"), dtype=float)
    with table(f"{ms_path}/SPECTRAL_WINDOW", readonly=True, ack=False) as spw:
        result["freq"] = float(spw.getcol("CHAN_FREQ")[0, 0])
    return result


def inspect(ms_path: str) -> None:
    """Print MS dimensions and visibility component amplitudes.

    Args:
        ms_path: Input Measurement Set path.

    Returns:
        None; writes tabular diagnostic output to stdout.
    """
    d = read_first(ms_path)
    print(f"Rows={d['n_rows']} integrations={d['n_times']} "
          f"antennas={len(d['ants'])} frequency={d['freq']/1e6:.6f} MHz")
    for name, vis in d["components"].items():
        amp = np.abs(vis)
        print(f"{name:18s} mean={amp.mean():.7g} "
              f"max={amp.max():.7g} nonzero={np.mean(amp>0):.6f}")


def read_saved_omm(path: str) -> dict:
    """Parse tab-sim's column-oriented used_orbits.json.

    Args:
        path: Saved orbital-record JSON path.

    Returns:
        First orbit record in Skyfield OMM-compatible dictionary form.
    """
    with open(path, encoding="utf-8") as handle:
        obj = json.load(handle)
    record = {
        key: value["0"] for key, value in obj.items()
        if isinstance(value, dict) and "0" in value
    }
    for name, default in {
        "CLASSIFICATION_TYPE": "U", "EPHEMERIS_TYPE": "0",
        "ELEMENT_SET_NO": "0", "REV_AT_EPOCH": "0",
        "MEAN_MOTION_DOT": 0.0, "MEAN_MOTION_DDOT": 0.0,
    }.items():
        record.setdefault(name, default)
    if "EPOCH" not in record or "NORAD_CAT_ID" not in record:
        raise ValueError("No usable OMM record found in saved-orbits file")
    return record


def phase_metrics(observed: np.ndarray, model: np.ndarray) -> tuple[float,float]:
    """Calculate weighted coherence and circular RMS residual phase.

    Args:
        observed: One complex visibility per baseline.
        model: Integrated model phasor per baseline.

    Returns:
        Tuple of (coherence, phase RMS in degrees).
    """
    weights = np.abs(observed)
    residual = observed * np.conj(model)
    unit = residual/np.maximum(np.abs(residual), 1e-30)
    result = np.sum(weights*unit)/np.sum(weights)
    aligned = np.angle(unit*np.exp(-1j*np.angle(result)))
    rms = np.sqrt(np.average(aligned*aligned, weights=weights))
    return float(np.abs(result)), float(np.rad2deg(rms))


def orbit_phase(ms_path: str, orbit_path: str, samples: int) -> None:
    """Compare tab-sim near-field phase prescription with simulated visibilities.

    Args:
        ms_path: Simulated Measurement Set directory.
        orbit_path: used_orbits.json from the exact same simulation.
        samples: Integration subsamples (33 in the Ghana test).

    Returns:
        None; prints predicted phase-coherence statistics.
    """
    if samples <= 0:
        raise ValueError("samples must be positive")
    d = read_first(ms_path)
    if "RFI_MODEL_DATA" not in d["components"]:
        raise ValueError("RFI_MODEL_DATA not found")
    vis = d["components"]["RFI_MODEL_DATA"]
    good = (~d["flag"]) & np.isfinite(vis) & (np.abs(vis) > 0)
    vis, a1, a2, uvw = vis[good], d["a1"][good], d["a2"][good], d["uvw"][good]
    offsets = ((np.arange(samples)+0.5)/samples - 0.5)*d["interval"]
    times = Time(d["time"]/86400 + offsets/86400, format="mjd", scale="utc")
    ts = load.timescale()
    satellite = EarthSatellite.from_omm(ts, read_saved_omm(orbit_path))
    ecef = d["ants"]
    locations = EarthLocation.from_geocentric(
        ecef[:, 0]*u.m, ecef[:, 1]*u.m, ecef[:, 2]*u.m
    )
    wave = c/d["freq"]
    with_w, without_w = [], []
    for tm in times:
        sf_tm = ts.from_datetime(tm.to_datetime(timezone=timezone.utc))
        sat_gcrs = satellite.at(sf_tm).position.km*1000.0
        ant_pos, _ = locations.get_gcrs_posvel(tm)
        ant_gcrs = ant_pos.xyz.to_value(u.m).T
        dist = np.linalg.norm(ant_gcrs-sat_gcrs[None, :], axis=1)
        delta = dist[a1]-dist[a2]
        without_w.append(np.exp(-2j*np.pi*delta/wave))
        with_w.append(np.exp(-2j*np.pi*(delta+uvw[:, 2])/wave))
    corrected = np.mean(with_w, axis=0)
    raw = np.mean(without_w, axis=0)
    print(f"Baseline count: {len(vis)}, dt={d['interval']:.6f} s, "
          f"subsamples={samples}")
    for name, model in (
        ("Near field + w", corrected), ("Near field (no w)", raw),
        ("Conjugate model", np.conj(corrected))
    ):
        coherence, rms_deg = phase_metrics(vis, model)
        print(f"{name:22s} coherence={coherence:.8f} "
              f"phase_RMS_deg={rms_deg:.5f}")
    LOGGER.info("Forward-model phase check complete")


def main() -> None:
    """Dispatch MS inspection or saved-orbit phase-model validation."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_inspect = sub.add_parser("inspect")
    p_inspect.add_argument("--ms", required=True)
    p_orbit = sub.add_parser("phase")
    p_orbit.add_argument("--ms", required=True)
    p_orbit.add_argument("--orbits", required=True)
    p_orbit.add_argument("--samples", type=int, default=33)
    args = parser.parse_args()
    if args.command == "inspect":
        inspect(args.ms)
    else:
        orbit_phase(args.ms, args.orbits, args.samples)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
