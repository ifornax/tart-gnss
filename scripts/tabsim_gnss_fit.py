#!/usr/bin/env python3
"""Fit topocentric GNSS direction in a tab-sim Measurement Set.

Inputs:
    MS path, snapshot index, nominal satellite range, fine-sample count,
    search half-width, random seed, and output JSON path.

Outputs:
    Orbit-blind Az/El estimate, fitted complex amplitude, profile
    least-squares residual, phase coherence, and optimiser statuses.

Logic:
    Construct GCRS antenna and ENU coordinates at each integration
    subsample, predict finite-range path differences with tab-sim's
    antenna-w correction, analytically profile a complex amplitude,
    and optimise Az/El. The orbit is NOT read by this script.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from casacore.tables import table
from scipy.constants import c
from scipy.optimize import differential_evolution, minimize

LOGGER = logging.getLogger(__name__)


def gcrs_positions(loc: EarthLocation, times: Time) -> np.ndarray:
    """Return site GCRS positions (metres).

    Args:
        loc: Scalar or vector EarthLocation.
        times: Array of observation times.

    Returns:
        Array (n_time, n_site, 3) in GCRS metres.
    """
    positions = []
    for t in times:
        xyz, _ = loc.get_gcrs_posvel(t)
        positions.append(np.asarray(xyz.xyz.to_value(u.m)).T.reshape(-1, 3))
    return np.stack(positions)


def geometry(ant_ecef: np.ndarray, times: Time):
    """Build GCRS antenna positions and local ENU basis vectors.

    Args:
        ant_ecef: Antenna ITRF/ECEF positions, shape (n_ant,3), metres.
        times: Fine integration timestamps.

    Returns:
        Tuple (antenna_GCRS, centre_GCRS, ENU_basis_GCRS, centre_location).
    """
    antennas = EarthLocation.from_geocentric(
        ant_ecef[:, 0] * u.m, ant_ecef[:, 1] * u.m, ant_ecef[:, 2] * u.m
    )
    centre_ecef = ant_ecef.mean(axis=0)
    centre_loc = EarthLocation.from_geocentric(*(centre_ecef * u.m))
    lat = centre_loc.lat.to_value(u.rad)
    lon = centre_loc.lon.to_value(u.rad)
    east = [-np.sin(lon), np.cos(lon), 0.0]
    north = [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)]
    up = [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]
    enu_ecef = np.array([east, north, up])
    antenna_gcrs = gcrs_positions(antennas, times)
    centre_gcrs = gcrs_positions(centre_loc, times)[:, 0, :]
    # Translate probe points in ECEF to recover the local frame in GCRS.
    probes_ecef = centre_ecef + 100.0 * enu_ecef
    probes = EarthLocation.from_geocentric(
        probes_ecef[:, 0]*u.m,
        probes_ecef[:, 1]*u.m,
        probes_ecef[:, 2]*u.m,
    )
    enu_gcrs = gcrs_positions(probes, times) - centre_gcrs[:, None, :]
    enu_gcrs /= np.linalg.norm(enu_gcrs, axis=-1, keepdims=True)
    return antenna_gcrs, centre_gcrs, enu_gcrs, centre_loc


def read_snapshot(ms_path: str, snapshot: int, column: str) -> dict:
    """Load one tab-sim snapshot and its metadata without reading orbital truth.

    Args:
        ms_path: Measurement Set directory.
        snapshot: Zero-based integration index.
        column: Visibility column name, DATA or RFI_MODEL_DATA.

    Returns:
        Dictionary of valid baseline data, pointing, frequency and geometry.
    """
    with table(ms_path, readonly=True, ack=False) as ms:
        if column not in ms.colnames():
            raise ValueError(f"Column {column} not present in MS")
        times = np.asarray(ms.getcol("TIME"))
        unique = np.unique(times)
        if not 0 <= snapshot < len(unique):
            raise ValueError(f"snapshot must be between 0 and {len(unique)-1}")
        rows = np.flatnonzero(times == unique[snapshot])
        vis = ms.getcol(column)[rows, 0, 0]
        flags = ms.getcol("FLAG")[rows, 0, 0]
        uvw = ms.getcol("UVW")[rows]
        a1 = ms.getcol("ANTENNA1")[rows]
        a2 = ms.getcol("ANTENNA2")[rows]
        dt = float(ms.getcol("INTERVAL")[rows[0]])
    with table(f"{ms_path}/ANTENNA", readonly=True, ack=False) as ant:
        xyz = np.asarray(ant.getcol("POSITION"), dtype=np.float64)
    with table(f"{ms_path}/SPECTRAL_WINDOW", readonly=True, ack=False) as spw:
        freq = float(spw.getcol("CHAN_FREQ")[0, 0])
    with table(f"{ms_path}/FIELD", readonly=True, ack=False) as field:
        ra, dec = np.asarray(field.getcol("PHASE_DIR")[0]).reshape(-1)[:2]
    valid = (~flags) & np.isfinite(vis) & np.all(np.isfinite(uvw), axis=1)
    if not valid.any():
        raise ValueError("No usable cross-correlations")
    return dict(
        vis=vis[valid], uvw=uvw[valid], a1=a1[valid], a2=a2[valid],
        antenna_ecef=xyz, freq=freq, ra=ra, dec=dec, dt=dt,
        ms_time=float(unique[snapshot]),
    )


def fit(ms_path: str, output: str, snapshot: int = 0,
        nominal_range_km: float = 21500.0, samples: int = 33,
        search_deg: float = 15.0, seed: int = 123,
        column: str = "DATA") -> dict:
    """Fit direction through an orbit-blind spherical-wave visibility model.

    Args:
        ms_path: Path to tab-sim MS.
        output: JSON result path.
        snapshot: Zero-based timestamp index.
        nominal_range_km: Fixed site-to-satellite range, kilometres.
        samples: Number of fine samples across the integration.
        search_deg: Angular half-width about the MS pointing.
        seed: Differential-evolution random seed.
        column: DATA (recommended), or RFI_MODEL_DATA for a noise-free test.

    Returns:
        Result dictionary, also written to output as JSON.
    """
    if nominal_range_km <= 0 or samples <= 0 or search_deg <= 0:
        raise ValueError("Range, samples, and search half-width must be positive")
    d = read_snapshot(ms_path, snapshot, column)
    offsets = ((np.arange(samples) + 0.5)/samples - 0.5)*d["dt"]
    fine_times = Time(
        d["ms_time"]/86400 + offsets/86400, format="mjd", scale="utc"
    )
    LOGGER.info("Building GCRS geometry: %d internal samples", samples)
    ant, centre, enu, loc = geometry(d["antenna_ecef"], fine_times)
    mid = Time(d["ms_time"]/86400, format="mjd", scale="utc")
    pointing = SkyCoord(
        ra=d["ra"]*u.rad, dec=d["dec"]*u.rad, frame="icrs"
    ).transform_to(AltAz(obstime=mid, location=loc, pressure=0*u.hPa))
    centre_az, centre_el = float(pointing.az.deg), float(pointing.alt.deg)
    baseline_w = d["uvw"][:, 2]
    wave = c/d["freq"]
    vis = d["vis"]
    norm_v = float(np.vdot(vis, vis).real)
    if norm_v <= 0:
        raise ValueError("Input visibility energy is zero")

    def predict(az_deg: float, el_deg: float):
        """Return profile objective, optimal complex amplitude and model."""
        az, el = np.deg2rad([az_deg, el_deg])
        east = np.sin(az)*np.cos(el)
        north = np.cos(az)*np.cos(el)
        up = np.sin(el)
        unit_gcrs = (
            east*enu[:, 0, :] + north*enu[:, 1, :] + up*enu[:, 2, :]
        )
        sat = centre + nominal_range_km*1000*unit_gcrs
        dist = np.linalg.norm(ant - sat[:, None, :], axis=-1)
        path = dist[:, d["a1"]] - dist[:, d["a2"]]
        phase = (-2*np.pi/wave)*(path + baseline_w[None, :])
        model = np.mean(np.exp(1j*phase), axis=0)
        denom = float(np.vdot(model, model).real)
        amplitude = np.vdot(model, vis)/denom
        residual = vis - amplitude*model
        score = float(np.vdot(residual, residual).real/norm_v)
        return score, amplitude, model

    def objective(theta: np.ndarray) -> float:
        """Return the profile residual for an Az/El trial."""
        return predict(float(theta[0]), float(theta[1]))[0]

    bounds = [
        (centre_az-search_deg, centre_az+search_deg),
        (max(0.0, centre_el-search_deg), min(90.0, centre_el+search_deg)),
    ]
    if bounds[1][0] >= bounds[1][1]:
        raise ValueError("Elevation search does not intersect [0,90] degrees")
    LOGGER.info("Searching Az [%.3f, %.3f], El [%.3f, %.3f]",
                bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
    global_opt = differential_evolution(
        objective, bounds=bounds, seed=seed, maxiter=100,
        popsize=16, tol=1e-8, polish=False
    )
    local_opt = minimize(
        objective, global_opt.x, method="Nelder-Mead",
        options={"xatol":1e-8, "fatol":1e-12, "maxiter":250}
    )
    theta = local_opt.x if local_opt.fun <= global_opt.fun else global_opt.x
    score, amplitude, model = predict(*theta)
    coherence = float(
        abs(np.vdot(model, vis)) /
        np.sqrt(np.vdot(model, model).real * norm_v)
    )
    result = {
        "ms": str(Path(ms_path).resolve()),
        "snapshot_index": snapshot,
        "timestamp_utc": mid.utc.isot,
        "nominal_range_km": nominal_range_km,
        "samples_per_integration": samples,
        "frequency_hz": d["freq"],
        "n_baselines": int(len(vis)),
        "search_centre_az_deg": centre_az,
        "search_centre_el_deg": centre_el,
        "search_half_width_deg": search_deg,
        "best_az_deg": float(theta[0] % 360),
        "best_el_deg": float(theta[1]),
        "fitted_complex_amplitude": [
            float(amplitude.real), float(amplitude.imag)
        ],
        "normalised_squared_residual": score,
        "phase_coherence": coherence,
        "global_optimiser_success": bool(global_opt.success),
        "local_optimiser_success": bool(local_opt.success),
        "visibility_column": column,
    }
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    for key, value in result.items():
        if key != "ms":
            print(f"{key}: {value}")
    LOGGER.info("Wrote %s", output)
    return result


def main() -> None:
    """Parse command-line options and perform the orbit-blind fit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ms_path")
    parser.add_argument("--output", default="direction_fit.json")
    parser.add_argument("--snapshot", type=int, default=0)
    parser.add_argument("--range-km", type=float, default=21500.0)
    parser.add_argument("--samples", type=int, default=33)
    parser.add_argument("--search-deg", type=float, default=15.0)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--column", choices=["DATA", "RFI_MODEL_DATA"],
                        default="DATA")
    args = parser.parse_args()
    fit(args.ms_path, args.output, args.snapshot, args.range_km,
        args.samples, args.search_deg, args.seed, args.column)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
