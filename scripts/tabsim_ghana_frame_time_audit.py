#!/usr/bin/env python3
"""Diagnose Ghana satellite direction bias from frame and time conventions.

Inputs:
    1SAT Measurement Set, matching used_orbits.json, snapshot index,
    internal samples and nominal range.

Outputs:
    JSON report and terminal comparisons of orbital directions, phase
    coherence, noise-free fitted direction and error relative to each truth.

Logic:
    Propagate saved OMM through all fine integration times; calculate orbital
    topocentric directions in the estimator's GCRS/ENU basis and in Skyfield.
    Compare the MS timestamp with an integration-time mean direction.
    Construct the moving-satellite near-field phasor model, and compare it
    with stationary-direction phasors and RFI_MODEL_DATA. The orbit is used
    only for evaluation, never to initialise or constrain the fit.
"""
from __future__ import annotations

import argparse
import json
import logging
from datetime import timezone
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.time import Time
from scipy.constants import c
from skyfield.api import EarthSatellite, load, wgs84

from tabsim_gnss_fit import geometry, read_snapshot
from tabsim_ghana_uncertainty import make_forward_model, fit_one
from tabsim_gnss_diagnostics import read_saved_omm

LOGGER = logging.getLogger(__name__)


def angles_from_basis(unit: np.ndarray, basis: np.ndarray) -> tuple[float, float]:
    """Calculate local Az/El from GCRS vectors and ENU basis.

    Args:
        unit: Satellite line-of-sight unit vector in GCRS, shape (3,).
        basis: ENU basis vectors as rows in GCRS, shape (3, 3).

    Returns:
        (azimuth, elevation) in degrees, azimuth east of north.
    """
    east, north, up = basis @ unit
    return (
        float(np.rad2deg(np.arctan2(east, north)) % 360.0),
        float(np.rad2deg(np.arcsin(np.clip(up, -1, 1)))),
    )


def offset_arcsec(reference: tuple[float, float],
                  comparison: tuple[float, float]) -> dict:
    """Express comparison-minus-reference in local tangent-plane arcseconds.

    Args:
        reference: (Az, El) degrees for chosen origin.
        comparison: (Az, El) degrees for compared measurement.

    Returns:
        East, North and radial differences in arcseconds.
    """
    delta_az = (comparison[0] - reference[0] + 180.0) % 360.0 - 180.0
    east = 3600.0 * delta_az * np.cos(np.deg2rad(reference[1]))
    north = 3600.0 * (comparison[1] - reference[1])
    return dict(east=float(east), north=float(north),
                radial=float(np.hypot(east, north)))


def phase_agreement(observed: np.ndarray, model: np.ndarray) -> dict:
    """Return profiled-complex-amplitude model error and phase coherence.

    Args:
        observed: Complex measured visibilities.
        model: Complex model vector (arbitrary overall scale).

    Returns:
        Normalised profiled residual and complex coherence.
    """
    model_power = float(np.vdot(model, model).real)
    vis_power = float(np.vdot(observed, observed).real)
    amplitude = np.vdot(model, observed) / model_power
    residual = observed - amplitude * model
    return dict(
        J=float(np.vdot(residual, residual).real / vis_power),
        coherence=float(abs(np.vdot(model, observed)) /
                        np.sqrt(model_power * vis_power)),
    )


def run(args: argparse.Namespace) -> dict:
    """Evaluate the orbit-frame and time-reference bias diagnostics.

    Args:
        args: CLI namespace with ms, orbits, snapshot, samples, range,
            search, grid spacing, and output path.

    Returns:
        JSON-serialisable diagnostic summary.
    """
    if args.samples < 1 or args.range_km <= 0:
        raise ValueError("samples and range must be positive")
    d = read_snapshot(args.ms, args.snapshot, "RFI_MODEL_DATA")
    timestamp = Time(d["ms_time"]/86400.0, format="mjd", scale="utc")
    offsets = ((np.arange(args.samples)+0.5)/args.samples - 0.5)*d["dt"]
    fine = Time(d["ms_time"]/86400.0 + offsets/86400.0,
                format="mjd", scale="utc")
    ant_gcrs, centre_gcrs, enu_gcrs, location = geometry(
        d["antenna_ecef"], fine
    )
    # The estimator's topocentric reference site is the antenna centroid.
    site_lat = float(location.lat.deg)
    site_lon = float(location.lon.deg)
    site_height = float(location.height.to_value(u.m))
    ts = load.timescale()
    orbit = EarthSatellite.from_omm(ts, read_saved_omm(args.orbits))
    observer = wgs84.latlon(latitude_degrees=site_lat,
                            longitude_degrees=site_lon,
                            elevation_m=site_height)

    samples_xyz = []
    skyfield_angles = []
    basis_angles = []
    unit_vectors = []
    distances = []
    for index, t in enumerate(fine):
        sf_t = ts.from_datetime(t.to_datetime(timezone=timezone.utc))
        position = orbit.at(sf_t).position.km*1000.0
        samples_xyz.append(position)
        line_of_sight = position-centre_gcrs[index]
        distances.append(float(np.linalg.norm(line_of_sight)))
        unit = line_of_sight/np.linalg.norm(line_of_sight)
        unit_vectors.append(unit)
        basis_angles.append(angles_from_basis(unit, enu_gcrs[index]))
        alt, az, _ = (orbit-observer).at(sf_t).altaz()
        skyfield_angles.append((float(az.degrees), float(alt.degrees)))

    mid = ts.from_datetime(timestamp.to_datetime(timezone=timezone.utc))
    position_mid = orbit.at(mid).position.km*1000.0
    # Compute a genuinely midpoint GCRS/ENU basis, not a nearest fine sample.
    _, mid_centre, mid_enu, _ = geometry(
        d["antenna_ecef"], Time([timestamp])
    )
    mid_los = position_mid-mid_centre[0]
    mid_unit = mid_los/np.linalg.norm(mid_los)
    basis_mid = angles_from_basis(mid_unit, mid_enu[0])
    alt_mid, az_mid, dist_mid = (orbit-observer).at(mid).altaz()
    skyfield_mid = (float(az_mid.degrees), float(alt_mid.degrees))
    middle_index = args.samples//2

    # Mean local ENU components: accounts for basis rotation over integration.
    mean_enu_components = np.mean(
        [enu_gcrs[k] @ unit_vectors[k] for k in range(args.samples)],
        axis=0
    )
    mean_enu_components /= np.linalg.norm(mean_enu_components)
    mean_basis = (
        float(np.rad2deg(np.arctan2(mean_enu_components[0],
                                    mean_enu_components[1])) % 360),
        float(np.rad2deg(np.arcsin(np.clip(mean_enu_components[2], -1, 1)))),
    )

    # Moving satellite, same path-length and w correction as tab-sim.
    positions = np.stack(samples_xyz)
    ranges = np.linalg.norm(
        ant_gcrs - positions[:, None, :], axis=-1
    )
    path = ranges[:, d["a1"]] - ranges[:, d["a2"]]
    wavelength = c/d["freq"]
    model_orbit = np.exp(
        -2j*np.pi/wavelength*(path+d["uvw"][None, :, 2])
    ).mean(axis=0)

    # Reuse exactly the same orbit-blind estimator as Monte Carlo.
    predict, pointing, _, _ = make_forward_model(
        d, args.samples, args.range_km
    )
    measured = np.asarray(d["vis"], dtype=complex)
    fit = fit_one(measured, predict, pointing,
                  args.search_deg, args.coarse_step_deg)
    fitted_angles = (fit["az_deg"], fit["el_deg"])
    model_fixed = predict(*fitted_angles)
    model_basis_mid = predict(*basis_mid)
    model_mean = predict(*mean_basis)

    report = {
        "utc": timestamp.utc.isot,
        "samples": args.samples,
        "integration_seconds": d["dt"],
        "nominal_range_km": args.range_km,
        "orbit_topocentric_range_km": float(dist_mid.km),
        "skyfield_mid_deg": list(skyfield_mid),
        "gcrs_enu_mid_deg": list(basis_mid),
        "gcrs_enu_time_mean_deg": list(mean_basis),
        "skyfield_minus_basis_mid_arcsec":
            offset_arcsec(basis_mid, skyfield_mid),
        "mean_minus_instantaneous_arcsec":
            offset_arcsec(basis_mid, mean_basis),
        "orbit_drift_start_to_end_arcsec":
            offset_arcsec(basis_angles[0], basis_angles[-1]),
        "fit_angles_deg": list(fitted_angles),
        "fit_minus_skyfield_arcsec":
            offset_arcsec(skyfield_mid, fitted_angles),
        "fit_minus_gcrs_mid_arcsec":
            offset_arcsec(basis_mid, fitted_angles),
        "fit_minus_gcrs_mean_arcsec":
            offset_arcsec(mean_basis, fitted_angles),
        "fit": fit,
        "phase_model_scores": {
            "moving_orbit": phase_agreement(measured, model_orbit),
            "stationary_fit": phase_agreement(measured, model_fixed),
            "stationary_midpoint_truth": phase_agreement(
                measured, model_basis_mid),
            "stationary_time_mean_truth": phase_agreement(
                measured, model_mean),
        },
        "notes": [
            "The orbital direction is kept out of fit_one().",
            "The mean direction uses time-averaged local ENU components.",
            "Models ignore slight time-dependent antenna beam amplitudes.",
            "Current MS w is held constant over the integration.",
        ],
    }
    print("\n=== GHANA ORBIT FRAME AND TIME DIAGNOSTIC ===")
    for label, pair in (
        ("Skyfield at MS timestamp", skyfield_mid),
        ("GCRS/ENU at MS timestamp", basis_mid),
        ("GCRS/ENU mean over integration", mean_basis),
        ("Noise-free fitted direction", fitted_angles),
    ):
        print(f"{label:37s} Az={pair[0]:.9f} El={pair[1]:.9f} deg")
    for key in ("skyfield_minus_basis_mid_arcsec",
                "mean_minus_instantaneous_arcsec",
                "orbit_drift_start_to_end_arcsec",
                "fit_minus_skyfield_arcsec",
                "fit_minus_gcrs_mid_arcsec",
                "fit_minus_gcrs_mean_arcsec"):
        delta = report[key]
        print(f"{key:37s} E={delta['east']:+.4f}\" "
              f"N={delta['north']:+.4f}\" |d|={delta['radial']:.4f}\"")
    print("\n=== FORWARD-MODEL COMPARISON ===")
    for name, metric in report["phase_model_scores"].items():
        print(f"{name:27s} J={metric['J']:.9g} "
              f"coherence={metric['coherence']:.12f}")
    print(f"\nLocal optimiser success: {fit['local_success']}; "
          f"on boundary: {fit['on_edge']}")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    LOGGER.info("Wrote %s", out)
    return report


def main() -> None:
    """Parse CLI inputs and run the deterministic direction-bias test."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ms", required=True)
    parser.add_argument("--orbits", required=True)
    parser.add_argument("--output",
        default="results/tabsim_validation/ghana_audit/frame_time_bias.json")
    parser.add_argument("--snapshot", type=int, default=0)
    parser.add_argument("--samples", type=int, default=33)
    parser.add_argument("--range-km", type=float, default=21500.0)
    parser.add_argument("--search-deg", type=float, default=3.0)
    parser.add_argument("--coarse-step-deg", type=float, default=0.5)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
