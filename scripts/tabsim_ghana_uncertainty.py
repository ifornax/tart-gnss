#!/usr/bin/env python3
"""Monte Carlo calibration of a Ghana TART GNSS direction estimator.

Inputs:
    A tab-sim Ghana Measurement Set, matching used_orbits.json, nominal range,
    list of per-visibility complex-amplitude SNRs, repetitions and random seed.

Outputs:
    One trial CSV, one per-SNR summary CSV, a covariance JSON for each SNR,
    a plot of localisation error versus SNR, and a run-metadata JSON.

Logic:
    Scale RFI_MODEL_DATA to a prescribed single-visibility SNR, generate fresh
    independent circular complex-Gaussian noise, and recover Az/El without
    orbit information using a spherical-wave, w-corrected forward model.
    Consult the saved orbit only AFTER fitting, to score positional errors.
    Estimate tangent-plane covariance from a training half and assess
    68% ellipse coverage on the independent held-out half.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
from datetime import timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from casacore.tables import table
from scipy.constants import c
from scipy.optimize import minimize
from scipy.stats import chi2
from skyfield.api import EarthSatellite, load, wgs84

from tabsim_gnss_fit import geometry, read_snapshot
from tabsim_gnss_diagnostics import read_saved_omm

LOGGER = logging.getLogger(__name__)


def observation_noise_sigma(ms_path: str) -> float:
    """Estimate receiver noise per real/imaginary component in MS units.

    Args:
        ms_path: Simulation Measurement Set with a NOISE_DATA column.

    Returns:
        Scalar standard deviation estimated from all unflagged NOISE_DATA.

    Logic:
        Under circular zero-mean complex Gaussian noise,
        E[abs(noise)**2] = 2*sigma**2.
    """
    with table(ms_path, readonly=True, ack=False) as ms:
        if "NOISE_DATA" not in ms.colnames():
            raise ValueError("NOISE_DATA required for empirical noise scale")
        noise = np.asarray(ms.getcol("NOISE_DATA")[:, 0, 0])
        flags = np.asarray(ms.getcol("FLAG")[:, 0, 0])
    valid = (~flags) & np.isfinite(noise)
    if not np.any(valid):
        raise ValueError("No valid NOISE_DATA samples")
    sigma = float(np.sqrt(np.mean(np.abs(noise[valid]) ** 2) / 2))
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("Invalid NOISE_DATA variance")
    return sigma


def make_forward_model(d: dict, samples: int, range_km: float):
    """Construct a reusable orbit-blind complex-visibility forward model.

    Args:
        d: Snapshot dictionary from tabsim_gnss_fit.read_snapshot.
        samples: Number of internal samples across each integration.
        range_km: Assumed fixed topocentric range in kilometres.

    Returns:
        Tuple (predict, pointing Az/El, centre EarthLocation, snapshot Time).
        predict(az,el) returns a complex model vector per baseline.
    """
    offsets = ((np.arange(samples) + 0.5) / samples - 0.5) * d["dt"]
    fine = Time(d["ms_time"] / 86400 + offsets / 86400, format="mjd", scale="utc")
    antenna_gcrs, centre_gcrs, enu, site = geometry(d["antenna_ecef"], fine)
    timestamp = Time(d["ms_time"] / 86400, format="mjd", scale="utc")
    pointing = SkyCoord(
        ra=d["ra"]*u.rad, dec=d["dec"]*u.rad, frame="icrs"
    ).transform_to(AltAz(obstime=timestamp, location=site, pressure=0*u.hPa))
    centre_az = float(pointing.az.deg)
    centre_el = float(pointing.alt.deg)
    wavelength = c / d["freq"]
    a1, a2 = d["a1"], d["a2"]
    w = d["uvw"][:, 2]

    def predict(az_deg: float, el_deg: float) -> np.ndarray:
        """Predict integrated baseline phasors for a direction and fixed range.

        Args:
            az_deg: Topocentric azimuth in degrees east of north.
            el_deg: Topocentric elevation in degrees.

        Returns:
            One integrated unit-amplitude complex model per baseline.
        """
        az, el = np.deg2rad([az_deg, el_deg])
        direction = (
            np.sin(az)*np.cos(el)*enu[:, 0, :]
            + np.cos(az)*np.cos(el)*enu[:, 1, :]
            + np.sin(el)*enu[:, 2, :]
        )
        satellite = centre_gcrs + range_km * 1000.0 * direction
        distances = np.linalg.norm(
            satellite[:, None, :] - antenna_gcrs, axis=-1
        )
        delta = distances[:, a1] - distances[:, a2]
        phase = (-2 * np.pi / wavelength) * (delta + w[None, :])
        return np.exp(1j * phase).mean(axis=0)

    return predict, (centre_az, centre_el), site, timestamp


def fit_one(vis: np.ndarray, predict, pointing: tuple[float, float],
            search_deg: float, coarse_step_deg: float) -> dict:
    """Recover Az/El from a complex visibility vector without orbital truth.

    Args:
        vis: Noise-perturbed baseline visibilities.
        predict: Orbit-blind finite-range model callable.
        pointing: MS phase-centre Az/El in degrees.
        search_deg: Half-width of bounded pointing-informed search.
        coarse_step_deg: Spacing of deterministic coarse-search grid.

    Returns:
        Dictionary with fitted Az/El, profile objective and optimiser flags.

    Logic:
        Search a grid centred on the MS pointing, select the highest-coherence
        grid point, and locally refine it with a penalised Nelder-Mead fit.
    """
    az0, el0 = pointing
    az_bounds = (az0-search_deg, az0+search_deg)
    el_bounds = (max(0.0, el0-search_deg), min(90.0, el0+search_deg))
    power = float(np.vdot(vis, vis).real)

    def objective(x: np.ndarray) -> float:
        """Return complex-amplitude-profiled normalised least-square loss."""
        az, el = float(x[0]), float(x[1])
        if not (az_bounds[0] <= az <= az_bounds[1]
                and el_bounds[0] <= el <= el_bounds[1]):
            return 2.0
        model = predict(az, el)
        denominator = float(np.vdot(model, model).real)
        amplitude = np.vdot(model, vis) / denominator
        residual = vis - amplitude * model
        return float(np.vdot(residual, residual).real / power)

    az_grid = np.linspace(*az_bounds, num=max(3, int(
        np.ceil(2*search_deg / coarse_step_deg))+1))
    el_grid = np.linspace(*el_bounds, num=max(3, int(
        np.ceil((el_bounds[1]-el_bounds[0])/coarse_step_deg))+1))
    grid = [np.array([a, e]) for a in az_grid for e in el_grid]
    start = min(grid, key=objective)
    result = minimize(
        objective, start, method="Nelder-Mead",
        options={"xatol": 2e-7, "fatol": 1e-12, "maxiter": 300}
    )
    x = result.x if result.fun <= objective(start) else start
    az_hat, el_hat = map(float, x)
    score = objective(x)
    on_edge = (min(abs(az_hat-az_bounds[0]), abs(az_hat-az_bounds[1]),
                   abs(el_hat-el_bounds[0]), abs(el_hat-el_bounds[1])) < 0.002)
    return {
        "az_deg": az_hat % 360.0,
        "el_deg": el_hat,
        "residual": score,
        "local_success": bool(result.success),
        "on_edge": bool(on_edge),
    }


def orbit_truth(omm_path: str, site, timestamp: Time) -> tuple[float, float]:
    """Compute independent topocentric orbit truth for the antenna-centroid site.

    Args:
        omm_path: Matching saved used_orbits.json file.
        site: Array-centroid EarthLocation from fitted MS antenna coordinates.
        timestamp: Snapshot UTC time.

    Returns:
        Azimuth and elevation in degrees. Not used by the fit objective.
    """
    ts = load.timescale()
    satellite = EarthSatellite.from_omm(ts, read_saved_omm(omm_path))
    observer = wgs84.latlon(
        latitude_degrees=float(site.lat.deg),
        longitude_degrees=float(site.lon.deg),
        elevation_m=float(site.height.to_value(u.m)),
    )
    moment = ts.from_datetime(timestamp.to_datetime(timezone=timezone.utc))
    altitude, azimuth, _ = (satellite - observer).at(moment).altaz()
    return float(azimuth.degrees), float(altitude.degrees)


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a table of dictionaries to CSV.

    Args:
        path: Destination CSV.
        rows: Nonempty list of homogeneously keyed dictionaries.

    Returns:
        None.
    """
    if not rows:
        raise ValueError("No rows to write")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyse_level(rows: list[dict], true_az: float, true_el: float) -> dict:
    """Estimate a tangent-plane covariance and held-out 68% ellipse coverage.

    Args:
        rows: Ordered trial records for one SNR level.
        true_az: Independent orbit azimuth in degrees.
        true_el: Independent orbit elevation in degrees.

    Returns:
        Summary with errors, sample covariance, bias and held-out coverage.

    Logic:
        First half of trials estimates covariance and bias. Remaining
        independent trials assess nominal 2-D Gaussian ellipse coverage.
        Both raw and bias-corrected coverage are reported.
    """
    errors = np.array([[r["east_arcsec"], r["north_arcsec"]] for r in rows])
    half = len(errors)//2
    if half < 3 or len(errors)-half < 3:
        raise ValueError("Require at least 6 trials per SNR")
    train, test = errors[:half], errors[half:]
    bias = train.mean(axis=0)
    covariance = np.cov(train, rowvar=False, ddof=1)
    eigenvalues = np.linalg.eigvalsh(covariance)
    if eigenvalues[0] <= 0:
        LOGGER.warning("Covariance is singular; coverage not defined")
        raw_cov = float("nan")
        corrected_cov = float("nan")
    else:
        threshold = float(chi2.ppf(0.68, df=2))
        precision = np.linalg.inv(covariance)
        m2_raw = np.einsum("ni,ij,nj->n", test, precision, test)
        centred = test-bias
        m2_corr = np.einsum("ni,ij,nj->n", centred, precision, centred)
        raw_cov = float(np.mean(m2_raw <= threshold))
        corrected_cov = float(np.mean(m2_corr <= threshold))
    norms = np.linalg.norm(errors, axis=1)
    return {
        "snr": rows[0]["snr"],
        "n_trials": len(rows),
        "fit_success_fraction": float(np.mean([r["local_success"] for r in rows])),
        "edge_fraction": float(np.mean([r["on_edge"] for r in rows])),
        "median_error_arcsec": float(np.median(norms)),
        "p68_error_arcsec": float(np.quantile(norms, 0.68)),
        "p95_error_arcsec": float(np.quantile(norms, 0.95)),
        "rmse_error_arcsec": float(np.sqrt(np.mean(norms**2))),
        "bias_east_arcsec": float(bias[0]),
        "bias_north_arcsec": float(bias[1]),
        "covariance_east_north_arcsec2": covariance.tolist(),
        "heldout_coverage68_raw": raw_cov,
        "heldout_coverage68_bias_corrected": corrected_cov,
        "ellipse68_chi2_threshold": float(chi2.ppf(0.68, df=2)),
        "truth_az_deg": true_az,
        "truth_el_deg": true_el,
    }


def draw_plot(path: Path, summaries: list[dict]) -> None:
    """Save localisation-error quantiles as a function of input SNR.

    Args:
        path: Output PNG filename.
        summaries: Per-SNR dictionaries from analyse_level.

    Returns:
        None; writes a PNG.
    """
    ordered = sorted(summaries, key=lambda r: r["snr"])
    xs = np.array([r["snr"] for r in ordered])
    median = np.array([r["median_error_arcsec"] for r in ordered])
    p68 = np.array([r["p68_error_arcsec"] for r in ordered])
    p95 = np.array([r["p95_error_arcsec"] for r in ordered])
    plt.figure(figsize=(7.5, 4.8))
    plt.plot(xs, median, marker="o", label="Median angular error")
    plt.plot(xs, p68, marker="s", label="68th percentile error")
    plt.plot(xs, p95, marker="^", label="95th percentile error")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Per-visibility complex-amplitude SNR")
    plt.ylabel("Position error relative to saved orbit (arcsec)")
    plt.title("Ghana TART: direction recovery versus injected SNR")
    plt.grid(alpha=0.25, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(args: argparse.Namespace) -> None:
    """Run independent noisy fits, score against orbit, and write deliverables.

    Args:
        args: Parsed CLI options (MS, OMM, SNRs, repetitions, model and paths).

    Returns:
        None. Saves plot, CSVs, covariance JSON and run-metadata JSON.
    """
    if args.reps < 6 or args.samples < 1 or args.range_km <= 0:
        raise ValueError("Require reps>=6, samples>=1 and range-km>0")
    if args.search_deg <= 0 or args.coarse_step_deg <= 0:
        raise ValueError("Search parameters must be positive")
    snrs = [float(v) for v in args.snrs.split(",")]
    if any(not np.isfinite(s) or s <= 0 for s in snrs):
        raise ValueError("SNRs must be positive finite numbers")

    d = read_snapshot(args.ms, args.snapshot, "RFI_MODEL_DATA")
    sigma = observation_noise_sigma(args.ms)
    signal = np.asarray(d["vis"], dtype=np.complex128)
    signal_rms = np.sqrt(np.mean(np.abs(signal)**2))
    predict, pointing, location, timestamp = make_forward_model(
        d, args.samples, args.range_km
    )
    LOGGER.info("Reference RMS satellite=%.6g, receiver sigma=%.6g",
                signal_rms, sigma)
    LOGGER.info("Pointing-informed search centre: az=%.5f el=%.5f",
                *pointing)

    rng = np.random.default_rng(args.seed)
    trial_rows = []
    grouping = []
    for snr in snrs:
        # E|noise|^2 = 2*sigma^2. Thus input per-visibility SNR
        # is RMS(signal) / RMS(complex noise).
        scale = snr*np.sqrt(2)*sigma/signal_rms
        level_rows = []
        for index in range(args.reps):
            noise = sigma*(rng.standard_normal(len(signal))
                           + 1j*rng.standard_normal(len(signal)))
            observation = signal*scale + noise
            fit = fit_one(
                observation, predict, pointing, args.search_deg,
                args.coarse_step_deg
            )
            fit.update({"snr": snr, "trial": index})
            level_rows.append(fit)
            if (index+1) % max(1, args.reps//5) == 0:
                LOGGER.info("SNR %g: %d/%d fits", snr, index+1, args.reps)
        grouping.append(level_rows)

    # Deliberately consult orbital truth ONLY after all fits complete.
    true_az, true_el = orbit_truth(args.orbits, location, timestamp)
    LOGGER.info("Independent truth: az=%.7f el=%.7f", true_az, true_el)
    summaries = []
    for rows in grouping:
        for r in rows:
            delta_az = (r["az_deg"]-true_az+180) % 360-180
            r["east_arcsec"] = delta_az*np.cos(np.deg2rad(true_el))*3600
            r["north_arcsec"] = (r["el_deg"]-true_el)*3600
            r["error_arcsec"] = float(np.hypot(
                r["east_arcsec"], r["north_arcsec"]
            ))
        summaries.append(analyse_level(rows, true_az, true_el))
        trial_rows.extend(rows)

    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    write_csv(root/"ghana_direction_trials.csv", trial_rows)
    scalars = [
        {k: v for k, v in s.items() if k != "covariance_east_north_arcsec2"}
        for s in summaries
    ]
    write_csv(root/"ghana_direction_snr_summary.csv", scalars)
    payload = {
        "frame": "local tangent plane at true direction: East, North",
        "units": "arcsec^2",
        "interpretation": (
            "Covariance of repeated Az/El estimation errors at fixed true "
            "direction, centred on TRAINING sample mean. This is NOT a "
            "calibrated single-observation posterior covariance."
        ),
        "snr_results": summaries,
    }
    (root/"ghana_direction_covariance.json").write_text(
        json.dumps(payload, indent=2)+"\n", encoding="utf-8"
    )
    metadata = {
        "ms": str(Path(args.ms).resolve()),
        "omm": str(Path(args.orbits).resolve()),
        "snapshot": args.snapshot,
        "timestamp_utc": timestamp.utc.isot,
        "seed": args.seed,
        "n_trials_per_snr": args.reps,
        "snrs": snrs,
        "sigma_real_imag": sigma,
        "snr_definition": (
            "RMS complex visibility signal / RMS complex receiver noise"
        ),
        "nominal_range_km": args.range_km,
        "samples": args.samples,
        "search_half_width_deg": args.search_deg,
        "coarse_step_deg": args.coarse_step_deg,
        "note": "Independent fits before orbit-truth evaluation; idealised simulation",
    }
    (root/"run_metadata.json").write_text(
        json.dumps(metadata, indent=2)+"\n", encoding="utf-8"
    )
    draw_plot(root/"ghana_localisation_error_vs_snr.png", summaries)
    for s in summaries:
        LOGGER.info(
            "SNR=%g median=%.3f arcsec p68=%.3f p95=%.3f "
            "heldout68(raw)=%.3f",
            s["snr"], s["median_error_arcsec"],
            s["p68_error_arcsec"], s["p95_error_arcsec"],
            s["heldout_coverage68_raw"],
        )
    LOGGER.info("Results written to %s", root)


def main() -> None:
    """Parse arguments and execute a reproducible Ghana uncertainty experiment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ms", required=True)
    parser.add_argument("--orbits", required=True)
    parser.add_argument("--snrs", default="3,5,10,20,50,100")
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--snapshot", type=int, default=0)
    parser.add_argument("--range-km", type=float, default=21500.0)
    parser.add_argument("--samples", type=int, default=33)
    parser.add_argument("--search-deg", type=float, default=3.0)
    parser.add_argument("--coarse-step-deg", type=float, default=1.0)
    parser.add_argument(
        "--output-dir", default="results/tabsim_validation/ghana_uncertainty"
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s"
    )
    main()
