"""``tart-localise``: Stage-3 command-line interface and orchestration.

Each station is processed independently (stations are never combined into one
interferometer, and their timestamps are never forced onto a common grid).

Exit status:

* ``0``: at least one usable scientific direction estimate was produced;
* ``3``: the run completed but produced no usable estimate;
* ``2``: invalid input or configuration (manifest, station, missing MS, ...);
* ``1``: an unexpected execution failure.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from tart_gnss import __version__
from tart_gnss.stage3 import loader as loader_mod
from tart_gnss.stage3 import output as output_mod
from tart_gnss.stage3 import quality as q
from tart_gnss.stage3.errors import (
    InsufficientBaselinesError,
    InvalidManifestError,
    InvalidVisibilityDataError,
    MissingMeasurementSetError,
    Stage3Error,
)
from tart_gnss.stage3.estimator import EstimatorOptions, SearchRegion, estimate_direction, resolution_rad
from tart_gnss.stage3.geometry import SPEED_OF_LIGHT_M_S, angular_separation_deg, azel_to_unit
from tart_gnss.stage3.imaging import image_direction
from tart_gnss.stage3.uncertainty import bootstrap_direction, estimate_direction_uncertainty
from tart_gnss.tle import format_utc

logger = logging.getLogger("tart_gnss.stage3")

EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_NO_ESTIMATES = 0, 1, 2, 3
METHODS = ("visibility", "image", "both")
#: Smallest topocentric range of any GNSS satellite (GLONASS altitude), used when no TLE range is available.
MIN_GNSS_RANGE_KM = 19100.0
#: Far-field phase error (rad) below which the plane-wave model is accepted.
NEAR_FIELD_PHASE_TOLERANCE_RAD = 0.01


def build_parser() -> argparse.ArgumentParser:
    """Argument parser for ``tart-localise``."""
    p = argparse.ArgumentParser(
        prog="tart-localise",
        description="Stage 3: estimate GNSS satellite directions from TART visibilities, per station.",
    )
    p.add_argument("--manifest", required=True, type=Path, help="Stage-2 manifest.json")
    sel = p.add_mutually_exclusive_group(required=True)
    sel.add_argument("--station", action="append", help="station_id or tart_id (repeatable)")
    sel.add_argument("--all-stations", action="store_true", help="process every station independently")
    p.add_argument("--output-dir", type=Path, default=Path("results/stage3"))
    p.add_argument("--method", choices=METHODS, default="visibility",
                   help="primary estimator (default: visibility); 'both' = visibility + image")
    p.add_argument("--image-validation", action="store_true",
                   help="also run the image-domain estimator as a cross-check")
    p.add_argument("--time-bin", type=float, default=1.0, help="time-bin width in seconds (default 1)")
    p.add_argument("--min-baselines", type=int, default=10, help="minimum valid baselines per bin (default 10)")
    p.add_argument("--search-radius-deg", type=float, default=10.0,
                   help="search-cap radius about the TLE prediction (default 10)")
    p.add_argument("--blind", action="store_true", help="search the whole visible hemisphere")
    p.add_argument("--robust-loss", choices=("linear", "huber", "soft_l1"), default="huber")
    p.add_argument("--bootstrap", type=int, default=0, help="baseline-bootstrap resamples per bin (default 0)")
    p.add_argument("--seed", type=int, default=20260625, help="random seed for the bootstrap")
    p.add_argument("--min-snr", type=float, default=5.0)
    p.add_argument("--max-flagged-fraction", type=float, default=0.5)
    p.add_argument("--allow-invalid-stage2", action="store_true",
                   help="process individually valid stations of a Stage-2 run that is not three-station valid")
    p.add_argument("--data-root", type=Path, default=None,
                   help="directory the manifest's relative MS paths are relative to (Stage-2 working directory)")
    p.add_argument("--no-plots", action="store_true", help="skip the diagnostic figures")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    return p


@dataclass
class TleReference:
    """TLE-predicted topocentric geometry for one station (reference only)."""

    line1: str
    line2: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float

    def predict(self, unix_s: np.ndarray):
        """Azimuth, elevation (deg) and range (km) at Unix times."""
        from skyfield.api import EarthSatellite, wgs84

        from tart_gnss.propagation import get_timescale, topocentric

        ts = get_timescale()
        times = ts.from_datetimes([datetime.fromtimestamp(float(u), tz=timezone.utc) for u in np.atleast_1d(unix_s)])
        sat = EarthSatellite(self.line1, self.line2, "target", ts)
        site = wgs84.latlon(self.latitude_deg, self.longitude_deg, elevation_m=self.altitude_m)
        az, el, rng, _ = topocentric(sat, site, times)
        return np.atleast_1d(az), np.atleast_1d(el), np.atleast_1d(rng)


@dataclass
class StationResult:
    """Rows and summary for one processed station."""

    rows: List[Dict[str, object]]
    summary: Dict[str, object]
    diagnostics: List[str] = field(default_factory=list)


def _nan() -> float:
    return float("nan")


def _base_row(si: loader_mod.StationInput, norad: int, t_mjd_s: float, method: str) -> Dict[str, object]:
    unix = float(loader_mod.mjd_seconds_to_unix(t_mjd_s))
    row: Dict[str, object] = {c: _nan() for c in output_mod.COLUMNS}
    row.update({
        "time_utc": format_utc(datetime.fromtimestamp(unix, tz=timezone.utc), 6),
        "time_mjd": t_mjd_s / 86400.0,
        "station_id": si.station_id,
        "tart_id": si.tart_id,
        "norad_id": norad,
        "primary_method": method,
        "fit_success": False,
        "optimizer_status": "not run",
        "usable": False,
        "quality_flags": [],
        "stage2_three_station_valid": si.stage2_three_station_valid,
    })
    return row


def process_station(
    si: loader_mod.StationInput,
    dataset: loader_mod.VisibilityDataset,
    *,
    norad_id: int,
    tle: Optional[TleReference],
    method: str,
    image_validation: bool,
    bin_s: float,
    options: EstimatorOptions,
    thresholds: q.QualityThresholds,
    n_boot: int,
    rng: np.random.Generator,
) -> StationResult:
    """Estimate a direction in every non-empty time bin of one station."""
    run_image = image_validation or method in ("image", "both")
    primary = "image" if method == "image" else "visibility"
    calib_warn = bool(si.calibration_warnings)
    bins = list(loader_mod.time_bins(dataset, bin_s))
    t_unix_all = loader_mod.mjd_seconds_to_unix(dataset.time_mjd_seconds)
    logger.info("%s: %d rows, %s .. %s, %.6g-%.6g MHz, %d/%d antennas valid, %d time bins of %.3g s",
                si.station_id, dataset.n_rows,
                format_utc(datetime.fromtimestamp(float(t_unix_all.min()), tz=timezone.utc)),
                format_utc(datetime.fromtimestamp(float(t_unix_all.max()), tz=timezone.utc)),
                dataset.frequency_hz.min() / 1e6, dataset.frequency_hz.max() / 1e6,
                int((~dataset.antenna_flagged).sum()) - len(set(si.zero_gain_antennas)),
                dataset.antenna_enu_m.shape[0], len(bins), bin_s)
    selections = [loader_mod.select_samples(dataset, b.rows, si.zero_gain_antennas) for b in bins]
    mids = np.array([
        float(np.average(s.samples.time_mjd_seconds, weights=s.samples.weight)) if s.samples.n
        else float(np.mean(dataset.time_mjd_seconds[b.rows]))
        for b, s in zip(bins, selections)
    ])
    if tle is not None:
        az_t, el_t, rng_t = tle.predict(loader_mod.mjd_seconds_to_unix(mids))
    else:
        az_t = el_t = rng_t = np.full(len(bins), np.nan)
    n_valid_bl = max((s.n_baselines_used for s in selections), default=0)
    logger.info("%s: up to %d valid baselines per bin (of %d)", si.station_id, n_valid_bl,
                max((s.n_baselines_total for s in selections), default=0))

    rows = []
    for k, (b, s) in enumerate(zip(bins, selections)):
        row = _base_row(si, norad_id, mids[k], primary)
        row.update(n_times=s.n_times, n_samples=s.samples.n, n_baselines_total=s.n_baselines_total,
                   n_baselines_used=s.n_baselines_used, flagged_fraction=s.flagged_fraction,
                   az_tle_deg=float(az_t[k]), el_tle_deg=float(el_t[k]), range_tle_km=float(rng_t[k]))
        if options.blind or not np.isfinite(az_t[k]):
            region = SearchRegion.hemisphere()
        else:
            region = SearchRegion.around_azel(az_t[k], el_t[k], options.search_radius_deg)
        if s.n_baselines_used < thresholds.min_baselines or s.samples.n < 3:
            row["quality_flags"] = q.assess(thresholds=thresholds, n_baselines_used=s.n_baselines_used,
                                            flagged_fraction=s.flagged_fraction, fitted=False,
                                            calibration_warning=calib_warn)
            row["optimizer_status"] = "not run: too few baselines"
            rows.append(row)
            continue
        samples = s.samples
        res_deg = float(np.degrees(resolution_rad(samples)))
        img = image_direction(samples, region) if run_image else None
        fitted = True
        cov_finite = True
        zen = False
        sigma_pos = None
        if primary == "visibility":
            try:
                est = estimate_direction(samples, None if options.blind else az_t[k],
                                         None if options.blind else el_t[k], options=options,
                                         station_id=si.station_id, n_baselines_used=s.n_baselines_used)
            except InsufficientBaselinesError as exc:
                row["optimizer_status"] = f"not run: {exc}"
                row["quality_flags"] = q.assess(thresholds=thresholds, n_baselines_used=0,
                                                flagged_fraction=s.flagged_fraction, fitted=False,
                                                calibration_warning=calib_warn)
                rows.append(row)
                continue
            cov = estimate_direction_uncertainty(est, samples)
            est.covariance = cov
            cov_finite = bool(np.isfinite(cov.sigma_l) and np.isfinite(cov.sigma_m))
            zen = cov.zenith_az_uncertain
            if cov_finite:
                el_r = np.radians(est.el_deg)
                s_az = cov.sigma_az_deg if np.isfinite(cov.sigma_az_deg) else 0.0
                sigma_pos = float(np.hypot(s_az * np.cos(el_r), cov.sigma_el_deg))
            row.update(
                l=est.l, m=est.m, sigma_l=cov.sigma_l, sigma_m=cov.sigma_m, cov_lm=cov.cov_lm,
                az_obs_deg=est.az_deg, el_obs_deg=est.el_deg, sigma_az_deg=cov.sigma_az_deg,
                sigma_el_deg=cov.sigma_el_deg, cov_az_el=cov.cov_az_el,
                complex_amp_real=est.complex_amplitude.real, complex_amp_imag=est.complex_amplitude.imag,
                apparent_amplitude=abs(est.complex_amplitude), chi2=est.chi2, reduced_chi2=est.reduced_chi2,
                rms_complex_residual=est.rms_complex_residual, median_abs_residual=est.median_abs_residual,
                max_abs_residual=est.max_abs_residual, sigma_vis=est.sigma_vis,
                residual_tail_ratio=est.residual_tail_ratio, snr=est.snr, fit_success=est.success,
                optimizer_status=est.optimizer_status,
            )
            boundary, invalid_geom, success = est.boundary, est.invalid_geometry, est.success
            snr, tail, rchi2 = est.snr, est.residual_tail_ratio, est.reduced_chi2
            s_vis = est.unit_vector
            if n_boot > 0:
                boot = bootstrap_direction(est, samples, n_boot, rng)
                row.update(az_boot_p16=boot.az_p16, az_boot_p50=boot.az_p50, az_boot_p84=boot.az_p84,
                           el_boot_p16=boot.el_p16, el_boot_p50=boot.el_p50, el_boot_p84=boot.el_p84)
        else:
            row.update(l=img.l, m=img.m, az_obs_deg=img.az_deg, el_obs_deg=img.el_deg, snr=img.image_snr,
                       fit_success=img.success, optimizer_status="image peak (no covariance)")
            boundary, invalid_geom, success = img.boundary, img.el_deg <= options.min_elevation_deg, img.success
            snr, tail, rchi2 = img.image_snr, _nan(), _nan()
            s_vis = None
        sep_img = None
        if img is not None:
            row.update(az_image_deg=img.az_deg, el_image_deg=img.el_deg, image_peak_value=img.peak_value,
                       image_snr=img.image_snr)
            if s_vis is not None:
                sep_img = float(angular_separation_deg(s_vis, azel_to_unit(img.az_deg, img.el_deg)))
                row["vis_image_separation_deg"] = sep_img
        if np.isfinite(az_t[k]) and np.isfinite(row["az_obs_deg"]):
            row["tle_separation_deg"] = float(angular_separation_deg(
                azel_to_unit(row["az_obs_deg"], row["el_obs_deg"]), azel_to_unit(az_t[k], el_t[k])))
        flags = q.assess(
            thresholds=thresholds, n_baselines_used=s.n_baselines_used, flagged_fraction=s.flagged_fraction,
            fitted=fitted, fit_success=success, boundary=boundary, invalid_geometry=bool(invalid_geom),
            snr=snr, residual_tail_ratio=tail, reduced_chi2=rchi2, covariance_finite=cov_finite,
            zenith_az_uncertain=zen, vis_image_separation_deg=sep_img, resolution_deg=res_deg,
            sigma_position_deg=sigma_pos, calibration_warning=calib_warn,
        )
        row["quality_flags"] = flags
        row["usable"] = q.is_usable(flags)
        rows.append(row)

    summary = _station_summary(si, dataset, rows, rng_t, n_valid_bl)
    logger.info("%s: fit success %.0f%%, usable %d/%d, median SNR %.3g, median vis-image offset %.3g deg",
                si.station_id, 100 * summary["fit_success_fraction"], summary["n_usable"], len(rows),
                summary["median_snr"] if summary["median_snr"] is not None else float("nan"),
                summary["median_vis_image_separation_deg"]
                if summary["median_vis_image_separation_deg"] is not None else float("nan"))
    return StationResult(rows, summary)


def near_field_diagnostic(max_baseline_m: float, frequency_hz: float, range_km: float,
                          range_source: str) -> Dict[str, object]:
    """Fresnel far-field metric ``eta = B^2 / (lambda R)`` and the implied phase error.

    The plane-wave approximation neglects a curvature phase of at most
    ``pi B^2 / (lambda R) = pi eta`` across the array.
    """
    lam = SPEED_OF_LIGHT_M_S / frequency_hz
    eta = max_baseline_m ** 2 / (lam * range_km * 1e3)
    err = float(np.pi * eta)
    return {
        "max_baseline_m": max_baseline_m,
        "wavelength_m": lam,
        "range_km": range_km,
        "range_source": range_source,
        "eta_B2_over_lambda_R": eta,
        "max_curvature_phase_rad": err,
        "tolerance_rad": NEAR_FIELD_PHASE_TOLERANCE_RAD,
        "far_field_valid": bool(err < NEAR_FIELD_PHASE_TOLERANCE_RAD),
        "model": "plane wave (far field)",
    }


def _station_summary(si, dataset, rows, rng_t, n_valid_bl) -> Dict[str, object]:
    fitted = [r for r in rows if not str(r["optimizer_status"]).startswith("not run")]
    succ = [r for r in fitted if r["fit_success"]]
    usable = [r for r in rows if r["usable"]]
    snrs = [r["snr"] for r in fitted if np.isfinite(r["snr"])]
    seps = [r["vis_image_separation_deg"] for r in fitted if np.isfinite(r["vis_image_separation_deg"])]
    cross = dataset.antenna1 != dataset.antenna2
    max_b = float(np.linalg.norm(dataset.baseline_enu_m[cross], axis=1).max()) if cross.any() else 0.0
    finite_rng = np.asarray(rng_t, dtype=float)
    finite_rng = finite_rng[np.isfinite(finite_rng)]
    if finite_rng.size:
        nf = near_field_diagnostic(max_b, float(dataset.frequency_hz.max()), float(finite_rng.min()),
                                   "minimum TLE-predicted range during the observation")
    else:
        nf = near_field_diagnostic(max_b, float(dataset.frequency_hz.max()), MIN_GNSS_RANGE_KM,
                                   "lower bound for GNSS satellites (no TLE range)")
    flag_counts: Dict[str, int] = {}
    for r in rows:
        for f in r["quality_flags"]:
            flag_counts[f] = flag_counts.get(f, 0) + 1
    t_unix = loader_mod.mjd_seconds_to_unix(dataset.time_mjd_seconds)
    return {
        **si.to_dict(),
        "n_rows": dataset.n_rows,
        "first_time_utc_ms": format_utc(datetime.fromtimestamp(float(t_unix.min()), tz=timezone.utc), 6),
        "last_time_utc_ms": format_utc(datetime.fromtimestamp(float(t_unix.max()), tz=timezone.utc), 6),
        "frequency_hz_range": [float(dataset.frequency_hz.min()), float(dataset.frequency_hz.max())],
        "n_channels": int(dataset.frequency_hz.size),
        "n_antennas_ms": int(dataset.antenna_enu_m.shape[0]),
        "max_valid_baselines_per_bin": int(n_valid_bl),
        "uvw_check": dataset.uvw_check,
        "near_field": nf,
        "n_time_bins": len(rows),
        "n_fitted": len(fitted),
        "fit_success_fraction": len(succ) / len(fitted) if fitted else 0.0,
        "n_usable": len(usable),
        "median_snr": float(np.median(snrs)) if snrs else None,
        "median_vis_image_separation_deg": float(np.median(seps)) if seps else None,
        "quality_flag_counts": flag_counts,
    }


def _git_state() -> Dict[str, object]:
    here = Path(__file__).resolve().parent
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=here, capture_output=True, text=True,
                                check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=here,
                                    capture_output=True, text=True, check=True).stdout.strip())
        return {"git_commit": commit, "git_dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": None, "git_dirty": None}


def _software() -> Dict[str, object]:
    from tart_gnss.provenance import _version, software_versions

    out = software_versions()
    for dist in ("scipy", "pandas", "pyarrow", "matplotlib", "python-casacore"):
        out[dist] = _version(dist)
    return out


def _run_date(manifest: Dict[str, object], inputs: Sequence[loader_mod.StationInput]) -> str:
    q_start = (manifest.get("prediction") or {}).get("query_start_utc")
    text = q_start or (inputs[0].first_time_utc if inputs else None) or "unknown"
    return text[:10].replace("-", "")


def run(args: argparse.Namespace, command_line: Sequence[str], log_path_holder: List[Path]) -> int:
    """Execute a Stage-3 run; returns the exit status."""
    if args.time_bin <= 0 or args.min_baselines < 1 or args.bootstrap < 0:
        raise InvalidManifestError("--time-bin must be > 0, --min-baselines >= 1, --bootstrap >= 0")
    try:
        options = EstimatorOptions(robust_loss=args.robust_loss, search_radius_deg=args.search_radius_deg,
                                   blind=args.blind)
    except ValueError as exc:
        raise InvalidManifestError(str(exc)) from exc
    manifest, manifest_sha = loader_mod.load_stage2_manifest(args.manifest)
    norad = int(manifest["satellite"]["norad_id"])
    keys = loader_mod.resolve_station_keys(manifest, None if args.all_stations else args.station)
    inputs: List[loader_mod.StationInput] = []
    skipped: Dict[str, str] = {}
    for key in keys:
        try:
            inputs.append(loader_mod.station_input(manifest, key, args.manifest, args.allow_invalid_stage2,
                                                   args.data_root))
        except InvalidManifestError as exc:
            if args.all_stations and args.allow_invalid_stage2 and "not VALID" in str(exc):
                skipped[key] = str(exc)
                logger.warning("skipping %s", exc)
                continue
            raise
    if not inputs:
        raise InvalidManifestError("no station with a valid Stage-2 Measurement Set to process")

    out_dir = args.output_dir / _run_date(manifest, inputs) / f"NORAD_{norad}"
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)
    log_path_holder.append(out_dir / "logs" / "stage3.log")
    sat = manifest["satellite"]
    thresholds = q.QualityThresholds(min_baselines=args.min_baselines, max_flagged_fraction=args.max_flagged_fraction,
                                     min_snr=args.min_snr)
    rng = np.random.default_rng(args.seed)
    all_rows: List[Dict[str, object]] = []
    summaries: Dict[str, object] = {}
    for si in inputs:
        logger.info("station %s (%s): input MS %s", si.station_id, si.tart_id, si.ms_path)
        dataset = loader_mod.load_station_ms(si.ms_path, si.station_id, si.latitude_deg, si.longitude_deg)
        if si.n_antennas and dataset.antenna_enu_m.shape[0] != si.n_antennas:
            raise InvalidVisibilityDataError(
                f"{si.station_id}: MS has {dataset.antenna_enu_m.shape[0]} antennas, Stage 2 recorded {si.n_antennas}")
        tle = None
        if sat.get("tle_line1") and sat.get("tle_line2"):
            tle = TleReference(sat["tle_line1"], sat["tle_line2"], si.latitude_deg, si.longitude_deg, si.altitude_m)
        weights_cal = not bool(np.all(dataset.weight[~dataset.flag] == 1.0))
        thresholds.weights_calibrated = weights_cal
        result = process_station(si, dataset, norad_id=norad, tle=tle, method=args.method,
                                 image_validation=args.image_validation, bin_s=args.time_bin, options=options,
                                 thresholds=thresholds, n_boot=args.bootstrap, rng=rng)
        result.summary["weights_calibrated"] = weights_cal
        all_rows.extend(result.rows)
        summaries[si.station_id] = result.summary

    tables = output_mod.write_tables(all_rows, out_dir)
    logger.info("wrote %s%s", tables["csv"], f" and {tables['parquet']}" if tables["parquet"] else "")
    figures: Dict[str, List[str]] = {}
    if not args.no_plots:
        from tart_gnss.stage3.plotting import make_all_plots

        figures = make_all_plots(tables["parquet"] or tables["csv"], out_dir / "diagnostics")
        logger.info("wrote %d diagnostic figures to %s", sum(len(v) for v in figures.values()), out_dir / "diagnostics")
    n_usable = sum(1 for r in all_rows if r["usable"])
    tle_epoch = sat.get("tle_epoch")
    doc = {
        "schema_version": output_mod.STAGE3_SCHEMA_VERSION,
        "stage": 3,
        "input_manifest": str(args.manifest),
        "input_manifest_sha256": manifest_sha,
        "stage2_three_station_valid": bool(manifest.get("three_station_valid")),
        "stage2_validity_reasons": manifest.get("validity_reasons") or [],
        "norad_id": norad,
        "satellite": {k: sat.get(k) for k in ("satellite_name", "constellation", "signal_relevance",
                                               "tle_epoch", "tle_line1", "tle_line2", "tle_source")},
        "stations": summaries,
        "stations_skipped": skipped,
        "estimator": {
            "primary_method": "image" if args.method == "image" else "visibility",
            "method_option": args.method,
            "image_validation": bool(args.image_validation or args.method in ("image", "both")),
            "measurement_equation": "V_pq = C exp(-2 pi i (nu/c) (r_q - r_p) . s); MS UVW = r_p - r_q (tart2ms, zenith)",
            "model": "single far-field point source, complex nuisance amplitude, per-channel frequency",
            **options.to_dict(),
            "search_centre": "blind hemisphere" if args.blind else "TLE prediction at each bin time (search bound only)",
            "optimizer": "scipy.optimize.least_squares (trf, bounded) after a profiled-chi2 grid search",
            "uncertainty_method": "sigma^2 (J^T J)^-1 from the analytic residual Jacobian; "
                                  "sigma^2 = chi2/dof (linear) or (1.4826 MAD)^2 (robust losses)",
            "snr_definition": "|C_hat| sqrt(sum w / 2) / sigma_vis, sigma_vis = 1.4826 MAD of residual components",
            "bootstrap": {"n_resamples": args.bootstrap, "unit": "baseline (antenna pair)",
                          "random_seed": args.seed},
        },
        "time_binning": {"bin_s": args.time_bin, "edges": "from the first timestamp of each station",
                         "output_time": "weight-averaged time of the samples used", "empty_bins": "omitted"},
        "quality_thresholds": thresholds.to_dict(),
        "tle_reference": {
            "use": "search-cap centre and plotting reference only; never a prior or a measurement",
            "tle_epoch": tle_epoch,
        },
        "software": {**_software(), **_git_state()},
        "provenance": {
            "command_line": list(command_line),
            "run_timestamp_utc": format_utc(datetime.now(timezone.utc)),
            "working_directory": str(Path.cwd()),
            "input_hashes": {"stage2_manifest": manifest_sha},
        },
        "outputs": {
            "directions_csv": tables["csv"],
            "directions_parquet": tables["parquet"],
            "diagnostics": figures,
            "log": str(out_dir / "logs" / "stage3.log"),
            "n_rows": len(all_rows),
            "n_usable": n_usable,
        },
    }
    path = output_mod.write_manifest(doc, out_dir / "stage3_manifest.json")
    logger.info("wrote %s (%d rows, %d usable)", path, len(all_rows), n_usable)
    return EXIT_OK if n_usable else EXIT_NO_ESTIMATES


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``stage3.log``."""

    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.lines: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-localise``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        args = build_parser().parse_args(argv)
    except SystemExit as exc:
        return EXIT_USAGE if exc.code else EXIT_OK
    logging.Formatter.converter = time.gmtime
    fmt = logging.Formatter("%(asctime)sZ %(levelname)-7s %(name)s: %(message)s", "%Y-%m-%dT%H:%M:%S")
    root = logging.getLogger()
    if not root.handlers:
        console = logging.StreamHandler()
        console.setFormatter(fmt)
        root.addHandler(console)
    root.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    for noisy in ("matplotlib", "PIL"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    collector = _ListHandler()
    collector.setFormatter(fmt)
    root.addHandler(collector)
    log_path: List[Path] = []
    try:
        return run(args, ["tart-localise", *argv], log_path)
    except (InvalidManifestError, MissingMeasurementSetError, InvalidVisibilityDataError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (Stage3Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)
        if log_path:
            log_path[0].write_text("\n".join(collector.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
