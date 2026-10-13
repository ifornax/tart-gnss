"""``tart-triangulate``: Stage-4 multi-station 3-D localisation from Stage-3 directions.

Reads ``stage3_manifest.json`` and its direction table, applies the Stage-3
quality policy, matches the stations in time, converts each measured direction
to an ECEF sight line from the station's frozen coordinates, intersects the
lines (covariance-weighted by default) and writes positions, per-station
contributions, ``stage4_manifest.json``, diagnostic plots and a log under
``<output-dir>/<YYYYMMDD>/NORAD_<id>/``.

No TLE enters the computation: TLE information in the Stage-3 manifest is
copied to the Stage-4 manifest for later comparison (Stage 6) only.

Exit codes: 0 at least one valid 3-D estimate; 3 completed but none valid;
2 invalid input or configuration; 1 unexpected failure.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from tart_gnss import __version__
from tart_gnss.stage4 import quality as q
from tart_gnss.stage4.errors import (
    CovarianceError,
    InsufficientStationsError,
    InvalidStage3ManifestError,
    MissingDirectionDataError,
    Stage4Error,
    TimeMatchError,
    TriangulationError,
)
from tart_gnss.stage4.geometry import SPEED_OF_LIGHT_M_S, angular_sigma_deg, ecef_to_geodetic
from tart_gnss.stage4.loader import apply_quality_policy, load_stage3_run, sight_line, sight_vector
from tart_gnss.stage4.output import STAGE4_SCHEMA_VERSION, read_table, write_manifest, write_tables
from tart_gnss.stage4.timing import angular_rates_deg_s, match_station_epochs, motion_during_offset_deg, time_seconds
from tart_gnss.stage4.triangulation import TriangulationResult, triangulate_lines, triangulate_weighted
from tart_gnss.stage4.uncertainty import bootstrap_position, covariance_components

logger = logging.getLogger("tart_gnss.stage4")

EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_NO_VALID = 0, 1, 2, 3
METHODS = ("weighted", "unweighted")
DEFAULT_SEED = 20260625
MJD_UNIX_EPOCH = 40587.0


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``stage4.log``."""

    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.lines: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))


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
    for dist in ("pandas", "pyarrow", "matplotlib"):
        out[dist] = _version(dist)
    return out


def build_parser() -> argparse.ArgumentParser:
    """Argument parser for ``tart-triangulate``."""
    p = argparse.ArgumentParser(
        prog="tart-triangulate",
        description="Stage 4: triangulate the 3-D (ITRS/ECEF) satellite position from the per-station "
                    "directions measured by tart-localise. TLEs are never used in the fit.",
    )
    p.add_argument("--stage3-manifest", required=True, type=Path, help="stage3_manifest.json from tart-localise")
    p.add_argument("--output-dir", type=Path, default=Path("results/stage4"),
                   help="output root (default results/stage4)")
    p.add_argument("--time-match-tolerance", type=float, default=0.75,
                   help="max |t_station - t_reference| in seconds (default 0.75)")
    p.add_argument("--min-stations", type=int, default=3,
                   help="stations required for a valid 3-D estimate (default 3; 2-station solutions are "
                        "written as degraded diagnostics)")
    p.add_argument("--method", choices=METHODS, default="weighted",
                   help="weighted (Stage-3 covariance, default) or unweighted closed form")
    p.add_argument("--bootstrap", type=int, default=0,
                   help="parametric bootstrap resamples per epoch (default 0 = off)")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED, help=f"bootstrap random seed (default {DEFAULT_SEED})")
    p.add_argument("--stations", nargs="+", default=None, help="restrict to these station ids")
    p.add_argument("--quality-policy", choices=q.POLICIES, default="strict",
                   help="strict (default): use only Stage-3 'usable' rows; permissive: also rows flagged "
                        "BOUNDARY_SOLUTION/LOW_SNR/TLE_UNASSOCIATED (epoch marked STAGE3_QUALITY_WARNING)")
    p.add_argument("--reference-time", default=None,
                   help="triangulate only the reference epoch nearest this UTC time (ISO 8601)")
    p.add_argument("--reference-station", default=None,
                   help="station supplying the reference epochs (default: the one with most valid rows)")
    p.add_argument("--max-condition-number", type=float, default=1e8,
                   help="POOR_GEOMETRY above this condition number of sum P_i (default 1e8)")
    p.add_argument("--max-normalised-residual", type=float, default=5.0,
                   help="HIGH_LINE_RESIDUAL above this covariance-normalised line residual (default 5)")
    p.add_argument("--max-reduced-chi2", type=float, default=5.0, help="HIGH_CHI2 threshold (default 5)")
    p.add_argument("--time-mismatch-sigma", type=float, default=1.0,
                   help="TIME_MISMATCH when the target moves more than this many angular sigmas during a "
                        "station's time offset (default 1)")
    p.add_argument("--no-plots", action="store_true", help="skip diagnostic plots")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    return p


def _iso(t_s: float) -> str:
    return datetime.fromtimestamp(t_s, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _nan() -> float:
    return float("nan")


def _solve_epoch(r, s, cov3, method: str, notes: List[str], flags: List[str]) -> TriangulationResult:
    """Run the requested estimator, falling back to unweighted if covariances are unusable."""
    if method == "weighted":
        try:
            return triangulate_weighted(r, s, cov3)
        except CovarianceError as exc:
            flags.append(q.COVARIANCE_SINGULAR)
            notes.append(f"weighted fit impossible ({exc}); unweighted fallback without covariance")
            return triangulate_lines(r, s)
    try:
        return triangulate_lines(r, s, cov3)
    except CovarianceError as exc:
        flags.append(q.COVARIANCE_SINGULAR)
        notes.append(f"{exc}; no position covariance")
        return triangulate_lines(r, s)


def process_epoch(ep, kept, run, rates, args, thresholds, rng) -> tuple:
    """Triangulate one matched epoch; return ``(position_row, contribution_rows, notes)``."""
    flags: List[str] = []
    notes: List[str] = []
    members = []
    for sid in ep.station_ids:
        row = kept.loc[ep.rows[sid]]
        coords = run.stations[sid]
        try:
            s, cov3, src = sight_line(row, coords)
        except CovarianceError as exc:
            # Direction finite but covariance unusable: keep the line, mark the covariance.
            notes.append(f"{sid}: {exc}")
            s, cov3, src = sight_vector(row, coords), np.full((3, 3), np.nan), "no covariance"
        except TriangulationError as exc:
            notes.append(f"{sid}: {exc}")
            continue
        members.append((sid, row, coords, s, cov3, src))
        flags += q.input_flags(row["quality_flags"])
        if coords.calibration_warnings:
            flags.append(q.CALIBRATION_WARNING)

    n = len(members)
    ref_t = ep.reference_time_s
    pos: Dict[str, object] = {
        "time_utc": _iso(ref_t), "time_mjd": ref_t / 86400.0 + MJD_UNIX_EPOCH, "norad_id": run.norad_id,
        "reference_station": ep.reference_station, "method": args.method, "n_stations": n,
        "stations": ";".join(m[0] for m in members), "max_time_offset_s": ep.max_offset_s,
        "fit_success": False, "converged": False,
    }
    if n < thresholds.min_stations:
        flags.append(q.INSUFFICIENT_STATIONS)
    res = None
    if n >= 2:
        r = np.array([m[2].ecef_m for m in members])
        s = np.array([m[3] for m in members])
        cov3 = np.array([m[4] for m in members])
        if not np.all(np.isfinite(cov3)):
            flags.append(q.COVARIANCE_SINGULAR)
            res = triangulate_lines(r, s)
        else:
            res = _solve_epoch(r, s, cov3, args.method, notes, flags)
        notes += res.notes
    contribs = []
    sigmas_theta = []
    motions = []
    for i, (sid, row, coords, s_i, cov3_i, src) in enumerate(members):
        off = ep.offsets_s[sid]
        rate = rates.get((sid, ep.rows[sid]), _nan())
        motion = motion_during_offset_deg(rate, off) if np.isfinite(rate) else _nan()
        sig = angular_sigma_deg(cov3_i) if np.all(np.isfinite(cov3_i)) else _nan()
        st_flags = []
        if np.isfinite(motion) and np.isfinite(sig) and motion > thresholds.time_mismatch_sigma * sig:
            st_flags.append(q.TIME_MISMATCH)
        sigmas_theta.append(sig)
        motions.append(motion)
        rng_m = float(res.station_ranges_m[i]) if res is not None else _nan()
        resid = float(res.line_residuals_m[i]) if res is not None else _nan()
        nres = float(res.normalised_residuals[i]) if res is not None else _nan()
        if res is not None:
            if rng_m <= 0:
                st_flags.append(q.NEGATIVE_RANGE)
            if np.isfinite(nres):
                if nres > thresholds.max_normalised_residual:
                    st_flags.append(q.HIGH_LINE_RESIDUAL)
            elif resid > thresholds.max_line_residual_m:
                st_flags.append(q.HIGH_LINE_RESIDUAL)
        flags += st_flags
        contribs.append({
            "time_utc": pos["time_utc"], "reference_time_utc": pos["time_utc"], "station_id": sid, "station_time_utc": _iso(ref_t + off),
            "station_time_offset_s": off,
            "station_x_ecef_m": coords.ecef_m[0], "station_y_ecef_m": coords.ecef_m[1],
            "station_z_ecef_m": coords.ecef_m[2],
            "az_obs_deg": row["az_obs_deg"], "el_obs_deg": row["el_obs_deg"], "sigma_az_deg": row["sigma_az_deg"],
            "sigma_el_deg": row["sigma_el_deg"], "cov_az_el": row["cov_az_el"], "l": row["l"], "m": row["m"],
            "sigma_theta_deg": sig, "direction_source": src,
            "sight_x": s_i[0], "sight_y": s_i[1], "sight_z": s_i[2],
            "range_m": rng_m, "line_residual_m": resid, "normalised_residual": nres,
            "angular_rate_deg_s": rate, "motion_during_offset_deg": motion,
            "light_time_s": rng_m / SPEED_OF_LIGHT_M_S if np.isfinite(rng_m) and rng_m > 0 else _nan(),
            "stage3_quality_flags": list(row["quality_flags"]), "stage4_station_flags": q.ordered(st_flags),
        })
    pos["max_motion_during_offset_deg"] = float(np.nanmax(motions)) if np.any(np.isfinite(motions)) else _nan()
    if res is not None:
        x = res.position_ecef_m
        pos.update({"x_ecef_m": x[0], "y_ecef_m": x[1], "z_ecef_m": x[2]})
        pos.update(covariance_components(res.covariance_ecef_m2))
        resid = res.line_residuals_m
        pos.update({
            "rms_line_residual_m": float(np.sqrt(np.mean(resid ** 2))), "max_line_residual_m": float(resid.max()),
            "max_normalised_residual": float(np.nanmax(res.normalised_residuals))
            if np.any(np.isfinite(res.normalised_residuals)) else _nan(),
            "condition_number": res.condition_number, "matrix_rank": res.matrix_rank,
            "min_crossing_angle_deg": res.min_crossing_angle_deg, "max_crossing_angle_deg": res.max_crossing_angle_deg,
            "chi2": res.chi2 if res.chi2 is not None else _nan(), "dof": res.dof,
            "reduced_chi2": res.reduced_chi2 if res.reduced_chi2 is not None else _nan(),
            "mean_light_time_s": float(np.mean(res.station_ranges_m)) / SPEED_OF_LIGHT_M_S,
            "solver": res.solver, "n_iterations": res.n_iterations, "converged": res.converged,
            "fit_success": res.success,
        })
        if res.solver == "lstsq":
            notes.append("SVD least-squares fallback used for the point estimate")
        if not res.success:
            flags.append(q.TRIANGULATION_FAILED)
        if res.matrix_rank < 3 or not (res.condition_number <= thresholds.max_condition_number):
            flags.append(q.POOR_GEOMETRY)
        if not res.covariance_ok and q.COVARIANCE_SINGULAR not in flags:
            flags.append(q.COVARIANCE_SINGULAR)
        if res.reduced_chi2 is not None and res.reduced_chi2 > thresholds.max_reduced_chi2:
            flags.append(q.HIGH_CHI2)
        if np.all(np.isfinite(x)):
            lat, lon, alt = ecef_to_geodetic(x)
            pos.update({"latitude_deg": lat, "longitude_deg": lon, "altitude_m": alt})
        if args.bootstrap > 0 and res.success and np.all(np.isfinite(cov3)):
            bs = bootstrap_position(r, s, cov3, args.bootstrap, rng, args.method,
                                    linear_sigma_m=res.sigma_ecef_m if res.covariance_ok else None,
                                    point_estimate_m=x)
            pos.update(bs.percentiles)
            pos.update({"bootstrap_n_success": bs.n_success, "bootstrap_failure_fraction": bs.failure_fraction})
            # The median-shift limit is never below 3x the sampling noise of a median (1.2533 sigma/sqrt(n)).
            if (bs.failure_fraction > thresholds.bootstrap_max_failure_fraction
                    or bs.max_sigma_ratio > thresholds.bootstrap_max_sigma_ratio
                    or bs.max_median_shift_sigma > max(thresholds.bootstrap_max_median_shift_sigma,
                                                       3.0 * 1.2533 / np.sqrt(max(bs.n_success, 1)))):
                flags.append(q.BOOTSTRAP_UNSTABLE)
    else:
        for c in ("x_ecef_m", "y_ecef_m", "z_ecef_m", "sigma_x_m", "sigma_y_m", "sigma_z_m"):
            pos[c] = _nan()
    pos["quality_flags"] = q.ordered(flags)
    pos["valid"] = q.is_valid(pos["quality_flags"]) and bool(pos["fit_success"])
    if notes:
        logger.debug("%s: %s", pos["time_utc"], "; ".join(notes))
    return pos, contribs, notes


def _parse_time(text: str) -> float:
    import pandas as pd

    try:
        return float(time_seconds([pd.Timestamp(text)])[0])
    except (ValueError, TypeError) as exc:
        raise TimeMatchError(f"cannot parse --reference-time {text!r}") from exc


def _validate_args(args) -> None:
    if args.time_match_tolerance < 0:
        raise InvalidStage3ManifestError("--time-match-tolerance must be >= 0")
    if args.min_stations < 2:
        raise InvalidStage3ManifestError("--min-stations must be >= 2")
    if args.bootstrap < 0:
        raise InvalidStage3ManifestError("--bootstrap must be >= 0")
    if not args.max_condition_number > 1:
        raise InvalidStage3ManifestError("--max-condition-number must be > 1")


def _median(values) -> float:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else float("nan")


def run(args: argparse.Namespace, command_line: Sequence[str], log_path_holder: List[Path]) -> int:
    """Execute Stage 4 for parsed arguments; return the exit status."""
    git = _git_state()
    started = datetime.now(timezone.utc)
    _validate_args(args)
    thresholds = q.Stage4Thresholds(
        min_stations=args.min_stations, max_condition_number=args.max_condition_number,
        max_normalised_residual=args.max_normalised_residual, max_reduced_chi2=args.max_reduced_chi2,
        time_mismatch_sigma=args.time_mismatch_sigma,
    )
    logger.info("input Stage-3 manifest: %s", args.stage3_manifest)
    run3 = load_stage3_run(args.stage3_manifest, args.stations)
    if len(run3.stations) < 2:
        raise InsufficientStationsError(f"only {len(run3.stations)} station(s) selected; at least 2 needed")
    logger.info("stations: %d (%s); direction table %s", len(run3.stations), ", ".join(run3.stations),
                run3.directions_path)
    df = run3.directions
    kept, rejected = apply_quality_policy(df, args.quality_policy)
    frac_rej = 1.0 - len(kept) / len(df) if len(df) else float("nan")
    logger.info("Stage-3 rows: %d, kept %d under '%s' policy (%.1f%% rejected: %s)", len(df), len(kept),
                args.quality_policy, 100 * frac_rej if len(df) else 0.0,
                ", ".join(f"{k} {v}" for k, v in sorted(rejected.items())) or "none")

    tables = {sid: kept[kept["station_id"] == sid] for sid in run3.stations}
    tables = {k: v for k, v in tables.items() if len(v)}
    rates: Dict[tuple, float] = {}
    for sid, t in tables.items():
        vecs = []
        for _, row in t.iterrows():
            try:
                vecs.append(sight_line(row, run3.stations[sid])[0])
            except (TriangulationError, CovarianceError):
                vecs.append(np.full(3, np.nan))
        rr = angular_rates_deg_s(time_seconds(t["time_utc"]), np.array(vecs))
        rates.update({(sid, idx): float(v) for idx, v in zip(t.index, rr)})

    ref_time = _parse_time(args.reference_time) if args.reference_time else None
    epochs = match_station_epochs(tables, args.time_match_tolerance, args.reference_station, ref_time) \
        if tables else []
    logger.info("matched epochs: %d (reference station %s, tolerance %.3f s)", len(epochs),
                epochs[0].reference_station if epochs else "-", args.time_match_tolerance)

    rng = np.random.default_rng(args.seed)
    positions, contributions = [], []
    n_fallback = 0
    for ep in epochs:
        pos, con, notes = process_epoch(ep, kept, run3, rates, args, thresholds, rng)
        n_fallback += int(any("fallback" in n for n in notes))
        positions.append(pos)
        contributions += con

    first_t = df["time_utc"].min() if len(df) else None
    date = first_t.strftime("%Y%m%d") if first_t is not None else "unknown"
    out = args.output_dir / date / f"NORAD_{run3.norad_id}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    log_path_holder.append(out / "logs" / "stage4.log")
    paths = write_tables(positions, contributions, out)

    valid = [p for p in positions if p["valid"]]
    solved = [p for p in positions if p["n_stations"] >= 2]
    n3 = [p for p in positions if p["n_stations"] >= args.min_stations]
    logger.info("epochs with >= %d stations: %d; triangulation success fraction %.3f; valid 3-D estimates: %d",
                args.min_stations, len(n3),
                (sum(bool(p["fit_success"]) for p in solved) / len(solved)) if solved else 0.0, len(valid))
    logger.info("median |time offset| %.3f s; median min crossing angle %.3f deg; median condition number %.3g",
                _median([abs(c["station_time_offset_s"]) for c in contributions]),
                _median([p.get("min_crossing_angle_deg") for p in solved]),
                _median([p.get("condition_number") for p in solved]))
    logger.info("median line residual %.3g km; median position sigma (x, y, z) = (%.3g, %.3g, %.3g) km",
                _median([c["line_residual_m"] for c in contributions]) / 1e3,
                *(_median([p.get(f"sigma_{a}_m") for p in solved]) / 1e3 for a in "xyz"))

    plots: List[str] = []
    if not args.no_plots and positions:
        try:
            from tart_gnss.stage4.plotting import make_plots

            key = "positions_parquet" if paths["positions_parquet"] else "positions_csv"
            ckey = "station_contributions_parquet" if paths["station_contributions_parquet"] \
                else "station_contributions_csv"
            plots = make_plots(read_table(paths[key]), read_table(paths[ckey]), out / "diagnostics",
                               f"NORAD {run3.norad_id}")
        except ImportError as exc:
            logger.warning("plots skipped: %s", exc)

    station_doc = {}
    for sid, c in run3.stations.items():
        rows = df[df["station_id"] == sid]
        con = [x for x in contributions if x["station_id"] == sid]
        station_doc[sid] = {
            **c.to_dict(),
            "n_stage3_rows": int(len(rows)), "n_kept": int(len(tables.get(sid, []))),
            "n_epochs_contributed": len(con),
            "median_line_residual_m": _median([x["line_residual_m"] for x in con]),
            "median_normalised_residual": _median([x["normalised_residual"] for x in con]),
            "median_abs_time_offset_s": _median([abs(x["station_time_offset_s"]) for x in con]),
        }
    flag_counts: Dict[str, int] = {}
    for p in positions:
        for f in p["quality_flags"]:
            flag_counts[f] = flag_counts.get(f, 0) + 1
    m3 = run3.manifest
    doc = {
        "schema_version": STAGE4_SCHEMA_VERSION,
        "stage": 4,
        "input_stage3_manifest": str(args.stage3_manifest),
        "input_sha256": run3.sha256,
        "input_directions": str(run3.directions_path),
        "input_directions_sha256": run3.directions_sha256,
        "norad_id": run3.norad_id,
        "coordinate_frame": "ITRS/ECEF",
        "units": {"position": "m", "covariance": "m^2", "angles": "deg", "time": "UTC"},
        "time_matching": {
            "tolerance_s": args.time_match_tolerance,
            "reference_station": epochs[0].reference_station if epochs else args.reference_station,
            "reference_station_rule": "--reference-station" if args.reference_station
            else "station with the most valid Stage-3 rows (ties: alphabetical)",
            "reference_time": args.reference_time,
            "interpolation": "none (nearest valid estimate within tolerance; each estimate used at most once)",
            "n_epochs": len(epochs),
            "n_epochs_min_stations": len(n3),
            "median_abs_offset_s": _median([abs(c["station_time_offset_s"]) for c in contributions]),
            "max_abs_offset_s": float(max((abs(c["station_time_offset_s"]) for c in contributions), default=0.0)),
            "angular_rate_source": "each station's own Stage-3 track (linear fit over a 120 s sliding window); never a TLE",
            "light_time": "rho/c recorded (~70 ms for GNSS); not corrected (no retarded-time modelling)",
        },
        "triangulation": {
            "method": args.method,
            "baseline_estimator": "argmin sum |P_i (x - r_i)|^2, closed form (sum P_i) x = sum P_i r_i",
            "weighted_estimator": "iteratively reweighted: W_i = (rho_i^2 E_i^T Sigma_theta,i E_i)^-1, "
                                  "start from the unweighted solution, tol 1 mm, max 50 iterations",
            "direction_input": "Stage-3 (l, m) and cov(l, m) (Az/El + cov_az_el only if (l, m) are missing)",
            "solver": "numpy.linalg.solve; SVD lstsq if the condition number exceeds 1e12",
            "n_epochs_numerical_fallback": n_fallback,
            "tle_used": False,
        },
        "uncertainty": {
            "linear_covariance": "Sigma_x = N^-1 (weighted) or A^-1 (sum P_i rho_i^2 Sigma_theta,i P_i) A^-1 "
                                 "(unweighted); first order in the angular errors, valid while they are small "
                                 "compared with the crossing angles",
            "bootstrap": args.bootstrap,
            "bootstrap_model": "Gaussian redraw of each direction from its Stage-3 covariance (tangent plane)",
            "random_seed": args.seed,
        },
        "quality_policy": {
            "policy": args.quality_policy,
            "always_rejected": sorted(q.STAGE3_FATAL_FLAGS),
            "rejected_if_strict": sorted(q.STAGE3_STRICT_FLAGS),
            "advisory": sorted(q.STAGE3_ADVISORY_FLAGS) + [q.STAGE3_CALIBRATION_FLAG],
            "n_stage3_rows": int(len(df)), "n_kept": int(len(kept)), "rejected": rejected,
            "rejected_fraction": frac_rej,
        },
        "quality_thresholds": thresholds.to_dict(),
        "blocking_flags": sorted(q.BLOCKING_FLAGS),
        "summary": {
            "n_epochs": len(positions), "n_valid": len(valid),
            "n_with_min_stations": len(n3),
            "flag_counts": flag_counts,
        },
        "stations": station_doc,
        "satellite_reference": {
            "note": "copied from the Stage-3 manifest for later comparison (Stage 6); not used by Stage 4",
            "satellite": m3.get("satellite"),
        },
        "software": {**_software(), **git, "tart_gnss": __version__},
        "provenance": {
            "command_line": list(command_line),
            "run_timestamp_utc": started.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z",
            "working_directory": os.getcwd(),
            "stage3_run": {"command_line": (m3.get("provenance") or {}).get("command_line"),
                           "git_commit": (m3.get("software") or {}).get("git_commit")},
        },
        "outputs": {**paths, "diagnostics": plots, "log": str(out / "logs" / "stage4.log"),
                    "manifest": str(out / "stage4_manifest.json")},
    }
    write_manifest(doc, out / "stage4_manifest.json")
    logger.info("outputs: %s", out)
    for k, v in paths.items():
        if v:
            logger.info("  %s: %s", k, v)
    if not valid:
        logger.warning("no valid 3-D estimate (flags: %s)",
                       ", ".join(f"{k} {v}" for k, v in sorted(flag_counts.items())) or "none")
        return EXIT_NO_VALID
    return EXIT_OK


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-triangulate``."""
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
        return run(args, ["tart-triangulate", *argv], log_path)
    except (InvalidStage3ManifestError, MissingDirectionDataError, TimeMatchError, InsufficientStationsError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (Stage4Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)
        if log_path:
            log_path[0].write_text("\n".join(collector.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
