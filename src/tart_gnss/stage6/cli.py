"""``tart-validate``: Stage-6 comparison of the Stage-5 trajectory with a TLE prediction.

Reads ``stage5_manifest.json`` and its trajectory (smoothed GCRS states with
full 6x6 covariances), selects exactly one reference TLE (``--tle-file`` or
``--tle-source stage1|download``), propagates it with SGP4 at the exact
Stage-5 epochs, transforms it to GCRS, and writes ``observed - TLE``
residuals (Cartesian, RTN, angular, range, along-track time equivalent),
Stage-5 uncertainties rotated into RTN, measurement-conditional (and, with an
ensemble of >= 3 TLEs, reference-aware) Mahalanobis distances, descriptive
statistics, TLE-age status, ``summary.json``, ``stage6_manifest.json``,
diagnostic plots and a log under ``<output-dir>/<YYYYMMDD>/NORAD_<id>/``.

The TLE is never fitted, shifted or used to modify the Stage-5 trajectory.

Exit codes: 0 comparison valid for precision analysis (or descriptive with
``--allow-stale-reference`` when staleness is the only failed condition);
3 completed but not valid for precision analysis; 2 invalid input or
configuration (including TLE_NOT_FOUND / TLE_NORAD_MISMATCH); 1 unexpected
failure (including TLE_FRAME_FAILURE).
"""

from __future__ import annotations

import argparse
import json
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
from tart_gnss.stage6 import quality as q
from tart_gnss.stage6.errors import (
    FrameTransformError,
    InvalidStage5ManifestError,
    MissingTrajectoryError,
    ReferenceCovarianceError,
    Stage6Error,
    TLENoradMismatchError,
    TLEPropagationError,
    TLEReferenceError,
)
from tart_gnss.stage6.frames import FRAMES
from tart_gnss.stage6.loader import load_stage5_trajectory, row_state
from tart_gnss.stage5.output import covariance_fields
from tart_gnss.stage6.output import STAGE6_SCHEMA_VERSION, read_table, write_manifest, write_tables
from tart_gnss.stage6.propagation import propagate_tle
from tart_gnss.stage6.reference_uncertainty import LABEL, MIN_ENSEMBLE, empirical_tle_covariance
from tart_gnss.stage6.residuals import StateTrajectory, compute_residuals, rows as residual_rows
from tart_gnss.stage6.rtn import covariance_to_rtn, rtn_basis
from tart_gnss.stage6.statistics import describe, effective_sample_size, lag1_autocorrelation, trend
from tart_gnss.stage6.tle_reference import (
    SOURCES,
    download_tle_reference,
    load_ensemble,
    load_tle_reference,
    stage1_tle_reference,
)
from tart_gnss.tle import format_utc

logger = logging.getLogger("tart_gnss.stage6")

EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_NOT_PRECISION = 0, 1, 2, 3
MJD_UNIX_EPOCH = 40587.0
RTN = ("radial", "along_track", "cross_track")
LIGHT_TIME_NOTE = ("Stages 4-5 do not correct light time: a measured position refers to about t - 70 ms for GNSS, "
                   "so an along-track residual of about -0.27 km (-0.07 s equivalent) is expected from that alone")


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``stage6.log``."""

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
    for dist in ("scipy", "pandas", "pyarrow", "matplotlib"):
        out[dist] = _version(dist)
    return out


def build_parser() -> argparse.ArgumentParser:
    """Argument parser for ``tart-validate``."""
    p = argparse.ArgumentParser(
        prog="tart-validate",
        description="Stage 6: compare the TLE-independent Stage-5 trajectory with a TLE (SGP4) prediction; "
                    "residuals are observed - TLE. The TLE is a reference, not exact truth.",
    )
    p.add_argument("--stage5-manifest", required=True, type=Path, help="stage5_manifest.json from tart-track")
    p.add_argument("--output-dir", type=Path, default=Path("results/stage6"),
                   help="output root (default results/stage6)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--tle-file", type=Path, help="TLE file holding the reference element set")
    src.add_argument("--tle-source", choices=SOURCES,
                     help="stage1: the exact Stage-1 element set from the run provenance; download: explicit "
                          "CelesTrak download (current elements only), archived before use")
    p.add_argument("--tle-download-dir", type=Path, default=Path("data/tle"),
                   help="archive directory for --tle-source download (default data/tle)")
    p.add_argument("--tle-ensemble-files", nargs="+", type=Path, default=[],
                   help="TLE files for the empirical reference dispersion")
    p.add_argument("--tle-ensemble-dir", type=Path, default=None,
                   help="directory of *.tle/*.txt files for the empirical reference dispersion")
    p.add_argument("--tle-ensemble-window-days", type=float, default=7.0,
                   help="ensemble sets within this many days of the trajectory midpoint (default 7)")
    p.add_argument("--tle-age-warning-days", type=float, default=3.0, help="TLE_AGE_WARNING above (default 3)")
    p.add_argument("--tle-age-stale-days", type=float, default=14.0,
                   help="TLE_STALE above, and no precision validation (default 14)")
    p.add_argument("--quality-policy", choices=q.POLICIES, default="strict",
                   help="strict (default): valid Stage-5 rows only; permissive: any row with a finite state")
    p.add_argument("--position-residual-warning-km", type=float, default=None,
                   help="HIGH_POSITION_RESIDUAL above this |dr| (no flag if omitted)")
    p.add_argument("--velocity-residual-warning-m-s", type=float, default=None,
                   help="HIGH_VELOCITY_RESIDUAL above this |dv| (no flag if omitted)")
    p.add_argument("--allow-stale-reference", action="store_true",
                   help="exit 0 for a descriptive comparison when the TLE age is the only failed condition")
    p.add_argument("--no-plots", action="store_true", help="skip diagnostic plots")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    return p


def _iso(t_s: float) -> str:
    return datetime.fromtimestamp(t_s, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _validate_args(args) -> None:
    if not 0 <= args.tle_age_warning_days <= args.tle_age_stale_days:
        raise InvalidStage5ManifestError("need 0 <= --tle-age-warning-days <= --tle-age-stale-days")
    if args.tle_ensemble_window_days <= 0:
        raise InvalidStage5ManifestError("--tle-ensemble-window-days must be > 0")
    for name in ("position_residual_warning_km", "velocity_residual_warning_m_s"):
        v = getattr(args, name)
        if v is not None and not v > 0:
            raise InvalidStage5ManifestError(f"--{name.replace('_', '-')} must be > 0")


def _finite(a) -> bool:
    return bool(np.all(np.isfinite(np.asarray(a, dtype=float))))


def run(args: argparse.Namespace, command_line: Sequence[str], log_path_holder: List[Path]) -> int:
    """Execute Stage 6 for parsed arguments; return the exit status."""
    git = _git_state()
    started = datetime.now(timezone.utc)
    _validate_args(args)
    logger.info("input Stage-5 manifest: %s", args.stage5_manifest)
    s5 = load_stage5_trajectory(args.stage5_manifest)
    if s5.norad_id is None:
        raise InvalidStage5ManifestError("Stage-5 run records no NORAD ID")
    df = s5.table
    n_rows = len(df)
    mid = s5.midpoint_s
    logger.info("NORAD %d; Stage-5 rows %d (table %s)", s5.norad_id, n_rows, s5.trajectory_path)

    if args.tle_file is not None:
        tle = load_tle_reference(s5.norad_id, args.tle_file, mid)
    elif args.tle_source == "stage1":
        tle = stage1_tle_reference(s5.manifest_path, s5.manifest, s5.norad_id)
    else:
        tle = download_tle_reference(s5.norad_id, args.tle_download_dir, mid)
    logger.info("reference TLE: %s (%s); epoch %s; %s", tle.name, tle.source_kind, format_utc(tle.epoch_utc, 3),
                tle.selection_rule)

    # Stage-5 input policy.
    use, sources, reasons, states, covs = [], [], [], np.full((n_rows, 6), np.nan), np.full((n_rows, 6, 6), np.nan)
    for k, (_, row) in enumerate(df.iterrows()):
        xs, ps = row_state(row, "smoothed")
        xf, pf = row_state(row, "filtered") if "x_gcrs_filtered_m" in df.columns else (np.full(6, np.nan), None)
        ok, src, why = q.stage5_row_decision(row["quality_flags"], bool(row["valid"]), _finite(xs), _finite(xf),
                                             args.quality_policy)
        use.append(ok)
        sources.append(src)
        reasons.append(why)
        if ok:
            states[k], covs[k] = (xs, ps) if src == "smoothed" else (xf, pf)
    use = np.array(use, bool)
    idx = np.flatnonzero(use)
    logger.info("Stage-5 rows compared under '%s' policy: %d of %d", args.quality_policy, len(idx), n_rows)

    # TLE propagation at every Stage-5 epoch (exact epochs, no interpolation).
    tle_traj = propagate_tle(tle, s5.times_s) if n_rows else None
    ages = tle_traj.age_days if tle_traj is not None else np.zeros(0)
    if len(idx):
        a = ages[idx]
        logger.info("TLE age: median %.3f d, max |age| %.3f d; %s propagation", float(np.median(a)),
                    float(np.max(np.abs(a))), "backward" if np.all(a < 0) else
                    ("forward" if np.all(a >= 0) else "forward and backward"))

    # Optional ensemble.
    ens_info: Dict[str, object] = {"requested": bool(args.tle_ensemble_files or args.tle_ensemble_dir)}
    ref_cov = None
    run_flags: List[str] = []
    if ens_info["requested"]:
        members, info = load_ensemble(s5.norad_id, args.tle_ensemble_files,
                                      [args.tle_ensemble_dir] if args.tle_ensemble_dir else [], mid,
                                      args.tle_ensemble_window_days)
        ens_info.update(info)
        ens_info["primary_in_ensemble"] = any(m.line1 == tle.line1 and m.line2 == tle.line2 for m in members)
        if len(members) >= MIN_ENSEMBLE and len(idx):
            try:
                ref_cov = empirical_tle_covariance(members, s5.times_s[idx])
                ens_info.update(ref_cov.summary())
                d = tle_traj.state_gcrs[idx] - ref_cov.mean_state
                ens_info["primary_minus_ensemble_mean"] = {
                    "median_position_m": float(np.median(np.linalg.norm(d[:, :3], axis=1))),
                    "median_velocity_m_s": float(np.median(np.linalg.norm(d[:, 3:], axis=1))),
                    "note": "the primary TLE remains the reference; the ensemble mean is not used in its place"}
            except ReferenceCovarianceError as exc:
                logger.warning("ensemble unusable: %s", exc)
        else:
            run_flags.append(q.TLE_ENSEMBLE_SMALL)
        logger.info("TLE ensemble: %d distinct set(s) in the %.1f-day window%s", len(members),
                    args.tle_ensemble_window_days, "" if ref_cov is not None else
                    f" (fewer than {MIN_ENSEMBLE}: no reference-aware statistics)")
    if ref_cov is None:
        run_flags.append(q.REFERENCE_UNCERTAINTY_UNKNOWN)
    ens_info.update({"n_tle_ensemble": ref_cov.n if ref_cov is not None else 0, "label": LABEL,
                     "min_members_for_total_covariance": MIN_ENSEMBLE})

    # Residuals.
    res = None
    if len(idx):
        obs = StateTrajectory(times_s=s5.times_s[idx], states=states[idx], covariances=covs[idx])
        ref = StateTrajectory(times_s=s5.times_s[idx], states=tle_traj.state_gcrs[idx])
        res = compute_residuals(obs, ref, ref_cov.covariance if ref_cov is not None else None)
        for note in res.notes:
            logger.info("%s", note)

    # Run-level statistics and flags.
    stats: Dict[str, object] = {}
    if res is not None:
        sig_rtn = np.sqrt(np.clip(np.diagonal(res.cov_r_rtn, axis1=1, axis2=2), 0, None))
        t_rel = res.times_s
        for j, c in enumerate(RTN):
            tr = trend(t_rel, res.dr_rtn[:, j], sig_rtn[:, j] if np.all(np.isfinite(sig_rtn[:, j])) else None)
            rho = lag1_autocorrelation(res.dr_rtn[:, j])
            med_sig = float(np.nanmedian(sig_rtn[:, j])) if np.any(np.isfinite(sig_rtn[:, j])) else float("nan")
            stats[c] = {**describe(res.dr_rtn[:, j]), "trend": tr, "lag1_autocorrelation": rho,
                        "effective_sample_size": effective_sample_size(len(res.times_s), rho),
                        "median_sigma_obs_m": med_sig}
            stats[f"velocity_{c}"] = describe(res.dv_rtn[:, j])
            if np.isfinite(med_sig) and abs(tr["change_over_span"]) > 3 * med_sig:
                run_flags.append(q.STRONG_RESIDUAL_TREND)
            if len(res.times_s) >= 5 and np.isfinite(rho) and rho > 3 / np.sqrt(len(res.times_s)):
                run_flags.append(q.RESIDUAL_AUTOCORRELATED)
        stats["position_3d"] = describe(res.position_norm_m)
        stats["velocity_3d"] = describe(res.velocity_norm_m_s)
        stats["angular_deg"] = describe(res.angular_deg)
        stats["range_m"] = describe(res.range_m)
        stats["along_track_time_s"] = describe(res.along_track_time_s)
        stats["position_mahalanobis_obs_cov"] = describe(res.d2_pos_obs)
        stats["state_mahalanobis_obs_cov"] = describe(res.d2_state_obs)
        stats["position_mahalanobis_total_cov"] = describe(res.d2_pos_total)
        stats["state_mahalanobis_total_cov"] = describe(res.d2_state_total)
        stats["fraction_position_d2_above_14.16"] = float(np.mean(res.d2_pos_obs[np.isfinite(res.d2_pos_obs)]
                                                                  > q.CONDITIONAL_D2_LIMIT)) \
            if np.any(np.isfinite(res.d2_pos_obs)) else float("nan")

    max_age = float(np.max(np.abs(ages[idx]))) if len(idx) else float("nan")
    conditions = {
        "tle_found_and_norad_matches": True,
        "tle_age_within_stale_limit": bool(len(idx) and max_age <= args.tle_age_stale_days),
        "stage5_trajectory_valid": bool(len(idx)),
        "frame_transformation_ok": True,
        "no_state_covariance_failure": bool(res is not None and res.covariance_ok.all()),
    }
    precise = q.precision_validation(conditions)
    if not precise:
        run_flags.append(q.NOT_PRECISION_VALIDATION)
    logger.info("precision validation: %s%s", precise, "" if precise else " (failed: " + ", ".join(
        k for k, v in conditions.items() if not v) + ")")

    # Rows.
    resid_list = residual_rows(res) if res is not None else []
    pos_in_res = {int(k): j for j, k in enumerate(idx)}
    out_rows, tle_rows = [], []
    for k in range(n_rows):
        row = df.iloc[k]
        t = float(s5.times_s[k])
        flags: List[str] = []
        s5flags = list(row["quality_flags"])
        rec: Dict[str, object] = {
            "time_utc": _iso(t), "time_mjd": t / 86400.0 + MJD_UNIX_EPOCH, "norad_id": s5.norad_id,
            "tle_epoch_utc": format_utc(tle.epoch_utc, 6), "tle_age_days": float(ages[k]),
            "tle_propagation_direction": "forward" if ages[k] >= 0 else "backward",
            "n_tle_ensemble": ref_cov.n if ref_cov is not None else 0,
            "stage5_state_source": sources[k], "stage5_valid": bool(row["valid"]),
            "precision_validation": precise, "stage5_quality_flags": s5flags,
        }
        xt = tle_traj.state_gcrs[k]
        rec.update({"tle_x_gcrs_m": xt[0], "tle_y_gcrs_m": xt[1], "tle_z_gcrs_m": xt[2], "tle_vx_gcrs_m_s": xt[3],
                    "tle_vy_gcrs_m_s": xt[4], "tle_vz_gcrs_m_s": xt[5]})
        te = tle_traj.position_ecef_m[k]
        rec.update({"tle_x_ecef_m": te[0], "tle_y_ecef_m": te[1], "tle_z_ecef_m": te[2]})
        flags += q.age_flags(float(ages[k]), args.tle_age_warning_days, args.tle_age_stale_days)
        if k in pos_in_res:
            j = pos_in_res[k]
            xo = states[k]
            rec.update({"obs_x_gcrs_m": xo[0], "obs_y_gcrs_m": xo[1], "obs_z_gcrs_m": xo[2],
                        "obs_vx_gcrs_m_s": xo[3], "obs_vy_gcrs_m_s": xo[4], "obs_vz_gcrs_m_s": xo[5]})
            oe = tle_traj.matrices[k].T @ xo[:3]
            rec.update({"obs_x_ecef_m": oe[0], "obs_y_ecef_m": oe[1], "obs_z_ecef_m": oe[2]})
            rec.update(resid_list[j])
            if ref_cov is not None:
                cref = covariance_to_rtn(res.bases[j], ref_cov.covariance[j][:3, :3])
                for i, c in enumerate(RTN):
                    rec[f"sigma_{c}_ref_m"] = float(np.sqrt(max(cref[i, i], 0.0)))
                rec.update(covariance_fields(ref_cov.covariance[j], "tle_dispersion"))
                if ref_cov.position_dispersion_m[j] > np.sqrt(np.trace(covs[k][:3, :3])):
                    flags.append(q.TLE_REFERENCE_DISPERSION_HIGH)
            if not res.covariance_ok[j]:
                flags.append(q.INVALID_STAGE5_COVARIANCE)
            d2 = res.d2_pos_obs[j]
            if np.isfinite(d2) and d2 > q.CONDITIONAL_D2_LIMIT:
                flags.append(q.HIGH_CONDITIONAL_MAHALANOBIS)
            if args.position_residual_warning_km is not None and \
                    rec["position_residual_3d_m"] > 1e3 * args.position_residual_warning_km:
                flags.append(q.HIGH_POSITION_RESIDUAL)
            if args.velocity_residual_warning_m_s is not None and \
                    rec["velocity_residual_3d_m_s"] > args.velocity_residual_warning_m_s:
                flags.append(q.HIGH_VELOCITY_RESIDUAL)
            if set(s5flags) & q.STAGE5_ADVISORY:
                flags.append(q.STAGE5_QUALITY_WARNING)
            flags += run_flags
            rec["valid"] = bool(np.isfinite(rec["position_residual_3d_m"]))
        else:
            flags.append(q.INVALID_STAGE5_STATE)
            rec["stage5_rejection_reason"] = reasons[k]
            rec["valid"] = False
            if q.NOT_PRECISION_VALIDATION in run_flags:
                flags.append(q.NOT_PRECISION_VALIDATION)
        rec["quality_flags"] = q.ordered(flags)
        out_rows.append(rec)
        tr_row = {"time_utc": rec["time_utc"], "norad_id": s5.norad_id, "tle_epoch_utc": rec["tle_epoch_utc"],
                  "tle_age_days": rec["tle_age_days"], "x_gcrs_m": xt[0], "y_gcrs_m": xt[1], "z_gcrs_m": xt[2],
                  "vx_gcrs_m_s": xt[3], "vy_gcrs_m_s": xt[4], "vz_gcrs_m_s": xt[5],
                  "x_ecef_m": te[0], "y_ecef_m": te[1], "z_ecef_m": te[2]}
        try:
            qb = rtn_basis(xt[:3], xt[3:])
            for name, v in zip(("r_hat", "t_hat", "n_hat"), qb):
                tr_row.update({f"{name}_x": v[0], f"{name}_y": v[1], f"{name}_z": v[2]})
        except Stage6Error:
            pass
        tle_rows.append(tr_row)

    if res is not None:
        logger.info("3-D position residual: median %.4g km, RMS %.4g km; 3-D velocity residual median %.4g m/s",
                    stats["position_3d"]["median"] / 1e3, stats["position_3d"]["rms"] / 1e3,
                    stats["velocity_3d"]["median"])
        logger.info("median RTN residual (km): R %.4g, T %.4g, N %.4g; median conditional position D^2 %.4g",
                    *(stats[c]["median"] / 1e3 for c in RTN), stats["position_mahalanobis_obs_cov"]["median"])

    first_t = df["time_utc"].min() if n_rows else None
    date = first_t.strftime("%Y%m%d") if first_t is not None else "unknown"
    out = args.output_dir / date / f"NORAD_{s5.norad_id}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    log_path_holder.append(out / "logs" / "stage6.log")
    paths = write_tables(out_rows, tle_rows, out)

    flag_counts: Dict[str, int] = {}
    for r in out_rows:
        for f in r["quality_flags"]:
            flag_counts[f] = flag_counts.get(f, 0) + 1
    t_used = s5.times_s[idx] if len(idx) else np.zeros(0)
    summary = {
        "residual_definition": "observed_minus_tle",
        "n_epochs": n_rows, "n_compared": int(len(idx)), "n_valid": int(sum(r["valid"] for r in out_rows)),
        "time_span_utc": [_iso(t_used[0]), _iso(t_used[-1])] if len(t_used) else None,
        "time_span_s": float(np.ptp(t_used)) if len(t_used) else 0.0,
        "selected_tle": tle.to_dict(),
        "tle_epoch_utc": format_utc(tle.epoch_utc, 6),
        "median_tle_age_days": float(np.median(ages[idx])) if len(idx) else None,
        "max_abs_tle_age_days": max_age if len(idx) else None,
        "precision_validation": precise, "precision_validation_conditions": conditions,
        "statistics": stats,
        "mahalanobis_note": "squared distances per epoch; *_obs_cov conditional on the TLE being fixed; epochs "
                            "are correlated, so no chi-square over epochs is formed",
        "ensemble": ens_info,
        "run_flags": q.ordered(run_flags) if run_flags else [],
        "flag_counts": flag_counts,
        "light_time_note": LIGHT_TIME_NOTE,
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=2) + "\n", encoding="utf-8")

    plots: List[str] = []
    if not args.no_plots and res is not None:
        try:
            from tart_gnss.stage6.plotting import make_plots

            rk = "residuals_parquet" if paths["residuals_parquet"] else "residuals_csv"
            tk = "tle_trajectory_parquet" if paths["tle_trajectory_parquet"] else "tle_trajectory_csv"
            plots = make_plots(read_table(paths[rk]), read_table(paths[tk]), out / "diagnostics",
                               f"NORAD {s5.norad_id}")
        except ImportError as exc:
            logger.warning("plots skipped: %s", exc)

    from tart_gnss.stage4.loader import sha256_file

    unchanged = (sha256_file(s5.trajectory_path) == s5.trajectory_sha256
                 and sha256_file(s5.manifest_path) == s5.sha256)
    if not unchanged:
        logger.error("the Stage-5 trajectory changed during the run")
    import sgp4

    doc = {
        "schema_version": STAGE6_SCHEMA_VERSION,
        "stage": 6,
        "input_stage5_manifest": str(args.stage5_manifest),
        "input_stage5_sha256": s5.sha256,
        "input_trajectory": str(s5.trajectory_path),
        "input_trajectory_sha256": s5.trajectory_sha256,
        "stage5_inputs_unchanged": unchanged,
        "norad_id": s5.norad_id,
        "tle_reference": {**tle.to_dict(), "sgp4_implementation": f"sgp4 {sgp4.__version__} via Skyfield "
                          "EarthSatellite (WGS72)", "propagation": "direct SGP4 at every Stage-5 epoch (no "
                          "interpolation)", "age_convention": "age = t_obs - t_TLE epoch; > 0 forward, < 0 backward",
                          "age_thresholds_days": {"warning": args.tle_age_warning_days,
                                                  "stale": args.tle_age_stale_days}},
        "tle_ensemble": ens_info,
        "frames": FRAMES,
        "residual_definition": "observed_minus_tle",
        "rtn_convention": "RTN = RIC: R = r/|r|, N = (r x v)/|r x v|, T = N x R, from the TLE-predicted state; "
                          "velocity residuals are inertial dv projected on these axes",
        "quality_policy": {"policy": args.quality_policy, "stage5_blocking": sorted(q.STAGE5_BLOCKING),
                           "stage5_advisory": sorted(q.STAGE5_ADVISORY), "n_compared": int(len(idx)),
                           "n_rows": n_rows, "rejected": {r: reasons.count(r) for r in set(reasons) if r}},
        "thresholds": {"position_residual_warning_km": args.position_residual_warning_km,
                       "velocity_residual_warning_m_s": args.velocity_residual_warning_m_s,
                       "conditional_position_d2": q.CONDITIONAL_D2_LIMIT},
        "precision_validation": {"value": precise, "conditions": conditions,
                                 "allow_stale_reference": args.allow_stale_reference,
                                 "meaning": "the comparison meets the pipeline's minimum conditions for precision "
                                            "analysis; it does not mean the TLE is exact"},
        "not_done": ["TLE fitted or shifted", "time offset subtracted", "RTN bias subtracted",
                     "Stage-5 covariance rescaled", "Stage-5 trajectory modified", "anomaly attribution"],
        "light_time_note": LIGHT_TIME_NOTE,
        "summary_file": str(out / "summary.json"),
        "software": {**_software(), **git, "tart_gnss": __version__},
        "provenance": {"command_line": list(command_line),
                       "run_timestamp_utc": started.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z",
                       "working_directory": os.getcwd(),
                       "stage5_run": {"command_line": (s5.manifest.get("provenance") or {}).get("command_line"),
                                      "git_commit": (s5.manifest.get("software") or {}).get("git_commit")}},
        "outputs": {**paths, "summary": str(out / "summary.json"), "diagnostics": plots,
                    "log": str(out / "logs" / "stage6.log"), "manifest": str(out / "stage6_manifest.json")},
    }
    write_manifest(doc, out / "stage6_manifest.json")
    logger.info("outputs: %s", out)
    for k, v in paths.items():
        if v:
            logger.info("  %s: %s", k, v)
    if precise:
        return EXIT_OK
    only_stale = [k for k, v in conditions.items() if not v] == ["tle_age_within_stale_limit"]
    if args.allow_stale_reference and only_stale:
        logger.warning("stale reference allowed (--allow-stale-reference): descriptive comparison only")
        return EXIT_OK
    logger.warning("comparison completed but not valid for precision analysis")
    return EXIT_NOT_PRECISION


def _clean(obj):
    from tart_gnss.stage4.output import _clean as c

    return c(obj)


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-validate``."""
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
        return run(args, ["tart-validate", *argv], log_path)
    except (InvalidStage5ManifestError, MissingTrajectoryError, TLENoradMismatchError, TLEReferenceError) as exc:
        label = "TLE_NORAD_MISMATCH: " if isinstance(exc, TLENoradMismatchError) else ""
        logger.error("%s%s", label, exc)
        print(f"error: {label}{exc}", file=sys.stderr)
        return EXIT_USAGE
    except (TLEPropagationError, FrameTransformError) as exc:
        logger.error("TLE_FRAME_FAILURE: %s", exc)
        print(f"error: TLE_FRAME_FAILURE: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    except (Stage6Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)
        if log_path:
            log_path[0].write_text("\n".join(collector.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
