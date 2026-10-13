"""``tart-track``: Stage-5 TLE-independent trajectory estimation from Stage-4 positions.

Reads ``stage4_manifest.json`` and its position table, applies the Stage-4
quality policy, validates every covariance, rotates positions and full
covariances from ITRS/ECEF into GCRS, initialises a position/velocity state
from the first measurements only, runs a forward extended Kalman filter and a
backward RTS smoother under the chosen dynamics, and writes the trajectory,
innovations, ``stage5_manifest.json``, diagnostic plots and a log under
``<output-dir>/<YYYYMMDD>/NORAD_<id>/``.

No TLE enters: the TLE block that Stage 4 copied for Stage 6 is not read.

Exit codes: 0 successful scientific run (>= 3 assimilated measurements and at
least one valid trajectory row); 3 completed but insufficient valid trajectory
estimates (including filter divergence); 2 invalid input or configuration;
1 unexpected failure.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from tart_gnss import __version__
from tart_gnss.stage5 import quality as q
from tart_gnss.stage5.dynamics import (
    ATOL,
    DEFAULT_PROCESS_ACCEL_SIGMA,
    DYNAMICS,
    GM_EARTH_M3_S2,
    INTEGRATOR,
    J2_EARTH,
    R_EARTH_M,
    RTOL,
    acceleration,
)
from tart_gnss.stage5.errors import (
    FilterDivergenceError,
    FrameTransformError,
    InitialisationError,
    InvalidCovarianceError,
    InvalidStage4ManifestError,
    MissingPositionDataError,
    SmootherError,
    Stage5Error,
)
from tart_gnss.stage5.filter import ROBUST_UPDATES, FilterConfig, FilterResult, ekf_filter
from tart_gnss.stage5.frames import FRAME_DESCRIPTION, earth_pole_gcrs, gcrs_to_ecef_positions
from tart_gnss.stage5.initialisation import LAST_RESORT_VELOCITY_SIGMA_M_S
from tart_gnss.stage5.innovations import BURN_IN, dynamics_diagnostics, summarise, whitened
from tart_gnss.stage5.loader import MeasurementSet, build_measurements, load_stage4_positions
from tart_gnss.stage5.output import (
    STAGE5_SCHEMA_VERSION,
    covariance_fields,
    read_table,
    write_manifest,
    write_tables,
)
from tart_gnss.stage5.smoother import SmoothedResult, rts_smooth

logger = logging.getLogger("tart_gnss.stage5")

EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_NO_VALID = 0, 1, 2, 3
MJD_UNIX_EPOCH = 40587.0


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``stage5.log``."""

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
    """Argument parser for ``tart-track``."""
    p = argparse.ArgumentParser(
        prog="tart-track",
        description="Stage 5: estimate a continuous position/velocity trajectory (GCRS, full 6x6 covariance) "
                    "from the Stage-4 3-D positions of tart-triangulate with an extended Kalman filter and RTS "
                    "smoother. TLEs are never used.",
    )
    p.add_argument("--stage4-manifest", required=True, type=Path, help="stage4_manifest.json from tart-triangulate")
    p.add_argument("--output-dir", type=Path, default=Path("results/stage5"),
                   help="output root (default results/stage5)")
    p.add_argument("--dynamics", choices=DYNAMICS, default="two-body",
                   help="two-body (default), constant-velocity (diagnostic baseline) or j2")
    p.add_argument("--process-accel-sigma", type=float, default=None,
                   help="white process acceleration sigma, m/s^2 (default per model: "
                        + ", ".join(f"{k} {v:g}" for k, v in DEFAULT_PROCESS_ACCEL_SIGMA.items()) + ")")
    p.add_argument("--initial-velocity-window", type=float, default=30.0,
                   help="seconds of Stage-4 data used for the initial velocity fit (default 30)")
    p.add_argument("--innovation-gate-probability", type=float, default=0.9973,
                   help="chi^2_3 probability of the innovation gate (default 0.9973)")
    p.add_argument("--quality-policy", choices=q.POLICIES, default="strict",
                   help="strict (default): Stage-4 valid epochs only; permissive: also POOR_GEOMETRY and "
                        "STAGE3_QUALITY_WARNING epochs")
    p.add_argument("--max-gap-s", type=float, default=120.0,
                   help="LARGE_TIME_GAP when more than this many seconds pass without an assimilated "
                        "measurement (default 120)")
    p.add_argument("--robust-update", choices=ROBUST_UPDATES, default="none",
                   help="none (default): gated measurements are not assimilated; inflate: assimilated with "
                        "covariance inflated to the gate")
    p.add_argument("--no-smoother", action="store_true", help="skip the RTS smoother (filtered states only)")
    p.add_argument("--no-plots", action="store_true", help="skip diagnostic plots")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    return p


def _iso(t_s: float) -> str:
    return datetime.fromtimestamp(t_s, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _median(values) -> float:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else float("nan")


def _validate_args(args) -> None:
    if args.process_accel_sigma is not None and not args.process_accel_sigma >= 0:
        raise InvalidStage4ManifestError("--process-accel-sigma must be >= 0")
    if not 0 < args.innovation_gate_probability < 1:
        raise InvalidStage4ManifestError("--innovation-gate-probability must be in (0, 1)")
    if not args.initial_velocity_window > 0:
        raise InvalidStage4ManifestError("--initial-velocity-window must be > 0")
    if not args.max_gap_s > 0:
        raise InvalidStage4ManifestError("--max-gap-s must be > 0")


@dataclass
class TrackingResult:
    """Everything one tracking run produces (tables are lists of row dicts)."""

    measurements: MeasurementSet
    config: FilterConfig
    filtered: Optional[FilterResult]
    smoothed: Optional[SmoothedResult]
    first_step: Optional[int]
    summary: Dict[str, object]
    trajectory: List[Dict[str, object]]
    innovations: List[Dict[str, object]]
    run_flags: List[str]
    smoother_error: Optional[str]


def _divergence_steps(fr: FilterResult) -> set:
    """Steps inside runs of >= DIVERGENCE_RUN consecutive gated (attempted) measurements."""
    out, run = set(), []
    for k in np.flatnonzero(fr.attempted):
        if fr.outlier[k] and not fr.accepted[k]:
            run.append(int(k))
        else:
            if len(run) >= q.DIVERGENCE_RUN:
                out.update(run)
            run = []
    if len(run) >= q.DIVERGENCE_RUN:
        out.update(run)
    return out


def track(mset: MeasurementSet, config: FilterConfig, smoother: bool = True, norad_id=None,
          policy: str = "strict") -> TrackingResult:
    """Filter and smooth a measurement set; build the trajectory and innovation rows."""
    meas = mset.measurements
    n = len(meas)
    usable = [k for k, m in enumerate(meas) if m.used and m.finite]
    run_flags: List[str] = []
    fr = sm = None
    k0 = usable[0] if usable else None
    smoother_error = None
    summary: Dict[str, object] = {}
    if len(usable) < q.MIN_MEASUREMENTS:
        run_flags.append(q.INSUFFICIENT_MEASUREMENTS)
    else:
        if config.dynamics == "j2" and config.pole_gcrs is None:
            config.pole_gcrs = earth_pole_gcrs(meas[k0].time_s)
        fr = ekf_filter(meas[k0:], config)
        summary = summarise(fr)
        run_flags += summary["flags"]
        if smoother:
            try:
                sm = rts_smooth(fr)
            except SmootherError as exc:
                smoother_error = str(exc)
                run_flags.append(q.SMOOTHER_FAILED)
                logger.error("RTS smoother failed: %s", exc)
    diverged = _divergence_steps(fr) if fr is not None else set()
    last_acc = (k0 + int(np.flatnonzero(fr.accepted)[-1])) if fr is not None and fr.accepted.any() else None
    init_end = None
    if fr is not None:
        fit_times = sorted(m.time_s for m in meas[k0:] if m.used and m.finite)[:fr.initial.n_points]
        init_end = fit_times[-1] if fit_times else None
    white = whitened(fr.innovation, fr.innovation_cov) if fr is not None else None
    rows4 = mset.rows
    mats = mset.matrices
    trajectory, innovations = [], []
    for k, m in enumerate(meas):
        row4 = rows4.iloc[k]
        flags: List[str] = []
        rec: Dict[str, object] = {
            "time_utc": _iso(m.time_s), "time_mjd": m.time_s / 86400.0 + MJD_UNIX_EPOCH, "norad_id": norad_id,
            "dynamics": config.dynamics,
            "stage4_x_ecef_m": row4["x_ecef_m"], "stage4_y_ecef_m": row4["y_ecef_m"],
            "stage4_z_ecef_m": row4["z_ecef_m"],
            "stage4_x_gcrs_m": m.position_gcrs_m[0], "stage4_y_gcrs_m": m.position_gcrs_m[1],
            "stage4_z_gcrs_m": m.position_gcrs_m[2],
            "measurement_used": bool(m.used), "measurement_accepted": False, "rejection_reason": m.reason,
            "stage4_quality_flags": list(m.stage4_flags),
        }
        for a, i in zip("xyz", range(3)):
            c = m.covariance_gcrs_m2[i, i]
            rec[f"stage4_sigma_{a}_gcrs_m"] = float(np.sqrt(c)) if np.isfinite(c) and c >= 0 else float("nan")
        if m.covariance_problem:
            flags.append(q.INVALID_STAGE4_COVARIANCE)
        if not m.used:
            flags.append(q.MEASUREMENT_REJECTED if np.all(np.isfinite(row4[["x_ecef_m", "y_ecef_m", "z_ecef_m"]]
                                                                         .to_numpy(dtype=float)))
                         else q.NO_MEASUREMENT)
        if m.used and q.has_advisory(m.stage4_flags, policy):
            flags.append(q.STAGE4_QUALITY_WARNING)
        if fr is None or k < k0:
            flags.append(q.INSUFFICIENT_MEASUREMENTS if fr is None else q.NOT_INITIALISED)
        else:
            j = k - k0
            xf, pf = fr.x_filt[j], fr.p_filt[j]
            rec.update({
                "x_gcrs_filtered_m": xf[0], "y_gcrs_filtered_m": xf[1], "z_gcrs_filtered_m": xf[2],
                "vx_gcrs_filtered_m_s": xf[3], "vy_gcrs_filtered_m_s": xf[4], "vz_gcrs_filtered_m_s": xf[5],
                "sigma_x_filtered_m": np.sqrt(pf[0, 0]), "sigma_y_filtered_m": np.sqrt(pf[1, 1]),
                "sigma_z_filtered_m": np.sqrt(pf[2, 2]),
                "measurement_accepted": bool(fr.accepted[j]), "time_since_update_s": fr.gap_s[j],
                "nis": fr.nis[j],
            })
            rec.update(covariance_fields(pf, "filtered"))
            if sm is not None:
                xs, ps = sm.x_smooth[j], sm.p_smooth[j]
                ecef = gcrs_to_ecef_positions(xs[None, :3], mats[k][None])[0]
                ce = mats[k].T @ ps[:3, :3] @ mats[k]
                rec.update({
                    "x_gcrs_smoothed_m": xs[0], "y_gcrs_smoothed_m": xs[1], "z_gcrs_smoothed_m": xs[2],
                    "vx_gcrs_smoothed_m_s": xs[3], "vy_gcrs_smoothed_m_s": xs[4], "vz_gcrs_smoothed_m_s": xs[5],
                    "x_ecef_smoothed_m": ecef[0], "y_ecef_smoothed_m": ecef[1], "z_ecef_smoothed_m": ecef[2],
                    "speed_smoothed_m_s": float(np.linalg.norm(xs[3:])),
                    "sigma_x_smoothed_m": np.sqrt(ps[0, 0]), "sigma_y_smoothed_m": np.sqrt(ps[1, 1]),
                    "sigma_z_smoothed_m": np.sqrt(ps[2, 2]), "sigma_vx_smoothed_m_s": np.sqrt(ps[3, 3]),
                    "sigma_vy_smoothed_m_s": np.sqrt(ps[4, 4]), "sigma_vz_smoothed_m_s": np.sqrt(ps[5, 5]),
                    "acceleration_m_s2": float(np.linalg.norm(acceleration(xs[:3], config.dynamics,
                                                                           config.pole_gcrs))),
                    "radius_m": float(np.linalg.norm(xs[:3])),
                    "cov_ecef_smoothed_xx": ce[0, 0], "cov_ecef_smoothed_xy": ce[0, 1],
                    "cov_ecef_smoothed_xz": ce[0, 2], "cov_ecef_smoothed_yy": ce[1, 1],
                    "cov_ecef_smoothed_yz": ce[1, 2], "cov_ecef_smoothed_zz": ce[2, 2],
                })
                rec.update(covariance_fields(ps, "smoothed"))
            if fr.outlier[j]:
                flags.append(q.INNOVATION_OUTLIER)
            if j in diverged:
                flags.append(q.FILTER_DIVERGENCE)
            if fr.gap_s[j] > config.max_gap_s:
                flags.append(q.LARGE_TIME_GAP)
            if last_acc is not None and k > last_acc:
                flags.append(q.EXTRAPOLATED)
            if fr.initial.weak and init_end is not None and m.time_s <= init_end:
                flags.append(q.INITIAL_VELOCITY_WEAK)
        flags += run_flags
        rec["quality_flags"] = q.ordered(flags)
        rec["valid"] = q.is_valid(rec["quality_flags"])
        trajectory.append(rec)
        if fr is not None and k >= k0 and np.isfinite(fr.nis[k - k0]):
            j = k - k0
            nu, s, w = fr.innovation[j], fr.innovation_cov[j], white[j]
            innovations.append({
                "time_utc": rec["time_utc"], "innovation_x_m": nu[0], "innovation_y_m": nu[1],
                "innovation_z_m": nu[2], "innovation_cov_xx": s[0, 0], "innovation_cov_xy": s[0, 1],
                "innovation_cov_xz": s[0, 2], "innovation_cov_yy": s[1, 1], "innovation_cov_yz": s[1, 2],
                "innovation_cov_zz": s[2, 2], "whitened_x": w[0], "whitened_y": w[1], "whitened_z": w[2],
                "nis": fr.nis[j], "gate_threshold": fr.gate_threshold, "accepted": bool(fr.accepted[j]),
                "attempted": bool(fr.attempted[j]), "inflation": fr.inflation[j],
                "stage4_quality_flags": list(m.stage4_flags), "stage5_quality_flags": rec["quality_flags"],
            })
    return TrackingResult(measurements=mset, config=config, filtered=fr, smoothed=sm, first_step=k0,
                          summary=summary, trajectory=trajectory, innovations=innovations,
                          run_flags=q.ordered(run_flags) if run_flags else [], smoother_error=smoother_error)


def run(args: argparse.Namespace, command_line: Sequence[str], log_path_holder: List[Path]) -> int:
    """Execute Stage 5 for parsed arguments; return the exit status."""
    git = _git_state()
    started = datetime.now(timezone.utc)
    _validate_args(args)
    logger.info("input Stage-4 manifest: %s", args.stage4_manifest)
    inp = load_stage4_positions(args.stage4_manifest)
    mset = build_measurements(inp.positions, args.quality_policy)
    meas = mset.measurements
    usable = [m for m in meas if m.used and m.finite]
    times = np.array([m.time_s for m in usable])
    logger.info("Stage-4 epochs: %d (table %s); retained after '%s' quality policy: %d (rejected: %s)",
                mset.n_input, inp.positions_path, args.quality_policy, len(usable),
                ", ".join(f"{k} {v}" for k, v in sorted(mset.rejected.items())) or "none")
    if not mset.sorted_input:
        logger.warning("Stage-4 timestamps were not monotonic; rows sorted by time")
    for note in mset.covariance_notes:
        logger.warning("invalid Stage-4 covariance %s", note)
    dts = np.diff(times)
    if len(times):
        logger.info("time span %s .. %s (%.1f s); median cadence %.3f s; maximum gap %.3f s",
                    _iso(times[0]), _iso(times[-1]), times[-1] - times[0],
                    float(np.median(dts)) if dts.size else float("nan"),
                    float(np.max(dts)) if dts.size else float("nan"))
    config = FilterConfig(dynamics=args.dynamics, process_accel_sigma_m_s2=args.process_accel_sigma,
                          gate_probability=args.innovation_gate_probability, robust_update=args.robust_update,
                          max_gap_s=args.max_gap_s, initial_velocity_window_s=args.initial_velocity_window)
    logger.info("dynamics %s; process acceleration sigma %.3g m/s^2%s", config.dynamics, config.sigma_a,
                "" if args.process_accel_sigma is not None else " (model default)")
    filter_error = None
    try:
        res = track(mset, config, smoother=not args.no_smoother, norad_id=inp.norad_id, policy=args.quality_policy)
    except FilterDivergenceError as exc:
        filter_error = str(exc)
        logger.error("filter failed: %s", exc)
        res = None
    fr = res.filtered if res is not None else None
    if fr is not None:
        init = fr.initial
        logger.info("initial position (GCRS) %s m; initial velocity %s m/s; method %s from %d point(s)%s",
                    np.array2string(init.state[:3], precision=0), np.array2string(init.state[3:], precision=1),
                    init.method, init.n_points, " [INITIAL_VELOCITY_WEAK]" if init.weak else "")
        for note in init.notes:
            logger.info("initialisation: %s", note)
        s = res.summary
        logger.info("innovation gate chi2_3(%.4f) = %.3f; rejected innovations %d of %d attempted; "
                    "median NIS %.3f, mean NIS %.3f (expected 2.37, 3)", config.gate_probability,
                    fr.gate_threshold, s["n_attempted"] - s["n_accepted"], s["n_attempted"], s["median_nis"],
                    s["mean_nis"])
        logger.info("lag-1 innovation autocorrelation %.3f (limit %.3f); innovation RMS %.4g km",
                    s["autocorrelation"]["1"], s["autocorrelation_limit_lag1"], s["innovation_rms_m"] / 1e3)
        logger.info("filter success: True; smoother success: %s", (res.smoothed is not None) if not
                    args.no_smoother else "skipped")
    elif res is not None:
        logger.warning("fewer than %d usable measurements; no trajectory estimated", q.MIN_MEASUREMENTS)
    traj = res.trajectory if res is not None else []
    innov = res.innovations if res is not None else []
    med_f = _median([np.sqrt(r.get("sigma_x_filtered_m", np.nan) ** 2 + r.get("sigma_y_filtered_m", np.nan) ** 2
                             + r.get("sigma_z_filtered_m", np.nan) ** 2) for r in traj])
    med_s = _median([np.sqrt(r.get("sigma_x_smoothed_m", np.nan) ** 2 + r.get("sigma_y_smoothed_m", np.nan) ** 2
                             + r.get("sigma_z_smoothed_m", np.nan) ** 2) for r in traj])
    logger.info("median 3-D position sigma: filtered %.4g km, smoothed %.4g km", med_f / 1e3, med_s / 1e3)

    first_t = inp.positions["time_utc"].min() if len(inp.positions) else None
    date = first_t.strftime("%Y%m%d") if first_t is not None else "unknown"
    out = args.output_dir / date / f"NORAD_{inp.norad_id}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    log_path_holder.append(out / "logs" / "stage5.log")
    paths = write_tables(traj, innov, out)

    plots: List[str] = []
    if not args.no_plots and traj and fr is not None:
        try:
            from tart_gnss.stage5.plotting import make_plots

            key = "trajectory_parquet" if paths["trajectory_parquet"] else "trajectory_csv"
            ikey = "innovations_parquet" if paths["innovations_parquet"] else "innovations_csv"
            plots = make_plots(read_table(paths[key]), read_table(paths[ikey]), out / "diagnostics",
                               f"NORAD {inp.norad_id} ({config.dynamics})")
        except ImportError as exc:
            logger.warning("plots skipped: %s", exc)

    valid = [r for r in traj if r["valid"] and np.isfinite(r.get("x_gcrs_filtered_m", np.nan))]
    n_acc = int(fr.accepted.sum()) if fr is not None else 0
    flag_counts: Dict[str, int] = {}
    for r in traj:
        for f in r["quality_flags"]:
            flag_counts[f] = flag_counts.get(f, 0) + 1
    m4 = inp.manifest
    doc = {
        "schema_version": STAGE5_SCHEMA_VERSION,
        "stage": 5,
        "input_stage4_manifest": str(args.stage4_manifest),
        "input_sha256": inp.sha256,
        "input_positions": str(inp.positions_path),
        "input_positions_sha256": inp.positions_sha256,
        "norad_id": inp.norad_id,
        "state_frame": "GCRS",
        "measurement_frame": "ITRS/ECEF",
        "frame_transformation": FRAME_DESCRIPTION,
        "units": {"position": "m", "velocity": "m/s", "covariance": "SI (m^2, m^2/s, m^2/s^2)", "time": "UTC"},
        "tle_used": False,
        "state_model": {
            "dynamics": config.dynamics,
            "state": ["x", "y", "z", "vx", "vy", "vz"],
            "gm_earth_m3_s2": GM_EARTH_M3_S2,
            "j2": J2_EARTH if config.dynamics == "j2" else None,
            "earth_radius_m": R_EARTH_M if config.dynamics == "j2" else None,
            "constants_source": "IERS Conventions (2010) Table 1.1 (GM equals WGS84, NIMA TR8350.2)",
            "integrator": INTEGRATOR if config.dynamics != "constant-velocity" else "closed form",
            "rtol": RTOL if config.dynamics != "constant-velocity" else None,
            "atol": ATOL if config.dynamics != "constant-velocity" else None,
            "state_transition": "variational equations integrated with the state" if config.dynamics !=
            "constant-velocity" else "F = [[I, dt I], [0, I]]",
            "j2_pole_gcrs": [float(v) for v in config.pole_gcrs] if config.pole_gcrs is not None else None,
        },
        "filter": {
            "type": "extended Kalman filter, Joseph-form update, Cholesky solves",
            "measurement_model": "H = [I3 0]; R = full Stage-4 covariance rotated to GCRS",
            "initialisation": {
                **(fr.initial.to_dict() if fr is not None else {}),
                "window_s": args.initial_velocity_window,
                "hierarchy": "GLS fit with full covariances (>=3 points) -> finite difference (2) -> zero "
                             f"velocity with {LAST_RESORT_VELOCITY_SIGMA_M_S:g} m/s sigma (1, flagged)",
                "prior_inflation": config.prior_inflation,
                "prior": "P0 = prior_inflation x fit covariance (window data are assimilated again by the "
                         "filter; the inflation keeps the reused information below 1/prior_inflation)",
                "weak_velocity_sigma_m_s": config.weak_velocity_sigma_m_s,
                "tle_used": False,
            },
            "success": fr is not None,
            "error": filter_error,
            "n_steps": len(fr.times_s) if fr is not None else 0,
            "n_assimilated": n_acc,
            "min_measurements": q.MIN_MEASUREMENTS,
            "max_gap_s": args.max_gap_s,
            "robust_update": args.robust_update,
        },
        "smoother": {
            "type": "Rauch-Tung-Striebel (nonlinear predictions, flow Jacobians from the filter)",
            "enabled": not args.no_smoother,
            "success": res is not None and res.smoothed is not None,
            "error": res.smoother_error if res is not None else None,
        },
        "process_noise": {
            "model": "white acceleration, Q = sigma_a^2 [[dt^4/4 I, dt^3/2 I], [dt^3/2 I, dt^2 I]]",
            "process_accel_sigma_m_s2": config.sigma_a,
            "source": "--process-accel-sigma" if args.process_accel_sigma is not None else
            "model default (physics: size of the unmodelled accelerations; not tuned against TLEs)",
            "defaults_m_s2": DEFAULT_PROCESS_ACCEL_SIGMA,
        },
        "innovation_gate": {
            "probability": config.gate_probability,
            "chi2_dof": 3,
            "threshold": fr.gate_threshold if fr is not None else None,
            "gated_measurements": "kept in innovations/trajectory with INNOVATION_OUTLIER; "
                                  + ("assimilated with inflated covariance" if args.robust_update == "inflate"
                                     else "not assimilated"),
            "statistics_burn_in": BURN_IN,
        },
        "innovation_statistics": res.summary if res is not None else {},
        "dynamics_diagnostics": dynamics_diagnostics(res.summary) if res is not None and res.summary else {},
        "quality_policy": {
            "policy": args.quality_policy,
            "always_rejected": sorted(q.STAGE4_FATAL_FLAGS) + ["fit_success false", "non-finite position",
                                                                "invalid covariance", "duplicate timestamp"],
            "rejected_if_strict": sorted(q.STAGE4_STRICT_FLAGS),
            "advisory": sorted(q.STAGE4_ADVISORY_FLAGS),
            "covariance_checks": {"asymmetry_tol": q.COVARIANCE_ASYMMETRY_TOL,
                                  "max_condition_number": q.MAX_COVARIANCE_CONDITION},
            "n_stage4_epochs": mset.n_input, "n_used": len(usable), "rejected": mset.rejected,
            "input_time_ordered": mset.sorted_input,
        },
        "blocking_flags": sorted(q.BLOCKING_FLAGS),
        "summary": {
            "n_rows": len(traj), "n_valid": len(valid), "n_assimilated": n_acc,
            "run_flags": res.run_flags if res is not None else [], "flag_counts": flag_counts,
            "time_span_s": float(times[-1] - times[0]) if len(times) else 0.0,
            "median_cadence_s": float(np.median(dts)) if dts.size else None,
            "max_gap_s": float(np.max(dts)) if dts.size else None,
            "median_position_sigma_filtered_m": med_f, "median_position_sigma_smoothed_m": med_s,
        },
        "stage4_reference": {
            "note": "Stage-4 provenance; the satellite/TLE block of the Stage-4 manifest is not read by Stage 5",
            "command_line": (m4.get("provenance") or {}).get("command_line"),
            "git_commit": (m4.get("software") or {}).get("git_commit"),
            "input_stage3_manifest": m4.get("input_stage3_manifest"),
        },
        "software": {**_software(), **git, "tart_gnss": __version__},
        "provenance": {
            "command_line": list(command_line),
            "run_timestamp_utc": started.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z",
            "working_directory": os.getcwd(),
        },
        "outputs": {**paths, "diagnostics": plots, "log": str(out / "logs" / "stage5.log"),
                    "manifest": str(out / "stage5_manifest.json")},
    }
    write_manifest(doc, out / "stage5_manifest.json")
    logger.info("outputs: %s", out)
    for k, v in paths.items():
        if v:
            logger.info("  %s: %s", k, v)
    if fr is None or n_acc < q.MIN_MEASUREMENTS or not valid:
        logger.warning("insufficient valid trajectory estimates (assimilated %d, valid rows %d; flags: %s)",
                       n_acc, len(valid), ", ".join(f"{k} {v}" for k, v in sorted(flag_counts.items())) or "none")
        return EXIT_NO_VALID
    return EXIT_OK


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-track``."""
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
        return run(args, ["tart-track", *argv], log_path)
    except (InvalidStage4ManifestError, MissingPositionDataError, InvalidCovarianceError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (InitialisationError, FrameTransformError, Stage5Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)
        if log_path:
            log_path[0].write_text("\n".join(collector.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
