"""``tart-diagnose``: Stage-7 sensor-network diagnosis and systematic-effect evidence.

Reads ``stage6_manifest.json``, loads and verifies the Stage-5/4/3/2 chain,
aligns all outputs by time and station, recomputes leave-one-station-out and
pairwise solutions with the Stage-4 solver, compares every station direction
with the Stage-5 trajectory, tests calibration / timing / geometry /
trajectory-model correlations, and assigns evidence tags, evidence levels
(0-5) and a descriptive confidence. Writes ``diagnosis``,
``station_influence``, ``pairwise``, ``correlations.csv``, ``summary.json``,
``stage7_manifest.json``, plots and a log under
``<output-dir>/<YYYYMMDD>/NORAD_<id>/``. No upstream output is modified and
no physical cause is ever asserted.

Exit codes: 0 diagnosis completed; 3 completed but fewer than 3 diagnosable
epochs (insufficient data, e.g. an invalid upstream chain); 2 invalid input
or configuration (including a broken provenance chain); 1 unexpected failure.
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
from tart_gnss.stage7 import quality as q
from tart_gnss.stage7.alignment import DEFAULT_TOLERANCE_S, build_epoch_table, build_station_table
from tart_gnss.stage7.attribution import EvidenceRules, assign_evidence, ordered
from tart_gnss.stage7.correlations import FAMILY_TAG, compute_correlations, strongest, supported_tags
from tart_gnss.stage7.errors import (
    AlignmentError,
    InvalidStage6ManifestError,
    ProvenanceChainError,
    Stage7Error,
)
from tart_gnss.stage7.influence import (
    DEFAULT_MAX_CONDITION,
    compute_leave_one_out,
    compute_station_direction_residuals,
    station_bias,
)
from tart_gnss.stage7.loader import load_diagnostic_chain
from tart_gnss.stage7.output import STAGE7_SCHEMA_VERSION, read_table, write_json, write_manifest, write_tables
from tart_gnss.stage7.statistics import METHODS

logger = logging.getLogger("tart_gnss.stage7")
EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_INSUFFICIENT = 0, 1, 2, 3


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``stage7.log``."""

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
    """Argument parser for ``tart-diagnose``."""
    p = argparse.ArgumentParser(
        prog="tart-diagnose",
        description="Stage 7: diagnose where an observed - TLE discrepancy comes from (station, instrument, "
                    "timing, geometry, trajectory model, reference) with transparent evidence rules. Reports "
                    "evidence patterns, never physical causes.",
    )
    p.add_argument("--stage6-manifest", required=True, type=Path, help="stage6_manifest.json from tart-validate")
    p.add_argument("--output-dir", type=Path, default=Path("results/stage7"), help="output root (default results/stage7)")
    p.add_argument("--quality-policy", choices=q.POLICIES, default="strict",
                   help="strict (default): diagnose only valid, fully linked epochs; hash mismatches are fatal")
    p.add_argument("--leave-one-out", action=argparse.BooleanOptionalAction, default=True,
                   help="leave-one-station-out re-triangulation (default on)")
    p.add_argument("--pairwise", action=argparse.BooleanOptionalAction, default=True,
                   help="pairwise two-station diagnostics (default on)")
    p.add_argument("--correlation-method", choices=METHODS, default="spearman", help="default spearman")
    p.add_argument("--min-samples", type=int, default=8,
                   help="minimum samples for a correlation / sample support (default 8)")
    p.add_argument("--min-correlation", type=float, default=0.5,
                   help="minimum |time-detrended correlation| to support a tag (default 0.5)")
    p.add_argument("--anomaly-threshold", type=float, default=None,
                   help="Stage-6 squared-Mahalanobis anomaly threshold (default: chi2_3 at 1 - 0.0027/n_epochs)")
    p.add_argument("--min-persistent-epochs", type=int, default=5,
                   help="consecutive strong anomalous epochs needed for persistence / Level 5 (default 5)")
    p.add_argument("--alignment-tolerance-s", type=float, default=DEFAULT_TOLERANCE_S,
                   help=f"time tolerance for joining stage outputs (default {DEFAULT_TOLERANCE_S})")
    p.add_argument("--no-plots", action="store_true", help="skip diagnostic plots")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    return p


def _validate_args(args) -> None:
    if args.min_samples < 3:
        raise InvalidStage6ManifestError("--min-samples must be >= 3")
    if not 0 < args.min_correlation < 1:
        raise InvalidStage6ManifestError("--min-correlation must be in (0, 1)")
    if args.anomaly_threshold is not None and not args.anomaly_threshold > 0:
        raise InvalidStage6ManifestError("--anomaly-threshold must be > 0")
    if args.min_persistent_epochs < 1:
        raise InvalidStage6ManifestError("--min-persistent-epochs must be >= 1")
    if not args.alignment_tolerance_s > 0:
        raise InvalidStage6ManifestError("--alignment-tolerance-s must be > 0")


def _median(v) -> float:
    a = np.asarray([x for x in v if x is not None], dtype=float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else float("nan")


def _mad(v) -> float:
    a = np.asarray(v, dtype=float)
    a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a)))) if a.size else float("nan")


def predicted_geometry(ep, rows, stations):
    """Noise-free geometry of one epoch: condition number and minimum crossing angle of the sight lines from
    the contributing stations to the Stage-5 position (the measured lines carry measurement noise, which would
    correlate spuriously with the residuals)."""
    from tart_gnss.stage4.triangulation import geometry_metrics

    x = ep["stage5_state_gcrs"]
    sids = [k for k in rows["station_id"] if k in stations] if len(rows) else []
    if len(sids) < 2 or not np.all(np.isfinite(x)):
        return float("nan"), float("nan")
    r = ep["ecef_to_gcrs"].T @ x[:3]
    s = np.array([(r - stations[k].ecef_m) / np.linalg.norm(r - stations[k].ecef_m) for k in sids])
    cond, _, amin, _ = geometry_metrics(s)
    return cond, amin


def station_summaries(inp, stations_tbl, influence_df, picks_info, corr_tags) -> Dict[str, Dict[str, object]]:
    """One descriptive summary per station over the diagnosed epochs (no combined health score).

    ``usable_fraction`` is the fraction of the station's contributions with a
    direction residual against the Stage-5 trajectory.
    """
    out = {}
    m3 = inp.manifests.get(3) or {}
    bias = station_bias(stations_tbl)
    for sid in sorted(inp.stations):
        g = stations_tbl[stations_tbl["station_id"] == sid] if len(stations_tbl) else stations_tbl
        li = influence_df[influence_df["station_id"] == sid] if len(influence_df) else influence_df
        n_ant = ((m3.get("stations") or {}).get(sid) or {}).get("n_antennas")
        zg = g["zero_gain_count"].to_numpy(float) if len(g) else np.zeros(0)
        tags = []
        if sid in picks_info["station_specific"]:
            tags.append("STATION_SPECIFIC")
        for t, tests in corr_tags.items():
            if any(x["scope"] == "station" and x["station_id"] == sid for x in tests):
                tags.append(t)
        health = []
        n = int(len(g))
        usable = float(np.mean(np.isfinite(g["direction_residual_angle_deg"].to_numpy(float)))) if n else 0.0
        if n and usable < 0.5:
            health.append("LOW_USABLE_FRACTION")
        if not inp.calibration_metadata_available:
            health.append(q.CALIBRATION_METADATA_MISSING)
        if n and g["time_mismatch"].any():
            health.append("TIME_MISMATCH_SEEN")
        if inp.stations[sid].calibration_warnings:
            health.append("CALIBRATION_WARNING")
        if picks_info["station_picks"].get(sid):
            health.append("PICKED_BY_STATION_TEST")
        out[sid] = {
            "station_id": sid, "n_epochs": n, "usable_fraction": usable,
            "median_direction_residual_deg": _median(g["direction_residual_angle_deg"]) if n else float("nan"),
            "median_line_residual_m": _median(g["line_residual_m"]) if n else float("nan"),
            "median_leave_one_out_shift_m": _median(li["loo_shift_3d_m"]) if len(li) else float("nan"),
            "median_leave_one_out_mahalanobis": _median(li["loo_mahalanobis"]) if len(li) else float("nan"),
            "n_leave_one_out_epochs": int(len(li)),
            "zero_gain_fraction": float(np.nanmedian(zg) / n_ant) if n_ant and np.any(np.isfinite(zg)) else float("nan"),
            "median_flagged_fraction": _median(g["flagged_fraction"]) if n and "flagged_fraction" in g else float("nan"),
            "calibration_change_count": int(np.nanmax(g["gain_solution_index"])) if n and np.any(
                np.isfinite(g["gain_solution_index"].to_numpy(float))) else None,
            "median_snr": _median(g["snr"]) if n and "snr" in g else float("nan"),
            "median_image_visibility_offset_deg": _median(g["vis_image_separation_deg"]) if n and
            "vis_image_separation_deg" in g else float("nan"),
            "timing_offset_median_s": _median(g["station_time_offset_s"]) if n else float("nan"),
            "timing_offset_mad_s": _mad(g["station_time_offset_s"]) if n else float("nan"),
            "direction_bias": bias.get(sid, {}),
            "station_test_picks": int(picks_info["station_picks"].get(sid, 0)),
            "evidence_tags": tags, "health_flags": health,
        }
    return out


def run(args: argparse.Namespace, command_line: Sequence[str], log_path_holder: List[Path]) -> int:
    """Execute Stage 7 for parsed arguments; return the exit status."""
    import pandas as pd

    git = _git_state()
    started = datetime.now(timezone.utc)
    _validate_args(args)
    logger.info("input Stage-6 manifest: %s", args.stage6_manifest)
    inp = load_diagnostic_chain(args.stage6_manifest, args.quality_policy)
    for lk in inp.links:
        logger.info("provenance Stage %d: %s (%s)", lk.stage, lk.path or lk.recorded, lk.status)
    logger.info("NORAD %s; stations %d (%s)", inp.norad_id, len(inp.stations), ", ".join(sorted(inp.stations)))
    epochs = build_epoch_table(inp, args.alignment_tolerance_s)
    if not len(epochs):
        raise AlignmentError("the Stage-6 residual table is empty")
    epochs["diagnosed"] = [q.epoch_used(r, args.quality_policy) for _, r in epochs.iterrows()]
    st = build_station_table(inp, epochs, args.alignment_tolerance_s)
    if len(st):
        st = st.join(pd.DataFrame(compute_station_direction_residuals(epochs, st, inp.stations), index=st.index))
        st["abs_station_time_offset_s"] = st["station_time_offset_s"].abs()
    n_used = int(epochs["diagnosed"].sum())
    logger.info("epochs: %d Stage-6 rows; %d diagnosed under '%s' policy (links: Stage 5 %d, Stage 4 %d)",
                len(epochs), n_used, args.quality_policy, int(epochs["link_stage5"].sum()),
                int(epochs["link_stage4"].sum()))

    used_ep = epochs[epochs["diagnosed"]]
    st_used = st[st["epoch_index"].isin(used_ep.index)] if len(st) else st
    infl, pairs, repro = compute_leave_one_out(used_ep, st_used, inp.stations, args.leave_one_out, args.pairwise,
                                               DEFAULT_MAX_CONDITION)
    infl_df = pd.DataFrame(infl)
    pairs_df = pd.DataFrame(pairs)
    logger.info("leave-one-out solutions: %d; pairwise solutions: %d; full-solution reproduction vs Stage 4: "
                "max %.3g m", len(infl_df), len(pairs_df), repro)

    # Epoch / station quantities used by the correlations.
    if len(infl_df):
        st = st.merge(infl_df[["epoch_index", "station_id", "loo_shift_3d_m", "loo_mahalanobis", "loo_dof",
                               "n_stations_remaining",
                               "loo_x_ecef_m", "loo_y_ecef_m", "loo_z_ecef_m", "loo_shift_radial_m",
                               "loo_shift_along_track_m", "loo_shift_cross_track_m", "loo_condition_number",
                               "loo_min_crossing_angle_deg", "degraded_geometry", "loo_d2_tle"]],
                      on=["epoch_index", "station_id"], how="left")
    epochs["abs_dr_along_track_m"] = epochs["dr_along_track_m"].abs()
    pg = [predicted_geometry(ep, st[st["epoch_index"] == ei] if len(st) else st, inp.stations)
          for ei, ep in epochs.iterrows()]
    epochs["predicted_condition_number"] = [g[0] for g in pg]
    epochs["predicted_min_crossing_angle_deg"] = [g[1] for g in pg]
    epochs["predicted_inverse_min_crossing"] = 1.0 / np.sin(np.radians(epochs["predicted_min_crossing_angle_deg"]))
    ext = {"max_loo_shift_m": [], "max_abs_station_time_offset_s": [], "mean_station_time_offset_s": [],
           "max_station_line_residual_m": [], "max_station_direction_residual_deg": [],
           "max_leave_one_out_mahalanobis": [], "most_influential_station": [], "max_pair_separation_z": []}
    for ei in range(len(epochs)):
        g = st[st["epoch_index"] == ei] if len(st) else st
        li = infl_df[infl_df["epoch_index"] == ei] if len(infl_df) else infl_df
        pp = pairs_df[pairs_df["epoch_index"] == ei] if len(pairs_df) else pairs_df
        ext["max_loo_shift_m"].append(float(li["loo_shift_3d_m"].max()) if len(li) else np.nan)
        ext["max_leave_one_out_mahalanobis"].append(float(li["loo_mahalanobis"].max()) if len(li) else np.nan)
        if len(li):
            key = li["loo_mahalanobis"].fillna(-1) if li["loo_mahalanobis"].notna().any() else li["loo_shift_3d_m"]
            ext["most_influential_station"].append(li.loc[key.idxmax(), "station_id"])
        else:
            ext["most_influential_station"].append("")
        off = g["station_time_offset_s"].to_numpy(float) if len(g) else np.zeros(0)
        ext["max_abs_station_time_offset_s"].append(float(np.nanmax(np.abs(off))) if off.size else np.nan)
        ext["mean_station_time_offset_s"].append(float(np.nanmean(off)) if off.size else np.nan)
        ext["max_station_line_residual_m"].append(float(g["line_residual_m"].max()) if len(g) else np.nan)
        ext["max_station_direction_residual_deg"].append(float(g["direction_residual_angle_deg"].max())
                                                         if len(g) else np.nan)
        ext["max_pair_separation_z"].append(float(pp["pair_separation_z"].max()) if len(pp) else np.nan)
    for k, v in ext.items():
        epochs[k] = v

    corr = compute_correlations(st[st["epoch_index"].isin(used_ep.index)] if len(st) else st,
                                epochs[epochs["diagnosed"]], args.correlation_method, args.min_samples,
                                args.min_correlation)
    corr_tags = supported_tags(corr)
    for fam in FAMILY_TAG:
        s = strongest(corr, fam)
        if s:
            logger.info("strongest %s correlation: %s vs %s (%s) partial %.3f, Holm p %.3g, n %d -> %s", fam,
                        s["x_variable"], s["y_variable"], s["station_id"] or "network", s["partial_coefficient"],
                        s["p_value_holm"], s["n_samples"], s["interpretation_tag"])

    m5 = inp.manifests[5]
    s5_flags = list((m5.get("summary") or {}).get("run_flags") or [])
    s5_resets = int((m5.get("filter") or {}).get("n_resets") or 0)
    precision = bool(inp.summary6.get("precision_validation", inp.manifests[6].get("precision_validation", {})
                                      .get("value", False)))
    links_ok = all(lk.status in ("verified", "input") for lk in inp.links) and inp.link(2) is not None
    data_complete = links_ok and inp.calibration_metadata_available and (
        not len(st) or bool(st["link_stage3"].all()))
    rules = EvidenceRules(anomaly_threshold=args.anomaly_threshold, min_persistent_epochs=args.min_persistent_epochs,
                          min_samples=args.min_samples, min_correlation=args.min_correlation)
    diag, info = assign_evidence(epochs, st, infl_df, pairs_df, corr_tags, s5_flags, s5_resets, precision,
                                 inp.calibration_metadata_available, data_complete, rules)
    for k in ("statistical_anomaly", "strong", "persistent", "station_pick", "station_pick_method",
              "anomaly_statistic", "anomaly_threshold", "diagnostic_tags", "evidence_level", "diagnostic_confidence"):
        epochs[k] = [d[k] for d in diag]
    # Persistence metrics on the Stage-6 anomaly runs.
    from tart_gnss.stage7.statistics import runs

    rl, rd = np.zeros(len(epochs), int), np.zeros(len(epochs))
    for a, b, n, d in runs(epochs["statistical_anomaly"].to_numpy(bool), epochs["time_s"].to_numpy(float)):
        rl[a:b + 1], rd[a:b + 1] = n, d
    epochs["anomaly_run_length"], epochs["anomaly_duration_s"] = rl, rd

    qflags_run = []
    if not inp.calibration_metadata_available:
        qflags_run.append(q.CALIBRATION_METADATA_MISSING)
    if any(lk.status not in ("verified", "input") for lk in inp.links):
        qflags_run.append(q.PROVENANCE_UNVERIFIED)
    qf = []
    for ei, ep in epochs.iterrows():
        f = list(qflags_run)
        if not ep["diagnosed"]:
            f.append(q.NOT_DIAGNOSED)
        else:
            if len(st) and not st.loc[st["epoch_index"] == ei, "link_stage3"].all():
                f.append(q.STAGE3_LINK_MISSING)
            if not (len(infl_df) and (infl_df["epoch_index"] == ei).any()):
                f.append(q.SUBSET_TEST_UNAVAILABLE)
        qf.append(f)
    epochs["quality_flags"] = qf

    lv = pd.to_numeric(epochs.loc[epochs["diagnosed"], "evidence_level"], errors="coerce")
    overall_level = int(lv.max()) if len(lv) and lv.notna().any() else 0
    overall_tags = set()
    for t in epochs.loc[epochs["diagnosed"], "diagnostic_tags"]:
        overall_tags |= set(t)
    overall_tags |= set(corr_tags)
    if info["station_specific"]:
        overall_tags.add("STATION_SPECIFIC")
    if info["trajectory_model_run_level"]:
        overall_tags.add("TRAJECTORY_MODEL_CORRELATED")
    if not precision:
        overall_tags.add("REFERENCE_LIMITED")
    if info["n_anomalous"]:
        overall_tags.discard("NORMAL")
    elif n_used:
        overall_tags.add("NORMAL")
    at = epochs.loc[epochs["diagnosed"] & (pd.to_numeric(epochs["evidence_level"], errors="coerce") == overall_level),
                    "diagnostic_confidence"]
    overall_conf = at.mode().iloc[0] if len(at) else "LOW"
    infl_station = ""
    n_loo_epochs = int(infl_df["epoch_index"].nunique()) if len(infl_df) else 0
    if n_loo_epochs >= rules.station_min_epochs:  # never named from one or two epochs
        med = infl_df.groupby("station_id")["loo_mahalanobis"].median()
        if med.notna().any():
            infl_station = str(med.idxmax())
        else:
            infl_station = str(infl_df.groupby("station_id")["loo_shift_3d_m"].median().idxmax())
    st_diag = st[st["epoch_index"].isin(np.flatnonzero(epochs["diagnosed"].to_numpy(bool)))] if len(st) else st
    summaries = station_summaries(inp, st_diag, infl_df, info, corr_tags)
    logger.info("anomalous epochs %d of %d (threshold %.3g on %s); persistent intervals %d; most influential station "
                "%s", info["n_anomalous"], n_used, info["anomaly_threshold"], info["anomaly_statistic"],
                sum(1 for r in info["anomaly_intervals"] if r["n_epochs"] >= args.min_persistent_epochs),
                infl_station or f"undetermined ({n_loo_epochs} leave-one-out epoch(s); "
                                f"{rules.station_min_epochs} needed)")
    logger.info("reference quality: precision_validation %s; maximum evidence level %d (%s); tags: %s", precision,
                overall_level, overall_conf, ", ".join(ordered(overall_tags)) or "none")

    # Rows for the tables.
    diag_rows = []
    for ei, ep in epochs.iterrows():
        rec = {c: ep.get(c) for c in ep.index if not isinstance(ep.get(c), np.ndarray)}
        rec["time_utc"] = str(ep["time_utc"])
        diag_rows.append(rec)
    infl_rows = []
    if len(st):
        for _, r in st.iterrows():
            rec = r.to_dict()
            rec["time_utc"] = str(r["time_utc"])
            rec["timing_warning"] = bool(r.get("time_mismatch", False))
            infl_rows.append(rec)
    pair_rows = []
    for r in pairs:
        p = r["pair_position"]
        pair_rows.append({**{k: v for k, v in r.items() if k != "pair_position"}, "time_utc": str(r["time_utc"]),
                          "norad_id": inp.norad_id, "pair_x_ecef_m": p[0], "pair_y_ecef_m": p[1],
                          "pair_z_ecef_m": p[2]})

    first = epochs["time_utc"].min()
    date = first.strftime("%Y%m%d") if hasattr(first, "strftime") else "unknown"
    out = args.output_dir / date / f"NORAD_{inp.norad_id}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    log_path_holder.append(out / "logs" / "stage7.log")
    paths = write_tables(diag_rows, infl_rows, pair_rows, corr, out)

    def strongest_doc(fam):
        s = strongest(corr, fam)
        if not s:
            return None
        return {k: s[k] for k in ("scope", "station_id", "x_variable", "y_variable", "coefficient",
                                  "partial_coefficient", "p_value_holm", "n_samples", "n_effective",
                                  "interpretation_tag")}

    summary = {
        "norad_id": inp.norad_id,
        "overall_evidence_level": overall_level, "overall_diagnostic_tags": ordered(overall_tags),
        "overall_confidence": overall_conf,
        "evidence_level_meaning": "evidence pattern only; no physical cause is asserted",
        "n_epochs": int(len(epochs)), "n_diagnosed": n_used, "n_anomalous_epochs": info["n_anomalous"],
        "anomaly_threshold": info["anomaly_threshold"], "anomaly_statistic": info["anomaly_statistic"],
        "anomaly_intervals": info["anomaly_intervals"],
        "persistent_anomaly_intervals": [r for r in info["anomaly_intervals"]
                                         if r["n_epochs"] >= args.min_persistent_epochs],
        "unexplained_intervals": info["unexplained_intervals"],
        "most_influential_station": infl_station or None,
        "n_leave_one_out_epochs": n_loo_epochs,
        "station_test": {"picks": info["station_picks"], "n_epochs_with_station_test": info["n_epochs_with_station_test"],
                         "station_specific": info["station_specific"]},
        "station_summaries": summaries,
        "strongest_calibration_correlation": strongest_doc("calibration"),
        "strongest_timing_correlation": strongest_doc("timing"),
        "strongest_geometry_correlation": strongest_doc("geometry"),
        "strongest_trajectory_model_correlation": strongest_doc("trajectory"),
        "strongest_localisation_correlation": strongest_doc("localisation"),
        "supported_correlation_tags": sorted(corr_tags),
        "n_correlation_tests": int(len(corr)),
        "reference_quality": {"precision_validation": precision,
                              "reference_limited": not precision,
                              "stage6_run_flags": inp.summary6.get("run_flags")},
        "precision_validation": precision,
        "stage5_dynamics": {"run_flags": s5_flags, "n_resets": s5_resets},
        "full_solution_reproduction_max_m": repro,
        "quality_flags": qflags_run,
    }
    write_json(summary, out / "summary.json")

    plots: List[str] = []
    if not args.no_plots and n_used:
        try:
            from tart_gnss.stage7.plotting import make_plots

            plots = make_plots(read_table(paths["diagnosis_parquet"] or paths["diagnosis_csv"]),
                               read_table(paths["station_influence_parquet"] or paths["station_influence_csv"]),
                               corr, out / "diagnostics", f"NORAD {inp.norad_id}")
        except ImportError as exc:
            logger.warning("plots skipped: %s", exc)

    doc = {
        "schema_version": STAGE7_SCHEMA_VERSION,
        "stage": 7,
        "input_stage6_manifest": str(args.stage6_manifest),
        "input_hashes": {**{f"stage{lk.stage}_manifest": lk.sha256 for lk in inp.links},
                         **{f"stage6_{k}": v for k, v in inp.table_hashes.items()}},
        "provenance_chain": [lk.to_dict() for lk in inp.links],
        "alignment": {"keys": ["time_utc", "station_id", "norad_id"], "tolerance_s": args.alignment_tolerance_s,
                      "n_epochs": int(len(epochs)), "n_linked_stage5": int(epochs["link_stage5"].sum()),
                      "n_linked_stage4": int(epochs["link_stage4"].sum()),
                      "n_station_rows": int(len(st)),
                      "n_station_rows_linked_stage3": int(st["link_stage3"].sum()) if len(st) else 0,
                      "n_station_rows_linked_stage2": int(st["link_stage2"].sum()) if len(st) else 0,
                      "row_number_joins": False},
        "leave_one_out": {"enabled": args.leave_one_out, "n_solutions": int(len(infl_df)),
                          "solver": "tart_gnss.stage4.triangulation.triangulate_weighted (reused)",
                          "normalised_influence": "shift^T (C_loo - C_full)^-1 shift (nested estimates); NaN -> "
                                                  "physical shift only",
                          "degraded_geometry": f"fewer than 3 stations remaining or condition number > "
                                               f"{DEFAULT_MAX_CONDITION:g}",
                          "full_solution_reproduction_max_m": repro},
        "pairwise": {"enabled": args.pairwise, "n_solutions": int(len(pairs_df)),
                     "separation_z": "closest-approach separation / sqrt(n^T (rho1^2 C1 + rho2^2 C2) n)"},
        "correlation_analysis": {"method": args.correlation_method, "min_samples": args.min_samples,
                                 "min_correlation": args.min_correlation, "alpha": rules.alpha,
                                 "decision_coefficient": "rank correlation after removing a linear time trend",
                                 "p_value": "effective sample size n(1-r_x r_y)/(1+r_x r_y); Holm over all tests",
                                 "n_tests": int(len(corr)), "supported_tags": sorted(corr_tags)},
        "evidence_rules": rules.to_dict(n_used),
        "quality_policy": {"policy": args.quality_policy, "n_diagnosed": n_used,
                           "missing_evidence": "recorded as quality flags; never treated as evidence of absence"},
        "software": {**_software(), **git, "tart_gnss": __version__},
        "provenance": {"command_line": list(command_line),
                       "run_timestamp_utc": started.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z",
                       "working_directory": os.getcwd()},
        "outputs": {**paths, "summary": str(out / "summary.json"), "diagnostics": plots,
                    "log": str(out / "logs" / "stage7.log"), "manifest": str(out / "stage7_manifest.json")},
    }
    write_manifest(doc, out / "stage7_manifest.json")
    logger.info("outputs: %s", out)
    if n_used < q.MIN_DIAGNOSED_EPOCHS:
        logger.warning("insufficient data for diagnosis: %d diagnosable epoch(s)", n_used)
        return EXIT_INSUFFICIENT
    return EXIT_OK


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-diagnose``."""
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
        return run(args, ["tart-diagnose", *argv], log_path)
    except (InvalidStage6ManifestError, ProvenanceChainError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (Stage7Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)
        if log_path:
            log_path[0].write_text("\n".join(collector.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
