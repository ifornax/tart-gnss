"""``tart-campaign``: population statistics across completed Stage-7 runs.

Discovers or takes Stage-7 manifests, loads each run and the upstream products
it needs (through the recorded provenance; nothing is rerun), applies the
quality policy and filters (listing every excluded run with its reason),
splits included runs into the precision-validation and descriptive cohorts,
and writes run, epoch and station tables, reporting tables, ``summary.json``,
``campaign_manifest.json``, figures and a log under
``<output-dir>/<campaign-name>/``.

Exit codes: 0 campaign analysis completed with >= 3 runs in the analysis
cohort; 3 completed but insufficient for campaign inference (descriptive
outputs still written); 2 invalid input or configuration; 1 unexpected failure.
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
from tart_gnss.campaign import filters as F
from tart_gnss.campaign.accuracy import TERMS, epoch_residual_distribution, performance
from tart_gnss.campaign.aggregation import build_run_tables, run_id_for
from tart_gnss.campaign.anomaly_metrics import evidence_summary, rates_vs_conditions
from tart_gnss.campaign.bootstrap import bootstrap_run_level, median_of
from tart_gnss.campaign.constellation_metrics import (
    compare,
    elevation_summary,
    geometry_summary,
    group_summary,
    satellite_summary,
    snr_summary,
    tle_age_summary,
)
from tart_gnss.campaign.discovery import collect_manifests
from tart_gnss.campaign.errors import (
    CampaignError,
    CampaignSchemaError,
    InvalidStage7ManifestError,
    ManifestDiscoveryError,
)
from tart_gnss.campaign.loader import load_campaign_run
from tart_gnss.campaign.output import CAMPAIGN_SCHEMA_VERSION, write_json, write_table
from tart_gnss.campaign.reproducibility import repeatability
from tart_gnss.campaign.schema import EPOCH_COLUMNS, RUN_COLUMNS, STATION_COLUMNS
from tart_gnss.campaign.station_metrics import stability_highlights, summarize_stations
from tart_gnss.campaign.statistics import sample_tier

logger = logging.getLogger("tart_gnss.campaign")
EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_INSUFFICIENT = 0, 1, 2, 3
DEFAULT_SEED = 20260625
GROUP_FIELDS = ("constellation", "norad_id", "station_set", "cohort", "observation_date", "tle_source")


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``campaign.log``."""

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
    """Argument parser for ``tart-campaign``."""
    p = argparse.ArgumentParser(
        prog="tart-campaign",
        description="Campaign analysis: population statistics across completed Stage-7 runs (runs are the "
                    "independent units). Nothing is reprocessed; residuals are TLE-relative, not absolute accuracy.",
    )
    p.add_argument("--input-root", type=Path, action="append", default=[],
                   help="directory searched recursively for stage7_manifest.json (repeatable)")
    p.add_argument("--stage7-manifest", type=Path, action="append", default=[], help="a Stage-7 manifest (repeatable)")
    p.add_argument("--output-dir", type=Path, default=Path("results/campaign"), help="default results/campaign")
    p.add_argument("--campaign-name", default="campaign", help="output sub-directory and label (default campaign)")
    p.add_argument("--quality-policy", choices=F.POLICIES, default="strict",
                   help="strict (default): runs with broken provenance are excluded; permissive: kept and flagged")
    p.add_argument("--precision-only", action=argparse.BooleanOptionalAction, default=True,
                   help="scientific summaries use the precision-validation cohort (default on); "
                        "--no-precision-only analyses all included runs")
    p.add_argument("--precision-validation-only", action="store_true",
                   help="exclude runs that are not precision validations from the campaign altogether")
    p.add_argument("--min-epochs-per-run", type=int, default=3, help="minimum diagnosed epochs per run (default 3)")
    p.add_argument("--min-runs-per-group", type=int, default=3, help="minimum runs per compared group (default 3)")
    p.add_argument("--group-by", nargs="+", choices=GROUP_FIELDS, default=["constellation"],
                   help="extra grouped summary fields (default constellation)")
    p.add_argument("--reference-frame", choices=("rtn", "gcrs"), default="rtn",
                   help="frame of the residual components summarised (default rtn)")
    p.add_argument("--bootstrap", type=int, default=1000, help="run-level bootstrap resamples (default 1000)")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED, help=f"bootstrap seed (default {DEFAULT_SEED})")
    p.add_argument("--norad-id", type=int, nargs="+", default=[], help="keep only these NORAD IDs")
    p.add_argument("--constellation", nargs="+", default=[], help="keep only these constellations")
    p.add_argument("--station", nargs="+", default=[], help="keep only runs including one of these stations")
    p.add_argument("--start", default=None, help="keep runs starting at or after this UTC time")
    p.add_argument("--end", default=None, help="keep runs starting at or before this UTC time")
    p.add_argument("--min-elevation", type=float, default=None, help="minimum per-run median elevation (deg)")
    p.add_argument("--max-tle-age-days", type=float, default=None, help="maximum |TLE age| (days)")
    p.add_argument("--max-evidence-level", type=int, default=None, help="maximum overall evidence level")
    p.add_argument("--publication-plots", action="store_true", help="larger fonts and vector PDF + PNG figures")
    p.add_argument("--no-plots", action="store_true", help="skip figures")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    return p


def _validate(args) -> None:
    if not args.input_root and not args.stage7_manifest:
        raise ManifestDiscoveryError("give --input-root and/or --stage7-manifest")
    if args.bootstrap < 0 or args.min_epochs_per_run < 1 or args.min_runs_per_group < 2:
        raise ManifestDiscoveryError("--bootstrap >= 0, --min-epochs-per-run >= 1 and --min-runs-per-group >= 2")
    if not args.campaign_name or "/" in args.campaign_name:
        raise ManifestDiscoveryError("--campaign-name must be a plain name")


def balance(runs) -> Dict[str, object]:
    """Share of the largest satellite / constellation / station set / day / elevation range, and flags."""
    out: Dict[str, object] = {"flags": []}
    n = len(runs)
    if not n:
        return out
    r = runs.copy()
    r["elevation_range"] = (np.floor(r["median_elevation_deg"].astype(float) / 10) * 10).astype("Int64").astype(str)
    for key in ("norad_id", "constellation", "station_set", "observation_date", "elevation_range"):
        counts = r[key].astype(str).value_counts()
        out[key] = {"n_distinct": int(len(counts)), "largest": counts.index[0], "largest_share": float(counts.iloc[0] / n)}
    if n >= 3 and out["norad_id"]["largest_share"] > 0.5:
        out["flags"].append("DOMINATED_BY_ONE_SATELLITE")
    if n >= 3 and out["station_set"]["largest_share"] > 0.8:
        out["flags"].append("DOMINATED_BY_ONE_STATION_SET")
    age = r["median_tle_age_days"].astype(float).abs()
    bins = np.digitize(age, [1, 3, 7, 14])
    if n >= 3 and np.bincount(bins).max() / n > 0.8:
        out["flags"].append("TLE_AGE_IMBALANCE")
    return out


def run(args: argparse.Namespace, command_line: Sequence[str], log_path_holder: List[Path]) -> int:
    """Execute the campaign analysis; return the exit status."""
    import pandas as pd

    git = _git_state()
    started = datetime.now(timezone.utc)
    _validate(args)
    out = args.output_dir / args.campaign_name
    (out / "logs").mkdir(parents=True, exist_ok=True)
    log_path_holder.append(out / "logs" / "campaign.log")
    logger.info("campaign %s; input roots %s; explicit manifests %d", args.campaign_name,
                [str(x) for x in args.input_root], len(args.stage7_manifest))
    manifests = collect_manifests(args.input_root, args.stage7_manifest)
    logger.info("Stage-7 manifests discovered: %d", len(manifests))

    rows, run_errors, epoch_rows, station_rows, ids = [], [], {}, {}, set()
    for mp in manifests:
        try:
            run_obj = load_campaign_run(mp)
        except CampaignSchemaError as exc:
            run_errors.append({"run_id": None, "manifest_path": str(mp), "norad_id": None, "reason": "SCHEMA_MISMATCH",
                               "quality_flags": [], "detail": str(exc)})
            continue
        except InvalidStage7ManifestError as exc:
            run_errors.append({"run_id": None, "manifest_path": str(mp), "norad_id": None,
                               "reason": "MISSING_STAGE7_OUTPUT", "quality_flags": [], "detail": str(exc)})
            continue
        diag = run_obj.diagnosis
        dsel = diag[diag["diagnosed"].fillna(False).astype(bool)] if len(diag) else diag
        rid = run_id_for(run_obj, str(dsel["time_utc"].min()) if len(dsel) else (
            str(diag["time_utc"].min()) if len(diag) else None))
        base, k = rid, 2
        while rid in ids:
            rid = f"{base}-{k}"
            k += 1
        ids.add(rid)
        r, e, s = build_run_tables(run_obj, rid, args.reference_frame)
        r["provenance_problems"] = run_obj.provenance_problems
        rows.append(r)
        epoch_rows[rid], station_rows[rid] = e, s
    logger.info("manifests loaded: %d (%d unreadable or unsupported)", len(rows), len(run_errors))
    filt = F.RunFilters(norad_ids=args.norad_id, constellations=args.constellation, stations=args.station,
                        start=args.start, end=args.end, min_elevation_deg=args.min_elevation,
                        max_tle_age_days=args.max_tle_age_days, max_evidence_level=args.max_evidence_level,
                        precision_validation_only=args.precision_validation_only,
                        min_epochs_per_run=args.min_epochs_per_run)
    inc, exc = F.split_runs(rows, filt, args.quality_policy)
    excluded = exc + run_errors

    def frame(rs, cols):
        df = pd.DataFrame(rs)
        for c in cols:
            if c not in df.columns:
                df[c] = None
        return df[cols]

    runs = frame(inc, RUN_COLUMNS)
    inc_ids = list(runs["run_id"])
    epochs = frame([x for rid in inc_ids for x in epoch_rows[rid]], EPOCH_COLUMNS)
    stations = frame([x for rid in inc_ids for x in station_rows[rid]], STATION_COLUMNS)
    for col in ("precision_validation", "stage5_dynamics_mismatch"):
        runs[col] = runs[col].astype(bool)
    science = runs[runs["precision_validation"]] if args.precision_only else runs
    n_prec = int(runs["precision_validation"].sum()) if len(runs) else 0
    logger.info("runs included %d, excluded %d; precision-validation runs %d; analysis cohort %d (%s)", len(runs),
                len(excluded), n_prec, len(science), "precision only" if args.precision_only else "all included")
    if len(runs):
        logger.info("satellites %d; constellations %s; stations %s; span %s .. %s", runs["norad_id"].nunique(),
                    sorted(runs["constellation"].unique()), sorted(stations["station_id"].unique()) if len(stations)
                    else [], runs["start_utc"].min(), runs["end_utc"].max())
    logger.info("bootstrap: %d run-level resamples, seed %d", args.bootstrap, args.seed)

    # Analyses.
    perf = performance(science, args.bootstrap, args.seed)
    perf_all_unc = performance(runs, args.bootstrap, args.seed)
    st_sum = summarize_stations(stations, args.bootstrap, args.seed, args.min_runs_per_group)
    highlights = stability_highlights(st_sum, args.min_runs_per_group)
    const_tab = group_summary(runs, "constellation", science, args.min_runs_per_group) if len(runs) else pd.DataFrame()
    const_tests = {m: compare(src, "constellation", m, args.min_runs_per_group)
                   for m, src in (("median_snr", runs), ("median_sigma_theta_deg", runs),
                                  ("median_stage4_sigma_3d_m", runs), ("median_position_residual_3d_m", science))}
    sat_tab = satellite_summary(runs, science, args.min_runs_per_group) if len(runs) else pd.DataFrame()
    set_tab = group_summary(runs, "station_set", science, args.min_runs_per_group) if len(runs) else pd.DataFrame()
    geo_tab, geo_dep = geometry_summary(runs, science) if len(runs) else (pd.DataFrame(), {})
    age_tab, age_dep = tle_age_summary(runs) if len(runs) else (pd.DataFrame(), {})
    el_tab = elevation_summary(stations) if len(stations) else pd.DataFrame()
    snr_tab = snr_summary(stations) if len(stations) else pd.DataFrame()
    ev_tab = evidence_summary(runs, epochs, stations)
    rates = rates_vs_conditions(runs, stations)
    repeat = repeatability(science if len(science) >= 2 else runs, stations)
    grouped = {g: group_summary(runs, g, science, args.min_runs_per_group) for g in args.group_by} if len(runs) else {}
    bal = balance(runs)

    flags = list(bal.get("flags", []))
    tier = sample_tier(len(science))
    if len(science) < 3:
        flags.append("INSUFFICIENT_CAMPAIGN_SIZE")
    if len(runs) and n_prec / len(runs) < 0.5:
        flags.append("LOW_PRECISION_VALIDATION_FRACTION")
    if any(x.get("reason") == "BROKEN_PROVENANCE" for x in excluded) or any(
            "BROKEN_PROVENANCE" in (r.get("quality_flags") or []) for r in inc):
        flags.append("BROKEN_PROVENANCE_PRESENT")
    if len(runs) and (runs["constellation"] == "Unknown").any():
        flags.append("MISSING_CONSTELLATION_METADATA")
    if len(const_tab) and const_tab["under_sampled"].any():
        flags.append("INSUFFICIENT_GROUP_SIZE")
    flags = sorted(set(flags)) or ["GOOD"]

    def ci(col, src):
        if not len(src):
            return None
        b = bootstrap_run_level(src, median_of(col), args.bootstrap, args.seed)
        return {"median": b.estimate, "ci95": [b.ci_low, b.ci_high], "n_runs": b.n_runs}

    lv_run = ev_tab[(ev_tab["unit"] == "run") & ev_tab["quantity"].str.startswith("evidence_level")]
    tag_run = ev_tab[(ev_tab["unit"] == "run") & ~ev_tab["quantity"].str.startswith("evidence_level")]
    summary = {
        "campaign_name": args.campaign_name, "terminology": TERMS,
        "n_discovered_runs": len(manifests), "n_included_runs": int(len(runs)), "n_excluded_runs": len(excluded),
        "n_precision_validation_runs": n_prec, "n_descriptive_only_runs": int(len(runs) - n_prec),
        "analysis_cohort": "precision_validation" if args.precision_only else "all_included",
        "n_analysis_runs": int(len(science)), "sample_size_tier": tier,
        "n_satellites": int(runs["norad_id"].nunique()) if len(runs) else 0,
        "n_constellations": int(runs["constellation"].nunique()) if len(runs) else 0,
        "constellations": sorted(runs["constellation"].unique()) if len(runs) else [],
        "n_unique_stations": int(stations["station_id"].nunique()) if len(stations) else 0,
        "observation_start": runs["start_utc"].min() if len(runs) else None,
        "observation_end": runs["end_utc"].max() if len(runs) else None,
        "total_observing_duration_s": float(runs["duration_s"].sum()) if len(runs) else 0.0,
        "precision_validation_fraction": float(n_prec / len(runs)) if len(runs) else None,
        "median_stage3_direction_uncertainty_deg": ci("median_sigma_theta_deg", runs),
        "median_stage4_3d_uncertainty_m": ci("median_stage4_sigma_3d_m", runs),
        "median_stage4_sigma_major_m": ci("median_sigma_major_m", runs),
        "median_stage4_sigma_minor_m": ci("median_sigma_minor_m", runs),
        "median_stage5_smoothed_uncertainty_m": ci("median_stage5_sigma_3d_m", runs),
        "median_tle_relative_residual_3d_m": ci("median_position_residual_3d_m", science),
        "median_tle_relative_radial_m": ci("median_radial_m", science),
        "median_tle_relative_along_track_m": ci("median_along_track_m", science),
        "median_tle_relative_cross_track_m": ci("median_cross_track_m", science),
        "median_tle_relative_velocity_m_s": ci("median_velocity_residual_m_s", science),
        "reference_frame": args.reference_frame,
        "performance_analysis_cohort": perf, "measurement_uncertainty_all_included": {
            k: perf_all_unc[k] for k in ("direction", "localisation_3d", "trajectory")},
        "epoch_residual_distribution_display_only": epoch_residual_distribution(
            epochs[epochs["run_id"].isin(set(science["run_id"]))]) if len(epochs) else {},
        "evidence_level_fractions_by_run": {r["quantity"]: r["fraction"] for _, r in lv_run.iterrows()},
        "diagnostic_tag_fractions_by_run": {r["quantity"]: {"fraction": r["fraction"], "count": r["count"],
                                                            "ci95": [r["ci95_low"], r["ci95_high"]]}
                                            for _, r in tag_run.iterrows()},
        "station_stability": highlights, "most_stable_stations": highlights["most_stable"],
        "most_frequently_flagged_stations": highlights["most_frequently_flagged"],
        "constellation_comparisons": const_tests, "geometry_dependence": geo_dep, "tle_age_dependence": age_dep,
        "repeatability": repeat, "balance": bal, "campaign_quality_flags": flags,
        "inference_unit": "run (epochs within a pass are correlated and never treated as independent)",
        "no_physical_cause_rates": True,
    }
    paths = {}
    for name, df in (("campaign_runs", runs), ("campaign_epochs", epochs), ("campaign_stations", stations)):
        p = write_table(df, out, name)
        paths[f"{name}_parquet"], paths[f"{name}_csv"] = p["parquet"], p["csv"]
    exc_df = pd.DataFrame([{k: x.get(k) for k in ("run_id", "manifest_path", "norad_id", "reason", "quality_flags",
                                                   "detail")} for x in excluded],
                          columns=["run_id", "manifest_path", "norad_id", "reason", "quality_flags", "detail"])
    paths["excluded_runs_csv"] = write_table(exc_df, out, "excluded_runs", parquet=False)["csv"]
    tables = out / "tables"
    for name, df in (("station_summary", st_sum), ("constellation_summary", const_tab), ("satellite_summary", sat_tab),
                     ("evidence_summary", ev_tab), ("tle_age_summary", age_tab), ("geometry_summary", geo_tab),
                     ("elevation_summary", el_tab), ("snr_summary", snr_tab), ("station_set_summary", set_tab),
                     ("anomaly_rates", rates)):
        paths[f"{name}_csv"] = write_table(df if isinstance(df, pd.DataFrame) else pd.DataFrame(), tables, name,
                                           parquet=False)["csv"]
    for g, df in grouped.items():
        paths[f"grouped_by_{g}_csv"] = write_table(df, tables, f"grouped_by_{g}", parquet=False)["csv"]
    write_json(summary, out / "summary.json")

    plots: List[str] = []
    if not args.no_plots and len(runs):
        try:
            from tart_gnss.campaign.plotting import make_plots

            plots = make_plots(runs, epochs, stations, st_sum, science, repeat, out / "diagnostics",
                               args.campaign_name, args.publication_plots)
        except ImportError as exc:
            logger.warning("plots skipped: %s", exc)

    m = summary["median_tle_relative_residual_3d_m"]
    logger.info("median Stage-4 3-D uncertainty %s m; median TLE-relative 3-D residual (analysis cohort) %s m",
                summary["median_stage4_3d_uncertainty_m"], m)
    logger.info("campaign quality flags: %s; sample-size tier %s", ", ".join(flags), tier)
    doc = {
        "schema_version": CAMPAIGN_SCHEMA_VERSION, "analysis": "campaign", "campaign_name": args.campaign_name,
        "input_manifests": [str(p) for p in manifests],
        "input_hashes": {r["manifest_path"]: r["manifest_sha256"] for r in rows},
        "included_runs": [{"run_id": r["run_id"], "manifest_path": r["manifest_path"], "cohort": r["cohort"]}
                          for r in inc],
        "excluded_runs": [{k: x.get(k) for k in ("run_id", "manifest_path", "reason")} for x in excluded],
        "filters": filt.to_dict(),
        "grouping": {"group_by": args.group_by, "min_runs_per_group": args.min_runs_per_group,
                     "reference_frame": args.reference_frame},
        "bootstrap": {"n": args.bootstrap, "seed": args.seed, "unit": "run", "interval": "percentile, 95 %"},
        "quality_policy": {"policy": args.quality_policy, "precision_only": args.precision_only,
                           "precision_validation_only": args.precision_validation_only},
        "no_reprocessing": "completed Stage-7 runs and their recorded upstream products are read; no stage is run",
        "software": {**_software(), **git, "tart_gnss": __version__},
        "provenance": {"command_line": list(command_line),
                       "run_timestamp_utc": started.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z",
                       "working_directory": os.getcwd()},
        "outputs": {**paths, "summary": str(out / "summary.json"), "diagnostics": plots,
                    "log": str(out / "logs" / "campaign.log"), "manifest": str(out / "campaign_manifest.json")},
    }
    write_json(doc, out / "campaign_manifest.json")
    logger.info("outputs: %s", out)
    if len(science) < 3:
        logger.warning("INSUFFICIENT_CAMPAIGN_SIZE: %d run(s) in the analysis cohort; outputs are descriptive only",
                       len(science))
        return EXIT_INSUFFICIENT
    return EXIT_OK


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-campaign``."""
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
        return run(args, ["tart-campaign", *argv], log_path)
    except ManifestDiscoveryError as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (CampaignError, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)
        if log_path:
            log_path[0].write_text("\n".join(collector.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
