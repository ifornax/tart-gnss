"""``tart-archive-scan``: archive-aware candidate discovery.

Lists every required ``(station, UTC day)`` archive prefix once (or reuses
cached snapshots), builds conservative filename-based candidate coverage per
station, intersects the stations, and, with ``--visibility-results``,
matches the shared candidate segments with every Stage-1 window and ranks
the intersections. Results are ``ARCHIVE_COVERAGE_CANDIDATE`` evidence only:
nothing is downloaded, no sample timestamp is verified, and ``tart-retrieve``
is never run (a dry-run command is printed for manual use).

With ``--all-stations`` the command runs network-wide discovery instead
(:mod:`tart_gnss.archive_scan.network_run`): every eligible registry station,
one sweep for intervals with at least ``--min-stations`` distinct physical
sites, and ranked station triads with proposed Stage-1 commands.

Exit status:

* ``0``: complete inventory, scan finished (even with zero candidates);
* ``3``: scan finished but the inventory is incomplete (listing errors or
  unknown days); outputs are written with a warning;
* ``2``: invalid arguments, station selection, interval or Stage-1 input;
* ``1``: unexpected technical failure, or every listing failed.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from tart_gnss import __version__
from tart_gnss.archive_scan import COVERAGE_CONFIDENCE, EVIDENCE_CLASS, network
from tart_gnss.archive_scan import output as out
from tart_gnss.archive_scan.cache import CACHE_SCHEMA_VERSION, InventoryCache
from tart_gnss.archive_scan.errors import ArchiveScanError, ArchiveScanInputError
from tart_gnss.archive_scan.intervals import (
    CommonSegment,
    count_files,
    gap_groups,
    intersect_all,
    k_station_segments,
    merge_intervals,
    station_candidate_segments,
)
from tart_gnss.archive_scan.inventory import (
    STATUS_EMPTY,
    STATUS_ERROR,
    STATUS_FOUND,
    STATUS_UNKNOWN,
    build_inventory,
    required_days,
)
from tart_gnss.archive_scan.matching import IN_BAND, filter_windows, load_windows, match_windows
from tart_gnss.archive_scan.ranking import DEFAULT_DURATION_CAP_S, RANKING_RULE, rank_candidates
from tart_gnss.stage2 import archive as archive_mod
from tart_gnss.stage2 import selection as sel
from tart_gnss.stage2.archive import ArchiveBackend
from tart_gnss.stage2.cli import resolve_run_stations
from tart_gnss.stage2.errors import Stage1SelectionError, Stage2Error
from tart_gnss.stage2.synchronise import total_seconds
from tart_gnss.stations import (
    Station,
    StationConfigError,
    default_station_file,
    load_aliases,
    load_catalogue_metadata,
    load_stations,
    resolve_stations,
)
from tart_gnss.tle import format_utc, sha256_file

logger = logging.getLogger("tart_gnss.archive_scan")

EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_INCOMPLETE = 0, 1, 2, 3
MANIFEST_SCHEMA_VERSION = "1.0"
CACHE_POLICIES = ("refresh", "reuse", "offline")


def build_parser() -> argparse.ArgumentParser:
    """Create the ``tart-archive-scan`` argument parser."""
    p = argparse.ArgumentParser(
        prog="tart-archive-scan",
        description="List the TART visibility archive once per station/day, find possible multi-station "
                    "archive overlaps from filenames (candidates, NOT validated data), and match them with "
                    "every Stage-1 common-visibility window. Nothing is downloaded.",
    )
    who = p.add_mutually_exclusive_group()
    who.add_argument("--stations", nargs="+", metavar="ID",
                     help="station ids or aliases (default: the frozen Stage-1 stations); with "
                          "--visibility-results the set must match Stage 1")
    who.add_argument("--all-stations", action="store_true",
                     help="NETWORK mode: every registry station with an archive_name and valid coordinates; "
                          "finds intervals and station triads at >= --min-stations distinct physical sites")
    p.add_argument("--start", help="scan start, ISO-8601 UTC (default: Stage-1 start); half-open [start, end)")
    p.add_argument("--end", help="scan end, ISO-8601 UTC (default: Stage-1 end)")
    p.add_argument("--visibility-results", metavar="PATH",
                   help="Stage-1 detailed JSON (or its CSV with the JSON alongside); omit for an inventory-only scan")
    p.add_argument("--min-common-duration", type=float, default=300.0, metavar="SECONDS",
                   help="minimum possible overlap duration of a reported candidate (default 300)")
    p.add_argument("--signal-relevance", nargs="+", choices=sel.SIGNAL_RELEVANCE_VALUES, metavar="VALUE",
                   help="keep Stage-1 windows with these classifications (default with Stage-1 input: GPS_L1; "
                        f"choices: {', '.join(sel.SIGNAL_RELEVANCE_VALUES)})")
    p.add_argument("--constellations", nargs="+", metavar="NAME", help="additional constellation filter (default: all)")
    p.add_argument("--max-abs-tle-age-days", type=float, metavar="D",
                   help="exclude candidates whose |TLE age| at the overlap midpoint exceeds D (counted, not silent)")
    p.add_argument("--station-file",
                   help="station catalogue for archive names (default: packaged frozen registry snapshot)")
    p.add_argument("--aliases", help="station alias file (default: bundled station_aliases.yaml)")
    p.add_argument("--file-span-assumption-s", type=float, default=archive_mod.DEFAULT_FILE_SPAN_MARGIN_S,
                   metavar="S", help="ASSUMED maximum duration before a file's filename time, for candidate "
                                     "discovery only (default 120)")
    p.add_argument("--gap-tolerance-s", type=float, default=0.0, metavar="S",
                   help="gap allowed when grouping segments for the diagnostic summary; gaps never count as "
                        "covered (default 0)")
    p.add_argument("--duration-cap-s", type=float, default=DEFAULT_DURATION_CAP_S, metavar="S",
                   help="ranking caps the overlap duration at this value so long windows cannot dominate (default 600)")
    p.add_argument("--cache-dir", default="data/archive_inventory", help="inventory cache (default data/archive_inventory)")
    p.add_argument("--cache-policy", choices=CACHE_POLICIES, default="refresh",
                   help="refresh: list live; reuse: valid cache younger than --cache-max-age-hours; "
                        "offline: cache only, missing days UNKNOWN (default refresh)")
    p.add_argument("--cache-max-age-hours", type=float, default=24.0, metavar="H", help="TTL for reuse (default 24)")
    p.add_argument("--output", help="main CSV: ranked candidates (default results/archive_candidates.csv), or "
                                    "ranked triads with --all-stations (default results/network_overlap.csv)")
    p.add_argument("--inventory-output", help="station/day CSV (default: alongside --output)")
    p.add_argument("--min-stations", type=int, metavar="N",
                   help="stations needed for an overlap (default: all selected); fewer gives DIAGNOSTIC rows only. "
                        "With --all-stations: distinct physical SITES per interval and station-set size (default 3)")
    net = p.add_argument_group("network mode (--all-stations only)")
    net.add_argument("--include-stations", nargs="+", metavar="ID", help="restrict eligible stations to these")
    net.add_argument("--exclude-stations", nargs="+", metavar="ID", help="exclude these stations")
    net.add_argument("--min-baseline-km", type=float, metavar="KM",
                     help="require EVERY pairwise geodesic baseline of a station set >= KM (default: no cut); a "
                          "site-diversity filter, not a triangulation-quality metric")
    net.add_argument("--site-tolerance-km", type=float, default=network.DEFAULT_SITE_TOLERANCE_KM, metavar="KM",
                     help="stations within KM (geodesic, single linkage) are one physical site (default 1)")
    net.add_argument("--max-triads", type=int, metavar="N",
                     help="export/print only the top N ranked station sets (counts and scanning are unaffected)")
    p.add_argument("--top", type=int, default=20, help="rows printed (all rows are written; default 20)")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


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


def _parse_time(text: str, name: str) -> datetime:
    try:
        return sel.parse_utc(text)
    except ValueError as exc:
        raise ArchiveScanInputError(f"{name}: {exc}") from exc


def _validate_numbers(args: argparse.Namespace) -> None:
    checks = [
        (args.min_common_duration >= 0, "--min-common-duration must be >= 0"),
        (args.file_span_assumption_s > 0, "--file-span-assumption-s must be > 0"),
        (args.gap_tolerance_s >= 0, "--gap-tolerance-s must be >= 0"),
        (args.duration_cap_s > 0, "--duration-cap-s must be > 0"),
        (args.cache_max_age_hours > 0, "--cache-max-age-hours must be > 0"),
        (args.top >= 0, "--top must be >= 0"),
        (args.max_abs_tle_age_days is None or args.max_abs_tle_age_days > 0, "--max-abs-tle-age-days must be > 0"),
        (args.min_baseline_km is None or args.min_baseline_km >= 0, "--min-baseline-km must be >= 0"),
        (args.site_tolerance_km >= 0, "--site-tolerance-km must be >= 0"),
        (args.max_triads is None or args.max_triads >= 1, "--max-triads must be >= 1"),
    ]
    for ok, msg in checks:
        if not ok:
            raise ArchiveScanInputError(msg)


def _reject_network_options(args: argparse.Namespace) -> None:
    """Network-only options are an input error without ``--all-stations``."""
    for opt, value in (("--include-stations", args.include_stations), ("--exclude-stations", args.exclude_stations),
                       ("--min-baseline-km", args.min_baseline_km), ("--max-triads", args.max_triads)):
        if value is not None:
            raise ArchiveScanInputError(f"{opt} needs --all-stations")
    if args.site_tolerance_km != network.DEFAULT_SITE_TOLERANCE_KM:
        raise ArchiveScanInputError("--site-tolerance-km needs --all-stations")


def _canonical(st: Station) -> str:
    return (st.tart_id or st.id).lower()


def resolve_scan_stations(args: argparse.Namespace, station_file: Path, aliases: Dict[str, str],
                          doc: Optional[dict]) -> Tuple[List[Station], Dict[str, Dict[str, object]]]:
    """Stations of the scan, with frozen Stage-1 coordinates when Stage-1 input is given.

    With Stage-1 input, the Stage-1 stations are resolved exactly as Stage 2
    does (:func:`tart_gnss.stage2.cli.resolve_run_stations`: identity and
    coordinates from Stage 1, archive names from the catalogue), and any
    ``--stations`` selection must be the same set of canonical TART ids.

    Raises
    ------
    ArchiveScanInputError
        For unknown stations or a selection that differs from Stage 1.
    """
    try:
        catalogue = load_stations(station_file)
        requested = resolve_stations(args.stations, catalogue, aliases) if args.stations else None
    except StationConfigError as exc:
        raise ArchiveScanInputError(str(exc)) from exc
    if doc is None:
        if not requested:
            raise ArchiveScanInputError("--stations is required without --visibility-results")
        return requested, {s.id: {"station_metadata_changed": False} for s in requested}
    infos = [sel.StationInfo(id=s["id"], name=s.get("name", s["id"]), latitude_deg=float(s["latitude_deg"]),
                             longitude_deg=float(s["longitude_deg"]), altitude_m=float(s["altitude_m"]),
                             tart_id=s.get("tart_id"), raw=s) for s in doc["stations"]]
    try:
        stations, checks = resolve_run_stations(SimpleNamespace(stations=infos), station_file, aliases)
    except Stage1SelectionError as exc:
        raise ArchiveScanInputError(str(exc)) from exc
    if requested is not None:
        want = sorted(_canonical(s) for s in requested)
        have = sorted(_canonical(s) for s in stations)
        if want != have:
            raise ArchiveScanInputError(
                f"--stations {want} do not match the frozen Stage-1 stations {have}; "
                "Stage-1 geometry is never changed by the scanner"
            )
    return stations, checks


def _scan_interval(args: argparse.Namespace, doc: Optional[dict]) -> Tuple[datetime, datetime]:
    if doc is None:
        if not args.start or not args.end:
            raise ArchiveScanInputError("--start and --end are required without --visibility-results")
        start, end = _parse_time(args.start, "--start"), _parse_time(args.end, "--end")
    else:
        s1_start, s1_end = sel.parse_utc(doc["run"]["start_utc"]), sel.parse_utc(doc["run"]["end_utc"])
        start = _parse_time(args.start, "--start") if args.start else s1_start
        end = _parse_time(args.end, "--end") if args.end else s1_end
        if start < s1_start or end > s1_end:
            raise ArchiveScanInputError(
                f"--start/--end ({format_utc(start)} -- {format_utc(end)}) may only restrict the Stage-1 "
                f"interval {format_utc(s1_start)} -- {format_utc(s1_end)}; run an inventory-only scan to go beyond it"
            )
    if end <= start:
        raise ArchiveScanInputError("--end must be after --start")
    return start, end


def build_segments(stations: Sequence[Station], file_times: Dict[str, List[datetime]], listings, days,
                   start: datetime, end: datetime, span_s: float, k: int,
                   ) -> Tuple[List[CommonSegment], Dict[str, list]]:
    """Per-station candidate segments and the shared segments covered by at least ``k`` stations.

    A segment's ``inventory_complete`` is False if any selected station's
    listing of a day that could hold files for it (see
    :func:`tart_gnss.archive_scan.inventory.required_days`) is not complete.
    """
    ids = [s.id for s in stations]
    per_station = {sid: station_candidate_segments(file_times[sid], span_s, (start, end)) for sid in ids}
    if k == len(ids):
        shared = [((a, b), tuple(sorted(ids))) for a, b in intersect_all([per_station[s] for s in ids])]
    else:
        shared = k_station_segments(per_station, k)
    scanned = set(days)
    segments = []
    for n, ((a, b), sids) in enumerate(shared, start=1):
        need = [d for d in required_days(a, b, span_s) if d in scanned]
        complete = all(listings[(sid, d)].complete for sid in ids for d in need)
        order = [s for s in ids if s in sids]
        segments.append(CommonSegment(
            segment_id=n, start_utc=a, end_utc=b, station_ids=tuple(order),
            file_counts={sid: count_files(file_times[sid], a, b, span_s) for sid in order},
            inventory_complete=complete, all_stations=len(sids) == len(ids),
        ))
    return segments, per_station


def classify(all_complete: bool, any_complete: bool, any_error: bool, stage1: bool,
             n_segments_min: int, n_candidates: int) -> Tuple[str, int]:
    """Exit classification and status (see the module docstring)."""
    if not any_complete:
        return "ARCHIVE_LISTING_FAILURE", EXIT_FAILURE if any_error else EXIT_INCOMPLETE
    code = EXIT_OK if all_complete else EXIT_INCOMPLETE
    if n_segments_min == 0:
        return ("NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS" if all_complete
                else "NO_CANDIDATE_OVERLAP_BUT_INCOMPLETE_INVENTORY"), code
    if not stage1:
        return ("COMMON_CANDIDATE_SEGMENTS_FOUND" if all_complete
                else "COMMON_CANDIDATE_SEGMENTS_FOUND_INVENTORY_INCOMPLETE"), code
    if n_candidates == 0:
        return ("NO_IN_BAND_VISIBILITY_INTERSECTIONS" if all_complete
                else "NO_IN_BAND_VISIBILITY_INTERSECTIONS_INVENTORY_INCOMPLETE"), code
    return ("CANDIDATES_FOUND" if all_complete else "CANDIDATES_FOUND_INVENTORY_INCOMPLETE"), code


def run(args: argparse.Namespace, argv: Sequence[str],
        archive: Optional[ArchiveBackend] = None,
        now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
        sleep: Callable[[float], None] = time.sleep) -> int:
    """Execute a scan for parsed arguments; returns the exit status.

    ``archive`` (a mock in tests) replaces the public S3 archive; ``now`` and
    ``sleep`` are injectable for cache-age and retry tests.
    """
    _validate_numbers(args)
    if args.cache_policy not in CACHE_POLICIES:
        raise ArchiveScanInputError(f"--cache-policy must be one of {CACHE_POLICIES}")
    if getattr(args, "all_stations", False):
        from tart_gnss.archive_scan.network_run import run_network

        return run_network(args, argv, archive=archive, now=now, sleep=sleep)
    _reject_network_options(args)
    if args.output is None:
        args.output = "results/archive_candidates.csv"
    run_ts = now()
    run_tag = run_ts.strftime("%Y%m%dT%H%M%SZ")

    doc = results_path = None
    stage1_sha = None
    if args.visibility_results:
        try:
            doc, results_path = sel.load_stage1_results(args.visibility_results)
        except Stage1SelectionError as exc:
            raise ArchiveScanInputError(str(exc)) from exc
        stage1_sha = sha256_file(results_path)
        logger.info("Stage-1 results %s (SHA-256 %s): %d windows", results_path, stage1_sha, len(doc["windows"]))
    station_file = Path(args.station_file) if args.station_file else default_station_file()
    try:
        aliases = load_aliases(args.aliases)
        registry_meta = load_catalogue_metadata(station_file)
    except StationConfigError as exc:
        raise ArchiveScanInputError(str(exc)) from exc
    stations, metadata_checks = resolve_scan_stations(args, station_file, aliases, doc)
    start, end = _scan_interval(args, doc)
    n = len(stations)
    k = n if args.min_stations is None else args.min_stations
    if not 2 <= k <= n:
        raise ArchiveScanInputError(f"--min-stations must be between 2 and {n}")
    relevance = args.signal_relevance
    if doc is not None and relevance is None:
        relevance = ["GPS_L1"]
    for st in stations:
        logger.info("Station %s: TART id %s, archive '%s'", st.id, st.tart_id, st.archive_name)
        if not st.archive_name:
            logger.error("Station %s has no archive_name; its days are listing errors", st.id)

    span = args.file_span_assumption_s
    days = required_days(start, end, span)
    host, bucket = archive_mod.ARCHIVE_HOST, archive_mod.ARCHIVE_BUCKET
    cache = InventoryCache(args.cache_dir)
    listings, stats = build_inventory(
        stations, days, (lambda: archive) if archive is not None else archive_mod.TartS3Archive,
        cache, args.cache_policy, args.cache_max_age_hours, host, bucket, __version__, now=now, sleep=sleep,
    )

    # Per-station filename times, deduplicated by full key across days.
    file_times: Dict[str, List[datetime]] = {}
    for st in stations:
        by_key = {o.key: o.filename_utc for d in days for o in listings[(st.id, d)].objects}
        file_times[st.id] = sorted(by_key.values())
    segments, per_station = build_segments(stations, file_times, listings, days, start, end, span, k)
    all_station_segments = [s for s in segments if s.all_stations]
    min_dur = args.min_common_duration
    n_seg_min = sum(1 for s in all_station_segments if s.duration_s >= min_dur)

    windows_total = windows_kept = 0
    exclusions: Dict[str, int] = {}
    ranked: List = []
    n_short = n_tle_excluded = 0
    if doc is not None:
        windows = load_windows(doc)
        windows_total = len(windows)
        kept, exclusions = filter_windows(windows, relevance, args.constellations)
        windows_kept = len(kept)
        candidates, n_short = match_windows(kept, segments, file_times, span, min_dur)
        if args.max_abs_tle_age_days is not None:
            limit = args.max_abs_tle_age_days
            n_tle_excluded = sum(1 for c in candidates if c.tle_age_days is not None and abs(c.tle_age_days) > limit)
            candidates = [c for c in candidates if c.tle_age_days is None or abs(c.tle_age_days) <= limit]
        ranked = rank_candidates(candidates, args.duration_cap_s)
        logger.info("Matched %d of %d Stage-1 windows (after filters) with %d shared segments: %d candidates "
                    "(%d shorter than %g s, %d excluded by TLE age)", windows_kept, windows_total, len(segments),
                    len(ranked), n_short, min_dur, n_tle_excluded)

    all_complete = all(L.complete for L in listings.values())
    any_complete = any(L.complete for L in listings.values())
    any_error = any(L.status == STATUS_ERROR for L in listings.values())
    full_candidates = [c for c in ranked if c.segment.all_stations]
    classification, code = classify(all_complete, any_complete, any_error, doc is not None, n_seg_min,
                                    len(full_candidates))

    # ---------------------------------------------------------------- outputs
    paths = out.output_paths(Path(args.output), Path(args.inventory_output) if args.inventory_output else None)
    ids = [s.id for s in stations]
    visibility_arg = str(args.visibility_results) if args.visibility_results else ""
    cand_rows = out.candidate_rows(ranked, ids, visibility_arg, n)
    out.write_csv(paths["candidates_csv"], out.candidate_columns(ids), cand_rows)
    out.write_csv(paths["candidate_stations_csv"], ["candidate_id", "ranking_order", "station_id",
                                                    "files_in_lookup_envelope"],
                  [{"candidate_id": c.candidate_id, "ranking_order": c.order, "station_id": sid,
                    "files_in_lookup_envelope": c.file_counts[sid]} for c in ranked for sid in c.segment.station_ids])
    matched: Dict[int, int] = {}
    for c in ranked:
        matched[c.segment.segment_id] = matched.get(c.segment.segment_id, 0) + 1
    seg_rows = out.segment_rows(segments, ids, n, matched)
    for r, s in zip(seg_rows, segments):
        r["meets_min_common_duration"] = s.duration_s >= min_dur
    out.write_csv(paths["common_segments_csv"], out.segment_columns(ids) + ["meets_min_common_duration"], seg_rows)
    out.write_csv(paths["station_days_csv"], out.STATION_DAY_COLUMNS, out.station_day_rows(listings, ids))

    scan_s = (end - start).total_seconds()
    summary_rows = []
    for st in stations:
        ls = [listings[(st.id, d)] for d in days]
        covered = total_seconds(per_station[st.id])
        times = file_times[st.id]
        summary_rows.append({
            "station_id": st.id, "archive_name": st.archive_name, "days_queried": len(ls),
            "days_found": sum(L.status == STATUS_FOUND for L in ls), "days_empty": sum(L.status == STATUS_EMPTY for L in ls),
            "days_error": sum(L.status == STATUS_ERROR for L in ls), "days_unknown": sum(L.status == STATUS_UNKNOWN for L in ls),
            "files": len(times), "total_bytes": sum(L.total_bytes for L in ls),
            "candidate_intervals": len(per_station[st.id]), "possible_covered_s": round(covered, 6),
            "scan_interval_s": scan_s, "possible_coverage_fraction": round(covered / scan_s, 9),
            "earliest_filename_utc": format_utc(times[0], 6) if times else None,
            "latest_filename_utc": format_utc(times[-1], 6) if times else None,
            "inventory_complete": all(L.complete for L in ls),
        })
    shared_s = total_seconds([(s.start_utc, s.end_utc) for s in all_station_segments])
    summary_rows.append({
        "station_id": "ALL_STATIONS", "archive_name": None, "days_queried": len(days),
        "days_found": None, "days_empty": None, "days_error": None, "days_unknown": None,
        "files": None, "total_bytes": None, "candidate_intervals": len(all_station_segments),
        "possible_covered_s": round(shared_s, 6), "scan_interval_s": scan_s,
        "possible_coverage_fraction": round(shared_s / scan_s, 9), "earliest_filename_utc": None,
        "latest_filename_utc": None, "inventory_complete": all_complete,
    })
    out.write_csv(paths["station_summary_csv"], out.STATION_SUMMARY_COLUMNS, summary_rows)

    in_band = [c for c in full_candidates if c.window.signal_relevance in IN_BAND]
    in_band_union_s = total_seconds(merge_intervals((c.overlap_start_utc, c.overlap_end_utc) for c in in_band))
    notes: List[str] = []
    if not all_complete:
        notes.append("inventory INCOMPLETE: listing errors or unknown days; absence of overlap is NOT established "
                     "where listings are missing")
    if k < n:
        notes.append(f"--min-stations {k}: rows covered by fewer than {n} stations are DIAGNOSTICS, "
                     "not valid retrieval inputs")
    if any(L.covers_full_day is False for L in listings.values()):
        notes.append("some day listings were made before the day had settled (partial days)")
    summary = {
        "interval": [format_utc(start), format_utc(end)], "stations": ids,
        "mode": "STAGE1_MATCHING" if doc is not None else "INVENTORY_ONLY",
        "min_stations": k, "all_station_mode": k == n, "min_common_duration_s": min_dur,
        "n_segments": len(all_station_segments), "n_segments_min_duration": n_seg_min,
        "n_diagnostic_segments": len(segments) - len(all_station_segments),
        "shared_possible_s": shared_s,
        "diagnostic_possible_s": total_seconds([(s.start_utc, s.end_utc) for s in segments if not s.all_stations]), "n_windows_total": windows_total, "n_windows_after_filters": windows_kept,
        "n_candidates": len(full_candidates), "n_diagnostic_candidates": len(ranked) - len(full_candidates),
        "n_in_band_candidates": len(in_band), "in_band_candidate_union_s": in_band_union_s,
        "n_intersections_below_min_duration": n_short, "n_excluded_tle_age": n_tle_excluded,
        "inventory_complete": all_complete, "classification": classification, "exit_code": code,
        "visibility_results": visibility_arg, "notes": notes,
    }
    write_candidates_json(paths["candidates_json"], summary, cand_rows)

    snapshots = sorted(
        ({"station_id": L.station_id, "day": L.day.isoformat(), "prefix": L.prefix, "status": L.status,
          "source": L.source, "checksum": L.checksum, "cache_path": L.cache_path} for L in listings.values()),
        key=lambda d: (ids.index(d["station_id"]), d["day"]))
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION, "analysis": "archive_scan", "run_tag": run_tag,
        "run_timestamp_utc": format_utc(run_ts, 3), "command_line": list(argv),
        "arguments": {k_: v for k_, v in vars(args).items()},
        "software": {"tart_gnss": __version__, "python": sys.version.split()[0], **_git_state()},
        "evidence_class": EVIDENCE_CLASS, "coverage_confidence": COVERAGE_CONFIDENCE,
        "evidence_classes_not_claimed": ["VALIDATED_THREE_STATION_DATA", "TIMESTAMP_OVERLAP_VERIFIED", "RADIO_DETECTED"],
        "archive": {"host": host, "bucket": bucket, "access": "anonymous, read-only listing (no downloads)",
                    "layout": "<archive_name>/vis/<year>/<month>/<day>/vis_<UTC>.hdf (no zero padding)"},
        "interval": {"start_utc": format_utc(start, 3), "end_utc": format_utc(end, 3), "convention": "half-open [start, end)",
                     "days_listed": [d.isoformat() for d in days]},
        "stations": [{**st.to_dict(), **metadata_checks.get(st.id, {})} for st in stations],
        "station_file": str(station_file), "station_file_sha256": sha256_file(station_file),
        "station_registry": registry_meta, "alias_file": args.aliases,
        "stage1": None if doc is None else {
            "path": str(results_path), "sha256": stage1_sha, "run": doc["run"],
            "tle_catalogue": doc.get("tle_catalogue"), "n_windows": windows_total,
        },
        "assumptions": {
            "file_span_assumption_s": span,
            "file_span_meaning": "a file named T is ASSUMED to hold samples in [T - span, T]; a search envelope, "
                                 "not a verified exposure time or a proven upper bound",
            "filename_time": "last-sample time in previously inspected files (Stage 2)",
            "gaps": "candidate intervals are unioned only where they overlap or touch; gaps are never covered",
            "possible_overlap_duration_s": "planning estimate under the span assumption; verified_overlap_duration_s "
                                           "is always null (no native timestamps were read)",
            "s3_last_modified": "recorded only; never used as observation time",
            "gap_tolerance_s": args.gap_tolerance_s,
        },
        "filters": {"signal_relevance": relevance, "constellations": args.constellations,
                    "max_abs_tle_age_days": args.max_abs_tle_age_days, "min_common_duration_s": min_dur,
                    "min_stations": k, "window_exclusions": exclusions,
                    "n_intersections_below_min_duration": n_short, "n_excluded_tle_age": n_tle_excluded},
        "ranking": {"rule": RANKING_RULE, "duration_cap_s": args.duration_cap_s, "in_band": list(IN_BAND)},
        "cache": {"dir": str(args.cache_dir), "policy": args.cache_policy, "max_age_hours": args.cache_max_age_hours,
                  "schema_version": CACHE_SCHEMA_VERSION},
        "listing": {**stats.to_dict(), "inventory_complete": all_complete,
                    "status_counts": {s: sum(L.status == s for L in listings.values())
                                      for s in (STATUS_FOUND, STATUS_EMPTY, STATUS_ERROR, STATUS_UNKNOWN)},
                    "snapshots": snapshots},
        "gap_groups": [{**g, "start": format_utc(g["start"], 6), "end": format_utc(g["end"], 6)}
                       for g in gap_groups([(s.start_utc, s.end_utc) for s in all_station_segments],
                                           args.gap_tolerance_s)] if args.gap_tolerance_s > 0 else None,
        "retrieval_handoff": {
            "command": "tart-retrieve --visibility-results <stage1> --rank <rank> --dry-run",
            "note": out.RETRIEVAL_NOTE,
            "not_emitted": "no --max-duration command is emitted: it centres on the predicted window midpoint",
        },
        "summary": summary,
        "outputs": {name: str(p) for name, p in paths.items()},
        "exit_classification": classification, "exit_code": code,
    }
    out.write_manifest(paths["manifest"], manifest, run_tag)
    print(out.format_report(summary, summary_rows, ranked, args.top))
    print(f"\nOutputs: {paths['candidates_csv']} (+ companions); manifest {paths['manifest']}")
    if code == EXIT_INCOMPLETE:
        logger.warning("Scan completed with an INCOMPLETE inventory (%s)", classification)
    return code


def write_candidates_json(path: Path, summary: Dict[str, object], rows: List[Dict[str, object]]) -> Path:
    """Candidates JSON: the run summary plus every candidate row (station counts as a map)."""
    import json

    items = []
    for r in rows:
        item = dict(r)
        item["station_file_counts"] = json.loads(item["station_file_counts"])
        items.append(item)
    return out.write_json(path, {"evidence_class": EVIDENCE_CLASS, "summary": summary, "candidates": items})


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``tart-archive-scan``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    logging.Formatter.converter = time.gmtime
    fmt = logging.Formatter("%(asctime)sZ %(levelname)-7s %(name)s: %(message)s", "%Y-%m-%dT%H:%M:%S")
    root = logging.getLogger()
    if not root.handlers:
        console = logging.StreamHandler()
        console.setFormatter(fmt)
        root.addHandler(console)
    root.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    for noisy in ("urllib3", "h5py"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    try:
        return run(args, ["tart-archive-scan", *argv])
    except ArchiveScanInputError as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (ArchiveScanError, Stage2Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE


if __name__ == "__main__":
    raise SystemExit(main())
