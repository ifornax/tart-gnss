"""``tart-archive-scan --all-stations``: network-wide archive overlap discovery.

Lists each eligible station's UTC-day prefixes once (shared cache with the
explicit-station mode), builds filename-based candidate coverage per
station, sweeps the whole network once for intervals with at least
``--min-stations`` distinct physical sites, and ranks the station sets
("triads") of distinct sites that share contiguous candidate coverage. See
:mod:`tart_gnss.archive_scan.network` for the algorithm.

Nothing is downloaded; no Stage-1 geometry is computed or reused (a Stage-1
JSON belongs to one frozen station set, so ``--visibility-results`` is
rejected here). Ready-to-copy ``tart-common-gnss`` commands are emitted for
the top sets as PROPOSED next actions, not completed checks.

Exit status (as the explicit mode): ``0`` complete (even with zero sets),
``3`` inventory incomplete, ``2`` invalid input, ``1`` failure or every
listing failed.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import logging
import shlex
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from tart_gnss import __version__
from tart_gnss.archive_scan import COVERAGE_CONFIDENCE, EVIDENCE_CLASS
from tart_gnss.archive_scan import network as net
from tart_gnss.archive_scan import output as out
from tart_gnss.archive_scan.cache import CACHE_SCHEMA_VERSION, InventoryCache
from tart_gnss.archive_scan.errors import ArchiveScanInputError
from tart_gnss.archive_scan.intervals import Interval, count_files, merge_intervals, station_candidate_segments
from tart_gnss.archive_scan.inventory import (
    STATUS_EMPTY,
    STATUS_ERROR,
    STATUS_FOUND,
    STATUS_UNKNOWN,
    build_inventory,
    required_days,
)
from tart_gnss.stage2 import archive as archive_mod
from tart_gnss.stage2.archive import ArchiveBackend
from tart_gnss.stage2.synchronise import total_seconds
from tart_gnss.stations import (
    StationConfigError,
    default_station_file,
    load_aliases,
    load_catalogue_metadata,
    resolve_stations,
)
from tart_gnss.tle import format_utc, sha256_file

logger = logging.getLogger("tart_gnss.archive_scan")

EXIT_OK, EXIT_FAILURE, EXIT_USAGE, EXIT_INCOMPLETE = 0, 1, 2, 3
NETWORK_MANIFEST_SCHEMA_VERSION = "1.0"
DEFAULT_NETWORK_OUTPUT = "results/network_overlap.csv"
DEFAULT_MIN_SITES = 3
HANDOFF_STATUS = "PROPOSED_NOT_RUN"

TRIAD_COLUMNS = [
    "rank", "triad_id", "station_ids", "site_group_ids", "n_stations", "possible_shared_duration_s",
    "longest_contiguous_duration_s", "n_contiguous_segments", "n_short_segments_excluded",
    "first_shared_start_utc", "last_shared_end_utc", "longest_segment_start_utc", "longest_segment_end_utc",
    "min_pairwise_baseline_km", "median_pairwise_baseline_km", "max_pairwise_baseline_km",
    "min_pairwise_chord_km", "pairwise_baselines_km", "inventory_complete", "coverage_class",
    "verified_shared_duration_s", "stage1_command", "archive_scan_command", "handoff_status",
]
TRIAD_SEGMENT_COLUMNS = [
    "rank", "triad_id", "segment_index", "start_utc", "end_utc", "possible_duration_s", "verified_duration_s",
    "station_ids", "file_counts_per_station", "min_station_files", "co_active_additional_stations",
    "inventory_complete", "inventory_status", "coverage_confidence", "evidence_class",
]
ACTIVE_SEGMENT_COLUMNS = [
    "segment_id", "start_utc", "end_utc", "possible_duration_s", "verified_duration_s", "n_active_stations",
    "n_distinct_sites", "active_station_ids", "site_group_ids", "active_station_ids_json",
    "meets_min_common_duration", "inventory_complete", "coverage_confidence",
]
NETWORK_STATION_SUMMARY_COLUMNS = [
    "station_id", "name", "tart_id", "latitude_deg", "longitude_deg", "altitude_m", "archive_name",
    "archive_verified", "site_group_id", "days_queried", "days_found", "days_empty", "days_error", "days_unknown",
    "first_day_with_files", "last_day_with_files", "files", "total_bytes", "candidate_intervals",
    "possible_covered_s", "possible_coverage_fraction", "earliest_filename_utc", "latest_filename_utc",
    "inventory_complete",
]
EXCLUDED_COLUMNS = ["station_id", "tart_id", "reason", "detail", "archive_name", "archive_verified"]


def network_output_paths(output: Path) -> Dict[str, Path]:
    """Output paths derived from the main CSV.

    For ``results/network_overlap.csv`` the companions are
    ``network_overlap.json``, ``network_active_segments.csv``,
    ``network_triad_segments.csv``, ``network_station_days.csv``,
    ``network_station_summary.csv``, ``network_excluded_stations.csv`` and
    ``network_archive_scan_manifest.json`` (a stem not ending in
    ``_overlap`` is used whole as the prefix).
    """
    output = Path(output)
    stem = output.stem
    base = stem[: -len("_overlap")] if stem.endswith("_overlap") else stem
    d = output.parent
    return {
        "triads_csv": output,
        "triads_json": output.with_suffix(".json"),
        "active_segments_csv": d / f"{base}_active_segments.csv",
        "triad_segments_csv": d / f"{base}_triad_segments.csv",
        "station_days_csv": d / f"{base}_station_days.csv",
        "station_summary_csv": d / f"{base}_station_summary.csv",
        "excluded_stations_csv": d / f"{base}_excluded_stations.csv",
        "manifest": d / f"{base}_archive_scan_manifest.json",
    }


def _t(dt: Optional[datetime], decimals: int = 6) -> Optional[str]:
    return format_utc(dt, decimals) if dt else None


def _whole_seconds(a: datetime, b: datetime) -> Tuple[datetime, datetime]:
    """Round an interval outward to whole UTC seconds (Stage-1 search boundaries)."""
    lo = a.replace(microsecond=0)
    hi = b if b.microsecond == 0 else b.replace(microsecond=0) + timedelta(seconds=1)
    return lo, hi


def stage1_command(ids: Sequence[str], start: datetime, end: datetime, station_file: Optional[str]) -> str:
    """Proposed ``tart-common-gnss`` run for one station set and candidate interval (NOT run)."""
    lo, hi = _whole_seconds(start, end)
    parts = ["tart-common-gnss", "--stations", *ids, "--start", format_utc(lo, 0), "--end", format_utc(hi, 0),
             "--min-elevation", "20", "--step", "10", "--tle-file", "<contemporaneous-gnss-tle-catalogue>"]
    if station_file:
        parts += ["--station-file", station_file]
    return " ".join(shlex.quote(p) if not p.startswith("<") else p for p in parts)


def archive_scan_command(ids: Sequence[str], station_file: Optional[str]) -> str:
    """Proposed explicit-station follow-up scan on the new triad-specific Stage-1 JSON (NOT run)."""
    parts = ["tart-archive-scan", "--stations", *ids, "--visibility-results", "<new-stage1-json>"]
    if station_file:
        parts += ["--station-file", station_file]
    return " ".join(shlex.quote(p) if not p.startswith("<") else p for p in parts)


def classify_network(n_eligible: int, all_complete: bool, any_complete: bool, any_error: bool,
                     n_before_baseline: int, n_after_baseline: int) -> Tuple[str, int]:
    """Exit classification of a network scan.

    Precedence: no eligible stations; every listing failed; sets found;
    sets found but all removed by the baseline cut; no qualifying overlap
    (complete inventory) versus undetermined (incomplete inventory).
    """
    if n_eligible == 0:
        return "NO_ELIGIBLE_ARCHIVE_STATIONS", EXIT_OK
    if not any_complete:
        return "ARCHIVE_LISTING_FAILURE", EXIT_FAILURE if any_error else EXIT_INCOMPLETE
    code = EXIT_OK if all_complete else EXIT_INCOMPLETE
    if n_after_baseline:
        return ("TRIADS_FOUND" if all_complete else "TRIADS_FOUND_INVENTORY_INCOMPLETE"), code
    if n_before_baseline:
        return "THREE_SITE_OVERLAP_BUT_BASELINE_FILTER_EXCLUDES_ALL", code
    if all_complete:
        return "NO_THREE_SITE_CANDIDATE_OVERLAP_WITH_COMPLETE_INVENTORY", code
    return "OVERLAP_UNDETERMINED_DUE_TO_INCOMPLETE_INVENTORY", code


def _resolve_list(ids: Optional[Sequence[str]], catalogue, aliases, opt: str):
    if not ids:
        return None
    try:
        return resolve_stations(ids, catalogue, aliases)
    except StationConfigError as exc:
        raise ArchiveScanInputError(f"{opt}: {exc}") from exc


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class _ActivityIndex:
    """Which stations are active somewhere in an interval (bisection over each station's segments)."""

    def __init__(self, station_segments: Mapping[str, Sequence[Interval]]):
        self._segs = {s: list(v) for s, v in station_segments.items()}
        self._ends = {s: [b for _, b in v] for s, v in self._segs.items()}

    def active(self, a: datetime, b: datetime) -> Tuple[str, ...]:
        hits = []
        for sid in sorted(self._segs):
            i = bisect.bisect_right(self._ends[sid], a)  # first segment ending after a
            if i < len(self._segs[sid]) and self._segs[sid][i][0] < b:
                hits.append(sid)
        return tuple(hits)


def run_network(args: argparse.Namespace, argv: Sequence[str],
                archive: Optional[ArchiveBackend] = None,
                now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                sleep: Callable[[float], None] = time.sleep) -> int:
    """Execute a network-wide scan for parsed arguments; returns the exit status."""
    from tart_gnss.archive_scan.cli import _git_state, _parse_time

    run_ts = now()
    run_tag = run_ts.strftime("%Y%m%dT%H%M%SZ")
    if args.visibility_results:
        raise ArchiveScanInputError(
            "--visibility-results cannot be used with --all-stations: a Stage-1 JSON is frozen to one station "
            "set and cannot be applied to other triads. Run the network scan without it, then run tart-common-gnss "
            "and the explicit '--stations ... --visibility-results' scan for a chosen triad")
    if args.inventory_output:
        raise ArchiveScanInputError("--inventory-output is not used with --all-stations (see --output)")
    if args.signal_relevance or args.constellations or args.max_abs_tle_age_days is not None:
        raise ArchiveScanInputError("Stage-1 filters (--signal-relevance, --constellations, --max-abs-tle-age-days) "
                                    "need --visibility-results, which --all-stations does not accept")
    if not args.start or not args.end:
        raise ArchiveScanInputError("--start and --end are required with --all-stations")
    start, end = _parse_time(args.start, "--start"), _parse_time(args.end, "--end")
    if end <= start:
        raise ArchiveScanInputError("--end must be after --start")
    min_sites = DEFAULT_MIN_SITES if args.min_stations is None else args.min_stations
    if min_sites < 2:
        raise ArchiveScanInputError("--min-stations must be >= 2")
    output = Path(args.output or DEFAULT_NETWORK_OUTPUT)

    station_file = Path(args.station_file) if args.station_file else default_station_file()
    try:
        aliases = load_aliases(args.aliases)
        registry_meta = load_catalogue_metadata(station_file)
        catalogue, invalid = net.load_registry_tolerant(station_file)
    except StationConfigError as exc:
        raise ArchiveScanInputError(str(exc)) from exc
    include = _resolve_list(args.include_stations, catalogue, aliases, "--include-stations")
    exclude = _resolve_list(args.exclude_stations, catalogue, aliases, "--exclude-stations")
    if include and exclude:
        both = sorted({s.id for s in include} & {s.id for s in exclude})
        if both:
            raise ArchiveScanInputError(f"stations both included and excluded: {', '.join(both)}")
    eligible, excluded = net.select_eligible(catalogue, include, exclude)
    excluded = invalid + excluded
    for e in excluded:
        logger.info("Excluded %s: %s %s", e.station_id, e.reason, e.detail)
    by_id = {s.id: s for s in eligible}
    baselines = net.Baselines(by_id)
    site_of = net.site_groups(eligible, baselines, args.site_tolerance_km)
    n_sites = len(set(site_of.values()))
    logger.info("Registry %s: %d valid stations, %d eligible archive stations at %d distinct sites "
                "(site tolerance %g km)", station_file, len(catalogue), len(eligible), n_sites, args.site_tolerance_km)
    if eligible and n_sites < min_sites:
        raise ArchiveScanInputError(
            f"--min-stations {min_sites} exceeds the {n_sites} distinct physical sites among the eligible stations "
            f"({', '.join(sorted(by_id))}); lower --min-stations or widen --include-stations")

    span = args.file_span_assumption_s
    days = required_days(start, end, span)
    host, bucket = archive_mod.ARCHIVE_HOST, archive_mod.ARCHIVE_BUCKET
    cache = InventoryCache(args.cache_dir)
    listings, stats = build_inventory(
        eligible, days, (lambda: archive) if archive is not None else archive_mod.TartS3Archive,
        cache, args.cache_policy, args.cache_max_age_hours, host, bucket, __version__, now=now, sleep=sleep,
    )
    ids = [s.id for s in eligible]
    file_times: Dict[str, List[datetime]] = {}
    for sid in ids:
        by_key = {o.key: o.filename_utc for d in days for o in listings[(sid, d)].objects}
        file_times[sid] = sorted(by_key.values())
    per_station = {sid: station_candidate_segments(file_times[sid], span, (start, end)) for sid in ids}
    scanned = set(days)

    def complete_for(sids: Sequence[str], a: datetime, b: datetime) -> bool:
        need = [d for d in required_days(a, b, span) if d in scanned]
        return all(listings[(s, d)].complete for s in sids for d in need)

    # ---------------------------------------------------------------- sweep and station sets
    active = net.network_active_segments(per_station, site_of, min_sites)
    for seg in active:
        seg.inventory_complete = complete_for(ids, seg.start_utc, seg.end_utc)
    index = _ActivityIndex(per_station)
    sets = net.candidate_station_sets(active, site_of, min_sites)
    triads_all = net.build_triads(
        sets, per_station, site_of, baselines, args.min_common_duration,
        file_counter=lambda s, a, b: count_files(file_times[s], a, b, span),
        completeness=complete_for, co_active_of=index.active)
    cut = args.min_baseline_km
    for tr in triads_all:
        tr.passes_baseline = cut is None or tr.baseline_stats()[0] >= cut
    rejected = [tr for tr in triads_all if not tr.passes_baseline]
    ranked = net.rank_triads([tr for tr in triads_all if tr.passes_baseline])
    exported = ranked if args.max_triads is None else ranked[: args.max_triads]
    logger.info("Network sweep: %d active segments with >= %d sites; %d candidate station sets, %d with a "
                "contiguous segment >= %g s, %d after the baseline cut", len(active), min_sites, len(sets),
                len(triads_all), args.min_common_duration, len(ranked))

    all_complete = all(L.complete for L in listings.values())
    any_complete = any(L.complete for L in listings.values())
    any_error = any(L.status == STATUS_ERROR for L in listings.values())
    classification, code = classify_network(len(eligible), all_complete, any_complete, any_error,
                                            len(triads_all), len(ranked))

    # ---------------------------------------------------------------- tables
    paths = network_output_paths(output)
    sf_arg = args.station_file
    triad_rows, seg_rows = [], []
    for tr in exported:
        mn, med, mx = tr.baseline_stats()
        longest = tr.longest_segment
        triad_rows.append({
            "rank": tr.rank, "triad_id": tr.triad_id, "station_ids": ";".join(tr.station_ids),
            "site_group_ids": ";".join(tr.site_group_ids), "n_stations": len(tr.station_ids),
            "possible_shared_duration_s": round(tr.possible_shared_s, 6),
            "longest_contiguous_duration_s": round(tr.longest_s, 6), "n_contiguous_segments": len(tr.segments),
            "n_short_segments_excluded": tr.n_short_segments,
            "first_shared_start_utc": _t(tr.segments[0].start_utc), "last_shared_end_utc": _t(tr.segments[-1].end_utc),
            "longest_segment_start_utc": _t(longest.start_utc), "longest_segment_end_utc": _t(longest.end_utc),
            "min_pairwise_baseline_km": round(mn, 3), "median_pairwise_baseline_km": round(med, 3),
            "max_pairwise_baseline_km": round(mx, 3),
            "min_pairwise_chord_km": round(min(p[3] for p in tr.pair_baselines), 3),
            "pairwise_baselines_km": json.dumps({f"{a}|{b}": round(g, 3) for a, b, g, _ in tr.pair_baselines},
                                                sort_keys=True),
            "inventory_complete": tr.inventory_complete, "coverage_class": COVERAGE_CONFIDENCE,
            "verified_shared_duration_s": None,
            "stage1_command": stage1_command(tr.station_ids, longest.start_utc, longest.end_utc, sf_arg),
            "archive_scan_command": archive_scan_command(tr.station_ids, sf_arg),
            "handoff_status": HANDOFF_STATUS,
        })
        for s in tr.segments:
            seg_rows.append({
                "rank": tr.rank, "triad_id": tr.triad_id, "segment_index": s.segment_index,
                "start_utc": _t(s.start_utc), "end_utc": _t(s.end_utc), "possible_duration_s": round(s.duration_s, 6),
                "verified_duration_s": None, "station_ids": ";".join(tr.station_ids),
                "file_counts_per_station": json.dumps(s.file_counts, sort_keys=True),
                "min_station_files": min(s.file_counts.values()),
                "co_active_additional_stations": ";".join(s.co_active), "inventory_complete": s.inventory_complete,
                "inventory_status": "COMPLETE" if s.inventory_complete else "INVENTORY_INCOMPLETE",
                "coverage_confidence": COVERAGE_CONFIDENCE, "evidence_class": EVIDENCE_CLASS,
            })
    active_rows = [{
        "segment_id": a.segment_id, "start_utc": _t(a.start_utc), "end_utc": _t(a.end_utc),
        "possible_duration_s": round(a.duration_s, 6), "verified_duration_s": None,
        "n_active_stations": len(a.station_ids), "n_distinct_sites": a.n_sites,
        "active_station_ids": ";".join(a.station_ids), "site_group_ids": ";".join(a.site_group_ids),
        "active_station_ids_json": json.dumps(list(a.station_ids)),
        "meets_min_common_duration": a.duration_s >= args.min_common_duration,
        "inventory_complete": a.inventory_complete, "coverage_confidence": COVERAGE_CONFIDENCE,
    } for a in active]
    scan_s = (end - start).total_seconds()
    summary_rows = []
    for st in eligible:
        ls = [listings[(st.id, d)] for d in days]
        with_files = [L.day for L in ls if L.n_files]
        times = file_times[st.id]
        covered = total_seconds(per_station[st.id])
        summary_rows.append({
            "station_id": st.id, "name": st.name, "tart_id": st.tart_id, "latitude_deg": st.latitude_deg,
            "longitude_deg": st.longitude_deg, "altitude_m": st.altitude_m, "archive_name": st.archive_name,
            "archive_verified": st.archive_verified, "site_group_id": site_of[st.id], "days_queried": len(ls),
            "days_found": sum(L.status == STATUS_FOUND for L in ls), "days_empty": sum(L.status == STATUS_EMPTY for L in ls),
            "days_error": sum(L.status == STATUS_ERROR for L in ls),
            "days_unknown": sum(L.status == STATUS_UNKNOWN for L in ls),
            "first_day_with_files": with_files[0].isoformat() if with_files else None,
            "last_day_with_files": with_files[-1].isoformat() if with_files else None,
            "files": len(times), "total_bytes": sum(L.total_bytes for L in ls),
            "candidate_intervals": len(per_station[st.id]), "possible_covered_s": round(covered, 6),
            "possible_coverage_fraction": round(covered / scan_s, 9),
            "earliest_filename_utc": _t(times[0]) if times else None,
            "latest_filename_utc": _t(times[-1]) if times else None,
            "inventory_complete": all(L.complete for L in ls),
        })
    excluded_rows = [{"station_id": e.station_id, "tart_id": e.tart_id, "reason": e.reason, "detail": e.detail,
                      "archive_name": e.archive_name, "archive_verified": e.archive_verified} for e in excluded]

    active_union_s = total_seconds([(a.start_utc, a.end_utc) for a in active])
    triad_union_s = total_seconds(merge_intervals((s.start_utc, s.end_utc) for tr in ranked for s in tr.segments))
    notes: List[str] = [
        f"coverage is {COVERAGE_CONFIDENCE}: filename times under a {span:g} s file-span assumption; "
        "no HDF5 sample timestamp was read (verified durations are null)",
        "S3 absence means no archived file, not that a telescope was not observing",
    ]
    if not all_complete:
        notes.append("inventory INCOMPLETE: listing errors or unknown days; absence of overlap is NOT established "
                     "where listings are missing")
    if any(L.covers_full_day is False for L in listings.values()):
        notes.append("some day listings were made before the day had settled (partial days)")
    if cut is not None:
        notes.append(f"--min-baseline-km {cut:g} is a site-diversity cut, not a triangulation-quality metric; "
                     "geometry is evaluated by Stage 1/4")
    summary = {
        "mode": "NETWORK", "interval": [format_utc(start), format_utc(end)],
        "registry_stations": len(catalogue) + sum(e.reason == net.INVALID_COORDINATES for e in invalid),
        "eligible_stations": ids, "n_eligible_stations": len(ids), "n_distinct_sites": n_sites,
        "n_excluded_stations": len(excluded), "site_tolerance_km": args.site_tolerance_km,
        "min_stations": min_sites, "min_common_duration_s": args.min_common_duration, "min_baseline_km": cut,
        "prefixes_required": stats.prefixes_required, "live_listings": stats.live_listings,
        "live_list_calls": stats.live_list_calls, "cache_hits": stats.cache_hits, "listing_failures": stats.failures,
        "unknown_listings": stats.unknown,
        "n_active_segments": len(active), "active_possible_s": active_union_s,
        "n_active_segments_min_duration": sum(a.duration_s >= args.min_common_duration for a in active),
        "n_candidate_station_sets": len(sets), "n_triads_before_baseline": len(triads_all),
        "n_triads_rejected_by_baseline": len(rejected), "n_triads": len(ranked), "n_triads_exported": len(exported),
        "triad_union_possible_s": triad_union_s, "inventory_complete": all_complete,
        "classification": classification, "exit_code": code, "notes": notes,
    }

    out.write_csv(paths["triads_csv"], TRIAD_COLUMNS, triad_rows)
    out.write_csv(paths["triad_segments_csv"], TRIAD_SEGMENT_COLUMNS, seg_rows)
    out.write_csv(paths["active_segments_csv"], ACTIVE_SEGMENT_COLUMNS, active_rows)
    out.write_csv(paths["station_days_csv"], out.STATION_DAY_COLUMNS, out.station_day_rows(listings, ids))
    out.write_csv(paths["station_summary_csv"], NETWORK_STATION_SUMMARY_COLUMNS, summary_rows)
    out.write_csv(paths["excluded_stations_csv"], EXCLUDED_COLUMNS, excluded_rows)
    segs_by_triad: Dict[str, list] = {}
    for r in seg_rows:
        segs_by_triad.setdefault(r["triad_id"], []).append(
            {**r, "file_counts_per_station": json.loads(r["file_counts_per_station"])})
    out.write_json(paths["triads_json"], {
        "evidence_class": EVIDENCE_CLASS, "summary": summary,
        "triads": [{**r, "pairwise_baselines_km": json.loads(r["pairwise_baselines_km"]),
                    "segments": segs_by_triad.get(r["triad_id"], [])} for r in triad_rows],
        "baseline_rejected": [{"triad_id": tr.triad_id, "min_pairwise_baseline_km": round(tr.baseline_stats()[0], 3),
                               "possible_shared_duration_s": round(tr.possible_shared_s, 6)}
                              for tr in sorted(rejected, key=lambda t: t.station_ids)],
    })

    snapshots = sorted(
        ({"station_id": L.station_id, "day": L.day.isoformat(), "prefix": L.prefix, "status": L.status,
          "source": L.source, "checksum": L.checksum, "cache_path": L.cache_path} for L in listings.values()),
        key=lambda d: (d["station_id"], d["day"]))
    manifest = {
        "schema_version": NETWORK_MANIFEST_SCHEMA_VERSION, "analysis": "network_archive_discovery",
        "run_tag": run_tag, "run_timestamp_utc": format_utc(run_ts, 3), "command_line": list(argv),
        "arguments": dict(vars(args)),
        "software": {"tart_gnss": __version__, "python": sys.version.split()[0], **_git_state()},
        "evidence_class": EVIDENCE_CLASS, "coverage_confidence": COVERAGE_CONFIDENCE,
        "evidence_classes_not_claimed": ["VALIDATED_THREE_STATION_DATA", "TIMESTAMP_OVERLAP_VERIFIED", "RADIO_DETECTED"],
        "archive": {"host": host, "bucket": bucket, "access": "anonymous, read-only listing (no downloads)",
                    "layout": "<archive_name>/vis/<year>/<month>/<day>/vis_<UTC>.hdf (no zero padding)"},
        "interval": {"start_utc": format_utc(start, 3), "end_utc": format_utc(end, 3), "convention": "half-open [start, end)",
                     "days_listed": [d.isoformat() for d in days]},
        "station_file": str(station_file), "station_file_sha256": sha256_file(station_file),
        "station_registry": registry_meta, "alias_file": args.aliases,
        "eligible_stations": [{**s.to_dict(), "site_group_id": site_of[s.id]} for s in eligible],
        "excluded_stations": excluded_rows,
        "site_grouping": {"method": "single linkage on geodesic separation", "tolerance_km": args.site_tolerance_km,
                          "groups": {g: sorted(s for s, gg in site_of.items() if gg == g)
                                     for g in sorted(set(site_of.values()))}},
        "baselines": {"method": net.BASELINE_METHOD, "min_baseline_km": cut, "policy": "all pairwise >= cut",
                      "pairs_km": {f"{a}|{b}": {"geodesic": round(baselines.geodesic_km(a, b), 3),
                                                 "chord": round(baselines.chord_km(a, b), 3)}
                                   for tr in triads_all for a, b, _, _ in tr.pair_baselines}},
        "assumptions": {
            "file_span_assumption_s": span,
            "file_span_meaning": "a file named T is ASSUMED to hold samples in [T - span, T]; a search envelope, "
                                 "not a verified exposure time or a proven upper bound",
            "filename_time": "last-sample time in previously inspected files (Stage 2)",
            "gaps": "candidate intervals are unioned only where they overlap or touch; gaps are never covered",
            "station_sets": "a set's segments are the intersection of its OWN members' segments; co-located "
                            "siblings never lend coverage",
            "s3_last_modified": "recorded only; never used as observation time",
            "registry": "registry presence and archive_verified are configuration, not evidence of data",
        },
        "quality_policies": {"min_stations_counts": "distinct physical sites", "min_common_duration_s":
                             args.min_common_duration, "min_duration_applies_to": "each contiguous station-set segment",
                             "inventory": "incomplete listings flag segments and give exit 3"},
        "ranking": {"version": net.NETWORK_RANKING_VERSION, "rule": net.NETWORK_RANKING_RULE,
                    "max_triads": args.max_triads, "max_triads_note": "limits exported rows only"},
        "cache": {"dir": str(args.cache_dir), "policy": args.cache_policy, "max_age_hours": args.cache_max_age_hours,
                  "schema_version": CACHE_SCHEMA_VERSION},
        "listing": {**stats.to_dict(), "inventory_complete": all_complete,
                    "status_counts": {s: sum(L.status == s for L in listings.values())
                                      for s in (STATUS_FOUND, STATUS_EMPTY, STATUS_ERROR, STATUS_UNKNOWN)},
                    "snapshots": snapshots},
        "handoff": {
            "status": HANDOFF_STATUS,
            "steps": ["tart-common-gnss for the chosen triad and segment (contemporaneous TLE catalogue)",
                      "tart-archive-scan --stations <triad> --visibility-results <new Stage-1 JSON>",
                      "tart-retrieve --visibility-results <new Stage-1 JSON> --rank <rank> --dry-run",
                      "download and Stage-2 validation only for viable candidates; then Stages 3-7"],
            "warning": "never reuse a Stage-1 JSON of another station set; GNSS visibility is not asserted until "
                       "Stage 1 runs",
        },
        "summary": summary,
        "outputs": {name: str(p) for name, p in paths.items()},
        "output_sha256": {name: _sha256(p) for name, p in paths.items() if name != "manifest"},
        "exit_classification": classification, "exit_code": code,
    }
    out.write_manifest(paths["manifest"], manifest, run_tag)
    print(format_network_report(summary, summary_rows, exported, args.top, sf_arg))
    print(f"\nOutputs: {paths['triads_csv']} (+ companions); manifest {paths['manifest']}")
    if code == EXIT_INCOMPLETE:
        logger.warning("Network scan completed with an INCOMPLETE inventory (%s)", classification)
    return code


def _hms(seconds: float) -> str:
    s = int(round(seconds))
    return f"{s // 3600:d}:{s % 3600 // 60:02d}:{s % 60:02d}"


def format_network_report(summary: Mapping[str, object], station_rows: Sequence[Mapping[str, object]],
                          triads: Sequence[net.Triad], top: int = 10, station_file: Optional[str] = None) -> str:
    """Human-readable terminal report of a network scan (numbers come from ``summary``)."""
    k = summary["min_stations"]
    lines = [
        "TART NETWORK ARCHIVE DISCOVERY",
        f"Interval: {summary['interval'][0]} -- {summary['interval'][1]} (half-open)",
        f"Registry stations: {summary['registry_stations']}",
        f"Eligible archive stations: {summary['n_eligible_stations']} (excluded: {summary['n_excluded_stations']})",
        f"Distinct geographic sites: {summary['n_distinct_sites']} (site tolerance {summary['site_tolerance_km']:g} km)",
        f"S3 station/day prefixes: {summary['prefixes_required']}  listed live: {summary['live_listings']}  "
        f"cache hits: {summary['cache_hits']}  failures: {summary['listing_failures']}  "
        f"unknown: {summary['unknown_listings']}",
        "",
    ]
    header = f"{'Station':<17}{'Site':<22}{'Days with data':>15}{'Files':>9}  {'First filename UTC':<21}{'Last filename UTC':<21}{'Inventory':>10}"
    lines += [header, "-" * len(header)]
    for r in station_rows:
        first = (r["earliest_filename_utc"] or "-")[:19]
        last = (r["latest_filename_utc"] or "-")[:19]
        lines.append(f"{r['station_id']:<17}{r['site_group_id']:<22}{r['days_found']:>15}{r['files']:>9}  "
                     f"{first:<21}{last:<21}{'complete' if r['inventory_complete'] else 'INCOMPLETE':>10}")
    lines += [
        "",
        f"Network active intervals with >= {k} sites: {summary['n_active_segments']} "
        f"({summary['active_possible_s'] / 3600:.2f} h possible, no double counting)",
        f"Station sets of {k} distinct sites with a contiguous segment >= {summary['min_common_duration_s']:g} s: "
        f"{summary['n_triads_before_baseline']}",
    ]
    if summary["min_baseline_km"] is not None:
        lines.append(f"After all-pairs baseline >= {summary['min_baseline_km']:g} km: {summary['n_triads']} "
                     f"({summary['n_triads_rejected_by_baseline']} rejected)")
    lines.append(f"Ranked eligible triads: {summary['n_triads']} (exported: {summary['n_triads_exported']})")
    if triads:
        lines += ["", f"  {'#':>3}  {'Triad':<40}{'Continuous':>11}{'Total':>11}{'Min baseline':>14}  Inventory"]
        for tr in triads[:top]:
            lines.append(f"  {tr.rank:>3}  {tr.triad_id:<40}{_hms(tr.longest_s):>11}{_hms(tr.possible_shared_s):>11}"
                         f"{tr.baseline_stats()[0]:>11.0f} km  {'complete' if tr.inventory_complete else 'INCOMPLETE'}")
        best = triads[0].longest_segment
        lines += ["", "Proposed next action (NOT run; GNSS visibility not yet evaluated):",
                  "  " + stage1_command(triads[0].station_ids, best.start_utc, best.end_utc, station_file)]
    lines += ["", f"Coverage: {COVERAGE_CONFIDENCE}; timestamps not native-validated.",
              f"Inventory completeness: {'COMPLETE' if summary['inventory_complete'] else 'INCOMPLETE'}",
              f"Exit classification: {summary['classification']}"]
    for note in summary.get("notes", []):
        lines.append(f"Note: {note}")
    return "\n".join(lines)
