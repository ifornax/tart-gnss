"""Scanner tables, JSON, manifest and the terminal report.

Station-specific columns (``files_<station>``) are generated from the
selected stations, never hard-coded; the same counts are also written as a
JSON map per row and as the long-form ``*_candidate_stations.csv`` table.
All tables are written with the standard-library :mod:`csv` module and
sorted deterministically.
"""

from __future__ import annotations

import csv
import json
import shlex
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from tart_gnss.archive_scan import COVERAGE_CONFIDENCE, EVIDENCE_CLASS
from tart_gnss.archive_scan.intervals import CommonSegment
from tart_gnss.archive_scan.inventory import DayListing
from tart_gnss.archive_scan.matching import Candidate
from tart_gnss.tle import format_utc

#: Note attached to every emitted retrieval command.
RETRIEVAL_NOTE = (
    "full-window Stage-2 dry run; tart-retrieve selects files around the predicted Stage-1 window "
    "(--max-duration centres on the window midpoint, not on this archive overlap), so compare the "
    "overlap_start_utc/overlap_end_utc interval with the dry-run file listing"
)


def output_paths(output: Path, inventory_output: Optional[Path] = None) -> Dict[str, Path]:
    """All output paths, derived from the candidates CSV path.

    For ``results/archive_candidates.csv`` the companions are
    ``results/archive_candidates.json``, ``results/archive_common_segments.csv``,
    ``results/archive_station_days.csv``, ``results/archive_station_summary.csv``,
    ``results/archive_candidate_stations.csv`` and
    ``results/archive_scan_manifest.json`` (a stem not ending in
    ``_candidates`` is used whole as the prefix).
    """
    output = Path(output)
    stem = output.stem
    base = stem[: -len("_candidates")] if stem.endswith("_candidates") else stem
    d = output.parent
    return {
        "candidates_csv": output,
        "candidates_json": output.with_suffix(".json"),
        "candidate_stations_csv": d / f"{base}_candidate_stations.csv",
        "common_segments_csv": d / f"{base}_common_segments.csv",
        "station_days_csv": Path(inventory_output) if inventory_output else d / f"{base}_station_days.csv",
        "station_summary_csv": d / f"{base}_station_summary.csv",
        "manifest": d / f"{base}_scan_manifest.json",
    }


def _t(dt: Optional[datetime], decimals: int = 6) -> Optional[str]:
    return format_utc(dt, decimals) if dt else None


def file_column(station_id: str) -> str:
    """Column name of a station's file count (``mg-tana`` -> ``files_mg_tana``)."""
    return "files_" + station_id.replace("-", "_").replace(" ", "_")


def retrieval_command(visibility_results: str, rank: int) -> str:
    """A ready-to-copy, safe ``tart-retrieve`` dry run for one Stage-1 rank (whole window)."""
    return f"tart-retrieve --visibility-results {shlex.quote(visibility_results)} --rank {int(rank)} --dry-run"


def coverage_status(segment: CommonSegment, n_selected: int) -> str:
    """``FILENAME_CANDIDATE``, ``INVENTORY_INCOMPLETE``, ``TWO_STATION_DIAGNOSTIC`` or ``PARTIAL_NETWORK_DIAGNOSTIC``."""
    if not segment.all_stations:
        return "TWO_STATION_DIAGNOSTIC" if len(segment.station_ids) == 2 else "PARTIAL_NETWORK_DIAGNOSTIC"
    if not segment.inventory_complete:
        return "INVENTORY_INCOMPLETE"
    return COVERAGE_CONFIDENCE


def candidate_rows(candidates: Sequence[Candidate], station_ids: Sequence[str],
                   visibility_results: str, n_selected: int) -> List[Dict[str, object]]:
    """One output row per ranked candidate."""
    rows = []
    for c in candidates:
        w = c.window
        row: Dict[str, object] = {
            "candidate_id": c.candidate_id,
            "ranking_order": c.order,
            "rank_stage1": w.rank,
            "norad_id": w.norad_id,
            "satellite_name": w.satellite_name,
            "constellation": w.constellation,
            "signal_relevance": w.signal_relevance,
            "signal_relevance_reason": w.signal_reason,
            "window_index": w.window_index,
            "n_windows_for_satellite": w.n_windows_for_satellite,
            "predicted_start_utc": _t(w.start_utc, 3),
            "predicted_end_utc": _t(w.end_utc, 3),
            "predicted_duration_s": round(w.duration_s, 3),
            "overlap_start_utc": _t(c.overlap_start_utc),
            "overlap_end_utc": _t(c.overlap_end_utc),
            "possible_overlap_duration_s": round(c.duration_s, 6),
            "verified_overlap_duration_s": None,
            "segment_id": c.segment.segment_id,
            "stations": ";".join(c.segment.station_ids),
            "n_stations": len(c.segment.station_ids),
        }
        for sid in station_ids:
            row[file_column(sid)] = c.file_counts.get(sid, 0)
        row.update({
            "station_file_counts": json.dumps(c.file_counts, sort_keys=True),
            "tle_epoch_utc": w.tle_epoch_text,
            "tle_age_days": None if c.tle_age_days is None else round(c.tle_age_days, 6),
            "stage1_min_mid_elevation_deg": None if w.min_mid_elevation_deg is None else round(w.min_mid_elevation_deg, 3),
            "inventory_complete": c.segment.inventory_complete,
            "coverage_status": coverage_status(c.segment, n_selected),
            "coverage_confidence": COVERAGE_CONFIDENCE,
            "evidence_class": EVIDENCE_CLASS,
            "flags": ";".join(c.flags),
            "ranking_score_or_order": c.order,
            "retrieval_command": retrieval_command(visibility_results, w.rank),
            "retrieval_note": RETRIEVAL_NOTE,
        })
        rows.append(row)
    return rows


def candidate_columns(station_ids: Sequence[str]) -> List[str]:
    """Header of the candidates CSV (also used when there are no rows)."""
    head = ["candidate_id", "ranking_order", "rank_stage1", "norad_id", "satellite_name", "constellation",
            "signal_relevance", "signal_relevance_reason", "window_index", "n_windows_for_satellite",
            "predicted_start_utc", "predicted_end_utc", "predicted_duration_s", "overlap_start_utc",
            "overlap_end_utc", "possible_overlap_duration_s", "verified_overlap_duration_s", "segment_id",
            "stations", "n_stations"]
    tail = ["station_file_counts", "tle_epoch_utc", "tle_age_days", "stage1_min_mid_elevation_deg",
            "inventory_complete", "coverage_status", "coverage_confidence", "evidence_class", "flags",
            "ranking_score_or_order", "retrieval_command", "retrieval_note"]
    return head + [file_column(s) for s in station_ids] + tail


def segment_rows(segments: Sequence[CommonSegment], station_ids: Sequence[str], n_selected: int,
                 matched: Mapping[int, int]) -> List[Dict[str, object]]:
    """One row per shared candidate segment."""
    rows = []
    for s in segments:
        row: Dict[str, object] = {
            "segment_id": s.segment_id, "start_utc": _t(s.start_utc), "end_utc": _t(s.end_utc),
            "possible_duration_s": round(s.duration_s, 6), "station_ids": ";".join(s.station_ids),
            "n_stations": len(s.station_ids),
        }
        for sid in station_ids:
            row[file_column(sid)] = s.file_counts.get(sid, 0)
        row.update({
            "file_counts_per_station": json.dumps(s.file_counts, sort_keys=True),
            "inventory_status": "COMPLETE" if s.inventory_complete else "INVENTORY_INCOMPLETE",
            "coverage_confidence": COVERAGE_CONFIDENCE,
            "coverage_status": coverage_status(s, n_selected),
            "n_matched_candidates": matched.get(s.segment_id, 0),
        })
        rows.append(row)
    return rows


def segment_columns(station_ids: Sequence[str]) -> List[str]:
    """Header of the common-segments CSV."""
    return (["segment_id", "start_utc", "end_utc", "possible_duration_s", "station_ids", "n_stations"]
            + [file_column(s) for s in station_ids]
            + ["file_counts_per_station", "inventory_status", "coverage_confidence", "coverage_status",
               "n_matched_candidates"])


STATION_DAY_COLUMNS = ["station_id", "archive_name", "day", "prefix", "status", "n_files", "total_bytes",
                       "earliest_filename_utc", "latest_filename_utc", "n_ignored_objects", "n_duplicate_keys",
                       "source", "listing_started_utc", "listing_completed_utc", "covers_full_day",
                       "attempts", "checksum", "cache_path", "error"]


def station_day_rows(listings: Mapping[Tuple[str, object], DayListing], station_ids: Sequence[str]) -> List[Dict[str, object]]:
    """One row per station/day, in station then day order."""
    order = {sid: i for i, sid in enumerate(station_ids)}
    rows = []
    for (sid, day), L in sorted(listings.items(), key=lambda kv: (order.get(kv[0][0], 99), kv[0][1])):
        rows.append({
            "station_id": sid, "archive_name": L.archive_name, "day": day.isoformat(), "prefix": L.prefix,
            "status": L.status, "n_files": L.n_files, "total_bytes": L.total_bytes,
            "earliest_filename_utc": _t(L.earliest_utc), "latest_filename_utc": _t(L.latest_utc),
            "n_ignored_objects": len(L.ignored_keys), "n_duplicate_keys": L.n_duplicate_keys, "source": L.source,
            "listing_started_utc": _t(L.listing_started_utc), "listing_completed_utc": _t(L.listing_completed_utc),
            "covers_full_day": L.covers_full_day, "attempts": L.attempts, "checksum": L.checksum,
            "cache_path": L.cache_path, "error": L.error,
        })
    return rows


STATION_SUMMARY_COLUMNS = ["station_id", "archive_name", "days_queried", "days_found", "days_empty", "days_error",
                           "days_unknown", "files", "total_bytes", "candidate_intervals", "possible_covered_s",
                           "scan_interval_s", "possible_coverage_fraction", "earliest_filename_utc",
                           "latest_filename_utc", "inventory_complete"]


def write_csv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> Path:
    """Write rows (``None`` -> empty cell, booleans as ``True``/``False``)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(columns), extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})
    return path


def write_json(path: Path, doc: object) -> Path:
    """Write indented JSON with a trailing newline."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2, default=str) + "\n", encoding="utf-8")
    return path


def write_manifest(path: Path, doc: Mapping[str, object], run_tag: str) -> Path:
    """Write the scan manifest; an existing one is renamed (``<stem>.<old run tag>.json``), never overwritten."""
    path = Path(path)
    if path.exists():
        try:
            old_tag = json.loads(path.read_text(encoding="utf-8"))["run_tag"]
        except (OSError, ValueError, KeyError, TypeError):
            old_tag = f"replaced_{run_tag}"
        archived = path.with_name(f"{path.stem}.{old_tag}.json")
        n = 1
        while archived.exists():
            archived = path.with_name(f"{path.stem}.{old_tag}_{n}.json")
            n += 1
        path.rename(archived)
    return write_json(path, doc)


def format_report(summary: Mapping[str, object], station_rows: Sequence[Mapping[str, object]],
                  candidates: Sequence[Candidate], top: int) -> str:
    """Human-readable terminal report (numbers come from ``summary``; nothing is hard-coded)."""
    lines = ["TART ARCHIVE SCAN", f"Interval: {summary['interval'][0]} -- {summary['interval'][1]} (half-open)",
             f"Stations: {', '.join(summary['stations'])}", f"Mode: {summary['mode']}", ""]
    header = f"{'Station':<16}{'Days queried':>13}{'Days with data':>16}{'Empty':>7}{'Files':>9}{'Listing errors':>16}{'Unknown':>9}"
    lines += [header, "-" * len(header)]
    for r in station_rows:
        if r["station_id"] == "ALL_STATIONS":
            continue
        lines.append(f"{r['station_id']:<16}{r['days_queried']:>13}{r['days_found']:>16}{r['days_empty']:>7}"
                     f"{r['files']:>9}{r['days_error']:>16}{r['days_unknown']:>9}")
    lines.append("")
    k = summary["min_stations"]
    label = "all-station" if summary["all_station_mode"] else f"{k}-of-{len(summary['stations'])}-station DIAGNOSTIC"
    lines += [
        f"Possible {label} candidate segments: {summary['n_segments']}",
        f"Candidate segments >= {summary['min_common_duration_s']:g} s: {summary['n_segments_min_duration']}",
        f"Possible shared candidate time (no double counting): {summary['shared_possible_s'] / 3600:.2f} h",
    ]
    if summary["mode"] == "STAGE1_MATCHING":
        lines += [f"Stage-1 windows: {summary['n_windows_total']} total, {summary['n_windows_after_filters']} after filters",
                  f"Satellite x segment intersections >= {summary['min_common_duration_s']:g} s: {summary['n_candidates']}"
                  + (f" ({summary['n_excluded_tle_age']} more excluded by TLE age)" if summary.get("n_excluded_tle_age") else ""),
                  f"In-band satellite intersections: {summary['n_in_band_candidates']}"]
    lines += [f"Inventory completeness: {'COMPLETE' if summary['inventory_complete'] else 'INCOMPLETE'}",
              f"Coverage class: {COVERAGE_CONFIDENCE} (NOT NATIVE VALIDATED; no radio detection implied)",
              f"Exit classification: {summary['classification']}", ""]
    if candidates:
        lines.append(f"Top candidates (showing {min(top, len(candidates))} of {len(candidates)}; all are in the CSV):")
        lines.append(f"  {'#':>3} {'rank':>5} {'NORAD':>6} {'relevance':<11} {'overlap start -- end (UTC)':<43}"
                     f"{'possible s':>11} {'TLE age d':>10}  status")
        for c in candidates[:top]:
            age = "n/a" if c.tle_age_days is None else f"{c.tle_age_days:+.2f}"
            lines.append(f"  {c.order:>3} {c.window.rank:>5} {c.window.norad_id:>6} {c.window.signal_relevance:<11} "
                         f"{format_utc(c.overlap_start_utc, 0)} -- {format_utc(c.overlap_end_utc, 0)}  "
                         f"{c.duration_s:>10.0f} {age:>10}  {coverage_status(c.segment, len(summary['stations']))}")
        lines += ["", "Retrieval (manual; whole Stage-1 window, dry run):",
                  f"  {retrieval_command(summary['visibility_results'], candidates[0].window.rank)}",
                  "  Note: tart-retrieve targets the predicted window, not the archive overlap; "
                  "--max-duration centres on the window midpoint and can miss a short overlap."]
    for note in summary.get("notes", []):
        lines.append(f"Note: {note}")
    return "\n".join(lines)
