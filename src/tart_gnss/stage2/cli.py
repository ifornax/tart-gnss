"""``retrieve-tart``: Stage-2 command-line interface and orchestration.

Exit status:

* ``0``: a validated three-station dataset was produced;
* ``3``: the run completed but the result is NOT a valid three-station
  observation (missing station, no overlap, data not retrieved, ...);
* ``2``: invalid input or selection;
* ``1``: an unexpected failure (e.g. the archive is unreachable for all stations).
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from tart_gnss import __version__
from tart_gnss.stage2 import archive as archive_mod
from tart_gnss.stage2 import download as download_mod
from tart_gnss.stage2 import manifest as manifest_mod
from tart_gnss.stage2 import measurement_set as ms_mod
from tart_gnss.stage2 import native as native_mod
from tart_gnss.stage2 import selection as sel
from tart_gnss.stage2.errors import (
    ArchiveUnavailableError,
    DownloadError,
    MeasurementSetConversionError,
    ObservationValidationError,
    Stage1SelectionError,
    Stage2Error,
)
from tart_gnss.stage2.synchronise import synchronise
from tart_gnss.stations import (
    Station,
    StationConfigError,
    default_station_file,
    load_aliases,
    load_stations,
    resolve_stations,
)
from tart_gnss.tle import format_utc, sha256_file

logger = logging.getLogger("tart_gnss.stage2")

DEFAULT_PADDING_S = 120.0
EXIT_VALID, EXIT_FAILURE, EXIT_USAGE, EXIT_NOT_VALID = 0, 1, 2, 3

#: Coordinate agreement (deg) required between Stage-1 results and the station file.
COORDINATE_MATCH_DEG = 1e-9


class _ListHandler(logging.Handler):
    """Keep formatted log records so they can be written to ``retrieval.log``."""

    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.lines: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))


@dataclass
class StationRun:
    """Mutable per-station state collected during a run."""

    station: Station
    stage1: sel.StationInfo
    query_status: str = "NOT_QUERIED"
    query_error: Optional[str] = None
    records: List[archive_mod.ObservationRecord] = field(default_factory=list)
    raw_files: List[Dict[str, object]] = field(default_factory=list)
    observations: List[native_mod.NativeObservation] = field(default_factory=list)
    coverage: Optional[native_mod.StationCoverage] = None
    ms: Optional[Dict[str, object]] = None
    ms_status: str = "NOT_REQUESTED"
    errors: List[str] = field(default_factory=list)


def build_parser() -> argparse.ArgumentParser:
    """Create the ``retrieve-tart`` argument parser."""
    p = argparse.ArgumentParser(
        prog="retrieve-tart",
        description="Stage 2: retrieve, validate and synchronise the TART observations of a "
                    "Stage-1 common GNSS visibility window, and convert each station to a Measurement Set.",
    )
    p.add_argument("--visibility-results", required=True, metavar="PATH",
                   help="Stage-1 detailed JSON (a Stage-1 CSV is accepted if its JSON sits next to it)")
    p.add_argument("--rank", type=int, help="Stage-1 rank to process")
    p.add_argument("--norad-id", type=int, help="select by NORAD catalogue number")
    p.add_argument("--window-index", type=int, help="1-based window of --norad-id, in time order")
    p.add_argument("--padding", type=float, default=DEFAULT_PADDING_S,
                   help="seconds added before and after the predicted window (default 120)")
    p.add_argument("--max-duration", type=float, default=None,
                   help="retrieve at most this many seconds, centred on the window midpoint "
                        "(for a first, short retrieval)")
    p.add_argument("--output-dir", default="data/observations", help="root output directory")
    p.add_argument("--station-file",
                   help="station catalogue for archive/API metadata (default: packaged frozen registry "
                        "snapshot); coordinates always come from the Stage-1 results")
    p.add_argument("--aliases", help="station alias file (default: bundled station_aliases.yaml)")
    p.add_argument("--download", action="store_true", help="download the native files")
    p.add_argument("--convert-ms", action="store_true", help="convert each station to a Measurement Set (needs --download)")
    p.add_argument("--keep-native", action="store_true",
                   help="keep native files after successful MS conversion (otherwise they are removed; "
                        "their hashes and provenance stay in the manifest)")
    p.add_argument("--allow-partial", action="store_true",
                   help="continue when a station has no data (the result stays marked invalid)")
    p.add_argument("--uncalibrated", action="store_true",
                   help="ask tart2ms NOT to apply the gains and phases stored in the native files")
    p.add_argument("--tart2ms", help="path to the tart2ms executable (default: auto-detect)")
    p.add_argument("--gap-factor", type=float, default=native_mod.DEFAULT_GAP_FACTOR,
                   help="a time step larger than this multiple of the median cadence counts as a gap (default 3)")
    p.add_argument("--dry-run", action="store_true",
                   help="select, resolve and query the archives only; download, convert and write nothing")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def resolve_run_stations(
    window: sel.Stage1Window, station_file: Path, aliases: Optional[Dict[str, str]] = None
) -> Tuple[List[Station], Dict[str, Dict[str, object]]]:
    """Match the Stage-1 stations against the current station catalogue.

    The catalogue supplies operational metadata only (archive name, API URL).
    Identity and coordinates are always taken from the Stage-1 results, so
    Stage-1 geometry is never silently recomputed with different coordinates
    and the station catalogue is never refreshed here. Stage-1 stations are
    looked up by ``tart_id`` first, then by Stage-1 id (aliases accepted, so
    results that use legacy names such as ``mauritius`` still resolve).

    If the catalogue's coordinates differ from the frozen Stage-1 values, a
    warning is logged and :func:`station_metadata_check` reports both.

    Returns
    -------
    stations, checks
        Stations to use (Stage-1 identity and coordinates, catalogue
        operational metadata) and Stage-1 id -> :func:`station_metadata_check`.

    Raises
    ------
    Stage1SelectionError
        If the number of stations is not 3 or a station is unknown.
    """
    if len(window.stations) != manifest_mod.N_STATIONS:
        raise Stage1SelectionError(
            f"Stage 2 needs exactly {manifest_mod.N_STATIONS} stations; the Stage-1 results have {len(window.stations)}"
        )
    try:
        catalogue = load_stations(station_file)
        found = []
        for s1 in window.stations:
            try:
                found += resolve_stations([s1.tart_id or s1.id], catalogue, aliases)
            except StationConfigError:
                if not s1.tart_id:
                    raise
                found += resolve_stations([s1.id], catalogue, aliases)
    except StationConfigError as exc:
        raise Stage1SelectionError(str(exc)) from exc
    stations = []
    checks: Dict[str, Dict[str, object]] = {}
    for cfg, s1 in zip(found, window.stations):
        checks[s1.id] = station_metadata_check(cfg, s1)
        if checks[s1.id]["station_metadata_changed"]:
            logger.warning(
                "Station '%s': coordinates in %s (%s, %s, %s) differ from the frozen Stage-1 values "
                "(%s, %s, %s); the Stage-1 values are used", s1.id, station_file, cfg.latitude_deg,
                cfg.longitude_deg, cfg.altitude_m, s1.latitude_deg, s1.longitude_deg, s1.altitude_m,
            )
        stations.append(replace(
            cfg, id=s1.id, name=s1.name, latitude_deg=s1.latitude_deg, longitude_deg=s1.longitude_deg,
            altitude_m=s1.altitude_m, tart_id=s1.tart_id or cfg.tart_id,
            source=str(s1.raw.get("source") or cfg.source),
        ))
    return stations, checks


def station_metadata_check(catalogue_station: Station, stage1: sel.StationInfo) -> Dict[str, object]:
    """Compare a catalogue station with the frozen Stage-1 station.

    Returns
    -------
    dict
        ``station_metadata_changed`` (true if any coordinate differs by more
        than :data:`COORDINATE_MATCH_DEG`) and, when changed, the catalogue's
        ``catalogue_coordinates`` and ``catalogue_source`` alongside the
        Stage-1 values kept in the manifest.
    """
    attrs = ("latitude_deg", "longitude_deg", "altitude_m")
    changed = any(abs(getattr(catalogue_station, a) - getattr(stage1, a)) > COORDINATE_MATCH_DEG for a in attrs)
    out: Dict[str, object] = {"station_metadata_changed": changed}
    if changed:
        out["catalogue_coordinates"] = {a: getattr(catalogue_station, a) for a in attrs}
        out["catalogue_source"] = catalogue_station.source
    return out


def run_directory(output_dir: Path, window: sel.Stage1Window) -> Path:
    """``<output>/<YYYYMMDD>/NORAD_<id>``, plus ``_w<k>`` if the satellite has several windows."""
    name = f"NORAD_{window.norad_id}"
    if window.n_windows_for_satellite > 1:
        name += f"_w{window.window_index}"
    return Path(output_dir) / f"{window.start_utc:%Y%m%d}" / name


def _write_stage1_window(run_dir: Path, window: sel.Stage1Window) -> Path:
    path = run_dir / "stage1_window.json"
    doc = {
        "stage1_result_file": str(window.results_path),
        "stage1_result_sha256": window.results_sha256,
        "rank": window.rank,
        "window_index": window.window_index,
        "run": window.run,
        "stations": [s.raw for s in window.stations],
        "window": window.window,
    }
    text = json.dumps(doc, indent=2) + "\n"
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))
        if old.get("window") != window.window or old.get("stations") != doc["stations"]:
            raise Stage1SelectionError(
                f"{path} holds a different Stage-1 window; use another --output-dir"
            )
        return path
    path.write_text(text, encoding="utf-8")
    return path


def _station_entry(sr: StationRun, q0: datetime, q1: datetime, window: sel.Stage1Window,
                   metadata_check: Dict[str, object]) -> Dict[str, object]:
    st = sr.station
    pred = window.window.get("stations", {}).get(st.id)
    return {
        "station_id": st.id,
        "tart_id": st.tart_id,
        "station_name": st.name,
        "latitude_deg": st.latitude_deg,
        "longitude_deg": st.longitude_deg,
        "altitude_m": st.altitude_m,
        "coordinate_source": st.source,
        **metadata_check,
        "archive_name": st.archive_name,
        "api_url": st.api_url,
        "archive_query_status": sr.query_status,
        "archive_query_error": sr.query_error,
        "archive_query_start": format_utc(q0, 3),
        "archive_query_end": format_utc(q1, 3),
        "archive_files_found": len(sr.records),
        "archive_records": [r.to_dict() for r in sr.records],
        "raw_files": sr.raw_files,
        "observation": {
            "files": [o.to_dict() for o in sr.observations],
            "coverage": sr.coverage.to_dict() if sr.coverage else None,
            "calibration_status": _calibration_status(sr),
            "calibration_source": "per-antenna 'gains'/'phases' datasets stored in each native HDF5 file "
                                  "(solutions in use on the telescope at recording time); not independently verified",
        },
        "measurement_set_status": sr.ms_status,
        "measurement_set": sr.ms,
        "predicted_geometry": pred,
        "errors": sr.errors,
    }


def _calibration_status(sr: StationRun) -> str:
    if not sr.observations:
        return "unknown (no validated native data)"
    if not all(o.has_gains and o.has_phases for o in sr.observations):
        return "gains/phases missing from some native files"
    n = len(sr.coverage.calibration_fingerprints) if sr.coverage else 0
    applied = sr.ms.get("calibration_applied_by_tart2ms") if sr.ms else None
    text = "native gains/phases present" + (" (constant)" if n == 1 else f" ({n} distinct solution sets)")
    if applied is True:
        text += "; applied by tart2ms (unverified)"
    elif applied is False:
        text += "; NOT applied (raw visibilities in MS)"
    return text


def format_summary(window: sel.Stage1Window, runs: Sequence[StationRun], sync: Dict[str, object],
                   valid: bool, reasons: Sequence[str], manifest_path: Optional[Path],
                   signal: str, dry_run: bool) -> str:
    """Human-readable terminal summary."""
    lines = ["TART GNSS — STAGE 2" + (" (DRY RUN)" if dry_run else ""), "",
             "Satellite:", f"  {window.satellite_name}",
             f"  NORAD: {window.norad_id}   Constellation: {window.constellation}   Signal: {signal}", "",
             "Predicted common window:",
             f"  {format_utc(window.start_utc, 0)} -- {format_utc(window.end_utc, 0)}"
             f"  ({window.duration_s / 60:.1f} min)", ""]
    header = f"{'Station':<34}{'Archive':<10}{'Files':>6}  {'Data interval (UTC)':<28}{'MS':<8}"
    lines += [header, "-" * len(header)]
    for sr in runs:
        if sr.coverage:
            interval = f"{sr.coverage.first_timestamp_utc:%H:%M:%S}-{sr.coverage.last_timestamp_utc:%H:%M:%S}"
        elif sr.records:
            interval = f"~{sr.records[0].filename_utc:%dT%H:%M}-{sr.records[-1].filename_utc:%dT%H:%M} (names)"
        else:
            interval = "-"
        status = {"FOUND": "OK", "NOT_FOUND": "MISSING", "ARCHIVE_ERROR": "ERROR"}.get(sr.query_status, sr.query_status)
        ms = {"VALID": "OK", "FAILED": "FAILED", "NOT_REQUESTED": "-", "SKIPPED": "skipped"}.get(sr.ms_status, sr.ms_status)
        lines.append(f"{sr.station.name[:33]:<34}{status:<10}{len(sr.records):>6}  {interval:<28}{ms:<8}")
    lines.append("")
    if sync.get("actual_common_start_utc"):
        lines += ["Actual three-station overlap:",
                  f"  start: {sync['actual_common_start_utc']}",
                  f"  end:   {sync['actual_common_end_utc']}",
                  f"  duration: {sync['actual_common_duration_seconds'] / 60:.2f} min",
                  f"  coverage of predicted window: {100 * sync['three_station_coverage_fraction']:.1f}%",
                  f"  coverage of retrieved interval: {100 * sync['three_station_retrieved_window_coverage_fraction']:.1f}%",
                  ""]
    lines.append(f"Three-station valid: {'YES' if valid else 'NO'}")
    for r in reasons:
        lines.append(f"  - {r}")
    if manifest_path:
        lines += ["", "Manifest:", f"  {manifest_path}"]
    return "\n".join(lines)


def run(args: argparse.Namespace, argv: Sequence[str],
        archive: Optional[archive_mod.ArchiveBackend] = None,
        converter: Callable[..., Dict[str, object]] = ms_mod.convert_to_ms,
        log_handler: Optional[_ListHandler] = None) -> int:
    """Execute Stage 2 for parsed arguments; returns the exit status.

    ``archive`` and ``converter`` can be injected (tests use mocks).
    """
    run_ts = datetime.now(timezone.utc)
    run_tag = run_ts.strftime("%Y%m%dT%H%M%SZ")
    if args.convert_ms and not args.download:
        raise Stage1SelectionError("--convert-ms requires --download")
    if args.padding < 0:
        raise Stage1SelectionError("--padding must be non-negative")

    doc, results_path = sel.load_stage1_results(args.visibility_results)
    window = sel.select_window(doc, results_path, args.rank, args.norad_id, args.window_index)
    signal, signal_note = sel.signal_relevance(window.constellation, window.satellite_name)
    logger.info("Selected rank %d: %s (NORAD %d, %s), window %d/%d, %s -- %s (%.1f min)",
                window.rank, window.satellite_name, window.norad_id, window.constellation,
                window.window_index, window.n_windows_for_satellite,
                format_utc(window.start_utc), format_utc(window.end_utc), window.duration_s / 60)
    logger.info("Signal relevance: %s (%s)", signal, signal_note)
    if signal not in ("GPS_L1", "GALILEO_E1", "BEIDOU_B1C"):
        logger.warning("%s is not known to transmit in the TART band; it may not be detectable in Stage 3",
                       window.satellite_name)
    if window.clipped:
        logger.warning("The Stage-1 window is clipped by the Stage-1 search interval; the real window may be longer")

    station_file = Path(args.station_file) if args.station_file else default_station_file()
    try:
        aliases = load_aliases(args.aliases)
    except StationConfigError as exc:
        raise Stage1SelectionError(str(exc)) from exc
    stations, metadata_checks = resolve_run_stations(window, station_file, aliases)
    for st in stations:
        logger.info("Station %s: TART id %s, archive '%s'", st.id, st.tart_id, st.archive_name)
    q0, q1 = sel.query_interval(window.start_utc, window.end_utc, args.padding, args.max_duration)
    logger.info("Archive query interval %s -- %s (padding %g s%s)", format_utc(q0), format_utc(q1), args.padding,
                f", limited to {args.max_duration:g} s" if args.max_duration else "")

    archive = archive or archive_mod.TartS3Archive()
    runs = [StationRun(station=st, stage1=s1) for st, s1 in zip(stations, window.stations)]
    for sr in runs:
        try:
            sr.records = archive_mod.find_observations(sr.station, q0, q1, archive)
            sr.query_status = "FOUND" if sr.records else "NOT_FOUND"
        except ArchiveUnavailableError as exc:
            sr.query_status, sr.query_error = "ARCHIVE_ERROR", str(exc)
            logger.error("%s: %s", sr.station.id, exc)
        logger.info("%-14s %s (%d files)", sr.station.id, sr.query_status, len(sr.records))
    if all(sr.query_status == "ARCHIVE_ERROR" for sr in runs):
        raise ArchiveUnavailableError("The TART archive could not be queried for any station")

    reasons: List[str] = []
    missing = [sr.station.id for sr in runs if sr.query_status != "FOUND"]
    if missing:
        reasons.append("no archive data for: " + ", ".join(missing))

    if args.dry_run:
        total = sum((r.size_bytes or 0) for sr in runs for r in sr.records)
        print(format_summary(window, runs, {}, False, reasons or ["dry run: no data retrieved or validated"],
                             None, signal, True))
        print(f"\nWould download {sum(len(sr.records) for sr in runs)} files "
              f"({total / 1e6:.1f} MB) into {run_directory(Path(args.output_dir), window)}")
        return EXIT_NOT_VALID

    run_dir = run_directory(Path(args.output_dir), window)
    run_dir.mkdir(parents=True, exist_ok=True)
    stage1_window_path = _write_stage1_window(run_dir, window)
    proceed = args.download and (not missing or args.allow_partial)
    if args.download and missing and not args.allow_partial:
        reasons.append("download skipped because a station has no data (use --allow-partial for diagnostics)")
    if not args.download:
        reasons.append("native data not retrieved (--download not given)")

    tart2ms_cmd = None
    if proceed and args.convert_ms:
        try:
            tart2ms_cmd = ms_mod.find_tart2ms(args.tart2ms) if converter is ms_mod.convert_to_ms else ["tart2ms"]
        except MeasurementSetConversionError as exc:
            raise Stage1SelectionError(str(exc)) from exc

    if proceed:
        for sr in runs:
            if not sr.records:
                continue
            raw_dir = run_dir / sr.station.id / "raw"
            for rec in sr.records:
                try:
                    raw = download_mod.download_observation(rec, raw_dir, archive)
                except DownloadError as exc:
                    sr.errors.append(str(exc))
                    logger.error("%s", exc)
                    continue
                entry = raw.to_dict()
                entry["retained"] = True
                try:
                    obs = native_mod.inspect_hdf5(raw.local_path, sr.station)
                except ObservationValidationError as exc:
                    entry["validation"] = f"INVALID: {exc}"
                    sr.errors.append(str(exc))
                    logger.error("%s", exc)
                    sr.raw_files.append(entry)
                    continue
                rec.start_utc, rec.end_utc = obs.first_timestamp_utc, obs.last_timestamp_utc
                if obs.last_timestamp_utc < q0 or obs.first_timestamp_utc > q1:
                    entry["validation"] = "valid; outside query interval (not used)"
                else:
                    entry["validation"] = "valid"
                    sr.observations.append(obs)
                sr.raw_files.append(entry)
            if sr.observations:
                sr.coverage = native_mod.summarise_coverage(sr.station.id, sr.observations, args.gap_factor)
                c = sr.coverage
                logger.info("%s: %d samples %s -- %s, cadence %.3f s, %d gaps, freq %s Hz, %s antennas",
                            sr.station.id, c.number_of_timestamps, format_utc(c.first_timestamp_utc),
                            format_utc(c.last_timestamp_utc), c.median_time_spacing_seconds or float("nan"),
                            len(c.gaps), c.frequencies_hz, c.n_antennas)
                if c.zero_gain_antennas:
                    logger.warning("%s: antennas %s have zero/non-finite gains; tart2ms flags their baselines",
                                   sr.station.id, c.zero_gain_antennas)
                if len(c.calibration_fingerprints) > 1:
                    logger.warning("%s: %d different gain/phase solution sets within the retrieved data",
                                   sr.station.id, len(c.calibration_fingerprints))
                if c.gaps:
                    logger.warning("%s: %d gaps longer than %.1f s", sr.station.id, len(c.gaps),
                                   c.gap_threshold_seconds or float("nan"))
                if len(c.frequencies_hz) > 1 or len(c.n_antennas) > 1:
                    sr.errors.append("inconsistent frequency or antenna count between files")
                if not all(o.frequency_consistent for o in sr.observations):
                    reasons.append(f"{sr.station.id}: observing frequency inconsistent with or missing from metadata")
            elif sr.records:
                reasons.append(f"{sr.station.id}: no valid native data in the query interval")

    with_data = {sr.station.id: sr.coverage for sr in runs if sr.coverage}
    predicted = (window.start_utc, window.end_utc)
    retrieved = (max(q0, window.start_utc), min(q1, window.end_utc))
    if retrieved[1] <= retrieved[0]:
        retrieved = predicted
    sync = synchronise(with_data, predicted, retrieved) if with_data else {
        "predicted_window_duration_s": window.duration_s, "stations": [],
        "station_coverage_fraction": {}, "actual_common_start_utc": None, "actual_common_end_utc": None,
        "actual_common_duration_seconds": 0.0, "three_station_coverage_fraction": 0.0,
        "station_pair_time_offset_statistics": {}, "no_common_time_reason": "no validated data",
    }
    if with_data and len(with_data) == manifest_mod.N_STATIONS:
        if sync.get("actual_common_start_utc"):
            logger.info("Actual three-station overlap %s -- %s (%.1f s); coverage of predicted window %.1f%%",
                        sync["actual_common_start_utc"], sync["actual_common_end_utc"],
                        sync["actual_common_duration_seconds"], 100 * sync["three_station_coverage_fraction"])
        else:
            reasons.append(str(sync["no_common_time_reason"]))
            logger.error("%s", sync["no_common_time_reason"])
    elif proceed:
        reasons.append(f"validated data from only {len(with_data)} of {manifest_mod.N_STATIONS} stations")

    if proceed and args.convert_ms:
        for sr in runs:
            if not sr.observations:
                sr.ms_status = "SKIPPED"
                continue
            ms_dir = run_dir / sr.station.id / "ms"
            ms_dir.mkdir(parents=True, exist_ok=True)
            files = [o.path for o in sorted(sr.observations, key=lambda o: o.first_timestamp_utc)]
            try:
                sr.ms = converter(files, ms_dir / f"{sr.station.id}.ms", sr.station,
                                  tart2ms_command=tart2ms_cmd, apply_calibration=not args.uncalibrated,
                                  expected_antennas=(sr.coverage.n_antennas[0] if sr.coverage and len(sr.coverage.n_antennas) == 1 else None),
                                  log_path=ms_dir / "tart2ms.log")
                sr.ms_status = "VALID"
            except MeasurementSetConversionError as exc:
                sr.ms_status = "FAILED"
                sr.errors.append(str(exc))
                reasons.append(f"{sr.station.id}: Measurement Set conversion failed")
                logger.error("%s", exc)
                continue
            if not args.keep_native:
                for entry in sr.raw_files:
                    p = Path(str(entry["local_path"]))
                    if p.exists():
                        p.unlink()
                        entry["retained"] = False
                logger.info("%s: native files removed after validated conversion (hashes kept)", sr.station.id)
    for sr in runs:
        if sr.errors and not any(r.startswith(sr.station.id) for r in reasons):
            reasons.append(f"{sr.station.id}: " + sr.errors[0])

    valid = not reasons and len(with_data) == manifest_mod.N_STATIONS and bool(sync.get("actual_common_start_utc"))
    if not valid and not reasons:
        reasons.append("three-station requirements not met")

    satellite = {
        "satellite_name": window.satellite_name, "norad_id": window.norad_id,
        "constellation": window.constellation, "signal_relevance": signal,
        "signal_relevance_note": signal_note,
        "tle_epoch": window.window["tle"]["epoch"], "tle_source": window.window["tle"]["source"],
        "tle_line1": window.window["tle"]["line1"], "tle_line2": window.window["tle"]["line2"],
    }
    prediction = {
        "predicted_common_start_utc": window.window["start_utc"],
        "predicted_common_end_utc": window.window["end_utc"],
        "predicted_duration_seconds": window.window["duration_seconds"],
        "predicted_window_clipped": window.clipped,
        "stage1_rank": window.rank, "stage1_window_index": window.window_index,
        "min_elevation_deg": window.run.get("min_elevation_deg"),
        "query_start_utc": format_utc(q0, 3), "query_end_utc": format_utc(q1, 3),
        "padding_seconds": args.padding, "max_duration_seconds": args.max_duration,
        "stage1_result_file": str(window.results_path),
        "stage1_result_sha256": window.results_sha256,
        "stage1_window_file": str(stage1_window_path),
        "satellite_position_itrs": window.window.get("satellite_position"),
        "line_of_sight_mid": window.window.get("line_of_sight_mid"),
    }
    station_entries = {sr.station.id: _station_entry(sr, q0, q1, window, metadata_checks[sr.station.id])
                       for sr in runs}
    provenance = {
        "run_tag": run_tag, "run_timestamp_utc": format_utc(run_ts, 3),
        "command_line": list(argv), "station_file": str(station_file),
        "station_file_sha256": sha256_file(station_file),
        "archive": {"host": archive_mod.ARCHIVE_HOST, "bucket": archive_mod.ARCHIVE_BUCKET,
                    "layout": "<archive_name>/vis/<year>/<month>/<day>/vis_<UTC>.hdf",
                    "file_span_margin_s": archive_mod.DEFAULT_FILE_SPAN_MARGIN_S},
        "gap_factor": args.gap_factor, "download": args.download, "convert_ms": args.convert_ms,
        "keep_native": args.keep_native, "allow_partial": args.allow_partial,
        "tart2ms_command": tart2ms_cmd, "run_directory": str(run_dir),
        "retrieval_log": str(run_dir / "retrieval.log"),
    }
    doc_m = manifest_mod.build_manifest(
        three_station_valid=valid, validity_reasons=reasons, satellite=satellite,
        prediction=prediction, stations=station_entries, synchronisation=sync,
        provenance=provenance, stage1_window=window.window,
    )
    manifest_path = manifest_mod.write_manifest(run_dir / "manifest.json", doc_m, run_tag)
    if log_handler is not None:
        with open(run_dir / "retrieval.log", "a", encoding="utf-8") as fh:
            fh.write(f"===== run {run_tag} =====\n" + "\n".join(log_handler.lines) + "\n")
    print(format_summary(window, runs, sync, valid, reasons, manifest_path, signal, False))
    return EXIT_VALID if valid else EXIT_NOT_VALID


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point for ``retrieve-tart``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    if args.rank is None and args.norad_id is None:
        print("error: select a window with --rank or --norad-id", file=sys.stderr)
        return EXIT_USAGE
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
    collector = _ListHandler()
    collector.setFormatter(fmt)
    root.addHandler(collector)
    try:
        return run(args, ["retrieve-tart", *argv], log_handler=collector)
    except Stage1SelectionError as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except (Stage2Error, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    finally:
        root.removeHandler(collector)


if __name__ == "__main__":
    raise SystemExit(main())
