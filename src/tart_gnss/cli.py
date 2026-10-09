"""Command-line interface: ``common-gnss`` / ``python common_gnss.py``.

Time convention: ``--start``/``--end`` are ISO-8601. A ``Z`` or explicit UTC
offset is honoured (and converted to UTC). A timezone-naive value is, by this
documented convention, interpreted as UTC and a log message says so.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence

from tart_gnss import __version__, output, propagation, provenance, search, tle
from tart_gnss.registry import RegistryError
from tart_gnss.stations import (
    StationConfigError,
    default_alias_file,
    default_station_file,
    load_aliases,
    load_catalogue_metadata,
    load_stations,
    resolve_stations,
)

logger = logging.getLogger("tart_gnss")

#: Number of stations required by the Stage-1 CLI (internals support any N).
STAGE1_N_STATIONS = 3

#: Default flag level for |t_obs - t_TLE| warnings, in days. This is a
#: configurable project convention, not an SGP4 accuracy bound.
DEFAULT_TLE_AGE_WARNING_DAYS = 3.0


class CLIError(Exception):
    """A user-input error reported with a clear message and exit status 2."""


def parse_utc(text: str, label: str = "time") -> datetime:
    """Parse an ISO-8601 time string into an aware UTC datetime.

    ``Z``/offsets are honoured; naive times are treated as UTC (documented
    CLI convention).

    Raises
    ------
    CLIError
        If the string is not valid ISO-8601.
    """
    raw = text.strip()
    candidate = raw[:-1] + "+00:00" if raw.endswith(("Z", "z")) else raw
    try:
        dt = datetime.fromisoformat(candidate)
    except ValueError as exc:
        raise CLIError(f"Invalid {label} '{text}': expected ISO-8601, e.g. 2026-10-01T00:00:00Z") from exc
    if dt.tzinfo is None:
        logger.info("%s '%s' has no timezone; interpreted as UTC (CLI convention)", label, text)
        dt = dt.replace(tzinfo=timezone.utc)
    elif dt.utcoffset() != timezone.utc.utcoffset(None):
        logger.info("%s '%s' converted to UTC", label, text)
    return dt.astimezone(timezone.utc)


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser."""
    p = argparse.ArgumentParser(
        prog="common-gnss",
        description=(
            "Find GNSS satellites simultaneously above a minimum elevation at three "
            "TART stations, with refined common-visibility windows ranked by duration."
        ),
        epilog="Times are UTC. Timezone-naive ISO-8601 inputs are interpreted as UTC.",
    )
    p.add_argument("--stations", nargs="+", required=True, metavar="ID",
                   help="exactly three stations: canonical TART ids (e.g. mu-udm) or aliases "
                        "(e.g. mauritius); see tart-stations-list")
    p.add_argument("--start", required=True, help="UTC ISO-8601 start, e.g. 2026-10-01T00:00:00Z")
    p.add_argument("--end", required=True, help="UTC ISO-8601 end")
    p.add_argument("--min-elevation", type=float, default=20.0, help="minimum elevation, deg (default 20)")
    p.add_argument("--step", type=float, default=10.0, help="coarse sampling interval, s (default 10)")
    p.add_argument("--tolerance", type=float, default=0.1,
                   help="window boundary refinement tolerance, s (default 0.1)")
    p.add_argument("--tle-file", help="local TLE catalogue to use exactly (recommended for reproducibility)")
    p.add_argument("--tle-dir", default="data/tle",
                   help="where downloaded TLE catalogues are archived (default data/tle)")
    p.add_argument("--constellations", nargs="+", metavar="NAME",
                   help=f"restrict to constellations: {' '.join(tle.CONSTELLATIONS)}")
    p.add_argument("--output", help="output CSV (default results/common_gnss_<YYYYMMDD>.csv)")
    p.add_argument("--station-file",
                   help="YAML/JSON station catalogue (default: packaged frozen TART registry snapshot)")
    p.add_argument("--aliases", help="station alias file (default: bundled station_aliases.yaml)")
    p.add_argument("--refresh-stations", action="store_true",
                   help="fetch the live TART registry first and use the refreshed catalogue "
                        "(written under --station-snapshot-dir); off by default for reproducibility")
    p.add_argument("--station-snapshot-dir", default="data/stations",
                   help="where --refresh-stations archives registry snapshots (default data/stations)")
    p.add_argument("--min-duration", type=float, default=0.0, help="minimum window duration, s")
    p.add_argument("--tle-age-warning-days", type=float, default=DEFAULT_TLE_AGE_WARNING_DAYS,
                   help="warn when |t_obs - t_TLE epoch| exceeds this, days (default 3; "
                        "a flag level, not an accuracy bound)")
    p.add_argument("--summary-rows", type=int, default=20, help="rows shown in the terminal summary")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def output_paths(csv_arg: Optional[str], start: datetime) -> dict:
    """Derive CSV, JSON and metadata paths from the ``--output`` argument."""
    csv_path = Path(csv_arg) if csv_arg else Path("results") / f"common_gnss_{start:%Y%m%d}.csv"
    if csv_path.suffix.lower() != ".csv":
        csv_path = csv_path.with_suffix(".csv")
    stem = csv_path.stem
    tag = stem[len("common_gnss_"):] if stem.startswith("common_gnss_") else stem
    return {
        "csv": csv_path,
        "json": csv_path.with_suffix(".json"),
        "metadata": csv_path.with_name(f"run_metadata_{tag}.json"),
    }


def refresh_stations(snapshot_dir: Path) -> Path:
    """Fetch the live registry and write a run-specific frozen catalogue.

    The catalogue is ``<snapshot_dir>/stations_registry_<UTC>.yaml`` and the
    raw snapshot is archived next to it, so the run stays reproducible.

    Raises
    ------
    CLIError
        If the registry cannot be fetched or parsed (no files are changed;
        the packaged frozen catalogue is named in the message).
    """
    from tart_gnss.registry_cli import refresh_station_catalogue
    from tart_gnss.registry import utc_stamp
    from tart_gnss.provenance import utc_now

    output = snapshot_dir / f"stations_registry_{utc_stamp(utc_now())}.yaml"
    try:
        result = refresh_station_catalogue(output=output, snapshot_dir=snapshot_dir,
                                           previous=default_station_file())
    except RegistryError as exc:
        raise CLIError(
            f"--refresh-stations failed: {exc}. Nothing was changed; the frozen catalogue "
            f"{default_station_file()} is still available (omit --refresh-stations to use it)"
        ) from exc
    if result.diff is not None and not result.diff.is_empty:
        logger.info("Registry changes relative to the packaged catalogue:\n%s", result.diff.format())
    return output


def run(args: argparse.Namespace, argv: Sequence[str]) -> int:
    """Execute a search for parsed arguments; returns the process exit code."""
    run_timestamp = provenance.utc_now()
    start = parse_utc(args.start, "--start")
    end = parse_utc(args.end, "--end")
    if end <= start:
        raise CLIError(f"--end ({args.end}) must be later than --start ({args.start})")
    if args.step <= 0:
        raise CLIError("--step must be positive")
    if args.tolerance <= 0:
        raise CLIError("--tolerance must be positive")
    if not -90 <= args.min_elevation <= 90:
        raise CLIError("--min-elevation must be within [-90, 90] degrees")
    if args.min_duration < 0:
        raise CLIError("--min-duration must be non-negative")
    if len(args.stations) != STAGE1_N_STATIONS:
        raise CLIError(f"Stage 1 requires exactly {STAGE1_N_STATIONS} stations; got {len(args.stations)}")

    if args.refresh_stations and args.station_file:
        raise CLIError("--refresh-stations and --station-file are mutually exclusive")
    if args.refresh_stations:
        station_file = refresh_stations(Path(args.station_snapshot_dir))
    else:
        station_file = Path(args.station_file) if args.station_file else default_station_file()
    alias_file = Path(args.aliases) if args.aliases else None
    try:
        aliases = load_aliases(alias_file)
        stations = resolve_stations(args.stations, load_stations(station_file), aliases)
        registry_meta = load_catalogue_metadata(station_file) or {}
    except StationConfigError as exc:
        raise CLIError(str(exc)) from exc
    logger.info("Station file: %s", station_file)
    if registry_meta:
        logger.info("Station registry: %s retrieved %s (sha256 %s)", registry_meta.get("source"),
                    registry_meta.get("retrieved_utc"), registry_meta.get("sha256"))
    for s in stations:
        logger.info(
            "Station %s (%s): lat %.7f deg, lon %.7f deg, alt %.2f m [source: %s]",
            s.id, s.name, s.latitude_deg, s.longitude_deg, s.altitude_m, s.source,
        )

    try:
        constellations = [tle.normalise_constellation(c) for c in args.constellations] if args.constellations else None
    except tle.TLEError as exc:
        raise CLIError(str(exc)) from exc

    if args.tle_file:
        tle_path = Path(args.tle_file)
    else:
        logger.info("No --tle-file given: downloading current GNSS TLEs from CelesTrak")
        tle_path = tle.download_gnss_catalogue(args.tle_dir, constellations or tle.CONSTELLATIONS)
    catalogue = tle.load_tle_file(tle_path)
    logger.info("TLE file: %s (sha256 %s)", catalogue.path, catalogue.sha256)
    logger.info("TLE source: %s; retrieved: %s", catalogue.source,
                tle.format_utc(catalogue.retrieved_utc) if catalogue.retrieved_utc else "unknown")
    records = tle.filter_by_constellation(catalogue.records, constellations)
    if not records:
        raise CLIError("No satellites left after constellation filtering")
    unknown = [r.name for r in records if r.constellation == tle.UNKNOWN_CONSTELLATION]
    if unknown:
        logger.warning("%d satellites have unknown constellation: %s", len(unknown), ", ".join(unknown[:10]))

    config = search.SearchConfig(
        start=start, end=end, min_elevation_deg=args.min_elevation, step_s=args.step,
        tolerance_s=args.tolerance, min_duration_s=args.min_duration,
    )
    logger.info("Interval %s .. %s; step %g s; min elevation %g deg; tolerance %g s",
                tle.format_utc(start), tle.format_utc(end), args.step, args.min_elevation, args.tolerance)
    age_summary = search.assess_tle_ages(records, config, args.tle_age_warning_days)

    result = search.find_common_windows(records, stations, config)
    ranked = output.rank_windows(result.windows)

    paths = output_paths(args.output, start)
    station_itrs = {s.id: propagation.station_itrs_km(s) for s in stations}
    output.write_csv(paths["csv"], ranked, len(stations), args.min_elevation)
    output.write_json(
        paths["json"],
        output.build_results_json(ranked, stations, station_itrs, config, catalogue, records),
    )
    meta = provenance.build_run_metadata(
        run_timestamp=run_timestamp, start=start, end=end, stations=stations,
        station_file=str(station_file), station_registry=registry_meta,
        alias_file=str(alias_file or default_alias_file()), min_elevation_deg=args.min_elevation, step_s=args.step,
        tolerance_s=args.tolerance, min_duration_s=args.min_duration, constellations=constellations,
        catalogue=catalogue, n_satellites_searched=len(records), n_windows=len(ranked),
        tle_age_summary=age_summary, command_line=list(argv),
        output_files={k: str(v) for k, v in paths.items()},
    )
    meta["failed_propagations"] = result.failed_satellites
    output.write_json(paths["metadata"], meta)
    logger.info("Outputs: %s, %s, %s", paths["csv"], paths["json"], paths["metadata"])

    print(output.format_summary(ranked, stations, config, args.summary_rows))
    print(f"\nOutputs:\n  {paths['csv']}\n  {paths['json']}\n  {paths['metadata']}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    """Console-script entry point. Returns 0 on success, 2 on input errors, 1 otherwise."""
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    logging.Formatter.converter = time.gmtime  # log timestamps in UTC
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)sZ %(levelname)-7s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    try:
        return run(args, ["common-gnss", *argv])
    except CLIError as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except (tle.TLEError, OSError) as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
