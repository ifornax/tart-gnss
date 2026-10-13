"""Command-line tools for the TART station registry.

* ``tart-stations-update``: fetch the live registry, archive the raw snapshot,
  and write a normalised station catalogue (optionally diffing it against the
  previous one). Fails if the registry cannot be retrieved; existing files are
  then left untouched.
* ``tart-stations-list``: show the stations of a (frozen) catalogue.

Exit status: 0 success, 1 registry/network failure, 2 input/configuration error.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import yaml

from tart_gnss import __version__, registry
from tart_gnss.stations import (
    Station,
    StationConfigError,
    default_alias_file,
    default_station_file,
    load_aliases,
    load_stations,
    normalise_station_id,
)

logger = logging.getLogger("tart_gnss")

#: Default location of registry snapshots and refreshed catalogues.
DEFAULT_SNAPSHOT_DIR = Path("data/stations")
DEFAULT_CATALOGUE_NAME = "stations_registry.yaml"


@dataclass
class RefreshResult:
    """Outcome of :func:`refresh_station_catalogue`."""

    download: registry.RegistryDownload
    parsed: registry.ParsedRegistry
    stations: Dict[str, Station]
    document: Dict[str, object]
    snapshot: Optional[registry.SnapshotFiles]
    catalogue_path: Optional[Path]
    diff: Optional[registry.RegistryDiff]
    diff_path: Optional[Path]


def refresh_station_catalogue(
    *,
    registry_url: str = registry.DEFAULT_REGISTRY_URL,
    output: Path,
    snapshot_dir: Path,
    previous: Optional[Path] = None,
    verify_archives: bool = False,
    strict: bool = False,
    dry_run: bool = False,
    timeout_s: float = 30.0,
    displacement_warning_m: float = registry.DEFAULT_DISPLACEMENT_WARNING_M,
) -> RefreshResult:
    """Fetch the registry and (unless ``dry_run``) write snapshot + catalogue.

    Parameters
    ----------
    registry_url
        Registry endpoint.
    output
        Catalogue to write. If it exists it is the default ``previous``.
    snapshot_dir
        Where raw snapshots, sidecars and diffs are archived.
    previous
        Catalogue to diff against and to carry archive verifications from.
    verify_archives
        Query the archive for every station (otherwise previous verifications
        are carried over and other stations stay unverified).
    strict
        Abort on any malformed registry entry.
    dry_run
        Fetch, parse and diff only; write nothing.

    Returns
    -------
    RefreshResult

    Raises
    ------
    RegistryError
        On fetch, schema or (strict) station errors. Nothing is written
        before the registry has been fetched and parsed successfully.
    StationConfigError
        If ``previous`` exists but cannot be loaded.
    """
    download = registry.download_registry(registry_url, timeout_s)
    retrieved = registry.utc_iso(download.retrieved_utc)
    parsed = registry.parse_registry(download.payload, retrieved, source=registry_url, strict=strict)

    previous = previous if previous is not None else (output if output.is_file() else None)
    old_stations = load_stations(previous) if previous is not None else None
    archive_meta = registry.archive_metadata_from_catalogue(old_stations) if old_stations else {}
    if verify_archives:
        for tart_id, status in registry.verify_archives(list(parsed.stations)).items():
            # A failed query (no timestamp) must not erase an earlier recorded verification.
            if status.archive_verified_utc is not None or tart_id not in archive_meta:
                archive_meta[tart_id] = status
    stations = registry.registry_to_stations(parsed.stations, archive_meta, download.sha256)
    diff = (registry.diff_catalogues(old_stations, stations, displacement_warning_m)
            if old_stations is not None else None)

    snapshot = catalogue_path = diff_path = None
    stem = f"registry_{registry.utc_stamp(download.retrieved_utc)}"
    snapshot_name = None if dry_run else f"{stem}.json"
    document = registry.build_catalogue_document(
        parsed, stations, source_url=registry_url, retrieved_utc=retrieved,
        sha256=download.sha256, snapshot_file=snapshot_name,
    )
    if diff is not None:
        document["metadata"]["previous_catalogue"] = str(previous)
        document["metadata"]["changes_from_previous"] = diff.to_dict()
    if not dry_run:
        snapshot = registry.write_snapshot(download, snapshot_dir, len(stations))
        if diff is not None:
            diff_path = snapshot_dir / f"{stem}.diff.json"
            diff_path.write_text(json.dumps({"previous_catalogue": str(previous), **diff.to_dict()},
                                            indent=2) + "\n", encoding="utf-8")
        catalogue_path = registry.write_catalogue(output, document)
    return RefreshResult(download, parsed, stations, document, snapshot, catalogue_path, diff, diff_path)


def _setup_logging(verbose: bool) -> None:
    logging.Formatter.converter = time.gmtime
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)sZ %(levelname)-7s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )


# --------------------------------------------------------------------------
# tart-stations-update
# --------------------------------------------------------------------------


def build_update_parser() -> argparse.ArgumentParser:
    """Argument parser for ``tart-stations-update``."""
    p = argparse.ArgumentParser(
        prog="tart-stations-update",
        description="Fetch the TART network registry, archive the raw snapshot and write a "
                    "normalised station catalogue.",
    )
    p.add_argument("--registry-url", default=registry.DEFAULT_REGISTRY_URL, help="registry endpoint")
    p.add_argument("--output", type=Path, default=DEFAULT_SNAPSHOT_DIR / DEFAULT_CATALOGUE_NAME,
                   help=f"catalogue to write (default {DEFAULT_SNAPSHOT_DIR / DEFAULT_CATALOGUE_NAME}); "
                        "to update the packaged default use src/tart_gnss/config/stations_registry.yaml")
    p.add_argument("--snapshot-dir", type=Path, default=DEFAULT_SNAPSHOT_DIR,
                   help=f"raw snapshot directory (default {DEFAULT_SNAPSHOT_DIR})")
    p.add_argument("--aliases", type=Path, help="alias file to check against the new catalogue "
                                                "(default: bundled station_aliases.yaml)")
    p.add_argument("--diff", nargs="?", const="", metavar="CATALOGUE", default=None,
                   help="always print the changes relative to CATALOGUE (default baseline: --output "
                        "if it exists, else the packaged catalogue); non-empty changes are always "
                        "printed and recorded")
    p.add_argument("--verify-archives", action="store_true",
                   help="check the visibility archive for a <tart_id>/vis/ prefix per station")
    p.add_argument("--strict", action="store_true", help="fail on any malformed registry entry")
    p.add_argument("--dry-run", action="store_true", help="fetch, parse and diff; write nothing")
    p.add_argument("--timeout", type=float, default=30.0, help="HTTP timeout, s (default 30)")
    p.add_argument("--displacement-warning-m", type=float, default=registry.DEFAULT_DISPLACEMENT_WARNING_M,
                   help="warn when a station moves more than this, m (a flag level, not a rejection rule)")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def _check_aliases(aliases: Dict[str, str], stations: Dict[str, Station]) -> None:
    known = {normalise_station_id(k) for k in stations}
    for alias, target in sorted(aliases.items()):
        if normalise_station_id(target) not in known:
            logger.warning("Alias '%s' -> '%s' has no station in the new catalogue", alias, target)


def main_update(argv: Optional[List[str]] = None) -> int:
    """Entry point of ``tart-stations-update``."""
    args = build_update_parser().parse_args(argv)
    _setup_logging(args.verbose)
    # Baseline for the diff and for carried-over archive verifications.
    if args.diff:
        previous: Optional[Path] = Path(args.diff)
    else:
        previous = None if args.output.is_file() else default_station_file()
    try:
        aliases = load_aliases(args.aliases)
        print("Fetching TART registry...")
        result = refresh_station_catalogue(
            registry_url=args.registry_url, output=args.output, snapshot_dir=args.snapshot_dir,
            previous=previous, verify_archives=args.verify_archives, strict=args.strict,
            dry_run=args.dry_run, timeout_s=args.timeout,
            displacement_warning_m=args.displacement_warning_m,
        )
    except registry.RegistryError as exc:
        logger.error("%s", exc)
        kept = args.output if args.output.is_file() else default_station_file()
        print(f"error: {exc}\nNo files were changed. The frozen catalogue {kept} remains available.",
              file=sys.stderr)
        return 1
    except StationConfigError as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return 2

    _check_aliases(aliases, result.stations)
    skipped = result.parsed.diagnostics
    print(f"Found {result.parsed.entry_count} registry entries.")
    print(f"Validated {len(result.stations)} stations.")
    print(f"{len(skipped)} entries skipped due to invalid/missing coordinates or duplicate ids.")
    for d in skipped:
        print(f"  - {d.tart_id or '#' + str(d.index)}: {d.reason}")
    if args.verify_archives:
        n_arch = sum(1 for s in result.stations.values() if s.archive_verified)
        print(f"Archive prefix verified for {n_arch} of {len(result.stations)} stations.")
    if result.diff is not None and (args.diff is not None or not result.diff.is_empty):
        print(f"Changes relative to {result.document['metadata'].get('previous_catalogue')}:")
        print(result.diff.format())
    if args.dry_run:
        print("Dry run: nothing written.")
        return 0
    print(f"Wrote registry snapshot {result.snapshot.raw} (sha256 {result.download.sha256})")
    if result.diff_path:
        print(f"Wrote change record {result.diff_path}")
    print(f"Wrote normalised station catalogue {result.catalogue_path}")
    return 0


# --------------------------------------------------------------------------
# tart-stations-list
# --------------------------------------------------------------------------


def build_list_parser() -> argparse.ArgumentParser:
    """Argument parser for ``tart-stations-list``."""
    p = argparse.ArgumentParser(prog="tart-stations-list",
                                description="List the stations of a TART station catalogue.")
    p.add_argument("--station-file", type=Path, help="catalogue (default: packaged frozen registry snapshot)")
    p.add_argument("--aliases", type=Path, help="alias file (default: bundled station_aliases.yaml)")
    p.add_argument("--archive-only", action="store_true",
                   help="only stations with an archive_name (verified archive prefix)")
    p.add_argument("--format", choices=("table", "json", "yaml"), default="table")
    p.add_argument("--verbose", "-v", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def station_rows(stations: Sequence[Station], aliases: Dict[str, str]) -> List[Dict[str, object]]:
    """One dict per station: its fields plus the aliases that point to it."""
    reverse: Dict[str, List[str]] = {}
    for alias, target in aliases.items():
        reverse.setdefault(normalise_station_id(target), []).append(alias)
    rows = []
    for st in stations:
        keys = {normalise_station_id(st.id)}
        if st.tart_id:
            keys.add(normalise_station_id(st.tart_id))
        names = sorted({a for k in keys for a in reverse.get(k, [])})
        rows.append({**st.to_dict(), "aliases": names})
    return rows


def format_table(rows: Sequence[Dict[str, object]]) -> str:
    """Fixed-width table of stations."""
    header = f"{'ID':<17} {'Alias':<26} {'Site':<44} {'Lat':>12} {'Lon':>12} {'Alt(m)':>8}  Archive"
    lines = [header, "-" * len(header)]
    for r in rows:
        lines.append(
            f"{r['id']:<17} {','.join(r['aliases']) or '-':<26} {str(r['name'])[:44]:<44} "
            f"{r['latitude_deg']:>12.7f} {r['longitude_deg']:>12.7f} {r['altitude_m']:>8.2f}  "
            f"{r['archive_name'] or '-'}"
        )
    return "\n".join(lines)


def main_list(argv: Optional[List[str]] = None) -> int:
    """Entry point of ``tart-stations-list``."""
    args = build_list_parser().parse_args(argv)
    _setup_logging(args.verbose)
    path = args.station_file or default_station_file()
    try:
        catalogue = load_stations(path)
        aliases = load_aliases(args.aliases if args.aliases else default_alias_file())
    except StationConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    stations = [s for s in catalogue.values() if s.archive_name or not args.archive_only]
    rows = station_rows(stations, aliases)
    if args.format == "json":
        print(json.dumps(rows, indent=2))
    elif args.format == "yaml":
        print(yaml.safe_dump(rows, sort_keys=False), end="")
    else:
        print(f"Catalogue: {path}")
        print(format_table(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main_update())
