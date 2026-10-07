"""Converting native TART files to Measurement Sets with ``tart2ms``, and validating them.

Each station is converted on its own, giving one MS per station. Stations are
never combined into a single interferometer.

``tart2ms`` runs as a subprocess with an argument list (no shell string). It is
called as::

    tart2ms --hdf <files...> --ms <out.ms>
            --skip-online-source-catalog --skip-celestial-source-catalog
            [--uncalibrated]

The source catalogues are skipped because Stage 2 needs no MODEL data, and
skipping them avoids network access during conversion.

Calibration
-----------
By default (``applycal``) ``tart2ms`` applies the per-antenna ``gains`` and
``phases`` stored in each HDF5 file, i.e. the solutions in use on the
telescope when the data were recorded. It flags baselines to antennas with
zero or non-finite gains. Those solutions are recorded, not verified, here.
With ``apply_calibration=False`` the raw visibilities are stored
(``--uncalibrated``).
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

from tart_gnss.stage2 import _ms_inspect
from tart_gnss.stage2.errors import MeasurementSetConversionError
from tart_gnss.stations import Station
from tart_gnss.tle import format_utc

logger = logging.getLogger(__name__)

TART2MS_FIXED_ARGS = ("--skip-online-source-catalog", "--skip-celestial-source-catalog")


def find_tart2ms(explicit: Optional[str] = None) -> List[str]:
    """Locate the ``tart2ms`` executable.

    Search order: ``explicit``, the current interpreter's ``bin`` directory, then ``PATH``.

    Raises
    ------
    MeasurementSetConversionError
        If ``tart2ms`` cannot be found.
    """
    candidates = [explicit] if explicit else []
    candidates += [str(Path(sys.executable).with_name("tart2ms")), shutil.which("tart2ms") or ""]
    for c in candidates:
        if c and Path(c).is_file():
            return [c]
    raise MeasurementSetConversionError(
        "tart2ms not found. Install it (pip install 'tart-gnss[ms]') or pass --tart2ms PATH"
    )


def _python_for(tart2ms_cmd: Sequence[str]) -> Optional[str]:
    """Return the interpreter in the shebang of the tart2ms script, if any."""
    try:
        first = Path(tart2ms_cmd[0]).read_text(errors="ignore").splitlines()[0]
    except (OSError, IndexError):
        return None
    return first[2:].strip().split()[0] if first.startswith("#!") else None


def inspect_ms(path, python: Optional[str] = None) -> Dict[str, object]:
    """Summarise a Measurement Set.

    Uses python-casacore in this process if it is importable. Otherwise the
    standalone inspector runs under ``python`` (e.g. the tart2ms interpreter).

    Raises
    ------
    MeasurementSetConversionError
        If the MS cannot be read.
    """
    try:
        import casacore.tables  # noqa: F401
        in_process = True
    except ImportError:
        in_process = False
    try:
        if in_process:
            return _ms_inspect.summarise(str(path))
        if not python:
            raise MeasurementSetConversionError("python-casacore unavailable: cannot validate the MS")
        proc = subprocess.run(
            [python, _ms_inspect.__file__, str(path)], capture_output=True, text=True, check=False
        )
        if proc.returncode != 0:
            raise MeasurementSetConversionError(f"MS inspection failed: {proc.stderr.strip()[-2000:]}")
        return json.loads(proc.stdout)
    except (RuntimeError, OSError, ValueError) as exc:
        raise MeasurementSetConversionError(f"Cannot read Measurement Set {path}: {exc}") from exc


def validate_ms(summary: Dict[str, object], expected_antennas: Optional[int] = None) -> Dict[str, object]:
    """Check an MS summary and return the manifest metadata.

    Requirements: non-zero rows; TIME, UVW and DATA columns; ANTENNA and
    SPECTRAL_WINDOW tables. If ``expected_antennas`` is given, the ANTENNA
    table must match it.

    Raises
    ------
    MeasurementSetConversionError
        Listing every failed requirement.
    """
    problems = []
    if not summary.get("number_of_rows"):
        problems.append("no rows")
    for key, label in (("has_time", "TIME column"), ("has_uvw", "UVW column"), ("has_data", "DATA column"),
                       ("has_antenna_table", "ANTENNA table"),
                       ("has_spectral_window_table", "SPECTRAL_WINDOW table")):
        if not summary.get(key):
            problems.append(f"missing {label}")
    if expected_antennas is not None and summary.get("number_of_antennas") not in (None, expected_antennas):
        problems.append(f"{summary.get('number_of_antennas')} antennas, expected {expected_antennas}")
    if problems:
        raise MeasurementSetConversionError(f"Invalid Measurement Set {summary.get('path')}: " + "; ".join(problems))
    to_utc = lambda x: format_utc(datetime.fromtimestamp(x, tz=timezone.utc), 6) if x is not None else None  # noqa: E731
    return {
        "path": summary["path"],
        "valid": True,
        "number_of_rows": summary["number_of_rows"],
        "number_of_times": summary.get("number_of_times"),
        "number_of_antennas": summary.get("number_of_antennas"),
        "number_of_channels": summary.get("number_of_channels"),
        "frequency_hz": summary.get("channel_frequencies_hz"),
        "total_bandwidth_hz": summary.get("total_bandwidth_hz"),
        "first_time_utc": to_utc(summary.get("first_time_unix")),
        "last_time_utc": to_utc(summary.get("last_time_unix")),
        "flagged_fraction": summary.get("flagged_fraction"),
        "max_baseline_m": summary.get("max_baseline_m"),
        "telescope_name": summary.get("telescope_name"),
        "columns": summary.get("columns"),
    }


def convert_to_ms(
    input_files: Sequence[str],
    output_ms,
    station: Station,
    tart2ms_command: Optional[Sequence[str]] = None,
    apply_calibration: bool = True,
    expected_antennas: Optional[int] = None,
    log_path=None,
    runner: Callable[..., subprocess.CompletedProcess] = subprocess.run,
    inspector: Callable[..., Dict[str, object]] = inspect_ms,
) -> Dict[str, object]:
    """Convert native TART visibility files into one station Measurement Set.

    Parameters
    ----------
    input_files
        Native HDF5 files of *one* station, in time order.
    output_ms
        Destination MS. It must not exist yet; it is never overwritten.
    station
        Station metadata (for logging and provenance).
    tart2ms_command
        Command prefix for tart2ms (default: :func:`find_tart2ms`).
    apply_calibration
        Apply the gains and phases stored in the files (tart2ms default).
    expected_antennas
        Antenna count the MS must contain, from the native files.
    log_path
        If given, tart2ms stdout and stderr are written here. tart2ms also
        writes its own ``tart2ms.<time>.log`` next to the MS.
    runner, inspector
        Injectable for testing.

    Returns
    -------
    dict
        Validated MS metadata plus the conversion command and calibration settings.

    Raises
    ------
    MeasurementSetConversionError
        If the output exists, tart2ms fails, or the MS fails validation.
    """
    output_ms = Path(output_ms).resolve()
    if not input_files:
        raise MeasurementSetConversionError(f"No input files for station '{station.id}'")
    if output_ms.exists():
        raise MeasurementSetConversionError(f"{output_ms} already exists; refusing to overwrite it")
    output_ms.parent.mkdir(parents=True, exist_ok=True)
    cmd_prefix = list(tart2ms_command) if tart2ms_command else find_tart2ms()
    args = [*cmd_prefix, "--hdf", *[str(Path(f).resolve()) for f in input_files], "--ms", str(output_ms),
            *TART2MS_FIXED_ARGS]
    if not apply_calibration:
        args.append("--uncalibrated")
    logger.info("tart2ms for %s: %d files -> %s", station.id, len(input_files), output_ms)
    started = datetime.now(timezone.utc)
    try:
        # tart2ms writes its own tart2ms.<time>.log into the working directory:
        # run it next to the MS so that log stays with the station's data.
        proc = runner(args, capture_output=True, text=True, check=False, cwd=str(output_ms.parent))
    except OSError as exc:
        raise MeasurementSetConversionError(f"Cannot run tart2ms: {exc}") from exc
    if log_path is not None:
        Path(log_path).write_text(
            f"$ {' '.join(args)}\n\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n",
            encoding="utf-8",
        )
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-5:]
        raise MeasurementSetConversionError(
            f"tart2ms failed for {station.id} (exit {proc.returncode}): " + " | ".join(tail)
        )
    if not output_ms.is_dir():
        raise MeasurementSetConversionError(f"tart2ms reported success but {output_ms} was not created")
    summary = inspector(output_ms, python=_python_for(cmd_prefix))
    meta = validate_ms(summary, expected_antennas)
    meta.update({
        "station_id": station.id,
        "conversion_command": args,
        "conversion_started_utc": format_utc(started, 3),
        "calibration_applied_by_tart2ms": apply_calibration,
        "conversion_log": str(log_path) if log_path else None,
    })
    logger.info("Validated %s: %d rows, %s antennas, first %s, last %s",
                output_ms, meta["number_of_rows"], meta["number_of_antennas"],
                meta["first_time_utc"], meta["last_time_utc"])
    return meta
