"""Stage-2 manifest validation, Measurement-Set loading, and time binning.

The Stage-2 ``manifest.json`` is the authoritative input. A station is only
processed when its Stage-2 Measurement Set is ``VALID``; a dataset whose
three-station validity is ``false`` is refused unless the caller explicitly
allows it (and Stage 3 then records ``stage2_three_station_valid: false``).

Measurement Sets are read with python-casacore (imported lazily). The ENU
baselines used by the forward model come from the ANTENNA table positions; the
``UVW`` column is read too and checked against the ``tart2ms`` convention
``UVW = r_ANTENNA1 - r_ANTENNA2`` with a zenith phase centre (``w ~ 0``). A
Measurement Set that breaks the convention is rejected rather than silently
reinterpreted.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage3.errors import (
    InvalidManifestError,
    InvalidVisibilityDataError,
    MissingMeasurementSetError,
)
from tart_gnss.stage3.geometry import ecef_to_enu, ecef_to_geodetic
from tart_gnss.stage3.visibility_model import VisibilitySamples
from tart_gnss.tle import sha256_file

logger = logging.getLogger(__name__)

#: Seconds from MJD 0 (1858-11-17) to the Unix epoch.
MJD_UNIX_OFFSET_S = 3506716800.0

#: Allowed mismatch between UVW and -(ENU baseline): absolute (m) + relative to the longest baseline.
UVW_CHECK_ABS_M = 0.02
UVW_CHECK_REL = 0.02


@dataclass
class StationInput:
    """Everything Stage 3 takes from the Stage-2 manifest for one station."""

    station_id: str
    tart_id: Optional[str]
    station_name: Optional[str]
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    ms_path: Path
    ms_path_in_manifest: str
    measurement_set: Dict[str, object]
    n_antennas: Optional[int]
    frequency_hz: List[float]
    first_time_utc: Optional[str]
    last_time_utc: Optional[str]
    calibration_status: Optional[str]
    calibration_applied: Optional[bool]
    calibration_fingerprints: List[str]
    zero_gain_antennas: List[int]
    stage2_errors: List[str]
    stage2_three_station_valid: bool
    calibration_warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable summary recorded in the Stage-3 manifest."""
        return {
            "station_id": self.station_id,
            "tart_id": self.tart_id,
            "station_name": self.station_name,
            "latitude_deg": self.latitude_deg,
            "longitude_deg": self.longitude_deg,
            "altitude_m": self.altitude_m,
            "measurement_set": str(self.ms_path),
            "measurement_set_in_manifest": self.ms_path_in_manifest,
            "n_antennas": self.n_antennas,
            "frequency_hz": self.frequency_hz,
            "first_time_utc": self.first_time_utc,
            "last_time_utc": self.last_time_utc,
            "calibration": {
                "status": self.calibration_status,
                "applied_by_tart2ms": self.calibration_applied,
                "n_solution_sets": len(self.calibration_fingerprints),
                "solution_fingerprints": self.calibration_fingerprints,
                "zero_gain_antennas": self.zero_gain_antennas,
                "warnings": self.calibration_warnings,
            },
            "stage2_errors": self.stage2_errors,
            "stage2_three_station_valid": self.stage2_three_station_valid,
        }


def load_stage2_manifest(path) -> Tuple[Dict[str, object], str]:
    """Read and structurally check a Stage-2 manifest; return it with its SHA-256.

    Raises
    ------
    InvalidManifestError
        If the file is unreadable, not a Stage-2 manifest, or lacks required sections.
    """
    path = Path(path)
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise InvalidManifestError(f"Cannot read Stage-2 manifest {path}: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("stage") != 2:
        raise InvalidManifestError(f"{path} is not a Stage-2 manifest (stage != 2)")
    missing = [k for k in ("three_station_valid", "satellite", "stations") if k not in doc]
    if missing:
        raise InvalidManifestError(f"{path}: missing {', '.join(missing)}")
    if not isinstance(doc["stations"], dict) or not doc["stations"]:
        raise InvalidManifestError(f"{path}: no stations")
    if doc["satellite"].get("norad_id") is None:
        raise InvalidManifestError(f"{path}: no target satellite (satellite.norad_id)")
    return doc, sha256_file(path)


def resolve_station_keys(manifest: Dict[str, object], requested: Optional[Sequence[str]]) -> List[str]:
    """Map requested names (station_id or tart_id, case-insensitive) to manifest keys.

    ``requested=None`` selects every station in the manifest.

    Raises
    ------
    InvalidManifestError
        If a name matches no station.
    """
    stations = manifest["stations"]
    if requested is None:
        return list(stations)
    lookup: Dict[str, str] = {}
    for key, st in stations.items():
        for alias in (key, st.get("station_id"), st.get("tart_id"), st.get("archive_name")):
            if alias:
                lookup[str(alias).lower()] = key
    out = []
    for name in requested:
        key = lookup.get(name.lower())
        if key is None:
            raise InvalidManifestError(
                f"Station '{name}' is not in the manifest (available: {', '.join(stations)})"
            )
        if key not in out:
            out.append(key)
    return out


def resolve_ms_path(
    ms_path: str,
    manifest_path: Path,
    run_directory: Optional[str] = None,
    data_root: Optional[Path] = None,
) -> Path:
    """Find a Measurement Set referenced (usually relatively) by a Stage-2 manifest.

    Candidates, in order: ``data_root / ms_path``; ``ms_path`` as given
    (absolute, or relative to the current directory); ``ms_path`` with the
    Stage-2 ``run_directory`` prefix replaced by the manifest's own directory;
    and the manifest directory's ancestors joined with ``ms_path``.

    Raises
    ------
    MissingMeasurementSetError
        If no candidate exists.
    """
    p = Path(ms_path)
    manifest_dir = manifest_path.resolve().parent
    candidates: List[Path] = []
    if data_root is not None:
        candidates.append(Path(data_root) / p)
    candidates.append(p)
    if run_directory and not p.is_absolute():
        try:
            candidates.append(manifest_dir / p.relative_to(run_directory))
        except ValueError:
            pass
    if not p.is_absolute():
        candidates += [parent / p for parent in manifest_dir.parents]
    for c in candidates:
        if c.is_dir():
            return c
    raise MissingMeasurementSetError(
        f"Measurement Set '{ms_path}' not found (tried: {', '.join(str(c) for c in candidates[:4])} ...). "
        "Use --data-root to point at the Stage-2 working directory."
    )


def station_input(
    manifest: Dict[str, object],
    key: str,
    manifest_path: Path,
    allow_invalid_stage2: bool = False,
    data_root: Optional[Path] = None,
    max_zero_gain_fraction: float = 0.25,
) -> StationInput:
    """Validate one manifest station and build its :class:`StationInput`.

    Raises
    ------
    InvalidManifestError
        If the three-station dataset is invalid and ``allow_invalid_stage2`` is
        false, or if this station's own Stage-2 Measurement Set is not valid.
    MissingMeasurementSetError
        If its Measurement Set cannot be found.
    """
    three_valid = bool(manifest.get("three_station_valid"))
    if not three_valid and not allow_invalid_stage2:
        reasons = "; ".join(manifest.get("validity_reasons") or []) or "no reason recorded"
        raise InvalidManifestError(
            f"Stage-2 dataset is not three-station valid ({reasons}). "
            "Pass --allow-invalid-stage2 to process individually valid stations."
        )
    st = manifest["stations"][key]
    station_id = st.get("station_id") or key
    if station_id != key:
        raise InvalidManifestError(f"Station key '{key}' disagrees with its station_id '{station_id}'")
    for coord in ("latitude_deg", "longitude_deg"):
        if st.get(coord) is None:
            raise InvalidManifestError(f"{key}: missing {coord}")
    ms = st.get("measurement_set") or {}
    if st.get("measurement_set_status") != "VALID" or not ms.get("valid") or not ms.get("path"):
        raise InvalidManifestError(
            f"{key}: Stage-2 Measurement Set status is {st.get('measurement_set_status')!r}, not VALID"
        )
    if ms.get("station_id") not in (None, station_id):
        raise InvalidManifestError(f"{key}: Measurement Set belongs to station {ms.get('station_id')!r}")
    if not ms.get("frequency_hz"):
        raise InvalidManifestError(f"{key}: no observing frequency recorded for the Measurement Set")
    if not ms.get("first_time_utc") or not ms.get("last_time_utc"):
        raise InvalidManifestError(f"{key}: no time coverage recorded for the Measurement Set")
    run_dir = (manifest.get("provenance") or {}).get("run_directory")
    path = resolve_ms_path(ms["path"], manifest_path, run_dir, data_root)

    obs = st.get("observation") or {}
    cov = obs.get("coverage") or {}
    zero_gain = [int(a) for a in cov.get("zero_gain_antennas") or []]
    fingerprints = list(cov.get("calibration_fingerprints") or [])
    applied = ms.get("calibration_applied_by_tart2ms")
    n_ant = ms.get("number_of_antennas")
    warnings = []
    if applied is False:
        warnings.append("calibration not applied (uncalibrated visibilities)")
    if len(fingerprints) > 1:
        warnings.append(f"gain/phase solutions changed during the observation ({len(fingerprints)} solution sets)")
    if n_ant and len(zero_gain) / float(n_ant) > max_zero_gain_fraction:
        warnings.append(f"{len(zero_gain)}/{n_ant} antennas have zero gain (flagged)")
    if st.get("errors"):
        warnings.append("Stage-2 errors recorded for this station")
    return StationInput(
        station_id=station_id,
        tart_id=st.get("tart_id"),
        station_name=st.get("station_name"),
        latitude_deg=float(st["latitude_deg"]),
        longitude_deg=float(st["longitude_deg"]),
        altitude_m=float(st.get("altitude_m") or 0.0),
        ms_path=path,
        ms_path_in_manifest=str(ms["path"]),
        measurement_set=ms,
        n_antennas=n_ant,
        frequency_hz=[float(f) for f in ms.get("frequency_hz") or []],
        first_time_utc=ms.get("first_time_utc"),
        last_time_utc=ms.get("last_time_utc"),
        calibration_status=obs.get("calibration_status"),
        calibration_applied=applied,
        calibration_fingerprints=fingerprints,
        zero_gain_antennas=zero_gain,
        stage2_errors=[str(e) for e in st.get("errors") or []],
        stage2_three_station_valid=three_valid,
        calibration_warnings=warnings,
    )


@dataclass
class VisibilityDataset:
    """One station's visibilities in the internal Stage-3 representation.

    Row arrays have length ``R`` (MS rows); ``vis``, ``flag`` and ``weight``
    have shape ``(R, n_chan)`` (the first correlation is used).
    """

    time_mjd_seconds: np.ndarray
    antenna1: np.ndarray
    antenna2: np.ndarray
    uvw_m: np.ndarray
    baseline_enu_m: np.ndarray
    vis: np.ndarray
    flag: np.ndarray
    weight: np.ndarray
    frequency_hz: np.ndarray
    station_id: str
    antenna_enu_m: np.ndarray
    antenna_flagged: np.ndarray
    ms_path: Optional[str] = None
    uvw_check: Dict[str, object] = field(default_factory=dict)

    @property
    def n_rows(self) -> int:
        """Number of MS rows."""
        return int(self.time_mjd_seconds.size)


def dataset_from_arrays(
    *,
    station_id: str,
    time_mjd_seconds: np.ndarray,
    antenna1: np.ndarray,
    antenna2: np.ndarray,
    uvw_m: np.ndarray,
    data: np.ndarray,
    flag: np.ndarray,
    weight: np.ndarray,
    frequency_hz: np.ndarray,
    antenna_positions_ecef_m: np.ndarray,
    latitude_deg: Optional[float] = None,
    longitude_deg: Optional[float] = None,
    flag_row: Optional[np.ndarray] = None,
    antenna_flag_row: Optional[np.ndarray] = None,
    ms_path: Optional[str] = None,
    check_uvw: bool = True,
) -> VisibilityDataset:
    """Validate raw MS-like arrays and build a :class:`VisibilityDataset`.

    ``data``/``flag`` may be (R, chan, corr) or (R, chan); ``weight`` may be
    (R,), (R, corr) or (R, chan[, corr]). The ENU frame is centred on the mean
    antenna position, at the given geodetic coordinates or, if absent, those of
    that mean position.

    Raises
    ------
    InvalidVisibilityDataError
        On inconsistent shapes, missing data, or a UVW convention mismatch.
    """
    t = np.asarray(time_mjd_seconds, dtype=float)
    a1 = np.asarray(antenna1, dtype=int)
    a2 = np.asarray(antenna2, dtype=int)
    n_rows = t.size
    if n_rows == 0:
        raise InvalidVisibilityDataError(f"{station_id}: Measurement Set has no rows")
    data = np.asarray(data)
    flag = np.asarray(flag, dtype=bool)
    if data.ndim == 3:
        if data.shape[2] > 1:
            logger.warning("%s: %d correlations present; using the first", station_id, data.shape[2])
        data, flag = data[:, :, 0], flag[:, :, 0]
    if data.ndim != 2 or data.shape[0] != n_rows or flag.shape != data.shape:
        raise InvalidVisibilityDataError(f"{station_id}: inconsistent DATA/FLAG shapes {data.shape}/{flag.shape}")
    freq = np.atleast_1d(np.asarray(frequency_hz, dtype=float)).ravel()
    if freq.size != data.shape[1] or not np.all(np.isfinite(freq)) or np.any(freq <= 0):
        raise InvalidVisibilityDataError(f"{station_id}: {freq.size} channel frequencies for {data.shape[1]} channels")
    w = np.asarray(weight, dtype=float)
    if w.ndim == 1:
        w = np.repeat(w[:, None], data.shape[1], axis=1)
    elif w.ndim == 2 and w.shape[1] != data.shape[1]:
        w = np.repeat(w[:, :1], data.shape[1], axis=1)  # WEIGHT is per correlation, not channel
    elif w.ndim == 3:
        w = w[:, :, 0]
    if w.shape != data.shape:
        raise InvalidVisibilityDataError(f"{station_id}: WEIGHT shape {np.shape(weight)} incompatible with DATA")
    if flag_row is not None:
        flag = flag | np.asarray(flag_row, dtype=bool)[:, None]

    pos = np.asarray(antenna_positions_ecef_m, dtype=float)
    n_ant = pos.shape[0]
    if a1.max(initial=0) >= n_ant or a2.max(initial=0) >= n_ant or min(a1.min(), a2.min()) < 0:
        raise InvalidVisibilityDataError(f"{station_id}: antenna index out of range for {n_ant} antennas")
    if latitude_deg is None or longitude_deg is None:
        latitude_deg, longitude_deg, _ = ecef_to_geodetic(pos.mean(axis=0))
    enu = ecef_to_enu(pos, latitude_deg, longitude_deg)
    baseline = enu[a2] - enu[a1]
    uvw = np.asarray(uvw_m, dtype=float)
    check = check_uvw_convention(uvw, baseline, a1 != a2)
    if check_uvw and not check["consistent"]:
        raise InvalidVisibilityDataError(
            f"{station_id}: UVW does not follow the tart2ms zenith convention UVW = r_ANTENNA1 - r_ANTENNA2 "
            f"(max |UVW + b| = {check['max_abs_uvw_plus_b_m']:.3g} m, max |w| = {check['max_abs_w_m']:.3g} m). "
            "The visibility phase convention cannot be trusted."
        )
    ant_flag = (np.zeros(n_ant, dtype=bool) if antenna_flag_row is None
                else np.asarray(antenna_flag_row, dtype=bool))
    return VisibilityDataset(
        time_mjd_seconds=t, antenna1=a1, antenna2=a2, uvw_m=uvw, baseline_enu_m=baseline,
        vis=data.astype(complex), flag=flag, weight=w, frequency_hz=freq, station_id=station_id,
        antenna_enu_m=enu, antenna_flagged=ant_flag, ms_path=ms_path, uvw_check=check,
    )


def check_uvw_convention(uvw_m: np.ndarray, baseline_enu_m: np.ndarray, cross: np.ndarray) -> Dict[str, object]:
    """Compare stored UVW with ``-(r_q - r_p)`` in ENU (zenith phase centre).

    Rows with an all-zero UVW (not computed) are ignored.
    """
    sel = cross & np.any(uvw_m != 0, axis=1)
    if not np.any(sel):
        return {"consistent": False, "n_rows_checked": 0, "max_abs_uvw_plus_b_m": float("nan"),
                "max_abs_uvw_minus_b_m": float("nan"), "max_abs_w_m": float("nan"),
                "convention": "UVW = r_ANTENNA1 - r_ANTENNA2 (zenith phase centre)"}
    u, b = uvw_m[sel], baseline_enu_m[sel]
    longest = float(np.linalg.norm(b, axis=1).max())
    tol = UVW_CHECK_ABS_M + UVW_CHECK_REL * longest
    plus = float(np.abs(u + b).max())
    minus = float(np.abs(u - b).max())
    max_w = float(np.abs(u[:, 2]).max())
    return {
        "convention": "UVW = r_ANTENNA1 - r_ANTENNA2 (zenith phase centre); model baseline b = -UVW",
        "consistent": bool(plus <= tol and max_w <= tol),
        "n_rows_checked": int(sel.sum()),
        "tolerance_m": tol,
        "max_abs_uvw_plus_b_m": plus,
        "max_abs_uvw_minus_b_m": minus,
        "max_abs_w_m": max_w,
    }


def load_station_ms(
    ms_path,
    station_id: str,
    latitude_deg: Optional[float] = None,
    longitude_deg: Optional[float] = None,
) -> VisibilityDataset:
    """Load and validate one station Measurement Set (requires python-casacore).

    Reads TIME, DATA, UVW, ANTENNA1/2, FLAG, FLAG_ROW, WEIGHT (SIGMA is read
    when WEIGHT is absent: ``w = 1/sigma^2``) and the ANTENNA,
    SPECTRAL_WINDOW, FIELD and OBSERVATION subtables.

    Raises
    ------
    MissingMeasurementSetError
        If ``ms_path`` does not exist.
    InvalidVisibilityDataError
        If the MS cannot be read or is inconsistent.
    """
    ms_path = Path(ms_path)
    if not ms_path.is_dir():
        raise MissingMeasurementSetError(f"Measurement Set not found: {ms_path}")
    try:
        from casacore.tables import table
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise InvalidVisibilityDataError(
            "python-casacore is required to read Measurement Sets (pip install 'tart-gnss[stage3]')"
        ) from exc
    try:
        with table(str(ms_path), ack=False) as main:
            cols = set(main.colnames())
            missing = [c for c in ("TIME", "DATA", "UVW", "ANTENNA1", "ANTENNA2", "FLAG") if c not in cols]
            if missing:
                raise InvalidVisibilityDataError(f"{ms_path}: missing columns {', '.join(missing)}")
            arrays = {c: main.getcol(c) for c in ("TIME", "DATA", "UVW", "ANTENNA1", "ANTENNA2", "FLAG")}
            arrays["FLAG_ROW"] = main.getcol("FLAG_ROW") if "FLAG_ROW" in cols else None
            if "WEIGHT" in cols:
                arrays["WEIGHT"] = main.getcol("WEIGHT")
            elif "SIGMA" in cols:
                sigma = main.getcol("SIGMA")
                arrays["WEIGHT"] = np.where(sigma > 0, 1.0 / np.where(sigma > 0, sigma, 1.0) ** 2, 0.0)
            else:
                arrays["WEIGHT"] = np.ones(main.nrows())
            if "DATA_DESC_ID" in cols and len(np.unique(main.getcol("DATA_DESC_ID"))) > 1:
                raise InvalidVisibilityDataError(f"{ms_path}: multiple spectral windows are not supported")
            keywords = set(main.getkeywords())
        for sub in ("ANTENNA", "SPECTRAL_WINDOW"):
            if sub not in keywords:
                raise InvalidVisibilityDataError(f"{ms_path}: missing {sub} subtable")
        with table(str(ms_path / "ANTENNA"), ack=False) as ant:
            positions = ant.getcol("POSITION")
            ant_flag = ant.getcol("FLAG_ROW") if "FLAG_ROW" in ant.colnames() else None
        with table(str(ms_path / "SPECTRAL_WINDOW"), ack=False) as spw:
            if spw.nrows() != 1:
                raise InvalidVisibilityDataError(f"{ms_path}: expected one spectral window, found {spw.nrows()}")
            freq = np.asarray(spw.getcol("CHAN_FREQ"), dtype=float).ravel()
        for sub in ("FIELD", "OBSERVATION"):
            if sub not in keywords:
                logger.warning("%s: no %s subtable", ms_path, sub)
    except InvalidVisibilityDataError:
        raise
    except (RuntimeError, OSError, KeyError) as exc:
        raise InvalidVisibilityDataError(f"Cannot read Measurement Set {ms_path}: {exc}") from exc
    return dataset_from_arrays(
        station_id=station_id, time_mjd_seconds=arrays["TIME"], antenna1=arrays["ANTENNA1"],
        antenna2=arrays["ANTENNA2"], uvw_m=arrays["UVW"], data=arrays["DATA"], flag=arrays["FLAG"],
        weight=arrays["WEIGHT"], frequency_hz=freq, antenna_positions_ecef_m=positions,
        latitude_deg=latitude_deg, longitude_deg=longitude_deg, flag_row=arrays["FLAG_ROW"],
        antenna_flag_row=ant_flag, ms_path=str(ms_path),
    )


@dataclass
class TimeBin:
    """MS rows falling in one time bin."""

    index: int
    start_mjd_seconds: float
    end_mjd_seconds: float
    rows: np.ndarray


#: Bin-edge conventions accepted by :func:`time_bins`.
BIN_ORIGINS = ("utc", "first-sample")


def time_bins(dataset: VisibilityDataset, bin_s: float, origin: str = "utc") -> Iterator[TimeBin]:
    """Yield the non-empty time bins of width ``bin_s`` seconds.

    ``origin="utc"`` (default): edges lie on a fixed absolute grid, at whole
    multiples of ``bin_s`` seconds of UTC MJD time, so they fall at whole
    multiples of ``bin_s`` after UTC midnight whenever ``bin_s`` divides
    86400 s. Every station therefore shares the same bins, which is what lets
    Stage 4 match stations in time. ``origin="first-sample"``: edges start at
    this dataset's first timestamp (the pre-v0.4.1 behaviour).

    Empty bins are skipped, so no timestamp is ever produced where there is
    no data.
    """
    if not bin_s > 0:
        raise ValueError("time bin must be positive")
    if origin not in BIN_ORIGINS:
        raise ValueError(f"unknown bin origin {origin!r}")
    t = dataset.time_mjd_seconds
    if origin == "utc":
        t0 = float(np.floor(t.min() / bin_s) * bin_s)
    else:
        t0 = float(t.min())
    idx = np.floor((t - t0) / bin_s + 1e-9).astype(np.int64)
    order = np.argsort(idx, kind="stable")
    uniq, starts = np.unique(idx[order], return_index=True)
    bounds = list(starts) + [order.size]
    for i, b in enumerate(uniq):
        rows = order[bounds[i]:bounds[i + 1]]
        yield TimeBin(int(b), t0 + b * bin_s, t0 + (b + 1) * bin_s, np.sort(rows))


@dataclass
class SampleSelection:
    """Valid samples of a time bin plus the selection bookkeeping."""

    samples: VisibilitySamples
    n_baselines_total: int
    n_baselines_used: int
    flagged_fraction: float
    n_times: int


def select_samples(
    dataset: VisibilityDataset,
    rows: np.ndarray,
    invalid_antennas: Sequence[int] = (),
) -> SampleSelection:
    """Keep valid cross-correlation samples of ``rows``.

    Excluded: auto-correlations; flagged samples/rows; baselines touching a
    flagged antenna or an ``invalid_antennas`` entry (Stage-2 zero-gain
    antennas); non-finite visibilities; zero, negative or non-finite weights.

    ``flagged_fraction`` is the fraction of cross-correlation samples in the
    bin that were excluded for any of these reasons.
    """
    a1, a2 = dataset.antenna1[rows], dataset.antenna2[rows]
    cross = a1 != a2
    bad_ant = dataset.antenna_flagged.copy()
    for a in invalid_antennas:
        if 0 <= a < bad_ant.size:
            bad_ant[a] = True
    vis = dataset.vis[rows]
    w = dataset.weight[rows]
    ok = (~dataset.flag[rows]) & np.isfinite(vis) & np.isfinite(w) & (w > 0)
    ok &= (cross & ~bad_ant[a1] & ~bad_ant[a2])[:, None]
    n_chan = vis.shape[1]
    total_cross = int(cross.sum()) * n_chan
    bl_all = a1.astype(np.int64) * 100000 + a2
    n_bl_total = int(np.unique(bl_all[cross]).size)
    n_bl_used = int(np.unique(bl_all[ok.any(axis=1)]).size)
    r, c = np.nonzero(ok)
    rr = rows[r]
    samples = VisibilitySamples(
        baseline_enu_m=dataset.baseline_enu_m[rr],
        frequency_hz=dataset.frequency_hz[c],
        vis=vis[r, c],
        weight=w[r, c],
        antenna1=dataset.antenna1[rr],
        antenna2=dataset.antenna2[rr],
        time_mjd_seconds=dataset.time_mjd_seconds[rr],
    )
    flagged = 1.0 - (r.size / total_cross) if total_cross else 1.0
    return SampleSelection(samples, n_bl_total, n_bl_used, float(flagged),
                           int(np.unique(dataset.time_mjd_seconds[rows]).size))


def mjd_seconds_to_unix(t_mjd_s) -> np.ndarray:
    """Convert MS TIME (MJD seconds, UTC) to Unix seconds."""
    return np.asarray(t_mjd_s, dtype=float) - MJD_UNIX_OFFSET_S
