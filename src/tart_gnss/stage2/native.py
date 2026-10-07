"""Inspecting and validating native TART HDF5 visibility files.

The archive's HDF5 layout (checked against files from za-rhodes, mg-tana and
mu-udm, 2026) contains:

``config``
    a JSON string with ``name``, ``frequency`` (Hz), ``bandwidth`` (Hz),
    ``lat``/``lon``/``alt``, ``num_antenna``, ...
``timestamp``
    ``(n_t,)`` ISO-8601 UTC strings, about 1 s apart;
``vis``
    ``(n_t, n_baselines)`` complex visibilities;
``baselines``
    ``(n_baselines, 2)`` antenna index pairs;
``antenna_positions``
    ``(n_ant, 3)`` ENU positions, metres;
``gains``, ``phases``
    per-antenna calibration solutions recorded with the data;
``phase_elaz``
    the phase centre (elevation, azimuth), degrees.

No separate per-visibility flag table exists in these files. ``tart2ms``
derives flags from zero or non-finite gains, and that is reported here too.

Validation reads the data itself. Fatal inconsistencies raise
:class:`ObservationValidationError`; minor anomalies are recorded as warnings.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage2.errors import ObservationValidationError
from tart_gnss.stage2.selection import TART_CENTRE_FREQUENCY_HZ
from tart_gnss.stations import Station
from tart_gnss.tle import format_utc

logger = logging.getLogger(__name__)

#: Station-identity gross-error check: file coordinates further than this (deg)
#: from the configured station mean the data belong to another site.
IDENTITY_TOLERANCE_DEG = 0.01

#: A time step larger than this multiple of the median cadence counts as a gap.
DEFAULT_GAP_FACTOR = 3.0


@dataclass
class NativeObservation:
    """Metadata and validation results for one native TART file.

    Times are aware UTC datetimes; frequencies in Hz; positions in degrees/metres.
    """

    path: str
    station_id: str
    timestamps: List[datetime] = field(repr=False)
    n_timestamps: int
    first_timestamp_utc: datetime
    last_timestamp_utc: datetime
    median_spacing_s: Optional[float]
    telescope_name: Optional[str]
    n_antennas: Optional[int]
    n_baselines: int
    vis_shape: Tuple[int, ...]
    vis_finite_fraction: float
    vis_all_zero: bool
    frequency_hz: Optional[float]
    bandwidth_hz: Optional[float]
    frequency_consistent: Optional[bool]
    file_latitude_deg: Optional[float]
    file_longitude_deg: Optional[float]
    file_altitude_m: Optional[float]
    position_offset_deg: Optional[float]
    has_gains: bool
    has_phases: bool
    zero_gain_antennas: List[int]
    calibration_fingerprint: Optional[str]
    phase_centre_el_az_deg: Optional[List[float]]
    flags_present: bool
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """JSON-serialisable summary (the full timestamp list is left out)."""
        d = asdict(self)
        d.pop("timestamps")
        d["first_timestamp_utc"] = format_utc(self.first_timestamp_utc, 6)
        d["last_timestamp_utc"] = format_utc(self.last_timestamp_utc, 6)
        d["vis_shape"] = list(self.vis_shape)
        return d


def _decode(value) -> str:
    return value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)


def _parse_timestamp(text: str) -> datetime:
    dt = datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    if dt.tzinfo is None:
        raise ObservationValidationError(f"Timezone-naive timestamp '{text}' in native file")
    return dt.astimezone(timezone.utc)


def inspect_hdf5(path, station: Station) -> NativeObservation:
    """Read and validate a native TART HDF5 visibility file.

    Parameters
    ----------
    path
        HDF5 file.
    station
        The station the file is expected to belong to.

    Returns
    -------
    NativeObservation

    Raises
    ------
    ObservationValidationError
        If the file is unreadable, timestamps are missing, unparsable or not
        increasing, visibilities are missing or inconsistent with the baselines,
        or the recorded position is grossly inconsistent with the station.
    """
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - environment issue
        raise ObservationValidationError("Reading native files needs h5py: pip install 'tart-gnss[stage2]'") from exc
    path = Path(path)
    warnings: List[str] = []
    try:
        handle = h5py.File(path, "r")
    except OSError as exc:
        raise ObservationValidationError(f"{path}: not a readable HDF5 file ({exc})") from exc
    with handle as f:
        for required in ("timestamp", "vis"):
            if required not in f:
                raise ObservationValidationError(f"{path}: missing '{required}' dataset")
        try:
            timestamps = [_parse_timestamp(_decode(t)) for t in f["timestamp"][()]]
        except ValueError as exc:
            raise ObservationValidationError(f"{path}: unparsable timestamps ({exc})") from exc
        if not timestamps:
            raise ObservationValidationError(f"{path}: no timestamps")
        secs = np.array([t.timestamp() for t in timestamps])
        if np.any(np.diff(secs) <= 0):
            raise ObservationValidationError(f"{path}: timestamps are not strictly increasing")
        vis = f["vis"][()]
        if vis.ndim != 2 or vis.shape[0] != len(timestamps):
            raise ObservationValidationError(
                f"{path}: vis shape {vis.shape} inconsistent with {len(timestamps)} timestamps"
            )
        n_baselines = int(vis.shape[1])
        config: Dict = {}
        if "config" in f:
            try:
                config = json.loads(_decode(f["config"][()][0] if f["config"].shape else f["config"][()]))
            except (ValueError, IndexError) as exc:
                warnings.append(f"config not parsable: {exc}")
        else:
            warnings.append("no 'config' dataset: frequency and position metadata unavailable")
        n_ant = config.get("num_antenna")
        if "antenna_positions" in f:
            n_pos = int(f["antenna_positions"].shape[0])
            if n_ant is not None and n_pos != n_ant:
                warnings.append(f"num_antenna={n_ant} but {n_pos} antenna positions")
            n_ant = n_ant if n_ant is not None else n_pos
        if "baselines" in f and int(f["baselines"].shape[0]) != n_baselines:
            raise ObservationValidationError(
                f"{path}: {f['baselines'].shape[0]} baselines but vis has {n_baselines} columns"
            )
        if n_ant is not None and n_baselines != n_ant * (n_ant - 1) // 2:
            warnings.append(f"{n_baselines} baselines is not n(n-1)/2 for {n_ant} antennas")
        gains = f["gains"][()] if "gains" in f else None
        phases = f["phases"][()] if "phases" in f else None
        phase_elaz = f["phase_elaz"][()].tolist() if "phase_elaz" in f else None
        flags_present = any(k in f for k in ("flags", "flag"))

    finite = np.isfinite(vis)
    finite_fraction = float(finite.mean()) if vis.size else 0.0
    if finite_fraction < 1.0:
        warnings.append(f"{100 * (1 - finite_fraction):.2f}% non-finite visibilities")
    all_zero = bool(vis.size == 0 or not np.any(vis[finite]))
    if all_zero:
        raise ObservationValidationError(f"{path}: visibilities are all zero or non-finite")

    freq = config.get("frequency")
    bw = config.get("bandwidth")
    if freq is None:
        warnings.append("observing frequency not recorded in file")
        freq_ok = None
    elif bw is None:
        warnings.append("bandwidth not recorded; frequency consistency not assessed")
        freq_ok = None
    else:
        freq_ok = bool(abs(float(freq) - TART_CENTRE_FREQUENCY_HZ) <= float(bw) / 2)
        if not freq_ok:
            warnings.append(
                f"recorded frequency {freq} Hz: 1575.42 MHz is outside the recorded band (bandwidth {bw} Hz)"
            )

    lat, lon, alt = config.get("lat"), config.get("lon"), config.get("alt")
    offset = None
    if lat is not None and lon is not None:
        offset = float(max(abs(lat - station.latitude_deg), abs(lon - station.longitude_deg)))
        if offset > IDENTITY_TOLERANCE_DEG:
            raise ObservationValidationError(
                f"{path}: file position ({lat}, {lon}) differs from station '{station.id}' "
                f"({station.latitude_deg}, {station.longitude_deg}) by {offset:.4f} deg"
            )
        if offset > 0:
            warnings.append(f"file position differs from configured station by {offset:.2e} deg")
    else:
        warnings.append("station position not recorded in file; identity checked by archive name only")

    zero_gain: List[int] = []
    fingerprint = None
    if gains is not None:
        g = np.asarray(gains, dtype=float)
        zero_gain = [int(i) for i in np.flatnonzero(~np.isfinite(g) | (g == 0))]
        digest = hashlib.sha256(np.ascontiguousarray(g).tobytes())
        if phases is not None:
            digest.update(np.ascontiguousarray(np.asarray(phases, dtype=float)).tobytes())
        fingerprint = digest.hexdigest()[:16]
    else:
        warnings.append("no gains in file")

    spacing = float(np.median(np.diff(secs))) if len(secs) > 1 else None
    obs = NativeObservation(
        path=str(path), station_id=station.id, timestamps=timestamps,
        n_timestamps=len(timestamps), first_timestamp_utc=timestamps[0],
        last_timestamp_utc=timestamps[-1], median_spacing_s=spacing,
        telescope_name=config.get("name"), n_antennas=n_ant, n_baselines=n_baselines,
        vis_shape=tuple(int(x) for x in vis.shape), vis_finite_fraction=finite_fraction,
        vis_all_zero=all_zero, frequency_hz=float(freq) if freq is not None else None,
        bandwidth_hz=float(bw) if bw is not None else None, frequency_consistent=freq_ok,
        file_latitude_deg=lat, file_longitude_deg=lon, file_altitude_m=alt,
        position_offset_deg=offset, has_gains=gains is not None, has_phases=phases is not None,
        zero_gain_antennas=zero_gain, calibration_fingerprint=fingerprint,
        phase_centre_el_az_deg=phase_elaz, flags_present=flags_present, warnings=warnings,
    )
    for w in warnings:
        logger.warning("%s: %s", path.name, w)
    return obs


@dataclass
class StationCoverage:
    """The real temporal coverage of one station, from all its native files.

    ``segments`` are ``(start, end)`` UTC pairs of continuous sampling. A step
    larger than ``gap_factor`` times the median cadence starts a new segment.
    """

    station_id: str
    timestamps: List[datetime] = field(repr=False)
    first_timestamp_utc: datetime
    last_timestamp_utc: datetime
    number_of_timestamps: int
    median_time_spacing_seconds: Optional[float]
    gap_threshold_seconds: Optional[float]
    gaps: List[Tuple[datetime, datetime]]
    segments: List[Tuple[datetime, datetime]]
    frequencies_hz: List[float]
    bandwidths_hz: List[float]
    n_antennas: List[int]
    calibration_fingerprints: List[str]
    zero_gain_antennas: List[int]
    duplicate_timestamps: int

    def to_dict(self) -> dict:
        """JSON-serialisable summary."""
        fmt = lambda pairs: [[format_utc(a, 6), format_utc(b, 6)] for a, b in pairs]  # noqa: E731
        return {
            "first_timestamp_utc": format_utc(self.first_timestamp_utc, 6),
            "last_timestamp_utc": format_utc(self.last_timestamp_utc, 6),
            "number_of_timestamps": self.number_of_timestamps,
            "median_time_spacing_seconds": self.median_time_spacing_seconds,
            "gap_threshold_seconds": self.gap_threshold_seconds,
            "n_gaps": len(self.gaps),
            "gaps": fmt(self.gaps),
            "segments": fmt(self.segments),
            "frequencies_hz": self.frequencies_hz,
            "bandwidths_hz": self.bandwidths_hz,
            "n_antennas": self.n_antennas,
            "calibration_fingerprints": self.calibration_fingerprints,
            "zero_gain_antennas": self.zero_gain_antennas,
            "duplicate_timestamps": self.duplicate_timestamps,
        }


def summarise_coverage(
    station_id: str, observations: Sequence[NativeObservation], gap_factor: float = DEFAULT_GAP_FACTOR
) -> StationCoverage:
    """Combine a station's files into one timeline and detect gaps.

    Parameters
    ----------
    observations
        Native observations of one station (any order).
    gap_factor
        A step larger than ``gap_factor x median cadence`` counts as a gap.

    Raises
    ------
    ObservationValidationError
        If ``observations`` is empty.
    """
    if not observations:
        raise ObservationValidationError(f"No native observations for station '{station_id}'")
    all_ts = sorted(t for o in observations for t in o.timestamps)
    unique: List[datetime] = []
    for t in all_ts:
        if not unique or t != unique[-1]:
            unique.append(t)
    duplicates = len(all_ts) - len(unique)
    secs = np.array([t.timestamp() for t in unique])
    steps = np.diff(secs)
    cadence = float(np.median(steps)) if steps.size else None
    threshold = gap_factor * cadence if cadence else None
    gaps: List[Tuple[datetime, datetime]] = []
    segments: List[Tuple[datetime, datetime]] = []
    seg_start = unique[0]
    for i, step in enumerate(steps):
        if threshold is not None and step > threshold:
            gaps.append((unique[i], unique[i + 1]))
            segments.append((seg_start, unique[i]))
            seg_start = unique[i + 1]
    segments.append((seg_start, unique[-1]))
    return StationCoverage(
        station_id=station_id, timestamps=unique, first_timestamp_utc=unique[0],
        last_timestamp_utc=unique[-1], number_of_timestamps=len(unique),
        median_time_spacing_seconds=cadence, gap_threshold_seconds=threshold,
        gaps=gaps, segments=segments,
        frequencies_hz=sorted({o.frequency_hz for o in observations if o.frequency_hz is not None}),
        bandwidths_hz=sorted({o.bandwidth_hz for o in observations if o.bandwidth_hz is not None}),
        n_antennas=sorted({o.n_antennas for o in observations if o.n_antennas is not None}),
        calibration_fingerprints=sorted({o.calibration_fingerprint for o in observations if o.calibration_fingerprint}),
        zero_gain_antennas=sorted({a for o in observations for a in o.zero_gain_antennas}),
        duplicate_timestamps=duplicates,
    )
