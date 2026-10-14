"""Synthetic diagnosis scenarios built through the real pipeline (Stages 3-6).

A scenario places a known truth (the PRN-26 TLE trajectory by the independent
TEME route of :mod:`tart_gnss.stage6.synthetic`, optionally time-shifted or
accelerated), observes it from the three example stations as Stage-3
directions with optional injected effects, and runs ``tart-triangulate``,
``tart-track`` and ``tart-validate`` on it, so every Stage-7 input carries
genuine provenance and hashes. A matching synthetic Stage-2 manifest supplies
per-file calibration metadata (zero-gain antennas, gain-solution
fingerprints). Nothing here is used by the Stage-7 pipeline itself.

Injectable effects (all optional):

``tau_s``              network-wide time shift of the truth (along-track offset)
``accel_m_s2``         unmodelled constant acceleration of the truth (GCRS vector)
``bias_deg``           ``{station: (dl, dm)}`` fixed direction bias
``calibration_bias``   ``(station, (dl, dm))`` bias applied only on epochs whose
                       ``flagged_fraction`` is high (``calibration_fraction`` of epochs)
``time_offsets_s``     ``{station: array}`` per-epoch sampling offsets (recorded)
``subsets``            ``{station: bool array}`` epochs at which the station observes
``random_metadata``    random-walk ``flagged_fraction``/``snr`` unrelated to the data
``tle_epoch_shift_d``  validate against an element set re-epoched by this many days
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage4.synthetic import EXAMPLE_STATIONS, direction_row, exact_lm, isotropic_lm_covariance
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices
from tart_gnss.stage6 import synthetic as syn6
from tart_gnss.tle import tle_epoch_utc

L1, L2 = syn6.EXAMPLE_TLE
#: Observations start 0.35 d after the TLE epoch (a fresh reference), when PRN 26 is 40-47 deg
#: above all three example stations.
T0 = tle_epoch_utc(L1, L2).timestamp() + 0.35 * 86400.0


@dataclass
class Scenario:
    """Parameters of one synthetic diagnosis scenario."""

    name: str
    n_epochs: int = 36
    cadence_s: float = 10.0
    sigma_deg: float = 0.005
    seed: int = 1
    tau_s: float = 0.0
    accel_m_s2: Optional[Sequence[float]] = None
    bias_deg: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    calibration_bias: Optional[Tuple[str, Tuple[float, float]]] = None
    calibration_fraction: float = 0.4
    time_offsets_s: Dict[str, np.ndarray] = field(default_factory=dict)
    subsets: Dict[str, np.ndarray] = field(default_factory=dict)
    random_metadata: bool = False
    tle_epoch_shift_d: float = 0.0
    with_stage2: bool = True
    min_stations: int = 3
    time_match_tolerance_s: float = 0.75
    stations: Dict[str, Tuple[float, float, float]] = field(default_factory=lambda: dict(EXAMPLE_STATIONS))


def epochs(sc: Scenario) -> np.ndarray:
    """Reference epochs (Unix s)."""
    return T0 + sc.cadence_s * np.arange(sc.n_epochs)


def truth_gcrs(sc: Scenario, t) -> np.ndarray:
    """True GCRS position at ``t``."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    x = syn6.independent_gcrs_states(L1, L2, t + sc.tau_s)[:, :3]
    if sc.accel_m_s2 is not None:
        x = x + 0.5 * np.outer((t - T0) ** 2, np.asarray(sc.accel_m_s2, float))
    return x


def truth_ecef(sc: Scenario, t: float) -> np.ndarray:
    """True ECEF position at one epoch."""
    m = itrs_to_gcrs_matrices([t])[0]
    return m.T @ truth_gcrs(sc, [t])[0]


def _random_walk(rng, n, lo, hi):
    w = np.cumsum(rng.normal(size=n))
    w = (w - w.min()) / (np.ptp(w) or 1.0)
    return lo + (hi - lo) * w


def stage3_rows(sc: Scenario, rng) -> list:
    """Stage-3 schema rows for every station and epoch, with injected effects and metadata."""
    t_ref = epochs(sc)
    n = len(t_ref)
    flagged_mask = rng.random(n) < sc.calibration_fraction
    rows = []
    for sid, (lat, lon, alt) in sc.stations.items():
        offs = np.asarray(sc.time_offsets_s.get(sid, np.zeros(n)), float)
        present = np.asarray(sc.subsets.get(sid, np.ones(n, bool)), bool)
        rw_flag = _random_walk(rng, n, 0.05, 0.6) if sc.random_metadata else None
        rw_snr = _random_walk(rng, n, 8.0, 30.0) if sc.random_metadata else None
        for k in range(n):
            if not present[k]:
                continue
            t = t_ref[k] + offs[k]
            l, m = exact_lm(truth_ecef(sc, t), lat, lon, alt)
            cov = isotropic_lm_covariance(l, m, np.radians(sc.sigma_deg))
            l, m = np.array([l, m]) + rng.multivariate_normal(np.zeros(2), cov)
            flagged = 0.1
            if sid in sc.bias_deg:
                l, m = l + np.radians(sc.bias_deg[sid][0]), m + np.radians(sc.bias_deg[sid][1])
            if sc.calibration_bias is not None and sc.calibration_bias[0] == sid:
                if flagged_mask[k]:
                    b = sc.calibration_bias[1]
                    l, m = l + np.radians(b[0]), m + np.radians(b[1])
                    flagged = 0.6
            if rw_flag is not None:
                flagged = float(rw_flag[k])
            row = direction_row(sid, float(t), float(l), float(m), cov)
            row.update({
                "flagged_fraction": flagged, "reduced_chi2": 1.0 + 0.1 * rng.random(),
                "vis_image_separation_deg": abs(rng.normal(0.0, 0.05)), "n_baselines_used": 250,
                "n_baselines_total": 276, "snr": float(rw_snr[k]) if rw_snr is not None else 20.0,
            })
            rows.append(row)
    return rows


def stage2_manifest(sc: Scenario, path: Path) -> Path:
    """Minimal Stage-2 manifest with per-file calibration metadata (one gain solution per station)."""
    t = epochs(sc)
    stations = {}
    for sid, (lat, lon, alt) in sc.stations.items():
        files = []
        for j, (a, b) in enumerate(((t[0] - 30, t[len(t) // 2]), (t[len(t) // 2] + 1e-3, t[-1] + 30))):
            files.append({
                "path": f"synthetic/{sid}/vis_{j}.hdf", "station_id": sid,
                "first_timestamp_utc": _iso(a), "last_timestamp_utc": _iso(b),
                "zero_gain_antennas": [], "calibration_fingerprint": f"{sid}-gains-0",
            })
        stations[sid] = {"station_id": sid, "latitude_deg": lat, "longitude_deg": lon, "altitude_m": alt,
                         "observation": {"files": files, "calibration_status": "synthetic"}}
    doc = {"stage": 2, "manifest_version": 1, "three_station_valid": True, "stations": stations,
           "synchronisation": {"station_pair_time_offset_statistics": {}}}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    return path


def _iso(t: float) -> str:
    return datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def build_chain(sc: Scenario, work: Path) -> Dict[str, object]:
    """Run Stages 3 (synthetic) -> 4 -> 5 -> 6; return paths and exit codes."""
    from tart_gnss.stage4 import cli as cli4
    from tart_gnss.stage4.synthetic import write_stage3_run
    from tart_gnss.stage5 import cli as cli5
    from tart_gnss.stage6 import cli as cli6

    work = Path(work)
    rng = np.random.default_rng(sc.seed)
    rows = stage3_rows(sc, rng)
    m3 = write_stage3_run(work / "stage3", rows, stations=sc.stations)
    if sc.with_stage2:
        m2 = stage2_manifest(sc, work / "stage2" / "manifest.json")
        doc = json.loads(m3.read_text())
        doc["input_manifest"] = str(m2)
        doc["input_manifest_sha256"] = hashlib.sha256(m2.read_bytes()).hexdigest()
        m3.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    c4 = cli4.main(["--stage3-manifest", str(m3), "--output-dir", str(work / "stage4"), "--no-plots",
                    "--min-stations", str(sc.min_stations), "--time-match-tolerance", str(sc.time_match_tolerance_s)])
    m4 = next((work / "stage4").glob("*/NORAD_*/stage4_manifest.json"))
    c5 = cli5.main(["--stage4-manifest", str(m4), "--output-dir", str(work / "stage5"), "--no-plots"])
    m5 = next((work / "stage5").glob("*/NORAD_*/stage5_manifest.json"))
    tle_sets = [(L1, L2)] if sc.tle_epoch_shift_d == 0 else [syn6.perturbed_tle(L1, L2,
                                                                               d_epoch_days=sc.tle_epoch_shift_d)]
    tle = syn6.write_tle_file(work / "reference.tle", tle_sets)
    c6 = cli6.main(["--stage5-manifest", str(m5), "--output-dir", str(work / "stage6"), "--no-plots",
                    "--tle-file", str(tle)])
    m6 = next((work / "stage6").glob("*/NORAD_*/stage6_manifest.json"))
    return {"stage3": m3, "stage4": m4, "stage5": m5, "stage6": m6, "codes": (c4, c5, c6)}
