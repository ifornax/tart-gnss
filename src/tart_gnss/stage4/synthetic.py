"""Synthetic Stage-3 outputs for Stage-4 tests and demonstrations.

Directions are computed from a known 3-D (ECEF) truth seen from each station,
perturbed by an isotropic angular error of ``sigma_deg`` drawn in the
``(l, m)`` plane with the matching covariance ``sigma^2 (J^T J)^-1`` (``J`` =
``d s / d(l, m)``), and written as a Stage-3 schema table plus a minimal
``stage3_manifest.json``. Nothing here is used by the Stage-4 pipeline itself.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage3.geometry import azel_jacobian, lm_to_azel
from tart_gnss.stage4.geometry import ecef_to_enu_vector, geodetic_to_ecef, lm_jacobian_enu

#: Real coordinates of the three stations of the 2026-06-25 run (Stage-2 manifest), for tests and demos.
EXAMPLE_STATIONS: Dict[str, Tuple[float, float, float]] = {
    "mauritius": (-20.2587508, 57.7591989, 20.0),
    "madagascar": (-18.8971446, 47.5551677, 1280.0),
    "south_africa": (-33.319413, 26.507688, 642.26),
}
#: A GNSS-like ECEF position above the three example stations (m) and a velocity (m/s).
EXAMPLE_POSITION_M = np.array([1.2e7, 9.0e6, -2.0e7])
EXAMPLE_VELOCITY_M_S = np.array([1500.0, -2500.0, 600.0])
T0_UNIX_S = 1782378840.0  # 2026-06-25T09:14:00Z


def linear_track(x0=EXAMPLE_POSITION_M, v=EXAMPLE_VELOCITY_M_S, t0: float = T0_UNIX_S) -> Callable:
    """``truth(t_unix_s) -> ECEF position`` for uniform motion."""
    x0, v = np.asarray(x0, dtype=float), np.asarray(v, dtype=float)
    return lambda t: x0 + v * (float(t) - t0)


def exact_lm(target_ecef_m, lat: float, lon: float, alt: float) -> Tuple[float, float]:
    """Direction cosines ``(l, m)`` of an ECEF point from a station (geometric, no light time)."""
    d = ecef_to_enu_vector(np.asarray(target_ecef_m, dtype=float) - geodetic_to_ecef(lat, lon, alt), lat, lon)
    d = d / np.linalg.norm(d)
    return float(d[0]), float(d[1])


def isotropic_lm_covariance(l: float, m: float, sigma_rad: float) -> np.ndarray:
    """``(l, m)`` covariance of an isotropic angular error ``sigma_rad``."""
    j = lm_jacobian_enu(l, m)
    return sigma_rad ** 2 * np.linalg.inv(j.T @ j)


def direction_row(station_id: str, t_unix_s: float, l: float, m: float, cov_lm: np.ndarray,
                  norad_id: int = 40534, flags: Sequence[str] = ()) -> Dict[str, object]:
    """One Stage-3 schema row for a measured ``(l, m)`` with covariance."""
    from datetime import datetime, timezone

    az, el = lm_to_azel(l, m)
    g = azel_jacobian(l, m)
    cae = np.degrees(1.0) ** 2 * (g @ cov_lm @ g.T)
    blocking = {"LOW_BASELINE_COUNT", "FIT_FAILED", "BOUNDARY_SOLUTION", "INVALID_GEOMETRY", "LOW_SNR",
                "TLE_UNASSOCIATED"}
    return {
        "time_utc": datetime.fromtimestamp(t_unix_s, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
        "time_mjd": t_unix_s / 86400.0 + 40587.0, "station_id": station_id, "tart_id": f"synthetic-{station_id}",
        "norad_id": norad_id, "primary_method": "visibility",
        "l": l, "m": m, "sigma_l": float(np.sqrt(cov_lm[0, 0])), "sigma_m": float(np.sqrt(cov_lm[1, 1])),
        "cov_lm": float(cov_lm[0, 1]),
        "az_obs_deg": float(az), "el_obs_deg": float(el), "sigma_az_deg": float(np.sqrt(cae[0, 0])),
        "sigma_el_deg": float(np.sqrt(cae[1, 1])), "cov_az_el": float(cae[0, 1]),
        "snr": 20.0, "fit_success": "FIT_FAILED" not in flags, "usable": not (set(flags) & blocking),
        "quality_flags": list(flags),
    }


def make_directions(
    truth: Callable,
    times_s: Dict[str, Iterable[float]],
    sigma_deg: float,
    rng: Optional[np.random.Generator],
    stations: Dict[str, Tuple[float, float, float]] = EXAMPLE_STATIONS,
    sigma_overrides: Optional[Dict[str, float]] = None,
    bias_deg: Optional[Dict[str, Tuple[float, float]]] = None,
) -> list:
    """Stage-3 rows for each station and time.

    ``rng=None`` gives noise-free directions (covariance still ``sigma_deg``).
    ``bias_deg[station] = (dl_deg, dm_deg)`` adds a fixed offset (a "bad station").
    """
    rows = []
    for sid, ts in times_s.items():
        lat, lon, alt = stations[sid]
        sig = np.radians((sigma_overrides or {}).get(sid, sigma_deg))
        for t in ts:
            l, m = exact_lm(truth(t), lat, lon, alt)
            cov = isotropic_lm_covariance(l, m, sig)
            if rng is not None:
                l, m = np.array([l, m]) + rng.multivariate_normal(np.zeros(2), cov)
            if bias_deg and sid in bias_deg:
                l, m = l + np.radians(bias_deg[sid][0]), m + np.radians(bias_deg[sid][1])
            rows.append(direction_row(sid, float(t), float(l), float(m), cov))
    return rows


def write_stage3_run(out_dir: Path, rows: list, stations: Dict[str, Tuple[float, float, float]] = EXAMPLE_STATIONS,
                     norad_id: int = 40534, calibration_warnings: Optional[Dict[str, list]] = None,
                     csv_only: bool = False) -> Path:
    """Write ``directions.{parquet,csv}`` and a minimal ``stage3_manifest.json``; return the manifest path."""
    from tart_gnss.stage3.output import write_manifest, write_tables

    out_dir = Path(out_dir)
    paths = write_tables(rows, out_dir)
    if csv_only and paths.get("parquet"):
        Path(paths["parquet"]).unlink()
        paths["parquet"] = None
    doc = {
        "schema_version": 1, "stage": 3, "input_manifest": None, "norad_id": norad_id,
        "satellite": {"satellite_name": "SYNTHETIC"},
        "stations": {
            sid: {"station_id": sid, "latitude_deg": lat, "longitude_deg": lon, "altitude_m": alt,
                  "calibration": {"warnings": list((calibration_warnings or {}).get(sid, []))}}
            for sid, (lat, lon, alt) in stations.items()
        },
        "software": {"git_commit": None}, "provenance": {"command_line": ["synthetic"]},
        "outputs": {"directions_csv": paths["csv"], "directions_parquet": paths["parquet"]},
    }
    return Path(write_manifest(doc, out_dir / "stage3_manifest.json"))


def write_json(path: Path, doc) -> None:
    """Small helper used by tests to rewrite a manifest."""
    Path(path).write_text(json.dumps(doc, indent=2), encoding="utf-8")
