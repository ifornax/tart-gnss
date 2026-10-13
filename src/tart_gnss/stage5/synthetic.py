"""Synthetic Stage-4 outputs for Stage-5 tests and demonstrations (no TLE anywhere).

Truth orbits are propagated *analytically* (Lagrange f and g functions with a
universal form of Kepler's equation), not with the integrator Stage 5 uses, so
the synthetic checks are not circular. Stage-4-like measurements are the true
ECEF position plus Gaussian noise drawn from an anisotropic covariance: large
along the line of sight from a reference station, small across it (the shape
of real Stage-4 covariances). Nothing here is used by the Stage-5 pipeline.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import numpy as np

from tart_gnss.stage5.dynamics import GM_EARTH_M3_S2
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices

#: 2026-06-25T09:14:00Z, the start of the real 2026-06-25 example.
T0_UNIX_S = 1782378840.0
#: Direction (ECEF) of a GNSS satellite above the three example stations (see stage4.synthetic).
EXAMPLE_DIRECTION_ECEF = np.array([1.2e7, 9.0e6, -2.0e7])
GNSS_RADIUS_M = 26_560e3
#: South Africa example station, reference for the line-of-sight covariance axis.
REFERENCE_STATION = (-33.319413, 26.507688, 642.26)


def kepler_propagate(state: np.ndarray, dt_s: float, mu: float = GM_EARTH_M3_S2) -> np.ndarray:
    """Exact two-body propagation of ``[r, v]`` by ``dt`` (elliptic orbits, any eccentricity < 1)."""
    r0, v0 = np.asarray(state[:3], float), np.asarray(state[3:], float)
    rn = np.linalg.norm(r0)
    a = 1.0 / (2.0 / rn - v0 @ v0 / mu)
    if a <= 0:
        raise ValueError("not an elliptic orbit")
    n = np.sqrt(mu / a ** 3)
    sig0 = (r0 @ v0) / np.sqrt(mu * a)
    c0 = 1.0 - rn / a
    m = n * dt_s
    de = m
    for _ in range(50):  # Newton on  M = dE + sig0 (1 - cos dE) - c0 sin dE
        f = de + sig0 * (1 - np.cos(de)) - c0 * np.sin(de) - m
        fp = 1 + sig0 * np.sin(de) - c0 * np.cos(de)
        step = f / fp
        de -= step
        if abs(step) < 1e-14:
            break
    f_ = 1 - a / rn * (1 - np.cos(de))
    g_ = dt_s - (de - np.sin(de)) / n
    r = f_ * r0 + g_ * v0
    r_n = np.linalg.norm(r)
    fd = -np.sqrt(mu * a) / (r_n * rn) * np.sin(de)
    gd = 1 - a / r_n * (1 - np.cos(de))
    return np.concatenate([r, fd * r0 + gd * v0])


def example_initial_state(t0_unix_s: float = T0_UNIX_S, radius_m: float = GNSS_RADIUS_M,
                          inclination_deg: float = 55.0) -> np.ndarray:
    """Circular GNSS-like GCRS state at ``t0`` above the example stations."""
    m = itrs_to_gcrs_matrices([t0_unix_s])[0]
    u = m @ (EXAMPLE_DIRECTION_ECEF / np.linalg.norm(EXAMPLE_DIRECTION_ECEF))
    z = np.array([0.0, 0.0, 1.0])
    zp = z - (z @ u) * u
    w = zp / np.linalg.norm(zp)
    p = np.cross(u, w)
    cb = np.clip(np.cos(np.radians(inclination_deg)) / np.linalg.norm(zp), -1, 1)
    h = cb * w + np.sqrt(1 - cb * cb) * p
    v_hat = np.cross(h, u)
    return np.concatenate([radius_m * u, np.sqrt(GM_EARTH_M3_S2 / radius_m) * v_hat])


def truth_states(state0: np.ndarray, t0_unix_s: float, times_unix_s: Iterable[float],
                 extra_accel_m_s2: Optional[np.ndarray] = None) -> np.ndarray:
    """True GCRS states at ``times``; ``extra_accel`` adds a constant unmodelled acceleration.

    Without ``extra_accel`` the propagation is the analytic Kepler solution;
    with it, the state is integrated numerically (``solve_ivp`` DOP853, tight
    tolerances) from ``t0``.
    """
    times = np.asarray(list(times_unix_s), dtype=float)
    if extra_accel_m_s2 is None:
        return np.array([kepler_propagate(state0, t - t0_unix_s) for t in times])
    from scipy.integrate import solve_ivp

    acc = np.asarray(extra_accel_m_s2, float)

    def rhs(_t, y):
        r = y[:3]
        return np.concatenate([y[3:], -GM_EARTH_M3_S2 * r / np.linalg.norm(r) ** 3 + acc])

    sol = solve_ivp(rhs, (0.0, float(times.max() - t0_unix_s)), np.asarray(state0, float), method="DOP853",
                    rtol=1e-12, atol=1e-6, t_eval=times - t0_unix_s, dense_output=False)
    return sol.y.T


def los_covariance(target_m: np.ndarray, station_m: np.ndarray, sigma_los_m: float,
                   sigma_transverse_m: float) -> np.ndarray:
    """Anisotropic 3x3 covariance: ``sigma_los`` along the station->target line, ``sigma_transverse`` across."""
    u = np.asarray(target_m, float) - np.asarray(station_m, float)
    u /= np.linalg.norm(u)
    return sigma_transverse_m ** 2 * np.eye(3) + (sigma_los_m ** 2 - sigma_transverse_m ** 2) * np.outer(u, u)


def stage4_rows(
    times_unix_s: Sequence[float],
    truth_gcrs: np.ndarray,
    sigma_los_m: float,
    sigma_transverse_m: float,
    rng: Optional[np.random.Generator],
    norad_id: int = 40534,
    outliers: Optional[Dict[int, np.ndarray]] = None,
    flags: Optional[Dict[int, Sequence[str]]] = None,
    station=REFERENCE_STATION,
) -> list:
    """Stage-4 schema rows: truth (rotated to ECEF) plus noise from the anisotropic covariance.

    ``outliers[k]`` adds an ECEF offset (m) to row ``k``; ``flags[k]`` sets its
    Stage-4 quality flags. ``rng=None`` gives noise-free positions.
    """
    from tart_gnss.stage4.geometry import ecef_to_geodetic, geodetic_to_ecef

    times = np.asarray(times_unix_s, float)
    mats = itrs_to_gcrs_matrices(times)
    st = geodetic_to_ecef(*station)
    rows = []
    for k, (t, xg, m) in enumerate(zip(times, truth_gcrs, mats)):
        x = m.T @ xg[:3]
        cov = los_covariance(x, st, sigma_los_m, sigma_transverse_m)
        y = x + (rng.multivariate_normal(np.zeros(3), cov) if rng is not None else 0.0)
        if outliers and k in outliers:
            y = y + np.asarray(outliers[k], float)
        lat, lon, alt = ecef_to_geodetic(y)
        f = list((flags or {}).get(k, []))
        rows.append({
            "time_utc": datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
            "time_mjd": t / 86400.0 + 40587.0, "norad_id": norad_id, "reference_station": "synthetic",
            "method": "weighted", "x_ecef_m": y[0], "y_ecef_m": y[1], "z_ecef_m": y[2],
            "cov_xx": cov[0, 0], "cov_xy": cov[0, 1], "cov_xz": cov[0, 2], "cov_yy": cov[1, 1],
            "cov_yz": cov[1, 2], "cov_zz": cov[2, 2],
            "sigma_x_m": np.sqrt(cov[0, 0]), "sigma_y_m": np.sqrt(cov[1, 1]), "sigma_z_m": np.sqrt(cov[2, 2]),
            "n_stations": 3, "stations": "a;b;c", "condition_number": 100.0, "matrix_rank": 3,
            "min_crossing_angle_deg": 5.0, "max_crossing_angle_deg": 10.0, "chi2": 3.0, "dof": 3,
            "reduced_chi2": 1.0, "fit_success": True, "converged": True, "valid": not (set(f) & {
                "TRIANGULATION_FAILED", "INSUFFICIENT_STATIONS", "NEGATIVE_RANGE", "COVARIANCE_SINGULAR",
                "POOR_GEOMETRY"}),
            "quality_flags": f or ["GOOD"], "latitude_deg": lat, "longitude_deg": lon, "altitude_m": alt,
        })
    return rows


def write_stage4_run(out_dir: Path, rows: list, norad_id: int = 40534, csv_only: bool = False,
                     with_tle_block: bool = False) -> Path:
    """Write ``positions.{csv,parquet}`` (+ empty contributions) and a minimal ``stage4_manifest.json``."""
    from tart_gnss.stage4.output import write_manifest, write_tables

    out_dir = Path(out_dir)
    paths = write_tables(rows, [], out_dir)
    if csv_only:
        for k in ("positions_parquet", "station_contributions_parquet"):
            if paths.get(k):
                Path(paths[k]).unlink()
                paths[k] = None
    doc = {
        "schema_version": 1, "stage": 4, "input_stage3_manifest": None, "norad_id": norad_id,
        "coordinate_frame": "ITRS/ECEF", "software": {"git_commit": None},
        "provenance": {"command_line": ["synthetic"]}, "outputs": paths,
    }
    if with_tle_block:
        doc["satellite_reference"] = {"note": "should never be read by Stage 5", "satellite": {"tle": "X"}}
    return Path(write_manifest(doc, out_dir / "stage4_manifest.json"))


def write_json(path: Path, doc) -> None:
    """Small helper used by tests to rewrite a manifest."""
    Path(path).write_text(json.dumps(doc, indent=2), encoding="utf-8")
