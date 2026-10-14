"""Synthetic inputs for Stage-6 tests and demonstrations.

* :func:`teme_states` / :func:`teme_to_gcrs` give the TLE-predicted GCRS
  state by an *independent* route: raw ``sgp4`` TEME output, rotated to the
  Earth-fixed frame with an explicit IAU-1982 GMST (Vallado), then to GCRS
  with the Stage-5 matrices. Skyfield's ``EarthSatellite`` is not used, so
  comparing the two routes tests the Stage-6 frame chain.
* :func:`perturbed_tle` writes valid TLE lines for nearby synthetic element
  sets (``sgp4init`` + ``sgp4.exporter``).
* :func:`write_stage5_run` writes a Stage-5-schema trajectory (with full 6x6
  covariances) and a minimal ``stage5_manifest.json``.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np

from tart_gnss.stage5.frames import itrs_to_gcrs_matrices, skyfield_times

#: The Stage-1 element set of GPS BIIF-9 (PRN 26) used throughout the 2026-06-25 examples.
EXAMPLE_TLE = (
    "1 40534U 15013A   26279.74180529 -.00000040  00000+0  00000+0 0  9991",
    "2 40534  53.1528 199.5676 0114192  41.6722 319.1990  2.00569358 84063",
)
OMEGA_TEME_RAD_S = 7.292115146706979e-5


def gmst82_rad(jd_ut1) -> np.ndarray:
    """Greenwich mean sidereal time, IAU 1982 (Vallado 2013, eq. 3-47), radians."""
    tut1 = (np.asarray(jd_ut1, dtype=float) - 2451545.0) / 36525.0
    sec = (67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 0.093104 * tut1 ** 2
           - 6.2e-6 * tut1 ** 3)
    return np.mod(np.radians(sec / 240.0), 2 * np.pi)


def teme_states(line1: str, line2: str, times_s) -> Tuple[np.ndarray, np.ndarray]:
    """Raw SGP4 TEME position (m) and velocity (m/s) at Unix seconds (UTC)."""
    from sgp4.api import Satrec

    t = np.atleast_1d(np.asarray(times_s, dtype=float))
    sat = Satrec.twoline2rv(line1, line2)
    jd = 2440587.5 + np.floor(t / 86400.0)
    fr = (t - np.floor(t / 86400.0) * 86400.0) / 86400.0
    err, r, v = sat.sgp4_array(jd, fr)
    if np.any(err):
        raise ValueError(f"sgp4 error codes {set(err.tolist())}")
    return r * 1e3, v * 1e3


def teme_to_gcrs(r_teme: np.ndarray, v_teme: np.ndarray, times_s) -> Tuple[np.ndarray, np.ndarray]:
    """TEME -> pseudo-Earth-fixed (GMST82, polar motion zero) -> GCRS (Stage-5 matrices)."""
    t = np.atleast_1d(np.asarray(times_s, dtype=float))
    th = gmst82_rad(skyfield_times(t).ut1)
    c, s = np.cos(th), np.sin(th)
    rot = np.zeros((len(t), 3, 3))
    rot[:, 0, 0], rot[:, 0, 1], rot[:, 1, 0], rot[:, 1, 1], rot[:, 2, 2] = c, s, -s, c, 1.0
    r_e = np.einsum("nij,nj->ni", rot, r_teme)
    w = np.array([0.0, 0.0, OMEGA_TEME_RAD_S])
    v_e = np.einsum("nij,nj->ni", rot, v_teme) - np.cross(w, r_e)
    m = itrs_to_gcrs_matrices(t)
    r_g = np.einsum("nij,nj->ni", m, r_e)
    v_g = np.einsum("nij,nj->ni", m, v_e + np.cross(w, r_e))
    return r_g, v_g


def independent_gcrs_states(line1: str, line2: str, times_s) -> np.ndarray:
    """``(N, 6)`` GCRS states by the independent route (not Skyfield's EarthSatellite)."""
    r, v = teme_states(line1, line2, times_s)
    rg, vg = teme_to_gcrs(r, v, times_s)
    return np.hstack([rg, vg])


def perturbed_tle(line1: str, line2: str, d_epoch_days: float = 0.0, d_mean_anomaly_deg: float = 0.0,
                  d_raan_deg: float = 0.0, d_inclination_deg: float = 0.0, d_mean_motion_frac: float = 0.0,
                  d_ecc: float = 0.0) -> Tuple[str, str]:
    """A nearby element set (valid lines with checksums), via ``sgp4init`` and ``sgp4.exporter``.

    A shifted epoch carries mean anomaly, argument of perigee and RAAN forward
    with SGP4's own secular rates, so an epoch-shifted set describes nearly
    the same trajectory; the ``d_*`` terms then add deliberate differences.
    """
    from sgp4 import exporter
    from sgp4.api import WGS72, Satrec

    s = Satrec.twoline2rv(line1, line2)
    n = Satrec()
    epoch = s.jdsatepoch + s.jdsatepochF - 2433281.5 + d_epoch_days
    dt_min = d_epoch_days * 1440.0  # carry the elements to the new epoch with SGP4's secular rates
    n.sgp4init(WGS72, "i", s.satnum, epoch, s.bstar, s.ndot, s.nddot, s.ecco + d_ecc,
               s.argpo + s.argpdot * dt_min, s.inclo + np.radians(d_inclination_deg),
               np.mod(s.mo + s.mdot * dt_min + np.radians(d_mean_anomaly_deg), 2 * np.pi),
               s.no_kozai * (1 + d_mean_motion_frac), np.mod(s.nodeo + s.nodedot * dt_min + np.radians(d_raan_deg),
                                                            2 * np.pi))
    n.intldesg, n.classification, n.elnum, n.revnum = s.intldesg, s.classification, s.elnum, s.revnum
    return exporter.export_tle(n)


def write_tle_file(path: Path, sets: Iterable[Tuple[str, str]], name: str = "GPS BIIF-9  (PRN 26)") -> Path:
    """Write three-line TLE sets to ``path``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{name}\n{a}\n{b}\n" for a, b in sets), encoding="utf-8")
    return path


def isotropic_state_covariance(sigma_pos_m: float, sigma_vel_m_s: float) -> np.ndarray:
    """Diagonal 6x6 covariance."""
    return np.diag([sigma_pos_m ** 2] * 3 + [sigma_vel_m_s ** 2] * 3)


def write_stage5_run(out_dir: Path, times_s: Sequence[float], states_gcrs: np.ndarray,
                     covariances: np.ndarray, norad_id: int = 40534, flags: Optional[dict] = None,
                     stage4_manifest: Optional[str] = None, valid: Optional[Sequence[bool]] = None) -> Path:
    """Stage-5-schema ``trajectory.{parquet,csv}`` and a minimal ``stage5_manifest.json``; return the manifest."""
    from tart_gnss.stage4.output import write_manifest
    from tart_gnss.stage5.frames import gcrs_to_ecef_positions
    from tart_gnss.stage5.output import covariance_fields, write_tables

    t = np.asarray(times_s, dtype=float)
    mats = itrs_to_gcrs_matrices(t)
    ecef = gcrs_to_ecef_positions(np.asarray(states_gcrs)[:, :3], mats)
    rows = []
    for k, (tk, x, p) in enumerate(zip(t, states_gcrs, covariances)):
        f = list((flags or {}).get(k, []))
        ok = True if valid is None else bool(valid[k])
        rows.append({
            "time_utc": datetime.fromtimestamp(tk, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
            "time_mjd": tk / 86400.0 + 40587.0, "norad_id": norad_id, "dynamics": "synthetic",
            "x_gcrs_smoothed_m": x[0], "y_gcrs_smoothed_m": x[1], "z_gcrs_smoothed_m": x[2],
            "vx_gcrs_smoothed_m_s": x[3], "vy_gcrs_smoothed_m_s": x[4], "vz_gcrs_smoothed_m_s": x[5],
            "x_gcrs_filtered_m": x[0], "y_gcrs_filtered_m": x[1], "z_gcrs_filtered_m": x[2],
            "vx_gcrs_filtered_m_s": x[3], "vy_gcrs_filtered_m_s": x[4], "vz_gcrs_filtered_m_s": x[5],
            "x_ecef_smoothed_m": ecef[k, 0], "y_ecef_smoothed_m": ecef[k, 1], "z_ecef_smoothed_m": ecef[k, 2],
            "valid": ok, "quality_flags": f or ["GOOD"], "stage4_quality_flags": [],
            **covariance_fields(p, "smoothed"), **covariance_fields(p, "filtered"),
        })
    out_dir = Path(out_dir)
    paths = write_tables(rows, [], out_dir)
    doc = {"schema_version": 1, "stage": 5, "norad_id": norad_id, "state_frame": "GCRS",
           "measurement_frame": "ITRS/ECEF", "input_stage4_manifest": stage4_manifest, "tle_used": False,
           "software": {"git_commit": None}, "provenance": {"command_line": ["synthetic"]}, "outputs": paths}
    return Path(write_manifest(doc, out_dir / "stage5_manifest.json"))
