#!/usr/bin/env python
"""Synthetic demonstration of Stage 6 (``tart-validate``).

Reference: the Stage-1 element set of GPS BIIF-9 (PRN 26, NORAD 40534),
epoch 2026-10-06T17:48Z; observations are placed 0.1 d after that epoch so
the reference is fresh.

Part A, known injected residuals
    A Stage-5-schema trajectory built from the TLE by an independent route
    (raw SGP4 TEME -> GMST82 -> GCRS, not Skyfield's EarthSatellite) plus a
    known RTN offset (+1000 m radial, +2000 m along-track, -500 m cross-track)
    and a known RTN velocity offset, with an anisotropic Stage-5-like
    covariance. Validated against the TLE alone and with a 5-member ensemble
    of nearby element sets (empirical reference dispersion).

Part B, end to end through Stage 5
    Stage-4-like noisy positions of the TLE trajectory shifted by tau = 2 s
    (a pure along-track displacement), tracked with ``tart-track`` and then
    validated: the along-track residual should be about |v| tau, the others
    about zero. A Monte Carlo with tau = 0 checks that the
    measurement-conditional position Mahalanobis distance follows chi^2_3 when
    the TLE equals the truth.

Usage::

    python scripts/stage6_synthetic_demo.py [--output-dir docs/examples/stage6] [--mc 100]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from tart_gnss.stage6 import cli as cli6
from tart_gnss.stage6 import synthetic as syn
from tart_gnss.stage6.output import read_table
from tart_gnss.stage6.rtn import rtn_bases
from tart_gnss.tle import tle_epoch_utc

L1, L2 = syn.EXAMPLE_TLE
T0 = tle_epoch_utc(L1, L2).timestamp() + 0.1 * 86400.0
DR_RTN = (1000.0, 2000.0, -500.0)
DV_RTN = (0.10, -0.05, 0.02)
SIGMA_RTN_M = (300.0, 800.0, 400.0)
SIGMA_V_RTN_M_S = (0.3, 0.2, 0.25)
TAU_S = 2.0


def inject(states, dr, dv):
    q = rtn_bases(states[:, :3], states[:, 3:])
    out = states.copy()
    out[:, :3] += np.einsum("nji,j->ni", q, np.asarray(dr))
    out[:, 3:] += np.einsum("nji,j->ni", q, np.asarray(dv))
    return out, q


def validate(manifest: Path, out: Path, tle: Path, *extra) -> tuple:
    code = cli6.main(["--stage5-manifest", str(manifest), "--output-dir", str(out), "--tle-file", str(tle), *extra])
    d = next(out.glob("*/NORAD_*"))
    return code, read_table(d / "residuals.parquet"), json.loads((d / "summary.json").read_text())


def fmt(v, unit, p=1):
    return "(" + ", ".join(f"{x:+.{p}f}" for x in v) + f") {unit}"


def part_a(out: Path, lines: list) -> None:
    t = T0 + np.arange(0.0, 1800.0, 30.0)
    truth = syn.independent_gcrs_states(L1, L2, t)
    obs, q = inject(truth, DR_RTN, DV_RTN)
    # Anisotropic Stage-5-like covariance, defined in RTN and rotated to GCRS.
    cov = np.zeros((len(t), 6, 6))
    for k in range(len(t)):
        cov[k, :3, :3] = q[k].T @ np.diag(np.square(SIGMA_RTN_M)) @ q[k]
        cov[k, 3:, 3:] = q[k].T @ np.diag(np.square(SIGMA_V_RTN_M_S)) @ q[k]
    m5 = syn.write_stage5_run(out / "part_a" / "stage5_input", t, obs, cov)
    tle = syn.write_tle_file(out / "reference.tle", [(L1, L2)])
    sets = [(L1, L2)] + [syn.perturbed_tle(L1, L2, d_epoch_days=e, d_mean_anomaly_deg=dm)
                         for e, dm in ((-3.0, 0.002), (-1.0, -0.001), (1.5, 0.0015), (3.0, -0.002))]
    ens = syn.write_tle_file(out / "part_a" / "ensemble" / "nearby_sets.tle", sets)
    c1, r1, s1 = validate(m5, out / "part_a" / "single_tle", tle)
    c2, r2, s2 = validate(m5, out / "part_a" / "with_ensemble", tle, "--tle-ensemble-dir", str(ens.parent), "--no-plots")
    med = lambda df, c: float(np.median(df[c]))  # noqa: E731
    lines += [
        "PART A: known injected residuals (observed built by an independent TEME -> GCRS route)",
        f"  TLE epoch {s1['tle_epoch_utc']}; median TLE age {s1['median_tle_age_days']:.3f} d (forward); "
        f"precision_validation {s1['precision_validation']}; exit codes {c1} (single TLE), {c2} (ensemble)",
        f"  injected RTN position offset           {fmt(DR_RTN, 'm')}",
        f"  recovered radial / along / cross       {fmt([med(r1, c) for c in ('dr_radial_m', 'dr_along_track_m', 'dr_cross_track_m')], 'm', 3)}",
        f"  injected RTN velocity offset           {fmt(DV_RTN, 'm/s', 3)}",
        f"  recovered RTN velocity residual        {fmt([med(r1, c) for c in ('dv_radial_m_s', 'dv_along_track_m_s', 'dv_cross_track_m_s')], 'm/s', 4)}",
        f"  3-D residual norm                      {med(r1, 'position_residual_3d_m'):.1f} m "
        f"(expected {np.linalg.norm(DR_RTN):.1f} m)",
        f"  geocentric angular residual            {med(r1, 'geocentric_angular_residual_deg') * 3600:.2f} arcsec",
        f"  equivalent along-track time offset     {med(r1, 'equivalent_along_track_time_offset_s'):.4f} s",
        f"  Stage-5 sigma R/T/N                    {fmt([med(r1, c) for c in ('sigma_radial_obs_m', 'sigma_along_track_obs_m', 'sigma_cross_track_obs_m')], 'm')}",
        f"  conditional position D^2 (TLE fixed)   {med(r1, 'position_mahalanobis_obs_cov'):.2f}  "
        f"(= (1000/300)^2 + (2000/800)^2 + (500/400)^2 = {(1000 / 300) ** 2 + (2000 / 800) ** 2 + (500 / 400) ** 2:.2f})",
        f"  conditional state D^2                  {med(r1, 'state_mahalanobis_obs_cov'):.2f}",
        f"  single TLE: total-covariance D^2 NaN: {bool(r1['position_mahalanobis_total_cov'].isna().all())}; flags "
        f"{sorted(set(f for fl in r1['quality_flags'] for f in fl))}",
        f"  ensemble: {s2['ensemble']['n_tle_ensemble']} sets (primary included: {s2['ensemble']['primary_in_ensemble']}); "
        f"empirical reference dispersion R/T/N "
        f"{fmt([med(r2, c) for c in ('sigma_radial_ref_m', 'sigma_along_track_ref_m', 'sigma_cross_track_ref_m')], 'm')}",
        f"  reference-aware position D^2           {med(r2, 'position_mahalanobis_total_cov'):.2f}; state D^2 "
        f"{med(r2, 'state_mahalanobis_total_cov'):.2f}",
        "",
    ]


def part_b(out: Path, n_mc: int, lines: list) -> None:
    import pandas as pd

    from tart_gnss.stage5 import cli as cli5
    from tart_gnss.stage5 import synthetic as syn5
    from tart_gnss.stage5.cli import track
    from tart_gnss.stage5.filter import FilterConfig
    from tart_gnss.stage5.loader import build_measurements
    from tart_gnss.stage6.propagation import propagate_tle
    from tart_gnss.stage6.residuals import StateTrajectory, compute_residuals
    from tart_gnss.stage6.tle_reference import TLEReference

    sig_los, sig_t = 3e3, 300.0
    t = T0 + np.arange(0.0, 1800.0, 10.0)
    rng = np.random.default_rng(6)
    shifted = syn.independent_gcrs_states(L1, L2, t + TAU_S)
    rows = syn5.stage4_rows(t, shifted, sig_los, sig_t, rng)
    m4 = syn5.write_stage4_run(out / "part_b" / "stage4_input", rows)
    c5 = cli5.main(["--stage4-manifest", str(m4), "--output-dir", str(out / "part_b" / "stage5")])
    m5 = next((out / "part_b" / "stage5").glob("*/NORAD_*/stage5_manifest.json"))
    c6, r, s = validate(m5, out / "part_b" / "stage6", out / "reference.tle")
    speed = float(np.median(np.linalg.norm(r[["tle_vx_gcrs_m_s", "tle_vy_gcrs_m_s", "tle_vz_gcrs_m_s"]], axis=1)))
    st = s["statistics"]
    lines += [
        "PART B: end to end, Stage-4-like positions -> tart-track -> tart-validate",
        f"  truth = TLE trajectory at t + {TAU_S:g} s; Stage-4 noise {sig_los / 1e3:g} km along the line of sight, "
        f"{sig_t:g} m across; {len(t)} epochs over 30 min; exit codes tart-track {c5}, tart-validate {c6}",
        f"  expected along-track residual |v| tau = {speed * TAU_S / 1e3:.3f} km; radial and cross-track ~0",
        f"  recovered median R / T / N             ({st['radial']['median'] / 1e3:+.3f}, "
        f"{st['along_track']['median'] / 1e3:+.3f}, {st['cross_track']['median'] / 1e3:+.3f}) km",
        f"  median Stage-5 sigma R / T / N         ({st['radial']['median_sigma_obs_m'] / 1e3:.3f}, "
        f"{st['along_track']['median_sigma_obs_m'] / 1e3:.3f}, {st['cross_track']['median_sigma_obs_m'] / 1e3:.3f}) km",
        f"  equivalent along-track time offset     median {st['along_track_time_s']['median']:.3f} s (injected {TAU_S:g} s)",
        f"  3-D residual median {st['position_3d']['median'] / 1e3:.3f} km; angular median "
        f"{st['angular_deg']['median'] * 3600:.2f} arcsec; velocity residual median {st['velocity_3d']['median']:.3f} m/s",
        f"  conditional position D^2 median {st['position_mahalanobis_obs_cov']['median']:.1f} (the time shift is "
        "a real discrepancy, so D^2 >> 3 is expected)",
        "",
    ]
    # Monte Carlo with tau = 0: the TLE is the truth, so D^2 should follow chi^2_3.
    from scipy.stats import chi2

    truth = syn.independent_gcrs_states(L1, L2, t)
    ref = TLEReference(40534, "PRN 26", L1, L2, tle_epoch_utc(L1, L2), "tle-file")
    tle_states = propagate_tle(ref, t).state_gcrs
    probe = [0, len(t) // 4, len(t) // 2, 3 * len(t) // 4, len(t) - 1]
    d2, restarts = [], 0
    mc_rng = np.random.default_rng(11)
    for _ in range(n_mc):
        df = pd.DataFrame(syn5.stage4_rows(t, truth, sig_los, sig_t, mc_rng))
        df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
        df["quality_flags"] = df["quality_flags"].map(list)
        res = track(build_measurements(df, "strict"), FilterConfig())
        restarts += int(res.summary["n_resets"] > 0)
        o = StateTrajectory(t[probe], res.smoothed.x_smooth[probe], res.smoothed.p_smooth[probe])
        d2 += list(compute_residuals(o, StateTrajectory(t[probe], tle_states[probe])).d2_pos_obs)
    d2 = np.array(d2)
    lines += [
        f"  Monte Carlo, tau = 0 ({n_mc} realisations x {len(probe)} epochs; {restarts} with a Stage-5 restart):",
        f"    conditional position D^2: mean {d2.mean():.2f}, median {np.median(d2):.2f} "
        f"(chi2_3: 3, {chi2.ppf(0.5, 3):.2f}); fraction <= 1/4/9: {np.mean(d2 <= 1):.3f} / {np.mean(d2 <= 4):.3f} / "
        f"{np.mean(d2 <= 9):.3f} (expected {chi2.cdf(1, 3):.3f} / {chi2.cdf(4, 3):.3f} / {chi2.cdf(9, 3):.3f})",
        "",
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", type=Path, default=Path("docs/examples/stage6"))
    ap.add_argument("--mc", type=int, default=100, help="Monte Carlo realisations for part B (default 100)")
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    lines = ["Stage-6 synthetic demonstration (scripts/stage6_synthetic_demo.py)",
             f"Reference TLE (NORAD 40534): {L1} / {L2}", "Residual convention: observed - TLE", ""]
    part_a(out, lines)
    part_b(out, args.mc, lines)
    (out / "demo_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
