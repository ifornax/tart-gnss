#!/usr/bin/env python
"""Synthetic demonstration of Stage 5 (``tart-track``). No TLE is used anywhere.

Truth: a circular GNSS-like orbit (radius 26 560 km, inclination 55 deg)
placed above the three example stations at 2026-06-25T09:14Z and propagated
*analytically* (Kepler), never with the Stage-5 integrator.

Part A, Stage-4-like positions
    30 min of positions at a jittered 10 s cadence with a 4-minute gap and
    one injected 150 km outlier. The noise has the anisotropic shape of real
    Stage-4 covariances: 30 km along the line of sight from South Africa,
    3 km across it. ``tart-track`` runs with two-body dynamics (default) and
    constant velocity (its default 1 m/s^2 and the two-body 1e-4 m/s^2).
    A Monte Carlo over 200 noise realisations compares the predicted state
    uncertainty with the empirical error (coverage and NEES).

Part B, end to end through the pipeline
    The same orbit observed as directions from the three stations
    (``tart_gnss.stage4.synthetic.make_directions``, 0.02 deg isotropic noise),
    triangulated with ``tart-triangulate`` and tracked with ``tart-track``.

Usage::

    python scripts/stage5_synthetic_demo.py [--output-dir docs/examples/stage5] [--mc 200]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from tart_gnss.stage5 import cli as cli5
from tart_gnss.stage5 import synthetic as syn
from tart_gnss.stage5.cli import track
from tart_gnss.stage5.filter import FilterConfig
from tart_gnss.stage5.frames import itrs_to_gcrs_matrices
from tart_gnss.stage5.loader import build_measurements
from tart_gnss.stage5.output import read_table, state_covariance

T0 = syn.T0_UNIX_S
CADENCE_S = 10.0
DURATION_S = 1800.0
GAP_S = (900.0, 1140.0)
OUTLIER = (40, np.array([0.0, 1.5e5, 0.0]))
SIGMA_LOS_M, SIGMA_T_M = 3e4, 3e3


def epochs(rng) -> np.ndarray:
    t = np.arange(0.0, DURATION_S, CADENCE_S) + rng.uniform(-2.0, 2.0, int(DURATION_S / CADENCE_S))
    t[0] = 0.0
    t = t[(t < GAP_S[0]) | (t > GAP_S[1])]
    return T0 + t


def frame(rows):
    import pandas as pd

    df = pd.DataFrame(rows)
    df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    df["quality_flags"] = df["quality_flags"].map(list)
    return df


def states(traj) -> np.ndarray:
    cols = [f"{c}_gcrs_smoothed_m" for c in "xyz"] + [f"v{c}_gcrs_smoothed_m_s" for c in "xyz"]
    return traj[cols].to_numpy(dtype=float)


def rmse(a, b) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def part_a(out: Path, x0: np.ndarray, n_mc: int, lines: list) -> dict:
    rng = np.random.default_rng(20260625)
    t = epochs(rng)
    truth = syn.truth_states(x0, T0, t)
    rows = syn.stage4_rows(t, truth, SIGMA_LOS_M, SIGMA_T_M, rng, outliers={OUTLIER[0]: OUTLIER[1]})
    manifest = syn.write_stage4_run(out / "stage4_input", rows)
    runs = {
        "two_body": ["--dynamics", "two-body"],
        "constant_velocity": ["--dynamics", "constant-velocity", "--no-plots"],
    }
    tables = {}
    for name, extra in runs.items():
        code = cli5.main(["--stage4-manifest", str(manifest), "--output-dir", str(out / name), *extra])
        d = next((out / name).glob("*/NORAD_*"))
        tables[name] = (code, read_table(d / "trajectory.parquet"), json.loads((d / "stage5_manifest.json").read_text()))
    code, traj, man = tables["two_body"]
    est = states(traj)
    s = man["innovation_statistics"]
    ps = [np.sqrt(np.trace(state_covariance(r, "smoothed")[:3, :3])) for _, r in traj.iterrows()]
    vs = np.sqrt(traj["sigma_vx_smoothed_m_s"] ** 2 + traj["sigma_vy_smoothed_m_s"] ** 2
                 + traj["sigma_vz_smoothed_m_s"] ** 2)
    mats = itrs_to_gcrs_matrices(t)
    stage4_err = [np.linalg.norm(m @ np.array([r["x_ecef_m"], r["y_ecef_m"], r["z_ecef_m"]]) - x[:3])
                  for m, r, x in zip(mats, rows, truth)]
    lines += [
        "PART A: Stage-4-like positions (no TLE anywhere)",
        f"  epochs {len(t)} over {DURATION_S / 60:.0f} min, jittered {CADENCE_S:.0f} s cadence, gap "
        f"{GAP_S[1] - GAP_S[0]:.0f} s, outlier at epoch {OUTLIER[0]} (+{np.linalg.norm(OUTLIER[1]) / 1e3:.0f} km)",
        f"  Stage-4 noise: sigma_LOS {SIGMA_LOS_M / 1e3:.0f} km, sigma_transverse {SIGMA_T_M / 1e3:.0f} km; "
        f"median single-epoch 3-D error {np.median(stage4_err) / 1e3:.1f} km",
        f"  tart-track two-body exit code {code}",
        f"  true initial position (GCRS, km)     {np.array2string(truth[0][:3] / 1e3, precision=3)}",
        f"  true initial velocity (GCRS, m/s)    {np.array2string(truth[0][3:], precision=2)}",
        f"  smoothed initial position (km)       {np.array2string(est[0][:3] / 1e3, precision=3)}",
        f"  smoothed initial velocity (m/s)      {np.array2string(est[0][3:], precision=2)}",
        f"  error at t0: position {np.linalg.norm(est[0][:3] - truth[0][:3]) / 1e3:.3f} km, velocity "
        f"{np.linalg.norm(est[0][3:] - truth[0][3:]):.2f} m/s",
        f"  position RMSE (all epochs) {rmse(est[:, :3], truth[:, :3]) / 1e3:.3f} km; velocity RMSE "
        f"{rmse(est[:, 3:], truth[:, 3:]):.2f} m/s",
        f"  predicted (smoothed) 3-D sigma: position median {np.median(ps) / 1e3:.3f} km (max "
        f"{np.max(ps) / 1e3:.3f} km, in the gap), velocity median {np.median(vs):.2f} m/s",
        f"  mean / median NIS {s['mean_nis']:.2f} / {s['median_nis']:.2f} (chi2_3: 3 / 2.37); innovation "
        f"rejection fraction {1 - s['accepted_fraction']:.3f} ({s['n_outliers']} of {s['n_attempted']}); "
        f"outlier epoch gated: {bool(traj['quality_flags'].iloc[OUTLIER[0]].count('INNOVATION_OUTLIER'))}",
        f"  lag-1 innovation autocorrelation {s['autocorrelation']['1']:.3f} (limit "
        f"{s['autocorrelation_limit_lag1']:.3f}); run flags {man['summary']['run_flags'] or 'none'}",
        "",
    ]
    # Dynamics comparison (TLE-independent diagnostics only).
    lines += ["  dynamics comparison (TLE-independent diagnostics):",
              f"    {'model':<34}{'sigma_a':>9}{'innov RMS km':>14}{'mean NIS':>10}{'median NIS':>11}"
              f"{'lag-1 acf':>10}{'accepted':>9}{'pos RMSE km':>12}"]
    df = frame(rows)
    for label, dyn, sa in (("two-body", "two-body", None), ("constant-velocity (default sigma_a)",
                                                              "constant-velocity", None),
                           ("constant-velocity (two-body sigma_a)", "constant-velocity", 1e-4)):
        res = track(build_measurements(df, "strict"), FilterConfig(dynamics=dyn, process_accel_sigma_m_s2=sa))
        s2 = res.summary
        lines.append(f"    {label:<34}{res.config.sigma_a:>9.0e}{s2['innovation_rms_m'] / 1e3:>14.2f}"
                     f"{s2['mean_nis']:>10.2f}{s2['median_nis']:>11.2f}{s2['autocorrelation']['1']:>10.3f}"
                     f"{s2['accepted_fraction']:>9.3f}{rmse(res.smoothed.x_smooth[:, :3], truth[:, :3]) / 1e3:>12.3f}")
    lines.append("")

    # Monte Carlo: predicted vs empirical error, coverage, NEES.
    from scipy.stats import chi2

    mc_rng = np.random.default_rng(7)
    dpos, dvel, nees, errs_p, errs_v, pred_p, pred_v = [], [], [], [], [], [], []
    probe = [0, len(t) // 4, len(t) // 2, len(t) - 1]
    for _ in range(n_mc):
        r_mc = syn.stage4_rows(t, truth, SIGMA_LOS_M, SIGMA_T_M, mc_rng)
        res = track(build_measurements(frame(r_mc), "strict"), FilterConfig())
        for k in probe:
            e = res.smoothed.x_smooth[k] - truth[k]
            p = res.smoothed.p_smooth[k]
            nees.append(e @ np.linalg.solve(p, e))
            dpos.append(e[:3] @ np.linalg.solve(p[:3, :3], e[:3]))
            dvel.append(e[3:] @ np.linalg.solve(p[3:, 3:], e[3:]))
            errs_p.append(e[:3] @ e[:3])
            errs_v.append(e[3:] @ e[3:])
            pred_p.append(np.trace(p[:3, :3]))
            pred_v.append(np.trace(p[3:, 3:]))
    dpos, dvel, nees = np.array(dpos), np.array(dvel), np.array(nees)
    lines += [f"  Monte Carlo ({n_mc} realisations, epochs {probe}; two-body, no outlier):",
              f"    3-D position error: empirical RMS {np.sqrt(np.mean(errs_p)) / 1e3:.3f} km vs predicted "
              f"{np.sqrt(np.mean(pred_p)) / 1e3:.3f} km",
              f"    3-D velocity error: empirical RMS {np.sqrt(np.mean(errs_v)):.3f} m/s vs predicted "
              f"{np.sqrt(np.mean(pred_v)):.3f} m/s",
              f"    mean NEES {nees.mean():.2f} (chi2_6: 6); per-run sd of the mean ~{np.sqrt(12 / n_mc):.2f}",
              "    coverage           1-sigma  2-sigma  3-sigma ellipsoid",
              f"    expected (chi2_3)  {chi2.cdf(1, 3):7.3f}  {chi2.cdf(4, 3):7.3f}  {chi2.cdf(9, 3):7.3f}",
              f"    position           {np.mean(dpos <= 1):7.3f}  {np.mean(dpos <= 4):7.3f}  {np.mean(dpos <= 9):7.3f}",
              f"    velocity           {np.mean(dvel <= 1):7.3f}  {np.mean(dvel <= 4):7.3f}  {np.mean(dvel <= 9):7.3f}",
              ""]
    return {"truth_initial_state_gcrs": truth[0].tolist(), "epochs_unix_s": t.tolist(),
            "outlier_epoch": OUTLIER[0], "outlier_offset_ecef_m": OUTLIER[1].tolist()}


def part_b(out: Path, x0: np.ndarray, lines: list) -> dict:
    from tart_gnss.stage4 import cli as cli4
    from tart_gnss.stage4.geometry import ecef_to_azel, geodetic_to_ecef
    from tart_gnss.stage4.synthetic import EXAMPLE_STATIONS, make_directions, write_stage3_run

    rng = np.random.default_rng(42)
    t = T0 + np.arange(0.0, 900.0, 10.0)

    def truth_ecef(tk):
        m = itrs_to_gcrs_matrices([tk])[0]
        return m.T @ syn.kepler_propagate(x0, tk - T0)[:3]

    min_el = min(ecef_to_azel(truth_ecef(tk), geodetic_to_ecef(*st), st[0], st[1])[1]
                 for st in EXAMPLE_STATIONS.values() for tk in t)
    rows = make_directions(truth_ecef, {sid: t for sid in EXAMPLE_STATIONS}, 0.02, rng)
    m3 = write_stage3_run(out / "end_to_end" / "stage3_input", rows)
    c4 = cli4.main(["--stage3-manifest", str(m3), "--output-dir", str(out / "end_to_end" / "stage4"), "--no-plots"])
    m4 = next((out / "end_to_end" / "stage4").glob("*/NORAD_*/stage4_manifest.json"))
    c5 = cli5.main(["--stage4-manifest", str(m4), "--output-dir", str(out / "end_to_end" / "stage5")])
    d5 = next((out / "end_to_end" / "stage5").glob("*/NORAD_*"))
    traj = read_table(d5 / "trajectory.parquet")
    man = json.loads((d5 / "stage5_manifest.json").read_text())
    truth = syn.truth_states(x0, T0, t)
    est = states(traj)
    pos4 = read_table(m4.parent / "positions.parquet")
    e4 = [np.linalg.norm(np.array([r.x_ecef_m, r.y_ecef_m, r.z_ecef_m]) - truth_ecef(tk))
          for tk, (_, r) in zip(t, pos4.iterrows())]
    s = man["innovation_statistics"]
    lines += [
        "PART B: end to end (directions -> tart-triangulate -> tart-track), 3 stations, 0.02 deg direction noise",
        f"  {len(t)} epochs over {(t[-1] - t[0]) / 60:.0f} min; minimum elevation {min_el:.1f} deg; exit codes "
        f"tart-triangulate {c4}, tart-track {c5}",
        f"  Stage-4 single-epoch 3-D error: median {np.median(e4) / 1e3:.1f} km",
        f"  Stage-5 smoothed: position RMSE {rmse(est[:, :3], truth[:, :3]) / 1e3:.2f} km, velocity RMSE "
        f"{rmse(est[:, 3:], truth[:, 3:]):.2f} m/s; median predicted 3-D sigma "
        f"{np.median(np.sqrt(traj.sigma_x_smoothed_m ** 2 + traj.sigma_y_smoothed_m ** 2 + traj.sigma_z_smoothed_m ** 2)) / 1e3:.2f} km",
        f"  mean / median NIS {s['mean_nis']:.2f} / {s['median_nis']:.2f}; lag-1 autocorrelation "
        f"{s['autocorrelation']['1']:.3f}; run flags {man['summary']['run_flags'] or 'none'}",
        "",
    ]
    return {"end_to_end_epochs_unix_s": t.tolist(), "direction_sigma_deg": 0.02}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", type=Path, default=Path("docs/examples/stage5"))
    ap.add_argument("--mc", type=int, default=200, help="Monte Carlo realisations (default 200)")
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    x0 = syn.example_initial_state()
    lines = ["Stage-5 synthetic demonstration (scripts/stage5_synthetic_demo.py)",
             "Truth: circular orbit r = 26 560 km, i = 55 deg, analytic Kepler propagation; no TLE.", ""]
    truth = {"t0_unix_s": T0, "initial_state_gcrs": x0.tolist(), "propagation": "analytic Kepler (f and g)"}
    truth.update(part_a(out, x0, args.mc, lines))
    truth.update(part_b(out, x0, lines))
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n", encoding="utf-8")
    (out / "demo_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
