"""Synthetic validation of Stage 5 (specification section 32)."""

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")

from tart_gnss.stage5 import synthetic as syn  # noqa: E402

from .helpers import run_cli, run_track, smoothed_states, times, truth  # noqa: E402


def _rmse(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def test_exact_two_body_trajectory_is_recovered():  # 32.1
    rng = np.random.default_rng(10)
    t = times(60)
    x = truth(t)
    res = run_track(syn.stage4_rows(t, x, 5e3, 1e3, rng))
    xs, ps = res.smoothed.x_smooth, res.smoothed.p_smooth
    pos_err = _rmse(xs[:, :3], x[:, :3])
    vel_err = _rmse(xs[:, 3:], x[:, 3:])
    sig_pos = np.median([np.sqrt(np.trace(p[:3, :3])) for p in ps])
    sig_vel = np.median([np.sqrt(np.trace(p[3:, 3:])) for p in ps])
    # Smoothing 60 measurements of 5 km (LOS) / 1 km beats a single Stage-4 epoch by a large factor.
    assert pos_err < 1.5e3 and pos_err < 3 * sig_pos
    assert vel_err < 10.0 and vel_err < 3 * sig_vel
    # The errors are consistent with the reported covariance (NEES of the full state, per epoch).
    nees = [(xs[k] - x[k]) @ np.linalg.solve(ps[k], xs[k] - x[k]) for k in range(len(t))]
    assert np.median(nees) < 20.0


def test_full_covariance_beats_diagonal_approximation():  # 32.2
    # Small transverse, large line-of-sight error: dropping the off-diagonal terms mis-weights the data
    # (less accurate) and mis-states the uncertainty (inconsistent); the full covariance is calibrated.
    t = times(30)
    x = truth(t)
    err, nees = {"full": [], "diag": []}, {"full": [], "diag": []}
    for seed in range(20):
        rows = syn.stage4_rows(t, x, 5e4, 1e3, np.random.default_rng(seed))
        runs = {"full": run_track(rows)}
        for r in rows:
            r["cov_xy"] = r["cov_xz"] = r["cov_yz"] = 0.0
        runs["diag"] = run_track(rows)
        for key, res in runs.items():
            e = res.smoothed.x_smooth[:, :3] - x[:, :3]
            err[key].append(np.sqrt(np.mean(np.sum(e ** 2, axis=1))))
            nees[key].append(np.mean([ek @ np.linalg.solve(p[:3, :3], ek) for ek, p in zip(e, res.smoothed.p_smooth)]))
    assert np.mean(err["full"]) < 0.9 * np.mean(err["diag"])
    assert 2.0 < np.mean(nees["full"]) < 4.2  # chi^2_3 mean 3
    assert np.mean(nees["diag"]) > 1.2 * np.mean(nees["full"])
    # The smoothed covariance stays anisotropic, its long axis along the station line of sight.
    from tart_gnss.stage4.geometry import geodetic_to_ecef

    p = runs["full"].smoothed.p_smooth[15][:3, :3]
    w, v = np.linalg.eigh(p)
    m = runs["full"].measurements.matrices[15]
    los = m @ (m.T @ x[15, :3] - geodetic_to_ecef(*syn.REFERENCE_STATION))
    assert w[-1] > 4 * w[0]
    assert abs(v[:, -1] @ los) / np.linalg.norm(los) > 0.95


def test_irregular_cadence():  # 32.3
    rng = np.random.default_rng(11)
    t = np.sort(syn.T0_UNIX_S + rng.uniform(0, 600, 50))
    x = truth(t)
    res = run_track(syn.stage4_rows(t, x, 5e3, 1e3, rng))
    assert np.allclose(np.diff(res.filtered.times_s), np.diff(t))
    assert _rmse(res.smoothed.x_smooth[:, :3], x[:, :3]) < 2e3
    assert 1.5 < res.summary["mean_nis"] < 4.5


def test_missing_epochs_increase_uncertainty(tmp_path):  # 32.4
    rng = np.random.default_rng(12)
    t = times(60)
    keep = np.r_[0:20, 40:60]  # 200 s with no measurement
    rows = syn.stage4_rows(t[keep], truth(t[keep]), 5e3, 1e3, rng)
    code, traj, _, _ = run_cli(tmp_path, rows, "--no-plots", "--max-gap-s", "60")
    assert code == 0
    sig = np.sqrt(traj["sigma_x_filtered_m"] ** 2 + traj["sigma_y_filtered_m"] ** 2 + traj["sigma_z_filtered_m"] ** 2)
    after_gap = 20
    assert "LARGE_TIME_GAP" in traj["quality_flags"].iloc[after_gap]
    assert traj["time_since_update_s"].iloc[after_gap] == pytest.approx(210.0)
    # The prediction across the gap is wider than any filtered state just before it.
    res = run_track(rows)
    p_pred = res.filtered.p_pred[after_gap]
    assert np.trace(p_pred[:3, :3]) > 4 * np.trace(res.filtered.p_filt[after_gap - 1][:3, :3])
    # Filtered uncertainty grows monotonically through a prediction-only stretch, and the smoothed
    # uncertainty inside the gap exceeds that of the same epoch with all data present.
    rows_all = syn.stage4_rows(t, truth(t), 5e3, 1e3, rng)
    full = run_track(rows_all)
    for k in range(20, 40):
        rows_all[k]["quality_flags"] = ["TRIANGULATION_FAILED"]
    gap = run_track(rows_all)
    tr_f = [np.trace(p[:3, :3]) for p in gap.filtered.p_filt[20:40]]
    assert np.all(np.diff(tr_f) > 0)
    s_gap = np.sqrt([np.trace(p[:3, :3]) for p in gap.smoothed.p_smooth])
    s_full = np.sqrt([np.trace(p[:3, :3]) for p in full.smoothed.p_smooth])
    assert s_gap[30] > 1.1 * s_full[30]
    assert np.all(sig.notna())


def test_outlier_is_detected_and_kept():  # 32.5
    rng = np.random.default_rng(13)
    t = times(40)
    rows = syn.stage4_rows(t, truth(t), 5e3, 1e3, rng, outliers={25: np.array([0.0, 8e4, 0.0])})
    res = run_track(rows)
    out = np.flatnonzero(res.filtered.outlier)
    assert list(out) == [25]
    assert not res.filtered.accepted[25]
    flags = res.trajectory[25]["quality_flags"]
    assert "INNOVATION_OUTLIER" in flags and res.trajectory[25]["valid"]
    assert any(r["time_utc"] == res.trajectory[25]["time_utc"] for r in res.innovations)


def test_process_noise_mismatch():  # 32.6
    # Truth has an unmodelled 5e-3 m/s^2 acceleration. With almost no process noise the filter is
    # inconsistent: measurements exceed the gate far more often than chance allows and the filter has to
    # restart; process noise of that size removes this. (A *constant* acceleration is not white, so its
    # innovations stay somewhat autocorrelated even then: DYNAMICS_MISMATCH may remain, as it should.)
    rng = np.random.default_rng(14)
    t = times(80, step=20.0)
    x = truth(t, extra_accel=np.array([1.0, -0.5, 0.5]) * 5e-3)
    rows = syn.stage4_rows(t, x, 100.0, 100.0, rng)
    tight = run_track(rows, sigma_a=1e-7).summary
    right = run_track(rows, sigma_a=5e-3).summary
    assert tight["n_resets"] >= 1 and tight["n_above_gate"] > tight["n_above_gate_limit"]
    assert {"HIGH_NIS", "DYNAMICS_MISMATCH"} <= set(tight["flags"])
    assert right["n_resets"] == 0 and right["n_above_gate"] == 0
    assert right["mean_nis"] < tight["mean_nis"]
    assert right["autocorrelation"]["1"] < tight["autocorrelation"]["1"]


def test_two_body_outperforms_constant_velocity():
    rng = np.random.default_rng(15)
    t = times(60, step=20.0)
    rows = syn.stage4_rows(t, truth(t), 1e3, 300.0, rng)
    tb = run_track(rows, sigma_a=1e-4)
    cv = run_track(rows, dynamics="constant-velocity", sigma_a=1e-4)
    assert tb.summary["innovation_rms_m"] < 0.5 * cv.summary["innovation_rms_m"]
    assert tb.summary["median_nis"] < cv.summary["median_nis"]
    assert tb.summary["accepted_fraction"] >= cv.summary["accepted_fraction"]
    assert abs(tb.summary["autocorrelation"]["1"]) < cv.summary["autocorrelation"]["1"]


def test_pipeline_runs_without_any_tle(tmp_path, monkeypatch):  # 32.7
    # Any TLE access fails loudly; the Stage-4 manifest even carries a TLE block, which must not be read.
    import skyfield.sgp4lib

    import tart_gnss.tle as tle

    def boom(*a, **k):
        raise AssertionError("TLE data used by Stage 5")

    for name in dir(tle):
        if callable(getattr(tle, name)) and not name.startswith("_") and name not in ("sha256_file", "format_utc"):
            obj = getattr(tle, name)
            if getattr(obj, "__module__", "") == tle.__name__:
                monkeypatch.setattr(tle, name, boom)
    monkeypatch.setattr(skyfield.sgp4lib.EarthSatellite, "__init__", boom)
    rng = np.random.default_rng(16)
    t = times(15)
    code, traj, _, man = run_cli(tmp_path, syn.stage4_rows(t, truth(t), 5e3, 1e3, rng), "--no-plots",
                                 with_tle_block=True)
    assert code == 0 and traj["valid"].all()
    assert man["tle_used"] is False
    assert "satellite" not in json.dumps(man["stage4_reference"]).lower().replace("satellite/tle", "")


import json  # noqa: E402
