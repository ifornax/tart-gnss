"""Synthetic validation (spec section 34) through the full tart-validate path."""

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("skyfield")

from tart_gnss.stage6 import synthetic as syn  # noqa: E402

from .helpers import L1, L2, TLE_EPOCH_S, inject_rtn, run_cli, stage5_run, times, tle_file, truth  # noqa: E402


def test_zero_residual(tmp_path):  # 34.1 (observed built by the independent TEME route)
    t = times(10)
    code, res, _, _ = run_cli(tmp_path, stage5_run(tmp_path, t, truth(t)), "--tle-file", str(tle_file(tmp_path)),
                              "--no-plots")
    assert code == 0
    assert res["position_residual_3d_m"].max() < 0.1 and res["velocity_residual_3d_m_s"].max() < 1e-3
    assert res["geocentric_angular_residual_deg"].max() < 1e-9 * 180


def test_known_rtn_and_velocity_offsets(tmp_path):  # 34.2 - 34.4
    t = times(10)
    obs = inject_rtn(truth(t), (1000.0, 2000.0, -500.0), (0.5, -0.2, 0.1))
    _, res, summ, _ = run_cli(tmp_path, stage5_run(tmp_path, t, obs), "--tle-file", str(tle_file(tmp_path)),
                              "--no-plots")
    assert np.allclose(res["dr_radial_m"], 1000.0, atol=0.1)
    assert np.allclose(res["dr_along_track_m"], 2000.0, atol=0.1)
    assert np.allclose(res["dr_cross_track_m"], -500.0, atol=0.1)
    assert np.allclose(res[["dv_radial_m_s", "dv_along_track_m_s", "dv_cross_track_m_s"]], [0.5, -0.2, 0.1],
                       atol=1e-3)
    assert np.allclose(res["position_residual_3d_m"], np.sqrt(1000 ** 2 + 2000 ** 2 + 500 ** 2), atol=0.1)
    assert summ["statistics"]["along_track"]["median"] == pytest.approx(2000.0, abs=0.1)
    # Cartesian residuals carry the observed - TLE sign.
    d = res[["obs_x_gcrs_m", "obs_y_gcrs_m", "obs_z_gcrs_m"]].to_numpy() - \
        res[["tle_x_gcrs_m", "tle_y_gcrs_m", "tle_z_gcrs_m"]].to_numpy()
    assert np.allclose(d, res[["dr_x_m", "dr_y_m", "dr_z_m"]].to_numpy(), atol=1e-6)


def test_covariance_in_rtn(tmp_path):  # 34.5
    t = times(4)
    x = truth(t)
    from tart_gnss.stage6.rtn import rtn_bases

    q = rtn_bases(x[:, :3], x[:, 3:])
    c_rtn = np.diag([100.0, 900.0, 400.0]) ** 2
    cov = np.zeros((4, 6, 6))
    for k in range(4):
        cov[k, :3, :3] = q[k].T @ c_rtn @ q[k]
        cov[k, 3:, 3:] = np.eye(3)
    m = syn.write_stage5_run(tmp_path / "s5", t, x, cov)
    _, res, _, _ = run_cli(tmp_path, m, "--tle-file", str(tle_file(tmp_path)), "--no-plots")
    assert np.allclose(res["sigma_radial_obs_m"], 100.0) and np.allclose(res["sigma_along_track_obs_m"], 900.0)
    assert np.allclose(res["sigma_cross_track_obs_m"], 400.0)
    assert np.allclose(res[["cov_radial_along_obs_m2", "cov_radial_cross_obs_m2", "cov_along_cross_obs_m2"]], 0.0,
                       atol=1e-3)


@pytest.mark.parametrize("offset_days,flags,precise", [
    (0.1, set(), True), (5.0, {"TLE_AGE_WARNING"}, True), (20.0, {"TLE_STALE"}, False),
    (-1.0, {"TLE_BACKPROPAGATED"}, True),
])
def test_tle_age_cases(tmp_path, offset_days, flags, precise):  # 34.6
    t = times(4, t0=TLE_EPOCH_S + offset_days * 86400)
    code, res, summ, _ = run_cli(tmp_path, stage5_run(tmp_path, t, truth(t)), "--tle-file",
                                 str(tle_file(tmp_path)), "--no-plots")
    got = set(res["quality_flags"].iloc[0]) & {"TLE_AGE_WARNING", "TLE_STALE", "TLE_BACKPROPAGATED"}
    assert got == flags and summ["precision_validation"] is precise and code == (0 if precise else 3)


def test_single_tle_has_no_reference_covariance(tmp_path):  # 34.7
    t = times(4)
    _, res, _, _ = run_cli(tmp_path, stage5_run(tmp_path, t, inject_rtn(truth(t), (3000.0, 0, 0))),
                           "--tle-file", str(tle_file(tmp_path)), "--no-plots")
    assert np.allclose(res["position_mahalanobis_obs_cov"], 9.0, rtol=1e-6)  # (3 km / 1 km)^2
    assert res["position_mahalanobis_total_cov"].isna().all()
    assert all("REFERENCE_UNCERTAINTY_UNKNOWN" in f for f in res["quality_flags"])


def test_ensemble_dispersion_and_total_covariance(tmp_path):  # 34.8
    t = times(6)
    sets = [(L1, L2)] + [syn.perturbed_tle(L1, L2, d_mean_anomaly_deg=d, d_epoch_days=e)
                         for d, e in ((0.004, -0.5), (-0.003, 0.7), (0.002, 1.5), (-0.001, -2.0))]
    ens = syn.write_tle_file(tmp_path / "ens" / "e.tle", sets)
    obs = inject_rtn(truth(t), (0.0, 3000.0, 0.0))
    _, res, summ, man = run_cli(tmp_path, stage5_run(tmp_path, t, obs), "--tle-file", str(tle_file(tmp_path)),
                                "--no-plots", "--tle-ensemble-dir", str(ens.parent))
    assert (res["n_tle_ensemble"] == 5).all() and man["tle_ensemble"]["label"] == "empirical reference dispersion"
    assert man["tle_ensemble"]["primary_in_ensemble"] is True
    assert res["position_mahalanobis_total_cov"].notna().all()
    # Adding the reference dispersion (mostly along-track) can only reduce the distance.
    assert (res["position_mahalanobis_total_cov"] < res["position_mahalanobis_obs_cov"]).all()
    assert (res["sigma_along_track_ref_m"] > res["sigma_radial_ref_m"]).all()
    assert not any("REFERENCE_UNCERTAINTY_UNKNOWN" in f for f in res["quality_flags"])
    assert summ["ensemble"]["n_tle_ensemble"] == 5
