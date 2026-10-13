"""End-to-end synthetic scenarios through ``tart-triangulate`` (spec section 25)."""

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from tart_gnss.stage4 import synthetic as syn  # noqa: E402

from .helpers import NAMIBIA, run_cli, truth_at  # noqa: E402

TRUTH = syn.linear_track()
T = [syn.T0_UNIX_S + k for k in range(20)]


def times(offsets=None, stations=syn.EXAMPLE_STATIONS):
    offsets = offsets or {}
    return {sid: [t + offsets.get(sid, 0.0) for t in T] for sid in stations}


def xyz(pos):
    return pos[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy(dtype=float)


@pytest.mark.parametrize("method", ["weighted", "unweighted"])
def test_noise_free_exact_recovery(tmp_path, method):
    rows = syn.make_directions(TRUTH, times(), 0.1, None)
    code, pos, con, man = run_cli(tmp_path, rows, "--method", method)
    assert code == 0
    assert len(pos) == 20 and pos["valid"].all() and (pos["n_stations"] == 3).all()
    err = np.linalg.norm(xyz(pos) - truth_at(pos, TRUTH), axis=1)
    assert err.max() < 1.0  # metres, at ~20 000 km range
    assert con["line_residual_m"].max() < 1.0
    assert all(f == [] for f in pos["quality_flags"])


def test_angular_noise_statistically_consistent(tmp_path):
    t = {sid: [syn.T0_UNIX_S + k for k in range(200)] for sid in syn.EXAMPLE_STATIONS}
    rows = syn.make_directions(TRUTH, t, 0.1, np.random.default_rng(11))
    code, pos, con, man = run_cli(tmp_path, rows)
    assert code == 0
    e = xyz(pos) - truth_at(pos, TRUTH)
    sig = pos[["sigma_x_m", "sigma_y_m", "sigma_z_m"]].to_numpy(dtype=float)
    pulls = e / sig
    assert np.all(np.abs(pulls.mean(axis=0)) < 0.25)
    assert np.all((0.85 < pulls.std(axis=0)) & (pulls.std(axis=0) < 1.15))
    cov = np.stack([pos[["cov_xx", "cov_xy", "cov_xz"]].to_numpy(float), pos[["cov_xy", "cov_yy", "cov_yz"]].to_numpy(float),
                    pos[["cov_xz", "cov_yz", "cov_zz"]].to_numpy(float)], axis=1)
    m2 = np.einsum("ni,ni->n", e, np.linalg.solve(cov, e[:, :, None])[:, :, 0])
    assert 2.5 < m2.mean() < 3.5
    assert 0.7 < np.nanmedian(pos["reduced_chi2"]) < 1.3


def test_poor_geometry_increases_condition_and_uncertainty(tmp_path):
    wide = syn.make_directions(TRUTH, times(), 0.1, None)
    lat, lon, alt = syn.EXAMPLE_STATIONS["mauritius"]
    close = {"a": (lat, lon, alt), "b": (lat + 0.5, lon, alt), "c": (lat, lon + 0.5, alt)}
    narrow = syn.make_directions(TRUTH, times(stations=close), 0.1, None, stations=close)
    _, pw, _, _ = run_cli(tmp_path / "w", wide)
    code, pn, _, _ = run_cli(tmp_path / "n", narrow, stations=close)
    assert pn["condition_number"].median() > 10 * pw["condition_number"].median()
    assert pn["min_crossing_angle_deg"].median() < 0.2 * pw["min_crossing_angle_deg"].median()
    for c in ("sigma_x_m", "sigma_y_m", "sigma_z_m"):
        assert pn[c].median() > 3 * pw[c].median()
    # A tiny network trips POOR_GEOMETRY at a stricter threshold, and those epochs are not valid.
    code, pp, _, _ = run_cli(tmp_path / "p", narrow, "--max-condition-number", "1e3", stations=close)
    assert code == 3 and all("POOR_GEOMETRY" in f for f in pp["quality_flags"]) and not pp["valid"].any()


def test_timestamp_offsets_and_motion_diagnostics(tmp_path):
    # 10 s cadence (as for 10 s Stage-3 bins), so the nearest sample really is 2-3 s away.
    long_t = {sid: [syn.T0_UNIX_S + 10.0 * k for k in range(60)] for sid in syn.EXAMPLE_STATIONS}
    off = {"madagascar": 2.0, "south_africa": -3.0}
    shifted = {sid: [t + off.get(sid, 0.0) for t in ts] for sid, ts in long_t.items()}
    rows = syn.make_directions(TRUTH, shifted, 0.002, None)
    code, pos, con, man = run_cli(tmp_path, rows, "--time-match-tolerance", "3.5",
                                  "--reference-station", "mauritius")
    full = pos[pos["n_stations"] == 3]
    assert len(full) >= 58
    assert np.allclose(full["max_time_offset_s"], 5.0, atol=1e-3)
    c = con[con["station_id"] == "south_africa"]
    assert np.allclose(c["station_time_offset_s"], -3.0, atol=1e-3)
    # Motion during the offset follows the true angular rate seen from that station.
    lat, lon, alt = syn.EXAMPLE_STATIONS["south_africa"]
    l0, m0 = syn.exact_lm(TRUTH(syn.T0_UNIX_S + 60), lat, lon, alt)
    l1, m1 = syn.exact_lm(TRUTH(syn.T0_UNIX_S + 61), lat, lon, alt)
    s0, s1 = [np.array([a, b, np.sqrt(1 - a * a - b * b)]) for a, b in ((l0, m0), (l1, m1))]
    true_rate = np.degrees(np.arccos(np.clip(s0 @ s1, -1, 1)))
    assert np.median(c["angular_rate_deg_s"]) == pytest.approx(true_rate, rel=0.05)
    assert np.median(c["motion_during_offset_deg"]) == pytest.approx(3.0 * true_rate, rel=0.05)
    # 3 s of motion is >> the 0.002 deg measurement error: TIME_MISMATCH on the station and the epoch.
    assert all("TIME_MISMATCH" in f for f in full["quality_flags"])
    assert all("TIME_MISMATCH" in f for f in c["stage4_station_flags"])
    assert full["valid"].all()  # a warning, not blocking
    assert 0.06 < full["mean_light_time_s"].median() < 0.075
    # Default tolerance: the offset stations do not match at all.
    code, pos2, _, _ = run_cli(tmp_path / "default", rows)
    assert code == 3 and (pos2["n_stations"] < 3).all()


def test_bad_station_identified_by_line_residuals(tmp_path):
    st = dict(syn.EXAMPLE_STATIONS, namibia=NAMIBIA)
    t = {sid: T for sid in st}
    for bad in st:
        rows = syn.make_directions(TRUTH, t, 0.05, np.random.default_rng(3), stations=st,
                                   bias_deg={bad: (0.8, -0.6)})
        code, pos, con, _ = run_cli(tmp_path / bad, rows, stations=st)
        worst = con.loc[con.groupby("time_utc")["normalised_residual"].idxmax(), "station_id"]
        assert (worst == bad).mean() >= 0.9, bad
        # The bad station is flagged at every epoch (a 20-sigma error also drags the fit off the others).
        mine = con[con["station_id"] == bad]
        assert all("HIGH_LINE_RESIDUAL" in f for f in mine["stage4_station_flags"])
        assert all("HIGH_CHI2" in f for f in pos["quality_flags"])


def test_two_station_case_is_degraded_diagnostic(tmp_path):
    two = {k: syn.EXAMPLE_STATIONS[k] for k in ("mauritius", "south_africa")}
    rows = syn.make_directions(TRUTH, times(stations=two), 0.1, None, stations=two)
    code, pos, con, man = run_cli(tmp_path, rows, stations=two)
    assert code == 3
    assert (pos["n_stations"] == 2).all() and not pos["valid"].any()
    assert all("INSUFFICIENT_STATIONS" in f for f in pos["quality_flags"])
    # The diagnostic solution is still computed (and exact here), with dof = 1.
    assert np.linalg.norm(xyz(pos) - truth_at(pos, TRUTH), axis=1).max() < 1.0
    assert (pos["dof"] == 1).all() and pos["fit_success"].all()
    code, pos2, _, _ = run_cli(tmp_path / "m2", rows, "--min-stations", "2", stations=two)
    assert code == 0 and pos2["valid"].all()


def test_negative_range_flagged(tmp_path):
    # Two stations looking away from each other: the lines' closest approach is behind both.
    st = {"west": (0.0, 20.0, 0.0), "east": (0.0, 30.0, 0.0)}
    rows = []
    for sid, az in (("west", 270.0), ("east", 90.0)):
        el = np.radians(20.0)
        l, m = np.sin(np.radians(az)) * np.cos(el), np.cos(np.radians(az)) * np.cos(el)
        for t in T[:3]:
            rows.append(syn.direction_row(sid, t, l, m, syn.isotropic_lm_covariance(l, m, np.radians(0.1))))
    code, pos, con, _ = run_cli(tmp_path, rows, "--min-stations", "2", stations=st)
    assert code == 3
    assert all("NEGATIVE_RANGE" in f for f in pos["quality_flags"]) and not pos["valid"].any()
    assert (con["range_m"] < 0).all()
    assert all("NEGATIVE_RANGE" in f for f in con["stage4_station_flags"])
