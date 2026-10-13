"""Recovery of injected synthetic sources (spec section 29)."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from tart_gnss.stage3.estimator import EstimatorOptions, SearchRegion, estimate_direction  # noqa: E402
from tart_gnss.stage3.geometry import angular_separation_deg, azel_to_unit  # noqa: E402
from tart_gnss.stage3.imaging import image_direction  # noqa: E402
from tart_gnss.stage3.quality import BOUNDARY_SOLUTION  # noqa: E402
from tart_gnss.stage3.uncertainty import estimate_direction_uncertainty  # noqa: E402

from .helpers import make_samples  # noqa: E402


def _err(est, az, el):
    return float(angular_separation_deg(est.unit_vector, azel_to_unit(az, el)))


def test_noise_free_recovery_to_numerical_precision(layout):
    sel, _ = make_samples(layout, 214.0, 43.0)
    est = estimate_direction(sel.samples, 216.0, 41.0)
    assert est.success and _err(est, 214.0, 43.0) < 1e-8
    assert np.isclose(est.complex_amplitude, 0.5 + 0.2j)


@pytest.mark.parametrize("az,el", [(90.0, 40.0), (270.0, 40.0), (0.0, 30.0), (180.0, 30.0)])
def test_no_east_west_or_north_south_flip(layout, az, el):
    """Blind search: a source East must not come back West (sign-convention test)."""
    sel, _ = make_samples(layout, az, el)
    est = estimate_direction(sel.samples, options=EstimatorOptions(blind=True))
    assert _err(est, az, el) < 1e-6


def test_gaussian_noise_pulls_are_consistent(layout):
    """Errors normalised by the reported sigma have unit scatter (honest uncertainties)."""
    pulls_l, pulls_m = [], []
    for seed in range(60):
        sel, _ = make_samples(layout, 200.0, 50.0, amp=0.6, noise=1.0, seed=seed)  # SNR ~ 7
        est = estimate_direction(sel.samples, 200.0, 50.0, options=EstimatorOptions(search_radius_deg=6.0,
                                                                                    robust_loss="linear"))
        cov = estimate_direction_uncertainty(est, sel.samples)
        true = azel_to_unit(200.0, 50.0)
        pulls_l.append((est.l - true[0]) / cov.sigma_l)
        pulls_m.append((est.m - true[1]) / cov.sigma_m)
    for pulls in (pulls_l, pulls_m):
        assert 0.7 < np.std(pulls) < 1.3
        assert abs(np.mean(pulls)) < 0.45


def test_huber_uncertainties_reasonable_with_gaussian_noise(layout):
    pulls = []
    for seed in range(40):
        sel, _ = make_samples(layout, 120.0, 60.0, amp=0.6, noise=1.0, seed=100 + seed)
        est = estimate_direction(sel.samples, 120.0, 60.0, options=EstimatorOptions(search_radius_deg=6.0))
        cov = estimate_direction_uncertainty(est, sel.samples)
        pulls.append((est.l - azel_to_unit(120.0, 60.0)[0]) / cov.sigma_l)
    assert 0.7 < np.std(pulls) < 1.4


def test_recovery_after_random_flagging(layout):
    sel, _ = make_samples(layout, 140.0, 35.0, noise=0.05, flag_fraction=0.4, seed=4)
    assert sel.flagged_fraction > 0.3
    est = estimate_direction(sel.samples, 145.0, 33.0)
    cov = estimate_direction_uncertainty(est, sel.samples)
    assert _err(est, 140.0, 35.0) < 5 * max(cov.sigma_l, cov.sigma_m) * 180 / np.pi / np.sin(np.radians(35)) + 0.05


def test_robust_loss_beats_least_squares_with_bad_baselines(layout):
    errs = {"linear": [], "huber": [], "soft_l1": []}
    for seed in range(12):
        sel, _ = make_samples(layout, 230.0, 45.0, amp=1.0, noise=0.1, seed=seed)
        s = sel.samples
        rng = np.random.default_rng(1000 + seed)
        bad = rng.choice(s.n, size=s.n // 10, replace=False)
        s.vis[bad] += 3.0 * np.exp(2j * np.pi * rng.uniform(size=bad.size))
        for loss in errs:
            est = estimate_direction(s, 230.0, 45.0, options=EstimatorOptions(robust_loss=loss))
            errs[loss].append(_err(est, 230.0, 45.0))
    assert np.median(errs["huber"]) < 0.5 * np.median(errs["linear"])
    assert np.median(errs["soft_l1"]) < 0.5 * np.median(errs["linear"])


def test_source_near_horizon_gives_valid_lm(layout):
    sel, _ = make_samples(layout, 75.0, 3.0, noise=0.02, seed=5)
    est = estimate_direction(sel.samples, 75.0, 6.0, options=EstimatorOptions(search_radius_deg=10.0))
    assert est.l ** 2 + est.m ** 2 <= 1.0
    assert est.el_deg >= 0.0
    assert _err(est, 75.0, 3.0) < 1.0


def test_zenith_source_flags_azimuth(layout):
    sel, _ = make_samples(layout, 0.0, 90.0, amp=0.3, noise=1.0, seed=6)
    est = estimate_direction(sel.samples, 10.0, 85.0)
    cov = estimate_direction_uncertainty(est, sel.samples)
    assert cov.zenith_az_uncertain
    assert np.isnan(cov.sigma_az_deg) and np.isfinite(cov.sigma_el_deg)
    assert np.isfinite(cov.sigma_l) and np.isfinite(cov.sigma_m)


def test_tle_offset_inside_radius_recovers_injected_position(layout):
    """MANDATORY: the fit follows the visibilities, not the TLE/search centre."""
    tle_az, tle_el = 210.0, 44.0
    true_az, true_el = 216.0, 48.5  # ~6 deg from the "TLE" direction
    offset = float(angular_separation_deg(azel_to_unit(true_az, true_el), azel_to_unit(tle_az, tle_el)))
    assert 5.0 < offset < 8.0
    sel, _ = make_samples(layout, true_az, true_el, noise=0.3, seed=7)
    est = estimate_direction(sel.samples, tle_az, tle_el, search_radius_deg=10.0)
    cov = estimate_direction_uncertainty(est, sel.samples)
    assert _err(est, true_az, true_el) < 0.2
    assert est.separation_from_centre_deg > offset - 0.2
    assert not est.boundary
    assert cov.sigma_az_deg < 1.0


def test_tle_offset_outside_radius_is_flagged_not_replaced(layout):
    tle_az, tle_el = 210.0, 44.0
    sel, _ = make_samples(layout, 230.0, 50.0, noise=0.05, seed=8)  # ~15 deg away
    est = estimate_direction(sel.samples, tle_az, tle_el, search_radius_deg=10.0)
    assert est.boundary
    assert est.separation_from_centre_deg > 5.0  # not collapsed to the TLE position
    from tart_gnss.stage3 import quality as q

    flags = q.assess(thresholds=q.QualityThresholds(), n_baselines_used=276, flagged_fraction=0.0, fitted=True,
                     fit_success=est.success, boundary=est.boundary, snr=est.snr,
                     residual_tail_ratio=est.residual_tail_ratio)
    assert BOUNDARY_SOLUTION in flags and not q.is_usable(flags)


def test_visibility_and_image_estimators_agree(layout):
    sel, _ = make_samples(layout, 160.0, 55.0, noise=0.2, seed=9)
    region = SearchRegion.around_azel(158.0, 53.0, 10.0)
    est = estimate_direction(sel.samples, 158.0, 53.0)
    img = image_direction(sel.samples, region)
    sep = float(angular_separation_deg(est.unit_vector, azel_to_unit(img.az_deg, img.el_deg)))
    assert sep < 0.5 and img.image_snr > 5


def test_multichannel_recovery(layout):
    sel, _ = make_samples(layout, 300.0, 40.0, noise=0.1, seed=10, frequency_hz=(1.570e9, 1.5754e9, 1.581e9))
    assert np.unique(sel.samples.frequency_hz).size == 3
    est = estimate_direction(sel.samples, 302.0, 41.0)
    assert _err(est, 300.0, 40.0) < 0.1


def test_brighter_unrelated_source_outside_cap_does_not_flag_genuine_target(layout):
    """Crowded-sky case: target inside the cap, a brighter satellite just outside it."""
    centre = (210.0, 44.0)
    sel, _ = make_samples(layout, 211.0, 45.0, amp=0.5, noise=0.05, seed=11,
                          extra_sources=[(228.0, 50.0, 1.5 + 0.0j)])  # ~13 deg from the centre
    est = estimate_direction(sel.samples, *centre, search_radius_deg=10.0)
    assert not est.boundary
    assert _err(est, 211.0, 45.0) < 1.0  # biased somewhat by the interferer's sidelobes, but on target
