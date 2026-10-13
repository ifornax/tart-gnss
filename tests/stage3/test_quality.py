"""Deterministic quality flags."""

from tart_gnss.stage3 import quality as q

T = q.QualityThresholds(min_baselines=10, max_flagged_fraction=0.5, min_snr=5.0)


def _assess(**kw):
    base = dict(thresholds=T, n_baselines_used=200, flagged_fraction=0.1, fitted=True, fit_success=True,
                snr=20.0, residual_tail_ratio=1.0, reduced_chi2=1.0)
    base.update(kw)
    return q.assess(**base)


def test_good_has_no_flags():
    assert _assess() == [] and q.is_usable([])


def test_each_flag_triggers():
    assert _assess(n_baselines_used=5, fitted=False) == [q.LOW_BASELINE_COUNT]
    assert q.HIGH_FLAGGED_FRACTION in _assess(flagged_fraction=0.7)
    assert q.FIT_FAILED in _assess(fit_success=False)
    assert q.FIT_FAILED in _assess(covariance_finite=False)
    assert q.BOUNDARY_SOLUTION in _assess(boundary=True)
    assert q.LOW_SNR in _assess(snr=3.0)
    assert q.HIGH_RESIDUAL in _assess(residual_tail_ratio=2.0)
    assert q.ZENITH_AZ_UNCERTAIN in _assess(zenith_az_uncertain=True)
    assert q.INVALID_GEOMETRY in _assess(invalid_geometry=True)
    assert q.CALIBRATION_WARNING in _assess(calibration_warning=True)
    assert q.TLE_UNASSOCIATED in _assess(tle_separation_deg=7.0)
    assert q.TLE_UNASSOCIATED not in _assess(tle_separation_deg=3.0)
    off = q.QualityThresholds(max_tle_separation_deg=None)
    assert q.TLE_UNASSOCIATED not in q.assess(thresholds=off, n_baselines_used=200, flagged_fraction=0.0,
                                              fitted=True, fit_success=True, snr=20.0, tle_separation_deg=40.0)


def test_reduced_chi2_only_used_with_calibrated_weights():
    assert q.HIGH_RESIDUAL not in _assess(reduced_chi2=50.0)
    cal = q.QualityThresholds(weights_calibrated=True)
    assert q.HIGH_RESIDUAL in q.assess(thresholds=cal, n_baselines_used=200, flagged_fraction=0.0, fitted=True,
                                       fit_success=True, snr=20.0, residual_tail_ratio=1.0, reduced_chi2=50.0)


def test_image_disagreement_threshold_uses_beam_and_sigma():
    assert q.IMAGE_DISAGREEMENT not in _assess(vis_image_separation_deg=1.0, resolution_deg=3.0)
    assert q.IMAGE_DISAGREEMENT in _assess(vis_image_separation_deg=2.0, resolution_deg=3.0, sigma_position_deg=0.2)
    assert q.IMAGE_DISAGREEMENT not in _assess(vis_image_separation_deg=2.0, resolution_deg=3.0,
                                               sigma_position_deg=1.0)


def test_multiple_flags_and_blocking():
    flags = _assess(boundary=True, snr=1.0, calibration_warning=True)
    assert flags == [q.BOUNDARY_SOLUTION, q.LOW_SNR, q.CALIBRATION_WARNING]
    assert not q.is_usable(flags)
    assert q.is_usable([q.CALIBRATION_WARNING, q.HIGH_RESIDUAL, q.IMAGE_DISAGREEMENT, q.HIGH_FLAGGED_FRACTION])
    assert not q.is_usable([q.LOW_SNR])
    assert not q.is_usable([q.TLE_UNASSOCIATED])
