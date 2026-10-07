from datetime import timedelta

import numpy as np
import pytest

from stage2_helpers import make_tart_hdf5, utc
from tart_gnss.stage2.errors import NoCommonTimeError
from tart_gnss.stage2.native import inspect_hdf5, summarise_coverage
from tart_gnss.stage2.synchronise import (common_interval, coverage_fraction, intersect_segments,
                                          nearest_separation, synchronise)
from tart_gnss.stations import Station


def hm(h, m):
    return utc(2026, 6, 25, h, m)


def test_spec_temporal_intersection():
    intervals = {"A": (hm(10, 0), hm(10, 20)), "B": (hm(10, 3), hm(10, 18)), "C": (hm(10, 1), hm(10, 25))}
    assert common_interval(intervals) == (hm(10, 3), hm(10, 18))


def test_no_intersection_fails():
    intervals = {"A": (hm(10, 0), hm(10, 5)), "B": (hm(10, 6), hm(10, 10)), "C": (hm(10, 0), hm(10, 30))}
    with pytest.raises(NoCommonTimeError):
        common_interval(intervals)


def test_touching_intervals_are_not_an_overlap():
    with pytest.raises(NoCommonTimeError):
        common_interval({"A": (hm(10, 0), hm(10, 5)), "B": (hm(10, 5), hm(10, 9))})


def test_coverage_fraction_and_segments():
    window = (hm(10, 0), hm(10, 20))
    segs = [(hm(9, 50), hm(10, 5)), (hm(10, 10), hm(10, 15))]
    assert coverage_fraction(segs, window) == pytest.approx(10 / 20)
    assert intersect_segments(segs, [(hm(10, 3), hm(10, 12))]) == [(hm(10, 3), hm(10, 5)), (hm(10, 10), hm(10, 12))]


def test_nearest_separation():
    a = [hm(10, 0) + timedelta(seconds=s) for s in (0, 1, 2)]
    b = [hm(10, 0) + timedelta(seconds=s) for s in (0.25, 1.25, 2.25)]
    assert np.allclose(nearest_separation(a, b), 0.25)


def _cov(tmp_path, sid, offset_s, n_files=3):
    st = Station(sid, sid, *{"alpha": (-20.0, 57.0, 0.0), "beta": (-19.0, 47.0, 1000.0),
                            "gamma": (-33.0, 26.0, 500.0)}[sid], "test")
    obs = [inspect_hdf5(make_tart_hdf5(tmp_path / f"{sid}{k}.hdf", sid,
                                       hm(10, 1) + timedelta(seconds=60 * k + offset_s)), st)
           for k in range(n_files)]
    return summarise_coverage(sid, obs)


def test_synchronise_three_stations(tmp_path):
    covs = {"alpha": _cov(tmp_path, "alpha", 0.0), "beta": _cov(tmp_path, "beta", 0.3),
            "gamma": _cov(tmp_path, "gamma", 10.0)}
    s = synchronise(covs, (hm(10, 0), hm(10, 4)))
    assert s["actual_common_start_utc"] == "2026-06-25T10:00:11.000000Z"
    assert s["actual_common_end_utc"] == "2026-06-25T10:03:00.000000Z"  # alpha ends first
    assert 0 < s["three_station_coverage_fraction"] < 1
    pair = s["station_pair_time_offset_statistics"]["alpha__beta"]
    assert pair["median_offset_s"] == pytest.approx(0.3, abs=1e-3)
    assert s["max_nearest_timestamp_separation_s"] <= 0.5


def test_synchronise_without_overlap(tmp_path):
    covs = {"alpha": _cov(tmp_path, "alpha", 0.0, 1), "beta": _cov(tmp_path, "beta", 600.0, 1)}
    s = synchronise(covs, (hm(10, 0), hm(10, 20)))
    assert s["actual_common_start_utc"] is None and s["no_common_time_reason"]
