import numpy as np

from tart_gnss.visibility import common_visibility, station_visibility


def test_spec_example():
    el = [
        [19, 21, 25, 18],  # station A
        [22, 23, 24, 25],  # station B
        [30, 30, 19, 30],  # station C
    ]
    assert common_visibility(el, 20.0).tolist() == [False, True, False, False]


def test_threshold_is_inclusive():
    assert common_visibility([[20.0], [20.0], [20.0]], 20.0).tolist() == [True]


def test_nan_is_not_visible():
    assert station_visibility([[np.nan, 25.0]], 20.0).tolist() == [[False, True]]


def test_generic_number_of_stations():
    el = np.full((5, 3), 30.0)
    el[4, 1] = 10.0
    assert common_visibility(el, 20.0).tolist() == [True, False, True]
    assert common_visibility([[25.0, 15.0]], 20.0).tolist() == [True, False]
