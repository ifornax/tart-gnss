import numpy as np
import pytest

from tart_gnss.windows import extract_windows, find_true_runs, refine_transition

F, T = False, True


def test_single_window():
    mask = [F, F, T, T, T, F]
    assert find_true_runs(mask) == [(2, 4)]
    windows = extract_windows(np.arange(6) * 10.0, mask)
    assert len(windows) == 1


def test_multiple_windows_kept_separate():
    mask = [F, T, T, F, F, T, T, F]
    windows = extract_windows(np.arange(8) * 10.0, mask)
    assert len(windows) == 2
    assert windows[0].end_s < windows[1].start_s


def test_no_windows():
    assert extract_windows(np.arange(4.0), [F, F, F, F]) == []


def test_clipped_flags():
    windows = extract_windows(np.arange(5.0), [T, T, F, T, T])
    assert windows[0].start_clipped and not windows[0].end_clipped
    assert windows[1].end_clipped and not windows[1].start_clipped


@pytest.mark.parametrize("crossing", [19 * 3600 + 2 * 60 + 21.37, 68545.0001, 68549.999])
def test_refined_crossing_within_bracket(crossing):
    t_false, t_true = 68540.0, 68550.0  # e.g. 19:02:20 FALSE, 19:02:30 TRUE

    def predicate(t):
        return t >= crossing

    refined = refine_transition(predicate, t_false, t_true, tolerance_s=0.01)
    assert t_false <= refined <= t_true
    assert abs(refined - crossing) <= 0.01
    assert predicate(refined)


def test_refined_window_boundaries():
    rise, set_ = 23.4, 71.8

    def predicate(t):
        return rise <= t <= set_

    offsets = np.arange(0.0, 101.0, 10.0)
    mask = [predicate(t) for t in offsets]
    (w,) = extract_windows(offsets, mask, predicate, tolerance_s=0.05)
    assert 20.0 <= w.start_s <= 30.0 and abs(w.start_s - rise) <= 0.05
    assert 70.0 <= w.end_s <= 80.0 and abs(w.end_s - set_) <= 0.05
    assert w.duration_s == pytest.approx(set_ - rise, abs=0.1)


def test_shape_mismatch():
    with pytest.raises(ValueError):
        extract_windows([0.0, 1.0], [True])
