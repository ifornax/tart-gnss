"""Extraction of continuous visibility windows from sampled boolean series.

The coarse sampling only *brackets* each FALSE->TRUE or TRUE->FALSE
transition; boundaries are then refined by bisection on the exact visibility
predicate. This module is independent of orbits and stations: it works on
elapsed-second offsets and a boolean predicate, which makes it testable with
synthetic data.

Limitation: excursions shorter than the coarse step (a satellite dipping
below/above the threshold and back between two samples) cannot be detected.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

DEFAULT_TOLERANCE_S = 0.1


@dataclass(frozen=True)
class Window:
    """A maximal continuous interval during which the predicate is True.

    Attributes
    ----------
    start_s, end_s
        Boundaries in elapsed seconds from the search start.
    start_clipped, end_clipped
        True if the boundary coincides with the search-interval edge, i.e.
        the real window may extend beyond it.
    """

    start_s: float
    end_s: float
    start_clipped: bool = False
    end_clipped: bool = False

    @property
    def duration_s(self) -> float:
        """Window length in seconds."""
        return self.end_s - self.start_s

    @property
    def mid_s(self) -> float:
        """Window midpoint in elapsed seconds."""
        return 0.5 * (self.start_s + self.end_s)


def find_true_runs(mask: Sequence[bool]) -> List[Tuple[int, int]]:
    """Locate contiguous runs of True samples.

    Parameters
    ----------
    mask
        One-dimensional boolean sequence.

    Returns
    -------
    list of (first_index, last_index)
        Inclusive sample-index ranges of each run, in time order.
    """
    m = np.asarray(mask, dtype=bool).astype(np.int8)
    if m.size == 0:
        return []
    edges = np.diff(np.concatenate(([0], m, [0])))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1
    return [(int(a), int(b)) for a, b in zip(starts, ends)]


def refine_transition(
    predicate: Callable[[float], bool],
    t_false: float,
    t_true: float,
    tolerance_s: float = DEFAULT_TOLERANCE_S,
    max_iter: int = 100,
) -> float:
    """Bisect a single predicate transition between two bracketing epochs.

    Parameters
    ----------
    predicate
        Function of elapsed seconds returning visibility.
    t_false
        Epoch at which the predicate is False.
    t_true
        Epoch at which the predicate is True (may be before or after
        ``t_false``).
    tolerance_s
        Stop when the bracket is narrower than this (seconds).

    Returns
    -------
    float
        The True-side end of the final bracket: an epoch at which the
        predicate holds, within ``tolerance_s`` of the transition and always
        inside the original bracket.
    """
    if tolerance_s <= 0:
        raise ValueError("tolerance_s must be positive")
    lo, hi = float(t_false), float(t_true)
    for _ in range(max_iter):
        if abs(hi - lo) <= tolerance_s:
            break
        mid = 0.5 * (lo + hi)
        if predicate(mid):
            hi = mid
        else:
            lo = mid
    return hi


def extract_windows(
    offsets_s: Sequence[float],
    mask: Sequence[bool],
    predicate: Optional[Callable[[float], bool]] = None,
    tolerance_s: float = DEFAULT_TOLERANCE_S,
) -> List[Window]:
    """Build refined continuous windows from a coarse boolean time series.

    Parameters
    ----------
    offsets_s
        Strictly increasing sample epochs (elapsed seconds).
    mask
        Predicate values at ``offsets_s``.
    predicate
        Exact predicate for boundary refinement. If ``None`` the coarse
        sample epochs of the first/last True samples are used.
    tolerance_s
        Boundary refinement tolerance, seconds.

    Returns
    -------
    list of Window
        One entry per contiguous True run; separate runs are never merged.
    """
    t = np.asarray(offsets_s, dtype=float)
    m = np.asarray(mask, dtype=bool)
    if t.shape != m.shape:
        raise ValueError("offsets_s and mask must have the same shape")
    if t.size > 1 and np.any(np.diff(t) <= 0):
        raise ValueError("offsets_s must be strictly increasing")
    out: List[Window] = []
    last = t.size - 1
    for i0, i1 in find_true_runs(m):
        start_clipped = i0 == 0
        end_clipped = i1 == last
        start = t[i0]
        end = t[i1]
        if predicate is not None:
            if not start_clipped:
                start = refine_transition(predicate, t[i0 - 1], t[i0], tolerance_s)
            if not end_clipped:
                end = refine_transition(predicate, t[i1 + 1], t[i1], tolerance_s)
        out.append(Window(float(start), float(end), start_clipped, end_clipped))
    return out
