"""Run-level (cluster) bootstrap confidence intervals.

Runs are resampled with replacement; the statistic is recomputed on the
resampled runs (and, for epoch- or station-level statistics, on all rows of
those runs). Percentile intervals at ``level`` (default 95 %). A fixed
``seed`` makes the intervals reproducible.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Optional

import numpy as np

from tart_gnss.campaign.errors import BootstrapError


@dataclass
class BootstrapResult:
    """Point estimate and percentile confidence interval."""

    estimate: float
    ci_low: float
    ci_high: float
    n_runs: int
    n_bootstrap: int
    level: float
    resampling_unit: str = "run"

    def to_dict(self):
        """JSON-serialisable form."""
        return asdict(self)


def bootstrap_run_level(run_table, statistic: Callable, n_bootstrap: int = 1000, seed: int = 20260625,
                        level: float = 0.95, rows=None, cluster: str = "run_id") -> BootstrapResult:
    """Bootstrap campaign statistics using runs as the resampling unit.

    ``statistic`` receives the resampled run table (or, if ``rows`` is given,
    the rows of ``rows`` belonging to the resampled runs, with duplicated runs
    repeated) and returns a float.
    """
    if n_bootstrap < 0:
        raise BootstrapError("n_bootstrap must be >= 0")
    import pandas as pd

    runs = list(run_table[cluster])
    n = len(runs)
    target = run_table if rows is None else rows
    est = float(statistic(target))
    nan = float("nan")
    if n < 2 or n_bootstrap == 0 or not np.isfinite(est):
        return BootstrapResult(est, nan, nan, n, n_bootstrap, level)
    rng = np.random.default_rng(seed)
    groups = {r: g for r, g in target.groupby(cluster)} if rows is not None else None
    vals = []
    for _ in range(n_bootstrap):
        pick = rng.integers(0, n, n)
        if rows is None:
            sample = run_table.iloc[pick]
        else:
            parts = [groups[runs[i]] for i in pick if runs[i] in groups]
            sample = pd.concat(parts) if parts else target.iloc[:0]
        v = statistic(sample)
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < max(10, n_bootstrap // 2):
        raise BootstrapError("too many undefined bootstrap replicates")
    a = (1 - level) / 2
    return BootstrapResult(est, float(np.quantile(vals, a)), float(np.quantile(vals, 1 - a)), n, n_bootstrap, level)


def median_of(column: str) -> Callable:
    """Statistic: median of a column's finite values (NaN if none)."""
    def f(df) -> float:
        v = np.asarray(df[column], float) if len(df) else np.zeros(0)
        v = v[np.isfinite(v)]
        return float(np.median(v)) if v.size else float("nan")

    return f


def fraction_of(predicate: Callable) -> Callable:
    """Statistic: fraction of rows for which ``predicate(row)`` is true."""
    def f(df) -> float:
        return float(np.mean([bool(predicate(r)) for _, r in df.iterrows()])) if len(df) else float("nan")

    return f


def optional(result: Optional[BootstrapResult]):
    """``to_dict`` or None."""
    return result.to_dict() if result is not None else None
