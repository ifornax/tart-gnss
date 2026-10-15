"""Population statistics with runs as the independent units.

Epochs within one pass are strongly correlated, so every inference here is
made on per-run aggregates (one value per run) or with a run-level
(cluster) bootstrap; epochs are never treated as independent samples for a
p-value. Group comparisons report effect sizes next to p-values and are only
computed when every group has at least ``min_runs_per_group`` runs.

P-values are permutation p-values, never large-sample approximations: with
a handful of runs (and tied per-run values) asymptotic Mann-Whitney,
Kruskal-Wallis or Spearman p-values can be orders of magnitude too small.
The label (or pairing) permutations are enumerated exactly when there are
at most :data:`EXACT_LIMIT` of them, otherwise :data:`MC_PERMUTATIONS`
random permutations with a fixed seed are used; ``p_method`` says which.
"""

from __future__ import annotations

from typing import Dict, Sequence

import numpy as np

MAD_TO_SIGMA = 1.4826
SAMPLE_TIERS = ((3, "descriptive_only"), (10, "preliminary"), (30, "basic_inference"))
EXACT_LIMIT = 100_000
MC_PERMUTATIONS = 10_000
PERMUTATION_SEED = 12345


def _permutation_p(values: np.ndarray, sizes, statistic) -> tuple:
    """Two-sided-by-construction permutation p of ``statistic(groups)`` (larger = more extreme)."""
    from itertools import combinations
    from math import comb

    n = len(values)
    obs = statistic(np.split(values, np.cumsum(sizes)[:-1]))
    if len(sizes) == 2 and comb(n, sizes[0]) <= EXACT_LIMIT:
        idx = np.arange(n)
        count = total = 0
        for c in combinations(range(n), sizes[0]):
            mask = np.zeros(n, bool)
            mask[list(c)] = True
            total += 1
            count += statistic([values[mask], values[idx[~mask]]]) >= obs - 1e-12
        return count / total, "exact permutation"
    rng = np.random.default_rng(PERMUTATION_SEED)
    count = 0
    for _ in range(MC_PERMUTATIONS):
        perm = rng.permutation(values)
        count += statistic(np.split(perm, np.cumsum(sizes)[:-1])) >= obs - 1e-12
    return (count + 1) / (MC_PERMUTATIONS + 1), f"Monte Carlo permutation ({MC_PERMUTATIONS})"


def describe(values) -> Dict[str, float]:
    """n, median, mean, RMS, MAD-based sigma, 16/84/95th percentiles, min and max of the finite values."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    nan = float("nan")
    if not v.size:
        return {"n": 0, "median": nan, "mean": nan, "rms": nan, "mad": nan, "p16": nan, "p84": nan, "p95": nan,
                "min": nan, "max": nan}
    med = float(np.median(v))
    return {"n": int(v.size), "median": med, "mean": float(v.mean()), "rms": float(np.sqrt(np.mean(v * v))),
            "mad": float(MAD_TO_SIGMA * np.median(np.abs(v - med))), "p16": float(np.percentile(v, 16)),
            "p84": float(np.percentile(v, 84)), "p95": float(np.percentile(v, 95)), "min": float(v.min()),
            "max": float(v.max())}


def sample_tier(n_runs: int) -> str:
    """Transparent sample-size guidance: descriptive_only < 3 <= preliminary < 10 <= basic_inference < 30 <= grouped_comparisons."""
    for limit, name in SAMPLE_TIERS:
        if n_runs < limit:
            return name
    return "grouped_comparisons"


def bin_edges_label(lo: float, hi: float, unit: str = "") -> str:
    """Readable bin label."""
    return f"{lo:g}-{hi:g}{unit}" if np.isfinite(hi) else f">{lo:g}{unit}"


def binned(df, x: str, ys: Sequence[str], edges: Sequence[float], unit: str = "", cluster: str = "run_id"):
    """Per-bin n (rows and runs) and descriptive statistics of each ``y`` (x in [lo, hi))."""
    import pandas as pd

    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = df[(df[x] >= lo) & (df[x] < hi)] if len(df) else df
        rec = {"bin": bin_edges_label(lo, hi, unit), "bin_low": lo, "bin_high": hi, "n_rows": int(len(sel)),
               "n_runs": int(sel[cluster].nunique()) if len(sel) and cluster in sel else 0}
        for y in ys:
            d = describe(sel[y]) if len(sel) and y in sel else describe([])
            rec.update({f"{y}_median": d["median"], f"{y}_rms": d["rms"], f"{y}_mad": d["mad"], f"{y}_p16": d["p16"],
                        f"{y}_p84": d["p84"]})
        rows.append(rec)
    return pd.DataFrame(rows)


def compare_groups(groups: Dict[str, Sequence[float]], min_runs_per_group: int = 3) -> Dict[str, object]:
    """Mann-Whitney U (2 groups) or Kruskal-Wallis (> 2) on per-run values, with effect sizes.

    Effect sizes: rank-biserial correlation (2 groups) or epsilon-squared
    (Kruskal-Wallis); plus the difference in medians. Groups smaller than
    ``min_runs_per_group`` make the test unavailable (``under_sampled``).
    """
    from scipy import stats

    clean = {k: np.asarray([x for x in v if np.isfinite(x)], float) for k, v in groups.items()}
    sizes = {k: int(len(v)) for k, v in clean.items()}
    out: Dict[str, object] = {"group_sizes": sizes, "test": None, "statistic": None, "p_value": None,
                              "effect_size": None, "effect_size_name": None,
                              "medians": {k: float(np.median(v)) if len(v) else None for k, v in clean.items()}}
    ok = {k: v for k, v in clean.items() if len(v) >= min_runs_per_group}
    out["under_sampled"] = sorted(k for k in clean if k not in ok)
    if len(ok) < 2 or np.ptp(np.concatenate(list(ok.values()))) == 0:
        out["test"] = "unavailable (fewer than two adequately sampled groups)" if len(ok) < 2 else "no variation"
        return out
    keys = sorted(ok)
    values = np.concatenate([ok[k] for k in keys])
    sizes = [len(ok[k]) for k in keys]
    if len(keys) == 2:
        a, b = ok[keys[0]], ok[keys[1]]
        u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
        centre = len(a) * len(b) / 2.0

        def stat(gs):
            return abs(stats.mannwhitneyu(gs[0], gs[1], alternative="two-sided")[0] - centre)

        p, how = _permutation_p(values, sizes, stat)
        out.update({"test": "Mann-Whitney U", "statistic": float(u), "p_value": float(p), "p_method": how,
                    "effect_size": float(1 - 2 * u / (len(a) * len(b))),
                    "effect_size_name": f"rank-biserial, positive when {keys[1]} > {keys[0]}",
                    "median_difference": float(np.median(b) - np.median(a))})
    else:
        h, _ = stats.kruskal(*[ok[k] for k in keys])
        p, how = _permutation_p(values, sizes, lambda gs: stats.kruskal(*gs)[0])
        n = sum(sizes)
        out.update({"test": "Kruskal-Wallis", "statistic": float(h), "p_value": float(p), "p_method": how,
                    "effect_size": float(max(h - len(keys) + 1, 0) / (n - len(keys))) if n > len(keys) else None,
                    "effect_size_name": "epsilon-squared"})
    return out


def dependence(x, y) -> Dict[str, object]:
    """Spearman correlation and Theil-Sen slope of per-run values (runs are the independent units)."""
    from scipy import stats

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    out = {"n_runs": int(ok.sum()), "spearman_r": None, "spearman_p": None, "p_method": None,
           "theil_sen_slope": None, "theil_sen_slope_low": None, "theil_sen_slope_high": None,
           "theil_sen_intercept": None, "note": None}
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return out
    r = float(stats.spearmanr(x, y)[0])
    p, how = spearman_permutation_p(x, y, r)
    slope, intercept, lo, hi = stats.theilslopes(y, x)
    out.update({"spearman_r": r, "spearman_p": p, "p_method": how, "theil_sen_slope": float(slope),
                "theil_sen_slope_low": float(lo), "theil_sen_slope_high": float(hi),
                "theil_sen_intercept": float(intercept),
                "note": "fewer than 5 runs: the slope interval is not meaningful" if len(x) < 5 else None})
    return out


def spearman_permutation_p(x, y, r_obs: float) -> tuple:
    """Two-sided permutation p of a Spearman coefficient (exact for n <= 8)."""
    from itertools import permutations

    from scipy import stats

    n = len(x)
    if n <= 8:
        rs = [stats.spearmanr(x, np.asarray(y)[list(perm)])[0] for perm in permutations(range(n))]
        rs = np.asarray(rs, float)
        return float(np.mean(np.abs(rs) >= abs(r_obs) - 1e-12)), "exact permutation"
    rng = np.random.default_rng(PERMUTATION_SEED)
    count = sum(abs(stats.spearmanr(x, rng.permutation(y))[0]) >= abs(r_obs) - 1e-12 for _ in range(MC_PERMUTATIONS))
    return float((count + 1) / (MC_PERMUTATIONS + 1)), f"Monte Carlo permutation ({MC_PERMUTATIONS})"


def wilson_interval(k: int, n: int, z: float = 1.96):
    """Wilson score interval of a proportion ``k / n`` (None for n = 0)."""
    if n <= 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(max(c - h, 0.0)), float(min(c + h, 1.0))
