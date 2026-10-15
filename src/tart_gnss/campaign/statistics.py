"""Population statistics with runs as the independent units.

Epochs within one pass are strongly correlated, so every inference here is
made on per-run aggregates (one value per run) or with a run-level
(cluster) bootstrap; epochs are never treated as independent samples for a
p-value. Group comparisons report effect sizes next to p-values and are only
computed when every group has at least ``min_runs_per_group`` runs.
"""

from __future__ import annotations

from typing import Dict, Sequence

import numpy as np

MAD_TO_SIGMA = 1.4826
SAMPLE_TIERS = ((3, "descriptive_only"), (10, "preliminary"), (30, "basic_inference"))


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
    if len(keys) == 2:
        a, b = ok[keys[0]], ok[keys[1]]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        out.update({"test": "Mann-Whitney U", "statistic": float(u), "p_value": float(p),
                    "effect_size": float(1 - 2 * u / (len(a) * len(b))),
                    "effect_size_name": f"rank-biserial, positive when {keys[1]} > {keys[0]}",
                    "median_difference": float(np.median(b) - np.median(a))})
    else:
        h, p = stats.kruskal(*[ok[k] for k in keys])
        n = sum(len(ok[k]) for k in keys)
        out.update({"test": "Kruskal-Wallis", "statistic": float(h), "p_value": float(p),
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
    out = {"n_runs": int(ok.sum()), "spearman_r": None, "spearman_p": None, "theil_sen_slope": None,
           "theil_sen_slope_low": None, "theil_sen_slope_high": None, "theil_sen_intercept": None}
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return out
    r, p = stats.spearmanr(x, y)
    slope, intercept, lo, hi = stats.theilslopes(y, x)
    out.update({"spearman_r": float(r), "spearman_p": float(p), "theil_sen_slope": float(slope),
                "theil_sen_slope_low": float(lo), "theil_sen_slope_high": float(hi),
                "theil_sen_intercept": float(intercept)})
    return out


def wilson_interval(k: int, n: int, z: float = 1.96):
    """Wilson score interval of a proportion ``k / n`` (None for n = 0)."""
    if n <= 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(max(c - h, 0.0)), float(min(c + h, 1.0))
