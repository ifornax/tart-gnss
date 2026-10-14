"""Diagnostic correlation analyses (calibration, localisation, timing, geometry, trajectory model).

Each test pairs an explanatory variable ``x`` from one *family* with a
residual measure ``y``, within one station (``scope = station``) or across
epochs (``scope = network``). A test *supports* its family's tag only if all
of these hold:

* ``n_samples >= min_samples``;
* ``x`` varies (at least 3 distinct values, or a two-valued indicator with at
  least 3 samples in each group, and, for geometry variables, a
  range of at least ``MIN_GEOMETRY_VARIATION`` of the median; otherwise
  ``INSUFFICIENT_VARIATION``);
* ``|partial_coefficient| >= min_correlation`` (effect size, after removing
  the common time trend);
* the Holm-adjusted p-value (over every test of the same family in the run,
  with the effective sample size) is below ``alpha``. Each family answers a
  separate question, so the multiple-testing correction is applied within
  it; with about ten tests per family at ``alpha = 0.01`` a single test needs
  p below about 1e-3. With dozens of tests per run some correlations will
  still look strong by chance, which is why effect size and sample count
  are required as well and why a supported correlation is evidence, not a
  conclusion.

P-values alone never promote a diagnosis. Expected relationships are
labelled ``EXPECTED``, never tags: poorer geometry -> larger Stage-4 sigma.
Geometry variables are computed from the *predicted* sight lines (stations to
the Stage-5 position), not the measured ones, whose noise would couple
spuriously to the residuals.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from tart_gnss.stage7.statistics import correlation, holm

DEFAULT_MIN_SAMPLES = 8
DEFAULT_MIN_CORRELATION = 0.5
DEFAULT_ALPHA = 0.01
MIN_GEOMETRY_VARIATION = 0.05

FAMILY_TAG = {
    "calibration": "CALIBRATION_CORRELATED",
    "localisation": "LOCALISATION_CORRELATED",
    "timing": "TIMING_CORRELATED",
    "geometry": "GEOMETRY_CORRELATED",
    "trajectory": "TRAJECTORY_MODEL_CORRELATED",
}

#: (family, x variable) per station row.
STATION_X: List[Tuple[str, str]] = [
    ("calibration", "zero_gain_count"), ("calibration", "gain_solution_index"), ("calibration", "flagged_fraction"),
    ("calibration", "n_baselines_used"),
    ("localisation", "snr"), ("localisation", "reduced_chi2"), ("localisation", "vis_image_separation_deg"),
    ("timing", "station_time_offset_s"), ("timing", "abs_station_time_offset_s"),
]
STATION_Y = ["direction_residual_angle_deg", "direction_residual_along_motion_deg", "line_residual_m",
             "loo_shift_3d_m"]
#: (family, x variable) per epoch.
EPOCH_X: List[Tuple[str, str]] = [
    ("geometry", "predicted_condition_number"), ("geometry", "predicted_inverse_min_crossing"),
    ("trajectory", "stage5_nis"),
    ("timing", "max_abs_station_time_offset_s"), ("timing", "mean_station_time_offset_s"),
]
EPOCH_Y = ["stage4_residual_3d_m", "position_residual_3d_m", "abs_dr_along_track_m", "max_loo_shift_m"]
EXPECTED = {("geometry", "stage4_sigma_3d_m")}


def compute_correlations(station_tbl, epoch_tbl, method: str = "spearman", min_samples: int = DEFAULT_MIN_SAMPLES,
                         min_correlation: float = DEFAULT_MIN_CORRELATION, alpha: float = DEFAULT_ALPHA):
    """Calculate diagnostic correlations with sample counts and significance metadata."""
    import pandas as pd

    tests: List[Dict[str, object]] = []
    if station_tbl is not None and len(station_tbl):
        for sid, g in station_tbl.groupby("station_id"):
            for fam, xv in STATION_X:
                if xv not in g.columns:
                    continue
                for yv in STATION_Y:
                    if yv not in g.columns:
                        continue
                    c = correlation(g[xv].to_numpy(float), g[yv].to_numpy(float), g["time_s"].to_numpy(float),
                                    method)
                    tests.append({"scope": "station", "station_id": sid, "family": fam, "x_variable": xv,
                                  "y_variable": yv, "method": method, **c})
    if epoch_tbl is not None and len(epoch_tbl):
        for fam, xv in EPOCH_X:
            if xv not in epoch_tbl.columns:
                continue
            for yv in EPOCH_Y + (["stage4_sigma_3d_m"] if fam == "geometry" else []):
                if yv not in epoch_tbl.columns:
                    continue
                xa = epoch_tbl[xv].to_numpy(float)
                c = correlation(xa, epoch_tbl[yv].to_numpy(float), epoch_tbl["time_s"].to_numpy(float), method)
                fin = xa[np.isfinite(xa)]
                if fam == "geometry" and fin.size and np.ptp(fin) < MIN_GEOMETRY_VARIATION * abs(np.median(fin)):
                    c["x_unique"], c["x_minority"] = 0, 0  # geometry effectively constant over the pass
                tests.append({"scope": "network", "station_id": "", "family": fam, "x_variable": xv,
                              "y_variable": yv, "method": method, **c})
    df = pd.DataFrame(tests)
    if not len(df):
        return pd.DataFrame(columns=["scope", "station_id", "family", "x_variable", "y_variable", "method",
                                     "coefficient", "partial_coefficient", "p_value", "p_value_holm",
                                     "n_samples", "n_effective", "interpretation_tag"])
    expected = df.apply(lambda r: (r["family"], r["y_variable"]) in EXPECTED, axis=1)
    p = df["p_value"].to_numpy(float).copy()
    p[expected.to_numpy()] = np.nan  # expected relations are not hypothesis tests
    df["p_value_holm"] = np.nan
    for fam in df["family"].unique():
        sel = (df["family"] == fam).to_numpy()
        df.loc[sel, "p_value_holm"] = holm(p[sel])
    tags = []
    for i, r in df.iterrows():
        if expected.iloc[i]:
            tags.append("EXPECTED")
        elif r["n_samples"] < min_samples:
            tags.append("INSUFFICIENT_SAMPLES")
        elif not (r["x_unique"] >= 3 or (r["x_unique"] == 2 and r.get("x_minority", 0) >= 3)) or not np.isfinite(
                r["partial_coefficient"]):
            tags.append("INSUFFICIENT_VARIATION")
        elif abs(r["partial_coefficient"]) >= min_correlation and r["p_value_holm"] < alpha:
            tags.append(FAMILY_TAG[r["family"]])
        else:
            tags.append("NOT_SUPPORTED")
    df["interpretation_tag"] = tags
    return df


def supported_tags(corr) -> Dict[str, List[Dict[str, object]]]:
    """Family tags supported by at least one test, with the supporting tests."""
    out: Dict[str, List[Dict[str, object]]] = {}
    if corr is None or not len(corr):
        return out
    for _, r in corr.iterrows():
        if r["interpretation_tag"] in FAMILY_TAG.values():
            out.setdefault(r["interpretation_tag"], []).append(r.to_dict())
    return out


def strongest(corr, family: str):
    """The test of ``family`` with the largest |partial coefficient| (None if none)."""
    if corr is None or not len(corr):
        return None
    g = corr[(corr["family"] == family) & corr["partial_coefficient"].notna()
             & ~corr["interpretation_tag"].isin(["EXPECTED", "INSUFFICIENT_SAMPLES"])]
    if not len(g):
        return None
    return g.loc[g["partial_coefficient"].abs().idxmax()].to_dict()
