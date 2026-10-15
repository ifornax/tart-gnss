"""Pass-to-pass repeatability for satellites observed more than once.

Per satellite with >= 2 passes (aggregated per pass, never per epoch):

* RTN: between-pass standard deviation and MAD of the per-pass median
  residual, against the within-pass median Stage-5 sigma of that component
  (``ratio`` > 1 means the passes scatter more than each pass's own
  uncertainty);
* direction-bias repeatability: per station, the between-pass scatter of the
  per-pass median direction residual;
* station-influence repeatability: between-pass scatter of the per-pass
  median leave-one-out D^2;
* evidence-tag repeatability: mean Jaccard similarity of the tag sets over
  all pairs of passes.
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict

import numpy as np

from tart_gnss.campaign.statistics import MAD_TO_SIGMA

COMPONENTS = (("radial", "median_radial_m", "median_sigma_radial_obs_m"),
              ("along_track", "median_along_track_m", "median_sigma_along_track_obs_m"),
              ("cross_track", "median_cross_track_m", "median_sigma_cross_track_obs_m"))


def _sd(v):
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(np.std(a, ddof=1)) if a.size > 1 else float("nan")


def _mad(v):
    a = np.asarray(v, float)
    a = a[np.isfinite(a)]
    return float(MAD_TO_SIGMA * np.median(np.abs(a - np.median(a)))) if a.size else float("nan")


def repeatability(runs, stations) -> Dict[str, object]:
    """Repeatability statistics per satellite (only satellites with >= 2 passes)."""
    out: Dict[str, object] = {}
    for nid, g in runs.groupby("norad_id"):
        if len(g) < 2:
            continue
        rec: Dict[str, object] = {"n_passes": int(len(g)), "run_ids": list(g["run_id"]),
                                  "cohorts": sorted(set(g["cohort"]))}
        for name, col, sig in COMPONENTS:
            sd, mad = _sd(g[col]), _mad(g[col])
            within = float(np.nanmedian(g[sig].astype(float))) if np.any(np.isfinite(g[sig].astype(float))) else \
                float("nan")
            rec[name] = {"between_pass_std_m": sd, "between_pass_mad_m": mad, "within_pass_median_sigma_m": within,
                         "ratio_between_to_within": sd / within if np.isfinite(sd) and within and np.isfinite(within)
                         else float("nan")}
        st = stations[stations["run_id"].isin(set(g["run_id"]))] if len(stations) else stations
        rec["direction_bias_between_pass_std_deg"] = {
            sid: _sd(h["median_direction_residual_deg"]) for sid, h in st.groupby("station_id")} if len(st) else {}
        rec["station_influence_between_pass_std"] = {
            sid: _sd(h["median_leave_one_out_mahalanobis"]) for sid, h in st.groupby("station_id")} if len(st) else {}
        tags = [set(t or []) for t in g["overall_diagnostic_tags"]]
        sims = [len(a & b) / len(a | b) if a | b else 1.0 for a, b in combinations(tags, 2)]
        rec["evidence_tag_jaccard_mean"] = float(np.mean(sims)) if sims else float("nan")
        rec["evidence_level_values"] = [int(x) for x in g["overall_evidence_level"]]
        out[str(nid)] = rec
    return out
