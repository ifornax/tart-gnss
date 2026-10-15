"""Evidence-level and diagnostic-tag population statistics.

Three separate denominators, never mixed:

* **run**: each run's overall Stage-7 evidence level and tags;
* **epoch**: each diagnosed epoch's level and tags;
* **station**: each (run, station) row, with the run's overall level and the
  station's own Stage-7 tags.

Fractions come with Wilson 95 % intervals at the run level (the independent
unit). Tags are evidence patterns: no physical-cause rates (manoeuvres, clock
or hardware failures, RFI events) are ever produced.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from tart_gnss.campaign.statistics import wilson_interval

TAGS = ("STATISTICAL_ANOMALY", "STATION_SPECIFIC", "CALIBRATION_CORRELATED", "TIMING_CORRELATED",
        "GEOMETRY_CORRELATED", "LOCALISATION_CORRELATED", "TRAJECTORY_MODEL_CORRELATED", "REFERENCE_LIMITED",
        "NETWORK_CONSISTENT", "UNEXPLAINED")


def level_fractions(levels) -> Dict[str, float]:
    """Fraction at each evidence level 0-5 (and n)."""
    v = np.array([float(x) for x in levels if x is not None and x == x])
    out = {"n": int(v.size)}
    for k in range(6):
        out[f"level_{k}"] = float(np.mean(v == k)) if v.size else float("nan")
    return out


def evidence_summary(runs, epochs, stations):
    """Level fractions by run, epoch and station, and tag frequencies (rows for ``evidence_summary.csv``)."""
    import pandas as pd

    rows: List[Dict[str, object]] = []
    for unit, df, lcol, tcol in (("run", runs, "overall_evidence_level", "overall_diagnostic_tags"),
                                 ("epoch", epochs, "evidence_level", "diagnostic_tags"),
                                 ("station", stations, "overall_evidence_level", "station_tags")):
        lf = level_fractions(df[lcol]) if len(df) else level_fractions([])
        for k in range(6):
            rows.append({"unit": unit, "quantity": f"evidence_level_{k}", "n": lf["n"],
                         "count": int(round(lf[f"level_{k}"] * lf["n"])) if lf["n"] else 0,
                         "fraction": lf[f"level_{k}"]})
        tags = TAGS if unit != "station" else ("STATION_SPECIFIC", "CALIBRATION_CORRELATED", "TIMING_CORRELATED",
                                               "LOCALISATION_CORRELATED")
        for t in tags:
            k = int(sum(t in (x or []) for x in df[tcol])) if len(df) else 0
            n = int(len(df))
            rows.append({"unit": unit, "quantity": t, "n": n, "count": k, "fraction": k / n if n else float("nan")})
    out = pd.DataFrame(rows)
    ci = [wilson_interval(int(r["count"]), int(r["n"])) if r["unit"] == "run" else (None, None)
          for _, r in out.iterrows()]
    out["ci95_low"] = [c[0] for c in ci]
    out["ci95_high"] = [c[1] for c in ci]
    out["ci_note"] = ["Wilson, runs independent" if r["unit"] == "run" else "no CI: rows within a run are correlated"
                      for _, r in out.iterrows()]
    return out


def rates_vs_conditions(runs, stations):
    """Fraction of runs carrying each tag versus station, constellation, crossing-angle bin and TLE-age bin."""
    import pandas as pd

    from tart_gnss.campaign.constellation_metrics import CROSSING_EDGES, TLE_AGE_EDGES

    rows = []

    def add(dim, key, g):
        n = len(g)
        for t in TAGS:
            k = int(sum(t in (x or []) for x in g["overall_diagnostic_tags"]))
            lo, hi = wilson_interval(k, n)
            rows.append({"condition": dim, "value": str(key), "tag": t, "n_runs": n, "count": k,
                         "fraction": k / n if n else float("nan"), "ci95_low": lo, "ci95_high": hi})

    if not len(runs):
        return pd.DataFrame()
    for c, g in runs.groupby("constellation"):
        add("constellation", c, g)
    for sid in sorted(stations["station_id"].unique()) if len(stations) else []:
        ids = set(stations.loc[stations["station_id"] == sid, "run_id"])
        add("station", sid, runs[runs["run_id"].isin(ids)])
    ang = runs["median_min_crossing_angle_deg"].astype(float)
    for lo, hi in zip(CROSSING_EDGES[:-1], CROSSING_EDGES[1:]):
        g = runs[(ang >= lo) & (ang < hi)]
        if len(g):
            add("min_crossing_angle_deg", f"{lo:g}-{hi:g}", g)
    age = runs["median_tle_age_days"].astype(float).abs()
    for lo, hi in zip(TLE_AGE_EDGES[:-1], TLE_AGE_EDGES[1:]):
        g = runs[(age >= lo) & (age < hi)]
        if len(g):
            add("abs_tle_age_days", f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}", g)
    return pd.DataFrame(rows)
