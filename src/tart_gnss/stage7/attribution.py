"""Evidence tags and evidence levels from transparent rules.

Statistical anomaly
    Stage-6 ``position_mahalanobis_total_cov`` when an empirical TLE reference
    covariance exists, else ``position_mahalanobis_obs_cov`` (conditional on
    the TLE being fixed). Default threshold: the chi^2_3 quantile at
    ``1 - 0.0027 / n`` (Bonferroni over the ``n`` diagnosed epochs), so a
    clean run of any length stays at Level 0 with probability ~0.997.
    ``--anomaly-threshold`` overrides it.

Per-epoch independent evidence
    Stage-6 residuals come from a smoothed track and are strongly
    correlated in time, so persistence and strength are judged on the
    epoch-by-epoch Stage-4 position minus the TLE, ``stage4_d2_tle``
    (chi^2_3 0.9973 = 14.16): *strong* when it exceeds that limit at the
    epoch, *persistent* when it does so in at least ``min_persistent_epochs``
    consecutive anomalous epochs.

Station test (epochs with >= 3 stations)
    * pairwise: station pairs whose sight lines miss each other by more than
      ``pair_z`` predicted sigmas are inconsistent; a station is picked if it
      belongs to every inconsistent pair and every pair without it is
      consistent;
    * reference: the Stage-4 solution disagrees with the TLE
      (``full_d2_tle`` > 14.16) and removing exactly one station makes the
      disagreement go away (``loo_d2_tle`` <= 14.16) - only when the
      reference is not limited.
    Both agree -> that station (cross-method agreement); one method -> that
    station; conflict -> ``AMBIGUOUS``. ``STATION_SPECIFIC`` is a run-level
    tag for a station picked at >= 3 epochs and at >= half of the epochs
    with >= 3 stations.

Tags (non-exclusive)
    ``STATISTICAL_ANOMALY``, ``STATION_SPECIFIC``, ``CALIBRATION_CORRELATED``,
    ``TIMING_CORRELATED``, ``LOCALISATION_CORRELATED``, ``GEOMETRY_CORRELATED``
    and ``TRAJECTORY_MODEL_CORRELATED`` (supported correlations, see
    :mod:`tart_gnss.stage7.correlations`, plus the epoch evidence below),
    ``REFERENCE_LIMITED``, ``NETWORK_CONSISTENT``, ``UNEXPLAINED``, ``NORMAL``.
    Epoch evidence: ``TIMING`` if a contributing station has ``TIME_MISMATCH``;
    ``LOCALISATION`` if Stage 4 flagged ``HIGH_LINE_RESIDUAL``/``HIGH_CHI2`` or a
    pair is inconsistent; ``TRAJECTORY_MODEL`` if Stage 5 raised
    ``DYNAMICS_MISMATCH``/``HIGH_NIS``/a filter restart for the run, or
    ``INNOVATION_OUTLIER``/``FILTER_DIVERGENCE`` at the epoch;
    ``REFERENCE_LIMITED`` if Stage 6 is not a precision validation or the
    epoch is ``TLE_STALE``/``TLE_REFERENCE_DISPERSION_HIGH``.

Evidence level (per anomalous epoch; 0 for a normal one)
    1  statistical anomaly
    2  + ``STATION_SPECIFIC`` (the station is picked at this epoch)
    3  + calibration / timing / localisation evidence
    4  ``NETWORK_CONSISTENT``: >= 3 stations, no leave-one-out shift
       significant (``loo_mahalanobis`` within the chi^2_dof 0.9973 quantile,
       ``dof = loo_dof``, for every station), no pair
       inconsistent, no station picked, no instrument / geometry /
       trajectory-model tag, reference not limited
    5  ``UNEXPLAINED``: Level 4, strong and persistent, calibration metadata
       available
    ``REFERENCE_LIMITED`` caps the level at 3; missing calibration metadata
    caps it at 4. Level 5 is an evidence level, never a physical cause.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List, Optional

import numpy as np

from tart_gnss.stage7.scoring import confidence
from tart_gnss.stage7.statistics import runs

D2_3_9973 = 14.156
INSTRUMENT_TAGS = ("CALIBRATION_CORRELATED", "TIMING_CORRELATED", "LOCALISATION_CORRELATED")
TAG_ORDER = ("NORMAL", "STATISTICAL_ANOMALY", "STATION_SPECIFIC", "CALIBRATION_CORRELATED", "TIMING_CORRELATED",
             "GEOMETRY_CORRELATED", "LOCALISATION_CORRELATED", "TRAJECTORY_MODEL_CORRELATED", "REFERENCE_LIMITED",
             "NETWORK_CONSISTENT", "UNEXPLAINED")


@dataclass
class EvidenceRules:
    """Thresholds of the evidence rules (recorded in the manifest)."""

    anomaly_threshold: Optional[float] = None
    anomaly_tail_probability: float = 0.0027
    stage4_threshold: float = D2_3_9973
    pair_z: float = 4.0
    min_persistent_epochs: int = 5
    min_samples: int = 8
    min_correlation: float = 0.5
    alpha: float = 0.01
    station_min_epochs: int = 3
    station_min_fraction: float = 0.5

    def threshold_for(self, n: int) -> float:
        """The Stage-6 anomaly threshold for ``n`` diagnosed epochs."""
        if self.anomaly_threshold is not None:
            return float(self.anomaly_threshold)
        from scipy.stats import chi2

        return float(chi2.ppf(1.0 - self.anomaly_tail_probability / max(n, 1), 3))

    def to_dict(self, n: int) -> Dict[str, object]:
        """JSON-serialisable rules including the threshold actually used."""
        return {**asdict(self), "anomaly_threshold_used": self.threshold_for(n),
                "anomaly_threshold_rule": "user value" if self.anomaly_threshold is not None else
                f"chi2_3 quantile at 1 - {self.anomaly_tail_probability} / n_epochs (Bonferroni)"}


def loo_limit(dof: int, probability: float = 0.9973) -> float:
    """chi^2_dof quantile used for a leave-one-out shift (inf for dof 0)."""
    if dof <= 0:
        return float("inf")
    from scipy.stats import chi2

    return float(chi2.ppf(probability, dof))


def ordered(tags) -> List[str]:
    """Unique tags in canonical order."""
    s = set(tags)
    return [t for t in TAG_ORDER if t in s] + sorted(s - set(TAG_ORDER))


def station_pick(pairs_ep, loo_ep, reference_ok: bool, rules: EvidenceRules):
    """``(station or None or 'AMBIGUOUS', method, agreement)`` for one epoch."""
    pick_pair = None
    if pairs_ep is not None and len(pairs_ep):
        bad = pairs_ep[pairs_ep["pair_separation_z"] > rules.pair_z]
        stations = set(pairs_ep["pair_station_1"]) | set(pairs_ep["pair_station_2"])
        if len(bad):
            cands = [s for s in stations if all(s in (a, b) for a, b in zip(bad["pair_station_1"],
                                                                          bad["pair_station_2"]))]
            good = [s for s in cands if all(
                z <= rules.pair_z for a, b, z in zip(pairs_ep["pair_station_1"], pairs_ep["pair_station_2"],
                                                     pairs_ep["pair_separation_z"]) if s not in (a, b))]
            if len(good) == 1:
                pick_pair = good[0]
            elif len(good) > 1:
                pick_pair = "AMBIGUOUS"
    pick_ref = None
    if reference_ok and loo_ep is not None and len(loo_ep):
        full = float(loo_ep["full_d2_tle"].iloc[0])
        if np.isfinite(full) and full > rules.stage4_threshold:
            gone = loo_ep[loo_ep["loo_d2_tle"] <= rules.stage4_threshold]
            if len(gone) == 1:
                pick_ref = gone["station_id"].iloc[0]
    named = [p for p in (pick_pair, pick_ref) if p not in (None, "AMBIGUOUS")]
    if pick_pair == "AMBIGUOUS" and pick_ref is not None:
        return pick_ref, "reference", False
    if len(named) == 2:
        return (named[0], "pairwise+reference", True) if named[0] == named[1] else ("AMBIGUOUS", "conflict", False)
    if len(named) == 1:
        return named[0], ("pairwise" if pick_pair == named[0] else "reference"), False
    return (pick_pair, "pairwise", False) if pick_pair == "AMBIGUOUS" else (None, "", False)


def assign_evidence(epochs, stations_tbl, influence, pairs, corr_tags: Dict[str, List[Dict[str, object]]],
                    stage5_run_flags: List[str], stage5_resets: int, precision_validation: bool,
                    calibration_available: bool, data_complete: bool, rules: EvidenceRules):
    """Assign diagnostic tags and evidence level per epoch; return ``(rows, run_info)``."""
    used = epochs["diagnosed"].to_numpy(bool)
    n_used = int(used.sum())
    thr = rules.threshold_for(n_used)
    stat_col = "position_mahalanobis_total_cov" if epochs["position_mahalanobis_total_cov"].notna().any() \
        else "position_mahalanobis_obs_cov"
    stat = epochs[stat_col].to_numpy(float)
    anomalous = used & np.isfinite(stat) & (stat > thr)
    strong = anomalous & (epochs["stage4_d2_tle"].to_numpy(float) > rules.stage4_threshold)
    t_s = epochs["time_s"].to_numpy(float)
    persistent = np.zeros(len(epochs), bool)
    for a, b, length, _ in runs(strong, t_s):
        if length >= rules.min_persistent_epochs:
            persistent[a:b + 1] = True
    traj_run = bool(set(stage5_run_flags) & {"DYNAMICS_MISMATCH", "HIGH_NIS"}) or stage5_resets > 0 \
        or "TRAJECTORY_MODEL_CORRELATED" in corr_tags
    geometry_run = "GEOMETRY_CORRELATED" in corr_tags
    net_tags = {t for t, tests in corr_tags.items() if any(x["scope"] == "network" for x in tests)}
    station_tags: Dict[str, set] = {}
    for t, tests in corr_tags.items():
        for x in tests:
            if x["scope"] == "station":
                station_tags.setdefault(x["station_id"], set()).add(t)

    # Station test per epoch.
    picks: Dict[int, tuple] = {}
    n3 = 0
    for ei in np.flatnonzero(used):
        loo_ep = influence[influence["epoch_index"] == ei] if len(influence) else None
        pairs_ep = pairs[pairs["epoch_index"] == ei] if len(pairs) else None
        if loo_ep is not None and len(loo_ep):
            n3 += 1
        ref_ok = precision_validation and "TLE_STALE" not in epochs.iloc[ei]["stage6_flags"]
        picks[int(ei)] = station_pick(pairs_ep, loo_ep, ref_ok, rules)
    counts: Dict[str, int] = {}
    for p in picks.values():
        if p[0] not in (None, "AMBIGUOUS"):
            counts[p[0]] = counts.get(p[0], 0) + 1
    specific = {s for s, c in counts.items()
                if c >= rules.station_min_epochs and c >= rules.station_min_fraction * max(n3, 1)}

    out = []
    for ei in range(len(epochs)):
        ep = epochs.iloc[ei]
        tags: List[str] = []
        comp: Dict[str, bool] = {}
        level = 0
        pick = picks.get(ei, (None, "", False))
        rec: Dict[str, object] = {"anomaly_statistic": stat_col, "anomaly_threshold": thr,
                                  "statistical_anomaly": bool(anomalous[ei]), "strong": bool(strong[ei]),
                                  "persistent": bool(persistent[ei]), "station_pick": pick[0] or "",
                                  "station_pick_method": pick[1]}
        if not used[ei]:
            rec.update({"diagnostic_tags": [], "evidence_level": np.nan, "diagnostic_confidence": "",
                        "confidence_components": {}})
            out.append(rec)
            continue
        loo_ep = influence[influence["epoch_index"] == ei] if len(influence) else influence
        pairs_ep = pairs[pairs["epoch_index"] == ei] if len(pairs) else pairs
        contributing = set(stations_tbl.loc[stations_tbl["epoch_index"] == ei, "station_id"]) \
            if len(stations_tbl) else set()
        reference_limited = (not precision_validation) or bool(
            set(ep["stage6_flags"]) & {"TLE_STALE", "TLE_REFERENCE_DISPERSION_HIGH"})
        pair_bad = bool(len(pairs_ep) and (pairs_ep["pair_separation_z"] > rules.pair_z).any())
        if anomalous[ei]:
            tags.append("STATISTICAL_ANOMALY")
            level = 1
            if pick[0] in specific:
                tags.append("STATION_SPECIFIC")
                level = 2
            ins = set(net_tags & set(INSTRUMENT_TAGS))
            for s in contributing:
                ins |= station_tags.get(s, set()) & set(INSTRUMENT_TAGS)
            if len(stations_tbl):
                rows = stations_tbl[stations_tbl["epoch_index"] == ei]
                if rows["time_mismatch"].any():
                    ins.add("TIMING_CORRELATED")
            if set(ep["stage4_flags"]) & {"HIGH_LINE_RESIDUAL", "HIGH_CHI2"} or pair_bad:
                ins.add("LOCALISATION_CORRELATED")
            tags += sorted(ins)
            if ins:
                level = max(level, 3)
            traj = traj_run or bool(set(ep["stage5_flags"]) & {"INNOVATION_OUTLIER", "FILTER_DIVERGENCE"})
            if traj:
                tags.append("TRAJECTORY_MODEL_CORRELATED")
            if geometry_run:
                tags.append("GEOMETRY_CORRELATED")
            if reference_limited:
                tags.append("REFERENCE_LIMITED")
            loo_clean = bool(len(loo_ep)) and all(
                np.isfinite(d) and d <= loo_limit(int(k)) for d, k in zip(loo_ep["loo_mahalanobis"], loo_ep["loo_dof"]))
            if (loo_clean and not pair_bad and pick[0] in (None, "") and not ins and not traj and not geometry_run
                    and not reference_limited):
                tags.append("NETWORK_CONSISTENT")
                level = 4
                if strong[ei] and persistent[ei] and calibration_available:
                    tags.append("UNEXPLAINED")
                    level = 5
            if reference_limited:
                level = min(level, 3)
            if not calibration_available:
                level = min(level, 4)
            comp = {
                "effect_size": bool(stat[ei] >= 2 * thr), "persistence": bool(persistent[ei]),
                "cross_method_agreement": bool(pick[2]) or (loo_clean and not pair_bad and len(pairs_ep) > 0),
            }
        else:
            tags.append("NORMAL")
            comp = {"effect_size": bool(np.isfinite(stat[ei]) and stat[ei] <= thr / 2), "persistence": True,
                    "cross_method_agreement": True}
        comp.update({"sample_support": n_used >= rules.min_samples, "reference_quality": bool(precision_validation),
                     "data_completeness": bool(data_complete)})
        rec.update({"diagnostic_tags": ordered(tags), "evidence_level": level,
                    "diagnostic_confidence": confidence(comp), "confidence_components": comp})
        out.append(rec)

    anomaly_runs = [{"start_utc": str(epochs.iloc[a]["time_utc"]), "end_utc": str(epochs.iloc[b]["time_utc"]),
                     "n_epochs": n, "duration_s": d} for a, b, n, d in runs(anomalous, t_s)]
    lvl = np.array([r["evidence_level"] if r["evidence_level"] == r["evidence_level"] else -1 for r in out])
    unexplained = [{"start_utc": str(epochs.iloc[a]["time_utc"]), "end_utc": str(epochs.iloc[b]["time_utc"]),
                    "n_epochs": n, "duration_s": d} for a, b, n, d in runs(lvl == 5, t_s)]
    longest = max((n for _, _, n, _ in runs(anomalous, t_s)), default=0)
    info = {
        "n_diagnosed": n_used, "anomaly_threshold": thr, "anomaly_statistic": stat_col,
        "n_anomalous": int(anomalous.sum()), "n_strong": int(strong.sum()), "n_persistent": int(persistent.sum()),
        "anomaly_intervals": anomaly_runs, "unexplained_intervals": unexplained,
        "longest_anomaly_run": longest,
        "station_picks": counts, "n_epochs_with_station_test": n3, "station_specific": sorted(specific),
        "trajectory_model_run_level": traj_run, "geometry_run_level": geometry_run,
    }
    return out, info
