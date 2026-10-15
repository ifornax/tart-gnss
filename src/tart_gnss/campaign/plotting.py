"""Campaign figures (matplotlib only); ``publication`` mode adds larger fonts and vector PDFs."""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

PLOT_NAMES = ("residual_distribution", "rtn_distribution", "residual_vs_elevation", "residual_vs_snr",
              "uncertainty_vs_crossing_angle", "station_stability", "constellation_comparison",
              "evidence_level_distribution", "tle_age_dependence", "repeatability")


def _plt(publication: bool):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 13 if publication else 9, "axes.titlesize": 14 if publication else 10,
                         "axes.labelsize": 13 if publication else 9, "legend.fontsize": 11 if publication else 8})
    return plt


def _save(fig, out: Path, name: str, plt, publication: bool) -> List[str]:
    fig.tight_layout()
    paths = [out / f"{name}.png"]
    fig.savefig(paths[0], dpi=200 if publication else 110)
    if publication:
        paths.append(out / f"{name}.pdf")
        fig.savefig(paths[1])
    plt.close(fig)
    return [str(p) for p in paths]


def _ecdf(ax, v, label):
    v = np.sort(np.asarray(v, float)[np.isfinite(np.asarray(v, float))])
    if v.size:
        ax.step(v, np.arange(1, v.size + 1) / v.size, where="post", label=f"{label} (n={v.size})")


def make_plots(runs, epochs, stations, station_summary, science, repeat, out_dir: Path, title: str,
               publication: bool = False) -> List[str]:
    """Write the ten campaign figures; return their paths."""
    plt = _plt(publication)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[str] = []
    sci_label = "precision cohort" if len(science) and science["precision_validation"].all() else "analysis cohort"

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    _ecdf(axes[0], science["median_position_residual_3d_m"] / 1e3, "per-run median")
    axes[0].set_xlabel("TLE-relative 3-D residual, per-run median (km)")
    axes[0].set_ylabel("cumulative fraction of runs")
    if axes[0].get_legend_handles_labels()[0]:
        axes[0].legend()
    sci_ep = epochs[epochs["run_id"].isin(set(science["run_id"]))] if len(epochs) else epochs
    if len(sci_ep):
        axes[1].hist(sci_ep["position_residual_3d_m"].astype(float) / 1e3, bins=40)
    axes[1].set_xlabel("TLE-relative 3-D residual, all epochs (km; display only, epochs correlated)")
    axes[1].set_ylabel("epochs")
    fig.suptitle(f"{title}: TLE-relative residual ({sci_label})")
    paths += _save(fig, out_dir, "residual_distribution", plt, publication)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (c, lab) in zip(axes, (("median_radial_m", "radial"), ("median_along_track_m", "along-track"),
                                   ("median_cross_track_m", "cross-track"))):
        v = science[c].astype(float) / 1e3
        ax.hist(v[np.isfinite(v)], bins=max(5, min(30, len(v))))
        ax.axvline(0, color="0.3", lw=1)
        ax.set_xlabel(f"{lab} residual, per-run median (km)")
    axes[0].set_ylabel("runs")
    fig.suptitle(f"{title}: RTN TLE-relative residuals ({sci_label})")
    paths += _save(fig, out_dir, "rtn_distribution", plt, publication)

    for name, x, xl in (("residual_vs_elevation", "median_elevation_deg", "station elevation (deg)"),
                        ("residual_vs_snr", "median_snr", "SNR")):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        axes[0].plot(stations[x].astype(float), stations["median_direction_residual_deg"].astype(float), "o", ms=4)
        axes[0].set_xlabel(f"per-run median {xl}")
        axes[0].set_ylabel("per-run median direction residual (deg)")
        axes[1].plot(stations[x].astype(float), stations["median_sigma_theta_deg"].astype(float), "o", ms=4)
        axes[1].set_xlabel(f"per-run median {xl}")
        axes[1].set_ylabel("per-run median direction sigma (deg)")
        for ax in axes:
            ax.grid(alpha=0.3)
        fig.suptitle(f"{title}: station direction performance vs {xl} (one point per run x station)")
        paths += _save(fig, out_dir, name, plt, publication)

    fig, ax = plt.subplots(figsize=(8, 5))
    for c, g in runs.groupby("constellation"):
        ax.semilogy(g["median_min_crossing_angle_deg"].astype(float), g["median_stage4_sigma_3d_m"].astype(float) / 1e3,
                    "o", label=c)
    ax.set_xlabel("per-run median minimum crossing angle (deg)")
    ax.set_ylabel("per-run median Stage-4 3-D sigma (km)")
    if ax.get_legend_handles_labels()[0]:
        ax.legend()
    ax.grid(alpha=0.3, which="both")
    ax.set_title(f"{title}: measurement uncertainty vs network geometry")
    paths += _save(fig, out_dir, "uncertainty_vs_crossing_angle", plt, publication)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    if len(station_summary):
        x = np.arange(len(station_summary))
        labels = list(station_summary["station_id"])
        axes[0].bar(x, station_summary["median_direction_residual_deg"],
                    yerr=station_summary["mad_direction_residual_deg"], capsize=3)
        axes[0].set_ylabel("median direction residual (deg) ± run-to-run MAD")
        axes[1].bar(x, station_summary["median_leave_one_out_mahalanobis"])
        axes[1].set_ylabel("median leave-one-out D²")
        axes[1].set_yscale("log")
        axes[2].bar(x, station_summary["fraction_runs_any_station_tag"])
        axes[2].set_ylabel("fraction of runs with a station tag")
        for ax in axes:
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=25)
            ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: station stability across runs")
    paths += _save(fig, out_dir, "station_stability", plt, publication)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    groups = sorted(runs["constellation"].unique())
    for ax, (c, lab, src) in zip(axes, (("median_snr", "per-run median SNR", runs),
                                        ("median_stage4_sigma_3d_m", "Stage-4 3-D sigma (m)", runs),
                                        ("median_position_residual_3d_m", "TLE-relative residual (m)", science))):
        data = [src.loc[src["constellation"] == g, c].astype(float).dropna() for g in groups]
        ax.boxplot([d if len(d) else [np.nan] for d in data])
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(groups, data)])
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: constellation comparison (per-run values; small groups are under-sampled)")
    paths += _save(fig, out_dir, "constellation_comparison", plt, publication)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    lv = np.arange(6)
    rl = runs["overall_evidence_level"].astype(float)
    axes[0].bar(lv, [np.mean(rl == k) if len(rl) else 0 for k in lv])
    axes[0].set_xlabel("overall evidence level (per run)")
    axes[0].set_ylabel("fraction of runs")
    el = epochs["evidence_level"].astype(float) if len(epochs) else np.zeros(0)
    axes[1].bar(lv, [np.mean(el == k) if len(el) else 0 for k in lv])
    axes[1].set_xlabel("evidence level (per epoch)")
    axes[1].set_ylabel("fraction of epochs")
    fig.suptitle(f"{title}: evidence levels (separate denominators)")
    paths += _save(fig, out_dir, "evidence_level_distribution", plt, publication)

    fig, ax = plt.subplots(figsize=(8, 5))
    for coh, g in runs.groupby("cohort"):
        ax.plot(g["median_tle_age_days"].astype(float).abs(), g["median_position_residual_3d_m"].astype(float) / 1e3,
                "o", label=coh)
    ax.set_xlabel("|TLE age| (days)")
    ax.set_ylabel("per-run median 3-D residual (km)")
    if ax.get_legend_handles_labels()[0]:
        ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}: residual vs TLE age (descriptive)")
    paths += _save(fig, out_dir, "tle_age_dependence", plt, publication)

    fig, ax = plt.subplots(figsize=(9, 5))
    sats = sorted(repeat)
    if sats:
        x = np.arange(len(sats))
        for k, (comp, lab) in enumerate((("radial", "R"), ("along_track", "T"), ("cross_track", "N"))):
            ax.bar(x + (k - 1) * 0.25, [repeat[s][comp]["between_pass_std_m"] / 1e3 for s in sats], width=0.25,
                   label=f"{lab} between-pass std")
            ax.plot(x + (k - 1) * 0.25, [repeat[s][comp]["within_pass_median_sigma_m"] / 1e3 for s in sats], "k_",
                    ms=14)
        ax.set_xticks(x)
        ax.set_xticklabels([f"NORAD {s}\n({repeat[s]['n_passes']} passes)" for s in sats])
        ax.legend()
    ax.set_ylabel("km (bars: between-pass std; ticks: within-pass sigma)")
    ax.set_title(f"{title}: pass-to-pass repeatability")
    ax.grid(alpha=0.3)
    paths += _save(fig, out_dir, "repeatability", plt, publication)
    return paths
