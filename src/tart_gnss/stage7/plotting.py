"""Stage-7 diagnostic plots (matplotlib, Agg backend)."""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

PLOT_NAMES = ("rtn_evidence", "station_influence", "station_direction_residuals", "calibration_correlations",
              "timing_vs_along_track", "geometry_vs_residual", "influence_heatmap")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save(fig, path: Path, plt) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)


def make_plots(diag, infl, corr, out_dir: Path, title: str) -> List[str]:
    """Write the seven diagnostic PNGs; return their paths."""
    plt = _plt()
    out_dir.mkdir(parents=True, exist_ok=True)
    d = diag[diag["diagnosed"].fillna(False).astype(bool)]
    if not len(d):
        return []
    t0 = d["time_utc"].min()
    tm = np.asarray((d["time_utc"] - t0).dt.total_seconds(), float) / 60.0
    xlabel = f"minutes after {t0:%Y-%m-%d %H:%M:%S} UTC"
    lvl = d["evidence_level"].astype(float).to_numpy()
    paths = []

    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, c in zip(axes, ("radial", "along_track", "cross_track")):
        sc = ax.scatter(tm, d[f"dr_{c}_m"].astype(float) / 1e3, c=lvl, cmap="viridis", vmin=0, vmax=5, s=14)
        ax.set_ylabel(f"d{c[0].upper()} (km)")
        ax.grid(alpha=0.3)
    fig.colorbar(sc, ax=axes, label="evidence level")
    axes[-1].set_xlabel(xlabel)
    fig.suptitle(f"{title}: observed - TLE RTN residuals coloured by evidence level")
    fig.savefig(out_dir / "rtn_evidence.png", dpi=110)
    plt.close(fig)
    paths.append(str(out_dir / "rtn_evidence.png"))

    stations = sorted(infl["station_id"].dropna().unique()) if len(infl) else []
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for i, sid in enumerate(stations):
        g = infl[infl["station_id"] == sid]
        axes[0].boxplot([g["loo_shift_3d_m"].dropna() / 1e3], positions=[i], widths=0.6)
        axes[1].boxplot([g["loo_mahalanobis"].dropna()], positions=[i], widths=0.6)
    for ax, lab in zip(axes, ("leave-one-out shift (km)", "leave-one-out D² (C_loo - C_full)")):
        ax.set_xticks(range(len(stations)))
        ax.set_xticklabels(stations, rotation=20)
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axes[1].axhline(14.156, color="C3", ls="--")
    fig.suptitle(f"{title}: station leave-one-out influence")
    paths.append(_save(fig, out_dir / "station_influence.png", plt))

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for sid in stations:
        g = infl[infl["station_id"] == sid]
        tg = np.asarray((g["time_utc"] - t0).dt.total_seconds(), float) / 60.0
        axes[0].plot(tg, g["direction_residual_angle_deg"], "o", ms=3, label=sid)
        axes[1].plot(tg, g["direction_residual_along_motion_deg"], "o", ms=3, label=sid)
    axes[0].set_ylabel("direction residual (deg)")
    axes[1].set_ylabel("along apparent motion (deg)")
    axes[1].set_xlabel(xlabel)
    axes[0].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: station direction vs Stage-5 trajectory (measured - predicted)")
    paths.append(_save(fig, out_dir / "station_direction_residuals.png", plt))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    c = corr[corr["family"].isin(["calibration", "localisation"])] if len(corr) else corr
    c = c[c["partial_coefficient"].notna()].copy() if len(c) else c
    if len(c):
        c["abs"] = c["partial_coefficient"].abs()
        c = c.sort_values("abs", ascending=False).head(10)
        labels = [f"{r.station_id}: {r.x_variable} / {r.y_variable}" for r in c.itertuples()]
        colors = ["C3" if t.endswith("_CORRELATED") else "0.6" for t in c["interpretation_tag"]]
        ax.barh(range(len(c)), c["partial_coefficient"], color=colors)
        ax.set_yticks(range(len(c)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
    ax.set_xlabel("time-detrended rank correlation (red: supported)")
    ax.set_xlim(-1, 1)
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}: strongest calibration / localisation correlations")
    paths.append(_save(fig, out_dir / "calibration_correlations.png", plt))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for sid in stations:
        g = infl[infl["station_id"] == sid]
        axes[0].plot(g["station_time_offset_s"], g["direction_residual_along_motion_deg"], "o", ms=3, label=sid)
    axes[0].set_xlabel("station time offset (s)")
    axes[0].set_ylabel("direction residual along motion (deg)")
    axes[0].legend(fontsize=8)
    axes[1].plot(d["mean_station_time_offset_s"].astype(float), d["dr_along_track_m"].astype(float) / 1e3, "o", ms=3)
    axes[1].set_xlabel("mean station time offset (s)")
    axes[1].set_ylabel("Stage-6 along-track residual (km)")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: timing vs along-track behaviour")
    paths.append(_save(fig, out_dir / "timing_vs_along_track.png", plt))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].semilogx(d["stage4_condition_number"].astype(float), d["stage4_residual_3d_m"].astype(float) / 1e3, "o",
                     ms=3)
    axes[0].set_xlabel("Stage-4 condition number")
    axes[0].set_ylabel("|Stage-4 position - TLE| (km)")
    axes[1].plot(d["stage4_min_crossing_angle_deg"].astype(float), d["position_residual_3d_m"].astype(float) / 1e3,
                 "o", ms=3)
    axes[1].set_xlabel("minimum crossing angle (deg)")
    axes[1].set_ylabel("Stage-6 3-D residual (km)")
    for ax in axes:
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(f"{title}: geometry dependence")
    paths.append(_save(fig, out_dir / "geometry_vs_residual.png", plt))

    fig, ax = plt.subplots(figsize=(10, 3 + 0.4 * max(len(stations), 1)))
    if len(infl) and infl["loo_mahalanobis"].notna().any():
        piv = infl.pivot_table(index="station_id", columns="time_utc", values="loo_mahalanobis")
        im = ax.imshow(np.log10(np.clip(piv.to_numpy(float), 1e-3, None)), aspect="auto", cmap="magma",
                       interpolation="nearest")
        fig.colorbar(im, ax=ax, label="log10 leave-one-out D²")
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels(piv.index)
        ax.set_xlabel("epoch index")
    else:
        ax.text(0.5, 0.5, "no leave-one-out epochs (fewer than 3 stations)", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_title(f"{title}: station-by-time influence")
    paths.append(_save(fig, out_dir / "influence_heatmap.png", plt))
    return paths
