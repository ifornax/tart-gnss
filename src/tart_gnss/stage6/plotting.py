"""Stage-6 diagnostic plots (matplotlib, Agg backend), drawn from the saved tables."""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

PLOT_NAMES = ("rtn_residuals", "position_residual_norm", "velocity_residuals", "angular_residual",
              "residual_vs_uncertainty", "mahalanobis", "tle_age", "trajectory_overlay")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _col(df, name) -> np.ndarray:
    return df[name].to_numpy(dtype=float)


def _save(fig, path: Path, plt) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)


def make_plots(residuals, tle_traj, out_dir: Path, title: str) -> List[str]:
    """Write the eight diagnostic PNGs; return their paths."""
    if len(residuals) == 0:
        return []
    plt = _plt()
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = residuals["time_utc"].min()
    tm = np.asarray((residuals["time_utc"] - t0).dt.total_seconds(), dtype=float) / 60.0
    xlabel = f"minutes after {t0:%Y-%m-%d %H:%M:%S} UTC"
    ok = residuals["valid"].fillna(False).to_numpy(dtype=bool)
    note = "observed - TLE"
    paths = []
    comps = (("radial", "R"), ("along_track", "T"), ("cross_track", "N"))

    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, (c, lab) in zip(axes, comps):
        d = _col(residuals, f"dr_{c}_m") / 1e3
        s = _col(residuals, f"sigma_{c}_obs_m") / 1e3
        ax.errorbar(tm[ok], d[ok], yerr=s[ok], fmt="o", ms=3, label=f"{lab}  ± Stage-5 1σ")
        if "sigma_" + c + "_ref_m" in residuals and np.any(np.isfinite(_col(residuals, f"sigma_{c}_ref_m"))):
            r = _col(residuals, f"sigma_{c}_ref_m") / 1e3
            ax.fill_between(tm, d - r, d + r, color="C1", alpha=0.25, label="TLE empirical dispersion")
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_ylabel(f"d{lab} (km)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel(xlabel)
    fig.suptitle(f"{title}: RTN position residuals ({note})")
    paths.append(_save(fig, out_dir / "rtn_residuals.png", plt))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(tm[ok], _col(residuals, "position_residual_3d_m")[ok] / 1e3, "o-", ms=3, label="|dr|")
    sig3 = np.sqrt(_col(residuals, "sigma_radial_obs_m") ** 2 + _col(residuals, "sigma_along_track_obs_m") ** 2
                   + _col(residuals, "sigma_cross_track_obs_m") ** 2) / 1e3
    ax.plot(tm[ok], sig3[ok], "--", color="0.4", label="Stage-5 3-D σ")
    ax.set_ylabel("km")
    ax.set_xlabel(xlabel)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}: 3-D position residual ({note})")
    paths.append(_save(fig, out_dir / "position_residual_norm.png", plt))

    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, (c, lab) in zip(axes, comps):
        d = _col(residuals, f"dv_{c}_m_s")
        s = _col(residuals, f"sigma_dv_{c}_obs_m_s")
        ax.errorbar(tm[ok], d[ok], yerr=s[ok], fmt="o", ms=3)
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_ylabel(f"dv{lab} (m/s)")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel(xlabel)
    fig.suptitle(f"{title}: inertial velocity residuals on the TLE RTN axes ({note})")
    paths.append(_save(fig, out_dir / "velocity_residuals.png", plt))

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    axes[0].plot(tm[ok], _col(residuals, "geocentric_angular_residual_deg")[ok], "o-", ms=3)
    axes[0].set_ylabel("geocentric angle (deg)")
    axes[1].plot(tm[ok], _col(residuals, "geocentric_range_residual_m")[ok] / 1e3, "o-", ms=3)
    axes[1].set_ylabel("geocentric range residual (km)")
    axes[1].set_xlabel(xlabel)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: angular and range residuals ({note})")
    paths.append(_save(fig, out_dir / "angular_residual.png", plt))

    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, (c, lab) in zip(axes, comps):
        d = _col(residuals, f"dr_{c}_m")
        s = _col(residuals, f"sigma_{c}_obs_m")
        ax.plot(tm[ok], (d / s)[ok], "o", ms=3)
        for k in (-3, 3):
            ax.axhline(k, color="C3", ls=":", lw=0.8)
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_ylabel(f"d{lab} / σ_obs")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel(xlabel)
    fig.suptitle(f"{title}: residual / Stage-5 sigma (conditional on the TLE)")
    paths.append(_save(fig, out_dir / "residual_vs_uncertainty.png", plt))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.semilogy(tm[ok], _col(residuals, "position_mahalanobis_obs_cov")[ok], "o", ms=3,
                label="position D² (Stage-5 cov only; TLE fixed)")
    tot = _col(residuals, "position_mahalanobis_total_cov")
    if np.any(np.isfinite(tot)):
        ax.semilogy(tm[ok], tot[ok], "s", ms=3, label="position D² (+ empirical TLE dispersion)")
    ax.axhline(14.156, color="C3", ls="--", label="χ²₃ 0.9973")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("squared Mahalanobis distance")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax.set_title(f"{title}: Mahalanobis distances (per epoch; epochs correlated)")
    paths.append(_save(fig, out_dir / "mahalanobis.png", plt))

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(tm, _col(residuals, "tle_age_days"), "-")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("TLE age (days; < 0 = back-propagated)")
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}: TLE age, epoch {residuals['tle_epoch_utc'].iloc[0]}")
    paths.append(_save(fig, out_dir / "tle_age.png", plt))

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, (i, j) in zip(axes, ((0, 1), (0, 2), (1, 2))):
        a, b = "xyz"[i], "xyz"[j]
        ax.plot(_col(tle_traj, f"{a}_gcrs_m") / 1e3, _col(tle_traj, f"{b}_gcrs_m") / 1e3, "-", color="C1",
                label="TLE (SGP4)")
        ax.plot(_col(residuals, f"obs_{a}_gcrs_m")[ok] / 1e3, _col(residuals, f"obs_{b}_gcrs_m")[ok] / 1e3, "o",
                ms=3, color="C0", label="observed (Stage 5)")
        ax.set_xlabel(f"{a} GCRS (km)")
        ax.set_ylabel(f"{b} GCRS (km)")
        ax.grid(alpha=0.3)
        ax.set_aspect("equal", adjustable="datalim")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"{title}: observed and TLE trajectories (GCRS)")
    paths.append(_save(fig, out_dir / "trajectory_overlay.png", plt))
    return paths
