"""Stage-5 diagnostic plots (matplotlib, Agg backend), drawn from the saved tables."""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

PLOT_NAMES = ("position_track", "velocity_track", "uncertainty_track", "innovations", "nis",
              "stage4_vs_smoothed", "covariance_eigenvalues")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _minutes(t, t0) -> np.ndarray:
    return np.asarray((t - t0).dt.total_seconds(), dtype=float) / 60.0


def _col(df, name) -> np.ndarray:
    return df[name].to_numpy(dtype=float)


def _save(fig, path: Path, plt) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)


def make_plots(trajectory, innovations, out_dir: Path, title: str) -> List[str]:
    """Write the seven diagnostic PNGs; return their paths."""
    from tart_gnss.stage5.output import state_covariance

    if len(trajectory) == 0:
        return []
    plt = _plt()
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = trajectory["time_utc"].min()
    tm = _minutes(trajectory["time_utc"], t0)
    xlabel = f"minutes after {t0:%Y-%m-%d %H:%M:%S} UTC"
    used = trajectory["measurement_accepted"].fillna(False).to_numpy(dtype=bool)
    paths = []

    # 1. Stage-4 measurements vs filtered and smoothed x, y, z (GCRS); residuals keep the km scale visible.
    fig, axes = plt.subplots(3, 2, figsize=(12, 8.5), sharex=True)
    for i, c in enumerate("xyz"):
        meas = _col(trajectory, f"stage4_{c}_gcrs_m") / 1e3
        sig = _col(trajectory, f"stage4_sigma_{c}_gcrs_m") / 1e3
        filt = _col(trajectory, f"{c}_gcrs_filtered_m") / 1e3
        smo = _col(trajectory, f"{c}_gcrs_smoothed_m") / 1e3
        ax = axes[i, 0]
        ax.errorbar(tm[~used], meas[~used], yerr=sig[~used], fmt="x", ms=4, color="0.6", label="Stage 4 (not used)")
        ax.errorbar(tm[used], meas[used], yerr=sig[used], fmt="o", ms=3, color="C0", alpha=0.6, label="Stage 4")
        ax.plot(tm, filt, "-", color="C1", label="filtered")
        ax.plot(tm, smo, "-", color="C3", label="smoothed")
        ax.set_ylabel(f"{c} GCRS (km)")
        ax = axes[i, 1]
        ref = smo if np.any(np.isfinite(smo)) else filt
        ax.errorbar(tm[used], (meas - ref)[used], yerr=sig[used], fmt="o", ms=3, color="C0", alpha=0.6)
        ax.plot(tm[~used], (meas - ref)[~used], "x", ms=4, color="0.6")
        ax.plot(tm, filt - ref, "-", color="C1")
        s = _col(trajectory, f"sigma_{c}_smoothed_m") / 1e3
        ax.fill_between(tm, -s, s, color="C3", alpha=0.25, label="smoothed 1-sigma")
        ax.axhline(0, color="C3", lw=1)
        ax.set_ylabel(f"{c} minus smoothed (km)")
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel(xlabel)
    for ax in axes.flat:
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: position track (GCRS)")
    paths.append(_save(fig, out_dir / "position_track.png", plt))

    # 2. Velocity and speed.
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, c in zip(axes.flat[:3], ("vx", "vy", "vz")):
        v = _col(trajectory, f"{c}_gcrs_smoothed_m_s") / 1e3
        s = _col(trajectory, f"sigma_{c}_smoothed_m_s") / 1e3
        ax.plot(tm, _col(trajectory, f"{c}_gcrs_filtered_m_s") / 1e3, "-", color="C1", label="filtered")
        ax.plot(tm, v, "-", color="C3", label="smoothed")
        ax.fill_between(tm, v - s, v + s, color="C3", alpha=0.25)
        ax.set_ylabel(f"{c} GCRS (km/s)")
    axes.flat[0].legend(fontsize=8)
    axes.flat[3].plot(tm, _col(trajectory, "speed_smoothed_m_s") / 1e3, "-", color="C3")
    axes.flat[3].set_ylabel("speed (km/s)")
    for ax in axes.flat:
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel(xlabel)
    fig.suptitle(f"{title}: velocity (GCRS)")
    paths.append(_save(fig, out_dir / "velocity_track.png", plt))

    # 3. Uncertainty vs time.
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for c, col in zip("xyz", ("C0", "C1", "C2")):
        axes[0].semilogy(tm, _col(trajectory, f"sigma_{c}_filtered_m") / 1e3, "--", color=col, label=f"{c} filtered")
        axes[0].semilogy(tm, _col(trajectory, f"sigma_{c}_smoothed_m") / 1e3, "-", color=col, label=f"{c} smoothed")
        axes[1].semilogy(tm, _col(trajectory, f"sigma_v{c}_smoothed_m_s"), "-", color=col, label=f"v{c} smoothed")
    axes[0].set_ylabel("position 1-sigma (km)")
    axes[1].set_ylabel("velocity 1-sigma (m/s)")
    axes[1].set_xlabel(xlabel)
    for ax in axes:
        ax.legend(fontsize=7, ncol=3)
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(f"{title}: state uncertainty (GCRS)")
    paths.append(_save(fig, out_dir / "uncertainty_track.png", plt))

    ti = _minutes(innovations["time_utc"], t0) if len(innovations) else np.zeros(0)
    acc = innovations["accepted"].to_numpy(dtype=bool) if len(innovations) else np.zeros(0, bool)
    # 4. Innovation components with +-sqrt(S).
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, c in zip(axes, "xyz"):
        if len(innovations):
            nu = _col(innovations, f"innovation_{c}_m") / 1e3
            s = np.sqrt(_col(innovations, f"innovation_cov_{c}{c}")) / 1e3
            ax.fill_between(ti, -s, s, color="C0", alpha=0.2, label="+-sqrt(S)")
            ax.plot(ti[acc], nu[acc], "o", ms=3, color="C0", label="assimilated")
            ax.plot(ti[~acc], nu[~acc], "x", ms=5, color="C3", label="not assimilated")
        ax.set_ylabel(f"innovation {c} (km)")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    axes[-1].set_xlabel(xlabel)
    fig.suptitle(f"{title}: innovations (GCRS)")
    paths.append(_save(fig, out_dir / "innovations.png", plt))

    # 5. NIS with the gate.
    fig, ax = plt.subplots(figsize=(9, 4.5))
    if len(innovations):
        nis = _col(innovations, "nis")
        ax.semilogy(ti[acc], nis[acc], "o", ms=3, label="assimilated")
        ax.semilogy(ti[~acc], nis[~acc], "x", ms=5, color="C3", label="not assimilated")
        ax.axhline(float(innovations["gate_threshold"].iloc[0]), color="C3", ls="--", label="gate")
        ax.axhline(3.0, color="0.4", ls=":", label="E[chi2_3] = 3")
    ax.set_ylabel("NIS")
    ax.set_xlabel(xlabel)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax.set_title(f"{title}: normalised innovation squared")
    paths.append(_save(fig, out_dir / "nis.png", plt))

    # 6. Stage-4 residuals about the smoothed trajectory (ECEF), raw and normalised by the Stage-4 sigma.
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for c in "xyz":
        d = _col(trajectory, f"stage4_{c}_ecef_m") - _col(trajectory, f"{c}_ecef_smoothed_m")
        axes[0].plot(tm[used], d[used] / 1e3, "o", ms=3, label=c)
        axes[0].plot(tm[~used], d[~used] / 1e3, "x", ms=4, color="0.6")
    gs = np.array([_col(trajectory, f"stage4_{c}_gcrs_m") - _col(trajectory, f"{c}_gcrs_smoothed_m") for c in "xyz"]).T
    norm = np.full(len(trajectory), np.nan)
    for i, row in trajectory.iterrows():
        if not np.all(np.isfinite(gs[i])):
            continue
        s4 = np.array([row[f"stage4_sigma_{c}_gcrs_m"] for c in "xyz"], dtype=float)
        if np.all(np.isfinite(s4)) and np.all(s4 > 0):
            norm[i] = np.sqrt(np.mean((gs[i] / s4) ** 2))
    axes[0].set_ylabel("Stage 4 - smoothed, ECEF (km)")
    axes[0].legend(fontsize=8)
    axes[1].plot(tm[used], norm[used], "o", ms=3)
    axes[1].axhline(1.0, color="0.4", ls=":")
    axes[1].set_ylabel("RMS residual / Stage-4 sigma")
    axes[1].set_xlabel(xlabel)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: Stage-4 positions about the smoothed trajectory")
    paths.append(_save(fig, out_dir / "stage4_vs_smoothed.png", plt))

    # 7. Covariance eigenvalues (position block, smoothed and filtered).
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for kind, ls in (("filtered", "--"), ("smoothed", "-")):
        ev_p, ev_v = [], []
        for _, row in trajectory.iterrows():
            try:
                p = state_covariance(row, kind)
            except (KeyError, ValueError):
                p = np.full((6, 6), np.nan)
            if np.all(np.isfinite(p)):
                ev_p.append(np.sqrt(np.clip(np.linalg.eigvalsh(p[:3, :3]), 0, None)) / 1e3)
                ev_v.append(np.sqrt(np.clip(np.linalg.eigvalsh(p[3:, 3:]), 0, None)))
            else:
                ev_p.append(np.full(3, np.nan))
                ev_v.append(np.full(3, np.nan))
        ev_p, ev_v = np.array(ev_p), np.array(ev_v)
        for i in range(3):
            axes[0].semilogy(tm, ev_p[:, i], ls, color=f"C{i}", label=f"{kind} axis {i + 1}")
            axes[1].semilogy(tm, ev_v[:, i], ls, color=f"C{i}")
    axes[0].set_ylabel("sqrt eigenvalue, position (km)")
    axes[1].set_ylabel("sqrt eigenvalue, velocity (m/s)")
    axes[1].set_xlabel(xlabel)
    axes[0].legend(fontsize=7, ncol=2)
    for ax in axes:
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(f"{title}: covariance principal axes")
    paths.append(_save(fig, out_dir / "covariance_eigenvalues.png", plt))
    return paths
