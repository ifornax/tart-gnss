"""Stage-4 diagnostic plots (matplotlib, Agg backend), drawn from the saved tables."""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

PLOT_NAMES = ("position_track_ecef", "geometry_condition", "station_residuals", "position_uncertainty",
              "timing_offsets")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _minutes(t, t0) -> np.ndarray:
    return np.asarray((t - t0).dt.total_seconds(), dtype=float) / 60.0


def make_plots(positions, contributions, out_dir: Path, title: str) -> List[str]:
    """Write the five diagnostic PNGs; return their paths (empty if nothing to plot)."""
    if len(positions) == 0:
        return []
    plt = _plt()
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = positions["time_utc"].min()
    tm = _minutes(positions["time_utc"], t0)
    xlabel = f"minutes after {t0:%Y-%m-%d %H:%M:%S} UTC"
    valid = positions["valid"].to_numpy(dtype=bool)
    paths = []

    # x(t), y(t), z(t) with 1-sigma bars, and the geodetic sub-point track.
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, c in zip(axes.flat[:3], "xyz"):
        v = positions[f"{c}_ecef_m"].to_numpy(dtype=float) / 1e3
        s = positions[f"sigma_{c}_m"].to_numpy(dtype=float) / 1e3
        ax.errorbar(tm[~valid], v[~valid], yerr=s[~valid], fmt="o", ms=3, color="0.6", alpha=0.6, label="not valid")
        ax.errorbar(tm[valid], v[valid], yerr=s[valid], fmt="o", ms=3, color="C0", label="valid")
        ax.set_ylabel(f"{c} ECEF (km)")
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.3)
    axes.flat[0].legend(fontsize=8)
    ax = axes.flat[3]
    ax.plot(positions["longitude_deg"][~valid], positions["latitude_deg"][~valid], "o", ms=3, color="0.6")
    ax.plot(positions["longitude_deg"][valid], positions["latitude_deg"][valid], "o", ms=3, color="C0")
    ax.set_xlabel("longitude (deg)")
    ax.set_ylabel("latitude (deg)")
    ax.set_title("geodetic sub-point of the solution", fontsize=9)
    ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: triangulated position (ITRS/ECEF)")
    fig.tight_layout()
    paths.append(_save(fig, out_dir / "position_track_ecef.png", plt))

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    axes[0].semilogy(tm, positions["condition_number"].astype(float), "o", ms=3)
    axes[0].set_ylabel("condition number of sum P_i")
    axes[1].plot(tm, positions["min_crossing_angle_deg"].astype(float), "o", ms=3, label="min")
    axes[1].plot(tm, positions["max_crossing_angle_deg"].astype(float), "s", ms=3, label="max")
    axes[1].set_ylabel("crossing angle (deg)")
    axes[1].set_xlabel(xlabel)
    axes[1].legend(fontsize=8)
    for a in axes:
        a.grid(alpha=0.3)
    fig.suptitle(f"{title}: sight-line geometry")
    fig.tight_layout()
    paths.append(_save(fig, out_dir / "geometry_condition.png", plt))

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    for sid, g in contributions.groupby("station_id"):
        tc = _minutes(g["time_utc"], t0)
        axes[0].plot(tc, g["line_residual_m"].astype(float) / 1e3, "o", ms=3, label=sid)
        axes[1].plot(tc, g["normalised_residual"].astype(float), "o", ms=3, label=sid)
    axes[0].set_ylabel("line residual (km)")
    axes[1].set_ylabel("normalised residual")
    axes[1].set_xlabel(xlabel)
    axes[0].legend(fontsize=8)
    for a in axes:
        a.grid(alpha=0.3)
    fig.suptitle(f"{title}: station line residuals")
    fig.tight_layout()
    paths.append(_save(fig, out_dir / "station_residuals.png", plt))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for c in "xyz":
        ax.plot(tm, positions[f"sigma_{c}_m"].astype(float) / 1e3, "o", ms=3, label=f"sigma_{c}")
    ax.set_ylabel("1-sigma (km)")
    ax.set_xlabel(xlabel)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}: position uncertainty (linear covariance)")
    fig.tight_layout()
    paths.append(_save(fig, out_dir / "position_uncertainty.png", plt))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for sid, g in contributions.groupby("station_id"):
        ax.plot(_minutes(g["time_utc"], t0), g["station_time_offset_s"].astype(float), "o", ms=3, label=sid)
    ax.set_ylabel("station time - reference time (s)")
    ax.set_xlabel(xlabel)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}: station timing offsets")
    fig.tight_layout()
    paths.append(_save(fig, out_dir / "timing_offsets.png", plt))
    return paths


def _save(fig, path: Path, plt) -> str:
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)
