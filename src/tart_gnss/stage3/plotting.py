"""Stage-3 diagnostic figures, generated from the saved direction table.

For each station:

* ``<station>_track.png``: azimuth and elevation vs time (visibility estimate
  with 1-sigma error bars, image estimate, TLE reference line);
* ``<station>_sky.png``: polar sky track (azimuth clockwise from North, radius
  = zenith angle);
* ``<station>_vis_image_offset.png``: visibility-image separation vs time;
* ``<station>_fit_quality.png``: SNR, reduced chi2, baseline count, flagged fraction.

Usable estimates are drawn filled, others (blocking quality flags) hollow.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np

from tart_gnss.stage3.geometry import wrap_azimuth_near
from tart_gnss.stage3.output import read_directions


def _pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _unwrap(az: np.ndarray, ref: float) -> np.ndarray:
    return wrap_azimuth_near(az, ref)


def make_station_plots(df, station_id: str, out_dir: Path) -> List[str]:
    """Write the four diagnostic figures for one station; return their paths."""
    plt = _pyplot()
    d = df[df["station_id"] == station_id].sort_values("time_utc")
    if d.empty:
        return []
    out_dir.mkdir(parents=True, exist_ok=True)
    t = d["time_utc"]
    usable = d["usable"].astype(bool).to_numpy()
    tle_ok = d["az_tle_deg"].notna().any()
    ref = float(np.nanmedian(d["az_tle_deg"])) if tle_ok else float(np.nanmedian(d["az_obs_deg"]))
    az = _unwrap(d["az_obs_deg"].to_numpy(float), ref)
    paths = []

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    for ax, val, sig, img, tle, lab in (
        (axes[0], az, d["sigma_az_deg"], _unwrap(d["az_image_deg"].to_numpy(float), ref),
         _unwrap(d["az_tle_deg"].to_numpy(float), ref), "Azimuth (deg)"),
        (axes[1], d["el_obs_deg"].to_numpy(float), d["sigma_el_deg"], d["el_image_deg"].to_numpy(float),
         d["el_tle_deg"].to_numpy(float), "Elevation (deg)"),
    ):
        sig = sig.to_numpy(float)
        ax.errorbar(t[usable], val[usable], yerr=sig[usable], fmt="o", ms=3, lw=0.8, color="C0",
                    label="visibility fit (usable)")
        if (~usable).any():
            ax.plot(t[~usable], val[~usable], "o", ms=3, mfc="none", color="C0", alpha=0.6,
                    label="visibility fit (flagged)")
        if np.isfinite(img).any():
            ax.plot(t, img, "x", ms=4, color="C1", label="image peak")
        if np.isfinite(tle).any():
            ax.plot(t, tle, "-", color="k", lw=1.2, label="TLE prediction (reference)")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc="best")
    axes[1].set_xlabel("Time (UTC)")
    fig.suptitle(f"{station_id}: measured direction vs time")
    fig.autofmt_xdate()
    paths.append(_save(fig, out_dir / f"{station_id}_track.png", plt))

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    za = 90.0 - d["el_obs_deg"].to_numpy(float)
    th = np.radians(d["az_obs_deg"].to_numpy(float))
    ax.plot(th[usable], za[usable], "o", ms=3, color="C0", label="visibility fit (usable)")
    if (~usable).any():
        ax.plot(th[~usable], za[~usable], "o", ms=3, mfc="none", color="C0", alpha=0.6, label="flagged")
    if d["az_image_deg"].notna().any():
        ax.plot(np.radians(d["az_image_deg"].to_numpy(float)), 90 - d["el_image_deg"].to_numpy(float),
                "x", ms=4, color="C1", label="image peak")
    if tle_ok:
        ax.plot(np.radians(d["az_tle_deg"].to_numpy(float)), 90 - d["el_tle_deg"].to_numpy(float),
                "-", color="k", label="TLE prediction (reference)")
    ax.set_rlim(0, 90)
    ax.set_rticks([15, 30, 45, 60, 75, 90])
    ax.set_yticklabels(["75", "60", "45", "30", "15", "0"], fontsize=7)
    ax.legend(fontsize=7, loc="lower left", bbox_to_anchor=(-0.1, -0.12))
    ax.set_title(f"{station_id}: sky track (radius = zenith angle)", fontsize=10)
    paths.append(_save(fig, out_dir / f"{station_id}_sky.png", plt))

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(t, d["vis_image_separation_deg"], "o-", ms=3, lw=0.6)
    ax.set_ylabel("Visibility-image separation (deg)")
    ax.set_xlabel("Time (UTC)")
    ax.grid(alpha=0.3)
    ax.set_title(f"{station_id}: visibility vs image estimator")
    fig.autofmt_xdate()
    paths.append(_save(fig, out_dir / f"{station_id}_vis_image_offset.png", plt))

    fig, axes = plt.subplots(4, 1, figsize=(9, 8), sharex=True)
    for ax, col, lab in zip(axes, ("snr", "reduced_chi2", "n_baselines_used", "flagged_fraction"),
                            ("SNR", "reduced chi2", "baselines used", "flagged fraction")):
        ax.plot(t, d[col], "o-", ms=2.5, lw=0.6)
        ax.set_ylabel(lab, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0].set_yscale("log")
    axes[-1].set_xlabel("Time (UTC)")
    fig.suptitle(f"{station_id}: fit quality")
    fig.autofmt_xdate()
    paths.append(_save(fig, out_dir / f"{station_id}_fit_quality.png", plt))
    return paths


def _save(fig, path: Path, plt) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)


def make_all_plots(table_path, out_dir: Path) -> Dict[str, List[str]]:
    """Read a saved Stage-3 table and plot every station in it."""
    df = read_directions(table_path)
    return {sid: make_station_plots(df, sid, out_dir) for sid in df["station_id"].unique()}
