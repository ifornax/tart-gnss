#!/usr/bin/env python3
"""Generate the figures of the Stage-1 technical report from a run's outputs.

Usage:
    python scripts/make_report_figures.py results/common_gnss_20261001.json docs/report/figures \
        [data/observations/20260625/NORAD_40534/manifest.json]

The optional Stage-2 manifest (with its native files kept, --keep-native)
adds the Stage-2 timing figure.

Requires the optional plotting dependency (``pip install .[plot]``).
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from tart_gnss import propagation  # noqa: E402
from tart_gnss.stations import Station  # noqa: E402
from tart_gnss.tle import TLERecord, tle_epoch_utc  # noqa: E402
from tart_gnss.visibility import common_visibility  # noqa: E402

# Reference categorical palette (validated, slots 1-4), text and grid inks.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CONSTELLATION_ORDER = ["GPS", "Galileo", "BeiDou", "GLONASS"]

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "lines.linewidth": 1.6, "savefig.dpi": 200,
})


def parse(t: str) -> datetime:
    """Parse an ISO UTC string into a naive datetime implicitly in UTC (for plotting)."""
    return datetime.fromisoformat(t.replace("Z", "+00:00")).replace(tzinfo=None)


def record_from_window(w: dict) -> TLERecord:
    l1, l2 = w["tle"]["line1"], w["tle"]["line2"]
    return TLERecord(w["satellite"], w["norad_id"], l1, l2, w["constellation"],
                     tle_epoch_utc(l1, l2), w["tle"]["source"])


def elevation_panels(doc: dict, out: Path) -> None:
    """Elevation vs time at each station for the top window of each constellation."""
    stations = [Station(**{k: s.get(k) for k in Station.__dataclass_fields__}) for s in doc["stations"]]
    start, end = parse(doc["run"]["start_utc"]), parse(doc["run"]["end_utc"])
    thr = doc["run"]["min_elevation_deg"]
    picks = []
    for c in CONSTELLATION_ORDER:
        ws = [w for w in doc["windows"] if w["constellation"] == c]
        if ws:
            picks.append(ws[0])
    offsets = np.arange(0, (end - start).total_seconds() + 1, 60.0)
    times = propagation.times_from_offsets(start.replace(tzinfo=timezone.utc), offsets)
    # naive datetimes in UTC so that matplotlib locators and formatters agree
    dts = [t.replace(tzinfo=None) for t in times.utc_datetime()]
    pos = [propagation.station_position(s) for s in stations]

    fig, axes = plt.subplots(len(picks), 1, figsize=(7.0, 1.9 * len(picks)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, w in zip(axes, picks):
        sat = propagation.make_satellite(record_from_window(w))
        el = propagation.elevations(sat, pos, times)
        common = common_visibility(el, thr)
        ax.fill_between(dts, -90, 90, where=common, color=SERIES[0], alpha=0.08, linewidth=0,
                        label="common window")
        for i, s in enumerate(stations):
            ax.plot(dts, el[i], color=SERIES[i], label=s.id)
        ax.axhline(thr, color=INK2, linestyle="--", linewidth=1.0)
        ax.set_ylim(0, 90)
        ax.set_ylabel("Elevation (deg)")
        ax.set_title(f"{w['satellite'].strip()}  (NORAD {w['norad_id']}, {w['constellation']})",
                     loc="left", fontsize=9, color=INK)
    axes[0].legend(ncol=4, loc="upper right", bbox_to_anchor=(1.0, 1.42), fontsize=8)
    axes[-1].xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 3)))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axes[-1].set_xlim(dts[0], dts[-1])
    axes[-1].set_xlabel(f"UTC on {start:%Y-%m-%d}")
    fig.tight_layout()
    fig.savefig(out / "elevation_examples.pdf")
    plt.close(fig)


def timeline(doc: dict, out: Path, n: int = 40) -> None:
    """Gantt chart of the n longest common windows, coloured by constellation."""
    ws = doc["windows"][:n]
    start, end = parse(doc["run"]["start_utc"]), parse(doc["run"]["end_utc"])
    fig, ax = plt.subplots(figsize=(7.0, 0.16 * len(ws) + 1.0))
    for i, w in enumerate(ws):
        a, b = parse(w["start_utc"]), parse(w["end_utc"])
        color = SERIES[CONSTELLATION_ORDER.index(w["constellation"])] if w["constellation"] in CONSTELLATION_ORDER else INK2
        ax.barh(i, mdates.date2num(b) - mdates.date2num(a), left=mdates.date2num(a),
                height=0.7, color=color, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(ws)))
    ax.set_yticklabels([f"{w['rank']}. {w['satellite'].strip()}" for w in ws], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(mdates.date2num(start), mdates.date2num(end))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 3)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel(f"UTC on {start:%Y-%m-%d}")
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=SERIES[i]) for i in range(4)]
    ax.legend(handles, CONSTELLATION_ORDER, ncol=4, loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "window_timeline.pdf")
    plt.close(fig)


def skyplot(doc: dict, out: Path) -> None:
    """Polar az/el track of the top-ranked non-geostationary window at each station."""
    stations = [Station(**{k: s.get(k) for k in Station.__dataclass_fields__}) for s in doc["stations"]]
    w = next(w for w in doc["windows"] if w["constellation"] in ("GPS", "Galileo", "GLONASS"))
    a, b = parse(w["start_utc"]), parse(w["end_utc"])
    run_start = parse(doc["run"]["start_utc"])
    offs = np.linspace((a - run_start).total_seconds(), (b - run_start).total_seconds(), 300)
    times = propagation.times_from_offsets(run_start.replace(tzinfo=timezone.utc), offs)
    sat = propagation.make_satellite(record_from_window(w))
    fig = plt.figure(figsize=(4.2, 4.2))
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    for i, s in enumerate(stations):
        az, el, _, _ = propagation.topocentric(sat, propagation.station_position(s), times)
        ax.plot(np.radians(az), 90 - el, color=SERIES[i], label=s.id)
        ax.plot(np.radians(az[0]), 90 - el[0], "o", color=SERIES[i], markersize=4,
                markeredgecolor="white", markeredgewidth=1)
    thr = doc["run"]["min_elevation_deg"]
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.full(200, 90 - thr), "--", color=INK2, linewidth=1)
    ax.set_rlim(0, 90)
    ax.set_rticks([30, 60, 90 - thr])
    ax.set_yticklabels(["60°", "30°", f"{thr:g}°"], fontsize=7, color=INK2)
    ax.set_rlabel_position(22.5)
    ax.set_title(f"{w['satellite'].strip()}, rank {w['rank']}", fontsize=9, color=INK)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "skyplot.pdf")
    plt.close(fig)


def stage2_timing(manifest_path: Path, out: Path) -> None:
    """Sample timelines of the three stations and nearest-sample offsets (Stage 2)."""
    import h5py

    from tart_gnss.stage2.synchronise import nearest_separation

    m = json.loads(manifest_path.read_text())
    ids = list(m["stations"])
    ts = {}
    for sid in ids:
        stamps = []
        for f in m["stations"][sid]["raw_files"]:
            if f.get("validation") == "valid" and Path(f["local_path"]).exists():
                with h5py.File(f["local_path"], "r") as h:
                    stamps += [datetime.fromisoformat(t.decode()) for t in h["timestamp"][()]]
        ts[sid] = sorted(stamps)
    sync = m["synchronisation"]
    c0, c1 = parse(sync["actual_common_start_utc"]), parse(sync["actual_common_end_utc"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw={"width_ratios": [1.6, 1]})
    ax1.axvspan(c0, c1, color=SERIES[0], alpha=0.08, linewidth=0, label="three-station overlap")
    for i, sid in enumerate(ids):
        x = [t.replace(tzinfo=None) for t in ts[sid]]
        ax1.plot(x, np.full(len(x), len(ids) - i), "|", color=SERIES[i], markersize=8, label=sid)
    ax1.set_yticks(range(1, len(ids) + 1))
    ax1.set_yticklabels(ids[::-1])
    ax1.set_ylim(0.4, len(ids) + 0.6)
    ax1.grid(axis="y", visible=False)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax1.set_xlabel(f"UTC on {c0:%Y-%m-%d}")
    ax1.set_title("(a) native samples (~1 s cadence)", loc="left", fontsize=9)
    from itertools import combinations
    for k, (a, b) in enumerate(combinations(ids, 2)):
        sel_a = [t for t in ts[a] if c0 <= t.replace(tzinfo=None) <= c1]
        dt = nearest_separation(sel_a, ts[b])
        ax2.hist(dt, bins=np.linspace(-0.05, 0.55, 61), color=SERIES[k], alpha=0.85, label=f"{b} - {a}")
    ax2.set_xlabel("nearest-sample offset (s)")
    ax2.set_ylabel("samples")
    ax2.set_title("(b) cross-station timing", loc="left", fontsize=9)
    ax2.legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(out / "stage2_timing.pdf")
    plt.close(fig)


def main() -> int:
    doc = json.loads(Path(sys.argv[1]).read_text())
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    elevation_panels(doc, out)
    timeline(doc, out)
    skyplot(doc, out)
    if len(sys.argv) > 3:
        stage2_timing(Path(sys.argv[3]), out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
