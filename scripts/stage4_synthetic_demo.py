#!/usr/bin/env python
"""Synthetic end-to-end demonstration of Stage 4 (``tart-triangulate``).

One 3-D truth: the TLE-predicted ITRS track of GPS PRN 26 (NORAD 40534) on
2026-06-25, displaced by a fixed 1000 km offset perpendicular to the mean
line of sight (so the displacement is well measured by the network). Three
synthetic TART-like stations at the real station coordinates observe this
truth: complex visibilities are simulated (with per-station timestamp
offsets), ``tart-localise`` measures Az/El at each station, and
``tart-triangulate`` recovers the 3-D position. The TLE is used only to build
the truth and as Stage 3's search-cap centre; Stage 4 never reads it.

Usage::

    python scripts/stage4_synthetic_demo.py [--output-dir docs/examples/stage4]
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from tart_gnss.stage3 import cli as cli3
from tart_gnss.stage3 import synthetic as syn3
from tart_gnss.stage3.loader import MJD_UNIX_OFFSET_S
from tart_gnss.stage4 import cli as cli4
from tart_gnss.stage4.geometry import ecef_to_azel, geodetic_to_ecef
from tart_gnss.stage4.output import read_table
from tart_gnss.stage4.timing import time_seconds

TLE = (
    "1 40534U 15013A   26279.74180529 -.00000040  00000+0  00000+0 0  9991",
    "2 40534  53.1528 199.5676 0114192  41.6722 319.1990  2.00569358 84063",
)
STATIONS = {
    "mauritius": ({"latitude_deg": -20.2587508, "longitude_deg": 57.7591989, "altitude_m": 20.0}, 0.00, 1),
    "madagascar": ({"latitude_deg": -18.8971446, "longitude_deg": 47.5551677, "altitude_m": 1280.0}, 0.25, 2),
    "south_africa": ({"latitude_deg": -33.319413, "longitude_deg": 26.507688, "altitude_m": 642.26}, -0.40, 3),
}
OFFSET_KM = 1000.0
N_TIMES = 30
CADENCE_S = 1.0247


def tle_ecef_m(unix_s: np.ndarray) -> np.ndarray:
    """TLE-predicted ITRS position (m) at Unix times (used only to build the truth)."""
    from skyfield.api import EarthSatellite
    from skyfield.framelib import itrs

    from tart_gnss.propagation import get_timescale

    ts = get_timescale()
    t = ts.from_datetimes([datetime.fromtimestamp(float(u), tz=timezone.utc) for u in np.atleast_1d(unix_s)])
    return EarthSatellite(*TLE, "target", ts).at(t).frame_xyz(itrs).m.T


def offset_vector_m(t0_unix: float) -> np.ndarray:
    """A fixed OFFSET_KM vector perpendicular to the mean line of sight at t0."""
    x = tle_ecef_m(np.array([t0_unix]))[0]
    los = np.mean([(x - geodetic_to_ecef(*site.values())) / np.linalg.norm(x - geodetic_to_ecef(*site.values()))
                   for site, _, _ in STATIONS.values()], axis=0)
    los /= np.linalg.norm(los)
    d = np.cross(los, [0.0, 0.0, 1.0])
    return OFFSET_KM * 1e3 * d / np.linalg.norm(d)


def build(work: Path):
    """Simulate visibilities from the 3-D truth; return (Stage-2 manifest path, truth function, offset)."""
    t0 = syn3.utc_to_mjd_seconds(2026, 6, 25, 9, 14, 0.0)
    t0_unix = t0 - MJD_UNIX_OFFSET_S
    off = offset_vector_m(t0_unix)
    entries = {}
    for sid, (site, dt, seed) in STATIONS.items():
        t = t0 + dt + np.arange(N_TIMES) * CADENCE_S
        x = tle_ecef_m(t - MJD_UNIX_OFFSET_S) + off
        r = geodetic_to_ecef(*site.values())
        azel = np.array([ecef_to_azel(xi, r, site["latitude_deg"], site["longitude_deg"])[:2] for xi in x])
        obs = syn3.simulate_observation(
            syn3.tart_like_layout(seed=seed), site["latitude_deg"], site["longitude_deg"], site["altitude_m"], t,
            azel[:, 0], azel[:, 1], amplitude=0.3 + 0.2j, noise_sigma=0.5, flag_fraction=0.05,
            rng=np.random.default_rng(seed))
        (work / sid / "ms").mkdir(parents=True)
        syn3.write_ms(obs, work / sid / "ms" / f"{sid}.ms")
        entries[sid] = {**site, "ms_path": f"{sid}/ms/{sid}.ms", "n_antennas": 24, "tart_id": f"synthetic-{sid}",
                        "first_time_utc": "2026-06-25T09:14:00Z", "last_time_utc": "2026-06-25T09:14:31Z"}
    manifest = syn3.write_manifest(work / "manifest.json", entries, norad_id=40534, tle_lines=TLE)
    return manifest, (lambda u: tle_ecef_m(np.array([u]))[0] + off), off


def report(out4: Path, truth, off) -> str:
    """Text report: truth vs recovered position, errors, sigmas, residuals, geometry."""
    pos = read_table(out4 / "positions.csv")
    con = read_table(out4 / "station_contributions.csv")
    t_unix = time_seconds(pos["time_utc"])
    lines = []
    errs, tle_sep, pulls, along_off, along_los = [], [], [], [], []
    u = off / np.linalg.norm(off)
    for i, r in pos.iterrows():
        x_true = truth(t_unix[i])
        x_hat = r[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy(dtype=float)
        sig = r[["sigma_x_m", "sigma_y_m", "sigma_z_m"]].to_numpy(dtype=float)
        errs.append(np.linalg.norm(x_hat - x_true))
        tle_sep.append(np.linalg.norm(x_hat - (x_true - off)))
        pulls.append((x_hat - x_true) / sig)
        along_off.append(float((x_hat - (x_true - off)) @ u))
        c_i = con[con["time_utc"] == r["time_utc"]]
        los = c_i[["sight_x", "sight_y", "sight_z"]].to_numpy(dtype=float).mean(axis=0)
        along_los.append(float((x_hat - x_true) @ (los / np.linalg.norm(los))))
    errs, tle_sep, pulls = np.array(errs), np.array(tle_sep), np.array(pulls)
    k = len(pos) // 2
    r = pos.iloc[k]
    x_true = truth(t_unix[k])
    x_hat = r[["x_ecef_m", "y_ecef_m", "z_ecef_m"]].to_numpy(dtype=float)
    c = con[con["time_utc"] == r["time_utc"]]
    lines.append(f"Stage-4 synthetic demo: truth = TLE track + {OFFSET_KM:.0f} km offset perpendicular to the mean LOS")
    lines.append(f"middle epoch {r['time_utc']:%Y-%m-%dT%H:%M:%S.%f}Z ({int(r['n_stations'])} stations, "
                 f"flags: {';'.join(r['quality_flags']) or 'GOOD'})")
    lines.append(f"{'':12s} {'x (km)':>14s} {'y (km)':>14s} {'z (km)':>14s}")
    lines.append(f"{'true':12s} " + " ".join(f"{v / 1e3:14.3f}" for v in x_true))
    lines.append(f"{'recovered':12s} " + " ".join(f"{v / 1e3:14.3f}" for v in x_hat))
    lines.append(f"{'error':12s} " + " ".join(f"{v / 1e3:14.3f}" for v in x_hat - x_true))
    lines.append(f"{'sigma':12s} " + " ".join(f"{r[f'sigma_{a}_m'] / 1e3:14.3f}" for a in 'xyz'))
    lines.append(f"3-D error {errs[k] / 1e3:.1f} km; distance from the TLE position {tle_sep[k] / 1e3:.1f} km "
                 f"(injected offset {OFFSET_KM:.0f} km)")
    lines.append(f"condition number {r['condition_number']:.1f}; crossing angles {r['min_crossing_angle_deg']:.2f}-"
                 f"{r['max_crossing_angle_deg']:.2f} deg; reduced chi2 {r['reduced_chi2']:.2f}; "
                 f"max time offset {r['max_time_offset_s']:.3f} s")
    lines.append("station        offset_s   range_km  line_resid_km  norm_resid  sigma_theta_deg")
    for _, s in c.iterrows():
        lines.append(f"{s['station_id']:13s} {s['station_time_offset_s']:8.3f} {s['range_m'] / 1e3:10.1f} "
                     f"{s['line_residual_m'] / 1e3:14.2f} {s['normalised_residual']:11.2f} {s['sigma_theta_deg']:16.3f}")
    lines.append(f"all {len(pos)} epochs: valid {int(pos['valid'].sum())}; median 3-D error {np.median(errs) / 1e3:.1f} km "
                 f"(max {errs.max() / 1e3:.1f}); median distance from TLE {np.median(tle_sep) / 1e3:.1f} km; "
                 f"pull std (x, y, z) = ({', '.join(f'{v:.2f}' for v in pulls.std(axis=0))})")
    along_off, along_los = np.array(along_off), np.array(along_los)
    lines.append(f"recovered displacement from the TLE along the injected offset direction: median "
                 f"{np.median(along_off) / 1e3:.1f} km (injected {OFFSET_KM:.0f} km; scatter "
                 f"{np.std(along_off) / 1e3:.1f} km); error along the mean line of sight: median "
                 f"{np.median(np.abs(along_los)) / 1e3:.1f} km, transverse: median "
                 f"{np.median(np.sqrt(np.maximum(errs ** 2 - along_los ** 2, 0))) / 1e3:.1f} km")
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output-dir", type=Path, default=Path("docs/examples/stage4"))
    args = p.parse_args()
    with tempfile.TemporaryDirectory(prefix="stage4_demo_") as tmp:
        work = Path(tmp)
        manifest, truth, off = build(work)
        s3_root = args.output_dir / "stage3_input"
        code3 = cli3.main(["--manifest", str(manifest), "--all-stations", "--time-bin", "1", "--no-plots",
                           "--output-dir", str(s3_root), "--data-root", str(work)])
    s3_manifest = s3_root / "20260625" / "NORAD_40534" / "stage3_manifest.json"
    code4 = cli4.main(["--stage3-manifest", str(s3_manifest), "--output-dir", str(args.output_dir),
                       "--bootstrap", "100"])
    out4 = args.output_dir / "20260625" / "NORAD_40534"
    text = report(out4, truth, off)
    text += f"\n\ntart-localise exit {code3}; tart-triangulate exit {code4}; outputs in {out4}\n"
    (args.output_dir / "demo_report.txt").write_text(text, encoding="utf-8")
    truth_doc = {"offset_ecef_m": off.tolist(), "construction": "TLE ITRS track (skyfield SGP4) + offset"}
    (args.output_dir / "truth.json").write_text(json.dumps(truth_doc, indent=2) + "\n", encoding="utf-8")
    print(text)
    return code4


if __name__ == "__main__":
    sys.exit(main())
