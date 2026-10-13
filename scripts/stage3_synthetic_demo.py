#!/usr/bin/env python
"""Synthetic end-to-end demonstration of Stage 3 (``tart-localise``).

Three synthetic TART-like stations (24 antennas each; real station
coordinates, SYNTHETIC visibilities) observe a point source that is
deliberately displaced from the real TLE track of GPS PRN 26 (NORAD 40534)
on 2026-06-25. The demo writes synthetic Measurement Sets and a Stage-2-style
manifest to a temporary directory, runs ``tart-localise`` on them, and
reports injected vs recovered directions. The recovered directions must
follow the injected source, not the TLE/search centre.

Usage::

    python scripts/stage3_synthetic_demo.py [--output-dir docs/examples/stage3]
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

from tart_gnss.stage3 import cli
from tart_gnss.stage3 import synthetic as syn
from tart_gnss.stage3.geometry import angular_separation_deg, azel_to_unit
from tart_gnss.stage3.loader import MJD_UNIX_OFFSET_S
from tart_gnss.stage3.output import read_directions

TLE = (
    "1 40534U 15013A   26279.74180529 -.00000040  00000+0  00000+0 0  9991",
    "2 40534  53.1528 199.5676 0114192  41.6722 319.1990  2.00569358 84063",
)
#: Real station coordinates (Stage-2 manifest); offsets (dAz, dEl) of the injected source from the TLE track.
STATIONS = {
    "mauritius": ({"latitude_deg": -20.2587508, "longitude_deg": 57.7591989, "altitude_m": 20.0}, (6.0, 2.0), 1),
    "madagascar": ({"latitude_deg": -18.8971446, "longitude_deg": 47.5551677, "altitude_m": 1280.0}, (-4.0, 3.0), 2),
    "south_africa": ({"latitude_deg": -33.319413, "longitude_deg": 26.507688, "altitude_m": 642.26}, (10.0, -4.0), 3),
}
N_TIMES = 20


def build(work: Path) -> Path:
    """Write the synthetic Measurement Sets and manifest; return the manifest path."""
    t = syn.utc_to_mjd_seconds(2026, 6, 25, 9, 14, 0.0) + np.arange(N_TIMES) * 1.0247
    entries = {}
    truth = {}
    for sid, (site, (daz, delv), seed) in STATIONS.items():
        tle = cli.TleReference(*TLE, site["latitude_deg"], site["longitude_deg"], site["altitude_m"])
        az_t, el_t, _ = tle.predict(t - MJD_UNIX_OFFSET_S)
        obs = syn.simulate_observation(
            syn.tart_like_layout(seed=seed), site["latitude_deg"], site["longitude_deg"], site["altitude_m"], t,
            az_t + daz, el_t + delv, amplitude=0.3 + 0.2j, noise_sigma=0.5, flag_fraction=0.05,
            rng=np.random.default_rng(seed))
        (work / sid / "ms").mkdir(parents=True)
        syn.write_ms(obs, work / sid / "ms" / f"{sid}.ms")
        entries[sid] = {**site, "ms_path": f"{sid}/ms/{sid}.ms", "n_antennas": 24, "tart_id": f"synthetic-{sid}",
                        "first_time_utc": "2026-06-25T09:14:00Z", "last_time_utc": "2026-06-25T09:14:20Z"}
        truth[sid] = {"az": (az_t + daz).tolist(), "el": (el_t + delv).tolist(), "t_mjd_s": t.tolist()}
    (work / "truth.json").write_text(json.dumps(truth))
    return syn.write_manifest(work / "manifest.json", entries, norad_id=40534, tle_lines=TLE)


def report(out: Path, work: Path) -> None:
    """Print injected vs recovered directions per station."""
    truth = json.loads((work / "truth.json").read_text())
    df = read_directions(out / "directions.csv")
    print(f"\n{'station':13s} {'injected az/el':>17s} {'recovered az/el':>17s} {'err':>7s} {'sigma az/el':>13s}"
          f" {'image az/el':>17s} {'vis-img':>8s} {'TLE sep':>8s} {'snr':>6s}  flags")
    for sid in STATIONS:
        d = df[df["station_id"] == sid].reset_index(drop=True)
        errs = []
        for _, r in d.iterrows():
            k = int(np.argmin(np.abs(np.array(truth[sid]["t_mjd_s"]) - r["time_mjd"] * 86400)))
            errs.append(float(angular_separation_deg(azel_to_unit(r["az_obs_deg"], r["el_obs_deg"]),
                                                     azel_to_unit(truth[sid]["az"][k], truth[sid]["el"][k]))))
        r = d.iloc[len(d) // 2]
        k = int(np.argmin(np.abs(np.array(truth[sid]["t_mjd_s"]) - r["time_mjd"] * 86400)))
        print(f"{sid:13s} {truth[sid]['az'][k]:8.3f} {truth[sid]['el'][k]:7.3f}  "
              f"{r['az_obs_deg']:8.3f} {r['el_obs_deg']:7.3f}  {errs[len(d) // 2]:6.3f}  "
              f"{r['sigma_az_deg']:5.3f} {r['sigma_el_deg']:5.3f}  {r['az_image_deg']:8.3f} {r['el_image_deg']:7.3f}  "
              f"{r['vis_image_separation_deg']:7.3f} {r['tle_separation_deg']:8.3f} {r['snr']:6.1f}  "
              f"{';'.join(r['quality_flags']) or 'GOOD'}")
        print(f"{'':13s} all {len(d)} bins: median error {np.median(errs):.3f} deg, max {np.max(errs):.3f} deg, "
              f"median sigma_el {np.median(d['sigma_el_deg']):.3f} deg, "
              f"min TLE separation {d['tle_separation_deg'].min():.2f} deg, usable {int(d['usable'].sum())}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output-dir", type=Path, default=Path("docs/examples/stage3"))
    args = p.parse_args()
    with tempfile.TemporaryDirectory(prefix="stage3_demo_") as tmp:
        work = Path(tmp)
        manifest = build(work)
        code = cli.main(["--manifest", str(manifest), "--all-stations", "--image-validation",
                         "--bootstrap", "20", "--time-bin", "1", "--output-dir", str(args.output_dir),
                         "--data-root", str(work)])
        out = args.output_dir / "20260625" / "NORAD_40534"
        report(out, work)
    print(f"\ntart-localise exit status {code}; outputs in {out}")
    return code


if __name__ == "__main__":
    sys.exit(main())
