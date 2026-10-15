"""Synthetic campaigns: many complete Stage 3-7 runs with known population properties.

Each run is a :class:`tart_gnss.stage7.synthetic.Scenario` built through the
real pipeline (``tart-triangulate``, ``tart-track``, ``tart-validate`` with a
TLE ensemble, ``tart-diagnose``). The default pool contains:

* three synthetic satellites: the PRN-26 example (GPS, NORAD 40534), and two
  relabelled neighbours 2 deg ahead / behind in mean anomaly: ``SYN-GAL``
  (Galileo, NORAD 90001, noisier: 0.008 deg, SNR 12) and ``SYN-BDS``
  (BeiDou, NORAD 90002, 0.006 deg, SNR 16); GPS runs use 0.005 deg, SNR 20;
* repeated passes on later sidereal days (the ground track repeats);
* an unstable station: Madagascar biased by 0.036 deg in three runs (without
  the TLE-age drift: a network offset together with a station bias can make
  Stage 7's TLE-based and pairwise station tests disagree, which Stage 7
  reports as ambiguous);
* varying geometry: two-station runs on each station pair (crossing angles
  about 2.8, 6.6 and 8.6 deg) and a four-station run with Namibia;
* varying TLE age: every GPS run carries an injected along-track drift of
  0.4 s (about 1.6 km) per day of TLE age (ages 0.35-20 d; ages > 14 d are
  stale, i.e. descriptive only);
* repeatability: three ``SYN-GAL`` passes with known along-track offsets of
  -1, 0, +1 s (between-pass standard deviation = |v| x 1 s).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from tart_gnss.stage4.synthetic import EXAMPLE_STATIONS
from tart_gnss.stage6 import synthetic as syn6
from tart_gnss.stage7.synthetic import L1, L2, SIDEREAL_DAY_S, Scenario, build_chain

NAMIBIA = (-22.612053, 17.856784, 1655.0)
GAL_TLE = syn6.perturbed_tle(L1, L2, d_mean_anomaly_deg=2.0, satnum=90001)
BDS_TLE = syn6.perturbed_tle(L1, L2, d_mean_anomaly_deg=-2.0, satnum=90002)
AGE_DRIFT_S_PER_DAY = 0.4
REPEAT_OFFSETS_S = (-1.0, 0.0, 1.0)
N_EPOCHS = 24


def _seed(name: str) -> int:
    return sum(ord(c) * (i + 1) for i, c in enumerate(name)) % 10007


def _gps(name, day, **kw):
    kw = {**_aged(day), **kw}
    return Scenario(name, n_epochs=N_EPOCHS, t0_offset_s=day * SIDEREAL_DAY_S, satellite_name="GPS BIIF-9 (PRN 26)",
                    constellation="GPS", snr=20.0, sigma_deg=0.005, seed=_seed(name), **kw)


def _gal(name, day, **kw):
    return Scenario(name, n_epochs=N_EPOCHS, t0_offset_s=day * SIDEREAL_DAY_S, tle=GAL_TLE, satellite_name="SYN-GAL",
                    constellation="Galileo", snr=12.0, sigma_deg=0.008, seed=_seed(name), **kw)


def _bds(name, day, **kw):
    return Scenario(name, n_epochs=N_EPOCHS, t0_offset_s=day * SIDEREAL_DAY_S, tle=BDS_TLE, satellite_name="SYN-BDS",
                    constellation="BeiDou", snr=16.0, sigma_deg=0.006, seed=_seed(name), **kw)


def _aged(day):
    """Observation ``day`` sidereal days after the TLE epoch + 0.35 d: age ~ day + 0.35, drift injected."""
    return {"tau_s": AGE_DRIFT_S_PER_DAY * (day + 0.35), "ensemble": day + 0.35 <= 14}


def default_pool() -> List[Tuple[str, Scenario]]:
    """The demonstration campaign (16 runs)."""
    def drop(sid):
        return {sid: np.zeros(N_EPOCHS, bool)}

    return [
        ("gps_clean_d0", _gps("gps_clean_d0", 0)),
        ("gps_clean_d1", _gps("gps_clean_d1", 1)),
        ("gps_clean_d2", _gps("gps_clean_d2", 2)),
        ("gps_unstable_d3", _gps("gps_unstable_d3", 3, bias_deg={"madagascar": (0.03, 0.02)}, tau_s=0.0)),
        ("gps_unstable_d4", _gps("gps_unstable_d4", 4, bias_deg={"madagascar": (0.03, 0.02)}, tau_s=0.0)),
        ("gps_unstable_d5", _gps("gps_unstable_d5", 5, bias_deg={"madagascar": (0.03, 0.02)}, tau_s=0.0)),
        ("gps_pair_mu_mg_d6", _gps("gps_pair_mu_mg_d6", 6, subsets=drop("south_africa"), min_stations=2)),
        ("gps_pair_mg_za_d8", _gps("gps_pair_mg_za_d8", 8, subsets=drop("mauritius"), min_stations=2)),
        ("gps_pair_mu_za_d10", _gps("gps_pair_mu_za_d10", 10, subsets=drop("madagascar"), min_stations=2)),
        ("gps_wide_d9", _gps("gps_wide_d9", 9, stations={**EXAMPLE_STATIONS, "namibia": NAMIBIA})),
        ("gps_aged_d7", _gps("gps_aged_d7", 7)),
        ("gps_aged_d12", _gps("gps_aged_d12", 12)),
        ("gps_aged_d20", _gps("gps_aged_d20", 20)),
        ("gal_repeat_m1", _gal("gal_repeat_m1", 0, tau_s=REPEAT_OFFSETS_S[0])),
        ("gal_repeat_0", _gal("gal_repeat_0", 1, tau_s=REPEAT_OFFSETS_S[1])),
        ("gal_repeat_p1", _gal("gal_repeat_p1", 2, tau_s=REPEAT_OFFSETS_S[2])),
        ("bds_clean_d0", _bds("bds_clean_d0", 0)),
    ]


def build_runs(work: Path, specs: List[Tuple[str, Scenario]], plots: bool = False) -> Dict[str, Path]:
    """Build each run through Stages 3-7; return ``name -> stage7_manifest.json``."""
    from tart_gnss.stage7 import cli as cli7

    out = {}
    for name, sc in specs:
        chain = build_chain(sc, Path(work) / name)
        cli7.main(["--stage6-manifest", str(chain["stage6"]), "--output-dir", str(Path(work) / name / "stage7"),
                   *([] if plots else ["--no-plots"])])
        out[name] = next((Path(work) / name / "stage7").glob("*/NORAD_*/stage7_manifest.json"))
    return out
