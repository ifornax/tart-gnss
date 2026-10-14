#!/usr/bin/env python
"""Synthetic end-to-end demonstration of Stage 7 (``tart-diagnose``).

Every scenario observes the PRN-26 TLE trajectory (or a known modification of
it) from the example stations, runs the real chain Stage 3 (synthetic
directions, with a matching synthetic Stage-2 manifest) -> ``tart-triangulate``
-> ``tart-track`` -> ``tart-validate`` -> ``tart-diagnose``, and compares the
diagnosis with the injected cause pattern. A false-positive control repeats a
clean network with random-walk metadata over many noise realisations.

Usage::

    python scripts/stage7_synthetic_demo.py [--output-dir docs/examples/stage7] [--fp-runs 20]
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np

from tart_gnss.stage4.synthetic import EXAMPLE_STATIONS
from tart_gnss.stage7 import cli as cli7
from tart_gnss.stage7.synthetic import Scenario, build_chain

NAMIBIA = (-22.612053, 17.856784, 1655.0)
BLOCKS = (np.arange(60) // 5) % 2 == 0


def scenarios():
    """(name, scenario, injected pattern, expectation check, expectation text)."""
    rng = np.random.default_rng(3)
    return [
        ("clean", Scenario("clean"), "none (consistent measurements, fresh TLE)",
         lambda s: s["overall_evidence_level"] == 0 and s["overall_diagnostic_tags"] == ["NORMAL"],
         "Level 0, NORMAL"),
        ("biased_station", Scenario("bias", bias_deg={"madagascar": (0.03, 0.02)}),
         "fixed 0.036 deg direction bias at madagascar",
         lambda s: "STATION_SPECIFIC" in s["overall_diagnostic_tags"] and s["overall_evidence_level"] >= 2
         and s["station_test"]["station_specific"] == ["madagascar"], "STATION_SPECIFIC (madagascar), Level >= 2"),
        ("calibration", Scenario("calibration", calibration_bias=("madagascar", (0.03, 0.02))),
         "0.036 deg bias at madagascar only on epochs with a high flagged fraction",
         lambda s: "CALIBRATION_CORRELATED" in s["overall_diagnostic_tags"] and s["overall_evidence_level"] >= 3,
         "CALIBRATION_CORRELATED, Level >= 3"),
        ("timing", Scenario("timing", time_offsets_s={"mauritius": rng.uniform(1.0, 5.0, 36)},
                            time_match_tolerance_s=6.0),
         "mauritius samples 1-5 s late (recorded), matched without interpolation",
         lambda s: "TIMING_CORRELATED" in s["overall_diagnostic_tags"], "TIMING_CORRELATED"),
        ("geometry", Scenario("geometry", n_epochs=60, seed=2, subsets={"south_africa": BLOCKS, "namibia": BLOCKS},
                              min_stations=2, stations={**EXAMPLE_STATIONS, "namibia": NAMIBIA}),
         "blocks of 4-station and narrow 2-station (2.8 deg) epochs, no bias",
         lambda s: "GEOMETRY_CORRELATED" in s["overall_diagnostic_tags"]
         and "NETWORK_CONSISTENT" not in s["overall_diagnostic_tags"], "GEOMETRY_CORRELATED, not NETWORK_CONSISTENT"),
        ("dynamics", Scenario("dynamics", accel_m_s2=(0.35, 0.35, 0.0)),
         "unmodelled 0.5 m/s^2 acceleration of the truth",
         lambda s: "TRAJECTORY_MODEL_CORRELATED" in s["overall_diagnostic_tags"] and s["overall_evidence_level"] < 4,
         "TRAJECTORY_MODEL_CORRELATED, Level < 4"),
        ("stale_reference", Scenario("stale", tau_s=10.0, tle_epoch_shift_d=20.0),
         "10 s along-track offset validated against a 20-day-old TLE",
         lambda s: "REFERENCE_LIMITED" in s["overall_diagnostic_tags"] and s["overall_evidence_level"] <= 3,
         "REFERENCE_LIMITED, Level <= 3"),
        ("network_unexplained", Scenario("network", tau_s=10.0),
         "10 s network-wide time shift of the truth (39 km along-track), clean instruments",
         lambda s: s["overall_evidence_level"] == 5 and {"NETWORK_CONSISTENT", "UNEXPLAINED"}
         <= set(s["overall_diagnostic_tags"]), "NETWORK_CONSISTENT + UNEXPLAINED, Level 5"),
        ("missing_stage2", Scenario("nostage2", tau_s=10.0, with_stage2=False),
         "as network_unexplained, but no Stage-2 calibration metadata",
         lambda s: s["overall_evidence_level"] == 4 and "CALIBRATION_METADATA_MISSING" in s["quality_flags"],
         "Level 4 (Level 5 blocked), CALIBRATION_METADATA_MISSING"),
    ]


def strongest(s):
    best = None
    for fam in ("calibration", "timing", "geometry", "trajectory_model", "localisation"):
        x = s.get(f"strongest_{fam}_correlation")
        if x and (best is None or abs(x["partial_coefficient"]) > abs(best[1]["partial_coefficient"])):
            best = (fam, x)
    if not best:
        return "none"
    fam, x = best
    return (f"{fam}: {x['x_variable']} vs {x['y_variable']} ({x['station_id'] or 'network'}) "
            f"r={x['partial_coefficient']:+.2f} -> {x['interpretation_tag']}")


PLOTTED = ("clean", "biased_station", "network_unexplained")


def run_one(sc, work: Path, plots: bool = False):
    chain = build_chain(sc, work)
    out = work / "stage7"
    code = cli7.main(["--stage6-manifest", str(chain["stage6"]), "--output-dir", str(out),
                      *([] if plots else ["--no-plots"])])
    d = next(out.glob("*/NORAD_*"))
    s = json.loads((d / "summary.json").read_text())
    import pandas as pd

    inf = pd.read_csv(d / "station_influence.csv")
    return code, s, inf


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", type=Path, default=Path("docs/examples/stage7"))
    ap.add_argument("--fp-runs", type=int, default=20, help="false-positive control realisations (default 20)")
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    lines = ["Stage-7 synthetic demonstration (scripts/stage7_synthetic_demo.py)",
             "Chain per scenario: synthetic Stage-3 directions (+ Stage-2 calibration manifest) -> tart-triangulate -> "
             "tart-track -> tart-validate (fresh TLE unless stated) -> tart-diagnose", ""]
    table = []
    for name, sc, injected, check, expected in scenarios():
        code, s, inf = run_one(sc, out / name, plots=name in PLOTTED)
        med_shift = inf.groupby("station_id")["loo_shift_3d_m"].median()
        ok = bool(check(s))
        table.append(ok)
        lines += [
            f"[{name}]",
            f"  injected cause pattern     {injected}",
            f"  expected                   {expected}",
            f"  recovered tags             {', '.join(s['overall_diagnostic_tags'])}",
            f"  evidence level             {s['overall_evidence_level']} (confidence {s['overall_confidence']}); "
            f"anomalous epochs {s['n_anomalous_epochs']}/{s['n_diagnosed']}; exit {code}",
            f"  most influential station   {s['most_influential_station']} (station test picks "
            f"{s['station_test']['picks'] or 'none'})",
            f"  leave-one-out shift        median per removed station: "
            + (", ".join(f"{k} {v / 1e3:.1f} km" for k, v in med_shift.items()) if len(med_shift) else "n/a"),
            f"  strongest correlation      {strongest(s)}",
            f"  matched expectation        {'YES' if ok else 'NO'}",
            "",
        ]
    # False-positive control: clean network, random-walk metadata.
    fp_tags, fp_level = 0, 0
    with tempfile.TemporaryDirectory() as tmp:
        for seed in range(args.fp_runs):
            _, s, _ = run_one(Scenario("fp", random_metadata=True, seed=1000 + seed), Path(tmp) / str(seed))
            if set(s["supported_correlation_tags"]) & {"CALIBRATION_CORRELATED", "LOCALISATION_CORRELATED"}:
                fp_tags += 1
            if s["overall_evidence_level"] > 0:
                fp_level += 1
    lines += [f"False-positive control: {args.fp_runs} clean networks with random-walk flagged-fraction / SNR "
              f"metadata: {fp_tags} with a supported calibration/localisation correlation, {fp_level} with an "
              f"evidence level above 0.",
              f"Scenarios matching expectation: {sum(table)} of {len(table)}.", ""]
    (out / "demo_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
