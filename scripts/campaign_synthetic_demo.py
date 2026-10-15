#!/usr/bin/env python
"""Synthetic campaign demonstration for ``tart-campaign``.

Builds the default 17-run synthetic pool (three satellites in three
constellations, repeated passes, an unstable station, varying geometry and
TLE age; see :mod:`tart_gnss.campaign.synthetic`) through the real Stage 3-7
chain, then runs ``tart-campaign`` on it and writes a short report.

Usage::

    python scripts/campaign_synthetic_demo.py [--output-dir docs/examples/campaign] [--bootstrap 1000]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from tart_gnss.campaign import cli
from tart_gnss.campaign.synthetic import AGE_DRIFT_S_PER_DAY, REPEAT_OFFSETS_S, build_runs, default_pool


def fmt_ci(d, scale=1.0, unit="", p=2):
    if not d:
        return "n/a"
    lo, hi = d["ci95"]
    return f"{d['median'] / scale:.{p}f} {unit} (95% CI {lo / scale:.{p}f}-{hi / scale:.{p}f}, {d['n_runs']} runs)"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", type=Path, default=Path("docs/examples/campaign"))
    ap.add_argument("--bootstrap", type=int, default=1000)
    args = ap.parse_args()
    runs_dir = args.output_dir / "runs"
    build_runs(runs_dir, default_pool())
    code = cli.main(["--input-root", str(runs_dir), "--output-dir", str(args.output_dir), "--campaign-name",
                     "synthetic", "--bootstrap", str(args.bootstrap), "--group-by", "constellation", "station_set"])
    d = args.output_dir / "synthetic"
    s = json.loads((d / "summary.json").read_text())
    st = pd.read_csv(d / "tables" / "station_summary.csv")
    co = pd.read_csv(d / "tables" / "constellation_summary.csv")
    g = s["geometry_dependence"]
    a = s["tle_age_dependence"]["residual_vs_abs_tle_age"]
    lines = [
        "Synthetic campaign demonstration (scripts/campaign_synthetic_demo.py)", f"tart-campaign exit code {code}", "",
        f"runs: {s['n_included_runs']} included of {s['n_discovered_runs']} discovered ({s['n_excluded_runs']} excluded); "
        f"satellites {s['n_satellites']}; constellations {', '.join(s['constellations'])}; stations "
        f"{s['n_unique_stations']}; sample-size tier {s['sample_size_tier']}",
        f"precision-validation fraction {s['precision_validation_fraction']:.3f} "
        f"({s['n_precision_validation_runs']} precision, {s['n_descriptive_only_runs']} descriptive only)",
        f"observing span {s['observation_start']} .. {s['observation_end']}; total {s['total_observing_duration_s']:.0f} s",
        f"campaign quality flags: {', '.join(s['campaign_quality_flags'])}", "",
        "Measurement uncertainty (all included runs; per-run medians, run-level bootstrap):",
        f"  Stage-3 direction sigma     {fmt_ci(s['median_stage3_direction_uncertainty_deg'], 1, 'deg', 4)}",
        f"  Stage-4 3-D sigma           {fmt_ci(s['median_stage4_3d_uncertainty_m'], 1e3, 'km')}",
        f"  Stage-4 major / minor axis  {fmt_ci(s['median_stage4_sigma_major_m'], 1e3, 'km')} / "
        f"{fmt_ci(s['median_stage4_sigma_minor_m'], 1e3, 'km')}",
        f"  Stage-5 smoothed 3-D sigma  {fmt_ci(s['median_stage5_smoothed_uncertainty_m'], 1e3, 'km')}", "",
        "TLE-relative residual (precision cohort; per-run medians):",
        f"  3-D                         {fmt_ci(s['median_tle_relative_residual_3d_m'], 1e3, 'km')}",
        f"  radial                      {fmt_ci(s['median_tle_relative_radial_m'], 1e3, 'km')}",
        f"  along-track                 {fmt_ci(s['median_tle_relative_along_track_m'], 1e3, 'km')}",
        f"  cross-track                 {fmt_ci(s['median_tle_relative_cross_track_m'], 1e3, 'km')}", "",
        "Station stability (per-run median direction residual across runs):",
        st[["station_id", "n_runs", "median_direction_residual_deg", "std_direction_residual_deg",
            "p84_direction_residual_deg", "fraction_runs_station_specific"]].round(4).to_string(index=False),
        f"  most stable: {s['most_stable_stations']}; most frequently flagged: {s['most_frequently_flagged_stations']}",
        "", "Constellation summary:",
        co[["constellation", "n_runs", "median_snr", "median_sigma_theta_deg", "median_stage4_sigma_3d_m",
            "tle_relative_median_position_residual_3d_m", "under_sampled"]].round(4).to_string(index=False),
        f"  SNR comparison: {s['constellation_comparisons']['median_snr']['test']}, p = "
        f"{s['constellation_comparisons']['median_snr']['p_value']:.2g}, "
        f"{s['constellation_comparisons']['median_snr']['effect_size_name']} = "
        f"{s['constellation_comparisons']['median_snr']['effect_size']:.2f}; under-sampled: "
        f"{s['constellation_comparisons']['median_snr']['under_sampled']}", "",
        "Geometry dependence:",
        f"  median Stage-4 3-D sigma by number of stations: "
        + ", ".join(f"{k}: {v / 1e3:.1f} km" for k, v in g["median_stage4_sigma_3d_m_by_n_stations"].items()),
        f"  two-station runs, log10 sigma vs crossing angle: Spearman "
        f"{g['log10_stage4_sigma_vs_crossing_angle_by_n_stations']['2']['spearman_r']:.2f}, Theil-Sen slope "
        f"{g['log10_stage4_sigma_vs_crossing_angle_by_n_stations']['2']['theil_sen_slope']:.3f} dex/deg", "",
        f"TLE-age dependence (injected {AGE_DRIFT_S_PER_DAY} s/day along-track drift on GPS runs = ~1.55 km/day):",
        f"  Spearman {a['spearman_r']:.2f} (p {a['spearman_p']:.2g}); Theil-Sen slope {a['theil_sen_slope'] / 1e3:.2f} "
        f"km/day (95% {a['theil_sen_slope_low'] / 1e3:.2f}-{a['theil_sen_slope_high'] / 1e3:.2f})", "",
        "Evidence-level fractions (per run): " + ", ".join(
            f"L{k[-1]} {v:.2f}" for k, v in s["evidence_level_fractions_by_run"].items()),
        "Diagnostic-tag fractions (per run, Wilson 95%): " + "; ".join(
            f"{k} {v['fraction']:.2f} [{v['ci95'][0]:.2f}, {v['ci95'][1]:.2f}]"
            for k, v in s["diagnostic_tag_fractions_by_run"].items() if v["count"]), "",
        f"Repeatability (SYN-GAL, passes offset by {REPEAT_OFFSETS_S} s; expected along-track std = |v| x 1 s ~ 3.87 km):",
    ]
    for nid, r in s["repeatability"].items():
        lines.append(f"  NORAD {nid}: {r['n_passes']} passes; between-pass std R/T/N = "
                     + ", ".join(f"{r[c]['between_pass_std_m'] / 1e3:.2f}" for c in ("radial", "along_track",
                                                                                   "cross_track"))
                     + " km; within-pass sigma R/T/N = "
                     + ", ".join(f"{r[c]['within_pass_median_sigma_m'] / 1e3:.2f}" for c in ("radial", "along_track",
                                                                                           "cross_track")) + " km")
    lines.append("  (NORAD 40534 passes differ by design: TLE-age drift, station bias, geometry; not a repeat test)")
    (args.output_dir / "demo_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
