#!/usr/bin/env python3
"""Re-score existing Ghana GNSS Monte Carlo fits against two time conventions.

Inputs:
    Trial CSV from tabsim_ghana_uncertainty.py; optional time-scale diagnostic
    JSON, otherwise explicitly supplied simulator and physical UTC Az/El.

Outputs:
    Rescored per-trial CSV, per-SNR summary CSV, covariance JSON, error-vs-SNR
    PNG, and a JSON comparison of reference directions.

Logic:
    Do NOT rerun visibility fits. Recalculate East/North angular errors for
    each saved Az/El against the chosen reference, preserving fit metadata.
    Estimate covariance and bias using the first half of each SNR group and
    evaluate raw/bias-corrected joint 68% coverage on held-out trials.
    The empirical covariance is conditional on the simulation and noise model.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

LOGGER = logging.getLogger(__name__)


def read_trials(filename: Path) -> list[dict]:
    """Read per-realisation Az/El fits with numeric and boolean conversion.

    Args:
        filename: Original Monte Carlo trial CSV.

    Returns:
        List of trial dictionaries preserving input columns.
    """
    with filename.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Trial CSV is empty")
    for row in rows:
        for key in ("az_deg", "el_deg", "snr", "residual"):
            row[key] = float(row[key])
        row["trial"] = int(row["trial"])
        for key in ("local_success", "on_edge"):
            row[key] = str(row[key]).strip().lower() == "true"
    return rows


def save_csv(path: Path, rows: list[dict]) -> None:
    """Write dictionaries as a CSV table.

    Args:
        path: Destination path.
        rows: Non-empty list of homogeneous dictionaries.

    Returns:
        None.
    """
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def errors(rows: list[dict], truth: tuple[float, float]) -> list[dict]:
    """Recalculate local tangent-plane errors relative to one direction.

    Args:
        rows: Saved direction-estimation trials.
        truth: (azimuth, elevation) degrees.

    Returns:
        New trial records with rescored East/North/radial errors in arcsec.
    """
    az_truth, el_truth = truth
    result = []
    for src in rows:
        row = src.copy()
        d_az = (row["az_deg"] - az_truth + 180) % 360 - 180
        row["east_arcsec"] = float(
            d_az * np.cos(np.deg2rad(el_truth)) * 3600
        )
        row["north_arcsec"] = float((row["el_deg"] - el_truth) * 3600)
        row["error_arcsec"] = float(np.hypot(
            row["east_arcsec"], row["north_arcsec"]
        ))
        result.append(row)
    return result


def summarise(group: list[dict], truth: tuple[float, float]) -> dict:
    """Calculate error quantiles, training covariance and held-out coverage.

    Args:
        group: Trials at a single SNR sorted by trial index.
        truth: Reference direction, Az/El degrees.

    Returns:
        Summary including covariance (East/North arcsec^2) and coverage.
    """
    group = sorted(group, key=lambda r: r["trial"])
    arr = np.array([[r["east_arcsec"], r["north_arcsec"]] for r in group])
    if len(arr) < 6:
        raise ValueError("At least six trials per SNR are necessary")
    n_train = len(arr)//2
    train, test = arr[:n_train], arr[n_train:]
    bias = train.mean(axis=0)
    covariance = np.cov(train, rowvar=False, ddof=1)
    norms = np.linalg.norm(arr, axis=1)
    threshold = float(chi2.ppf(0.68, 2))
    if np.linalg.eigvalsh(covariance)[0] > 0:
        inverse = np.linalg.inv(covariance)
        raw_distance = np.einsum("ni,ij,nj->n", test, inverse, test)
        corrected = test-bias
        centred_distance = np.einsum(
            "ni,ij,nj->n", corrected, inverse, corrected
        )
        coverage_raw = float(np.mean(raw_distance <= threshold))
        coverage_corrected = float(np.mean(centred_distance <= threshold))
    else:
        coverage_raw = float("nan")
        coverage_corrected = float("nan")

    # Convert tangent plane arcsec^2 to Az/El degrees^2.
    D_inverse = np.diag([
        1.0/(3600*np.cos(np.deg2rad(truth[1]))), 1.0/3600
    ])
    covariance_degrees = D_inverse @ covariance @ D_inverse
    return {
        "snr": float(group[0]["snr"]),
        "n_trials": len(group),
        "fit_success_fraction": float(np.mean(
            [r["local_success"] for r in group]
        )),
        "edge_fraction": float(np.mean([r["on_edge"] for r in group])),
        "median_error_arcsec": float(np.median(norms)),
        "p68_error_arcsec": float(np.quantile(norms, 0.68)),
        "p95_error_arcsec": float(np.quantile(norms, 0.95)),
        "rmse_error_arcsec": float(np.sqrt(np.mean(norms**2))),
        "bias_east_arcsec": float(bias[0]),
        "bias_north_arcsec": float(bias[1]),
        "covariance_east_north_arcsec2": covariance.tolist(),
        "covariance_azel_deg2": covariance_degrees.tolist(),
        "heldout_coverage68_raw": coverage_raw,
        "heldout_coverage68_bias_corrected": coverage_corrected,
        "truth_az_deg": truth[0],
        "truth_el_deg": truth[1],
    }


def plot_error(path: Path, summaries: list[dict], label: str) -> None:
    """Save logarithmic SNR-versus-error-quantile plot.

    Args:
        path: Output PNG.
        summaries: Per-SNR error statistics.
        label: Time-convention label for chart title.

    Returns:
        None.
    """
    ordered = sorted(summaries, key=lambda s: s["snr"])
    x = [r["snr"] for r in ordered]
    plt.figure(figsize=(8, 5))
    for key, label_text in (
        ("median_error_arcsec", "Median"),
        ("p68_error_arcsec", "68th percentile"),
        ("p95_error_arcsec", "95th percentile"),
    ):
        plt.plot(x, [r[key] for r in ordered], marker="o", label=label_text)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Per-visibility complex-amplitude SNR")
    plt.ylabel("Angular error (arcsec)")
    plt.title(f"Ghana GNSS direction error versus SNR — {label}")
    plt.grid(True, alpha=0.25, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(args: argparse.Namespace) -> None:
    """Reanalyse saved trials using physical-UTC and simulator truth separately.

    Args:
        args: CLI arguments with trials, reference angles and destination.

    Returns:
        None; writes both references' results without modifying input CSV.
    """
    if args.time_check:
        info = json.loads(Path(args.time_check).read_text(encoding="utf-8"))
        refs = {
            "simulator_ut1": (
                float(info["simulator_ut1"]["az_deg"]),
                float(info["simulator_ut1"]["el_deg"]),
            ),
            "physical_utc": (
                float(info["utc"]["az_deg"]), float(info["utc"]["el_deg"])
            ),
        }
    else:
        if any(x is None for x in (
            args.sim_az, args.sim_el, args.utc_az, args.utc_el
        )):
            raise ValueError(
                "Provide --time-check JSON, or all four Az/El values"
            )
        refs = {
            "simulator_ut1": (args.sim_az, args.sim_el),
            "physical_utc": (args.utc_az, args.utc_el),
        }

    original = read_trials(Path(args.trials))
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    comparison = {
        "reference_directions": {
            k: {"az_deg": v[0], "el_deg": v[1]}
            for k, v in refs.items()
        },
        "note": (
            "The same 1200 direction fits are rescored; no optimisation "
            "has been rerun. Covariances are training-sample scatter."
        ),
        "results": {},
    }

    for ref_name, truth in refs.items():
        output_dir = root/ref_name
        output_dir.mkdir(parents=True, exist_ok=True)
        revised = errors(original, truth)
        by_snr = {}
        for trial in revised:
            by_snr.setdefault(trial["snr"], []).append(trial)
        summaries = [
            summarise(by_snr[snr], truth) for snr in sorted(by_snr)
        ]
        save_csv(output_dir/"direction_trials_rescored.csv", revised)
        save_csv(output_dir/"direction_snr_summary.csv", [
            {k: v for k, v in s.items()
             if not k.startswith("covariance_")}
            for s in summaries
        ])
        (output_dir/"direction_covariance.json").write_text(
            json.dumps({
                "reference": ref_name,
                "frame": "local tangent plane East/North at reference direction",
                "covariance_units": "arcsec^2 and az/el deg^2",
                "snr_results": summaries,
                "caveat": (
                    "Empirical training-sample scatter, not a validated "
                    "single-observation posterior covariance"
                ),
            }, indent=2)+"\n", encoding="utf-8"
        )
        plot_error(output_dir/"localisation_error_vs_snr.png",
                   summaries, ref_name)
        comparison["results"][ref_name] = summaries
        for row in summaries:
            LOGGER.info(
                "%s SNR=%g median=%.3f arcsec coverage68=%.3f "
                "bias=(%.3f,%.3f) arcsec",
                ref_name, row["snr"], row["median_error_arcsec"],
                row["heldout_coverage68_raw"],
                row["bias_east_arcsec"], row["bias_north_arcsec"],
            )

    (root/"reference_comparison.json").write_text(
        json.dumps(comparison, indent=2)+"\n", encoding="utf-8"
    )
    LOGGER.info("Wrote rescored outputs under %s", root)


def main() -> None:
    """Parse command-line paths and directional reference values."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", required=True)
    parser.add_argument("--time-check", default=None)
    parser.add_argument("--sim-az", type=float, default=None)
    parser.add_argument("--sim-el", type=float, default=None)
    parser.add_argument("--utc-az", type=float, default=None)
    parser.add_argument("--utc-el", type=float, default=None)
    parser.add_argument("--output-dir", default="results/tabsim_validation/ghana_truth_corrected")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
