#!/usr/bin/env python3
"""Compare UTC and tab-sim UT1-labelled Julian dates in Ghana GNSS visibilities.

Inputs:
    Ghana 1SAT Measurement Set and matching used_orbits.json, snapshot index
    and number of simulation subsamples.

Outputs:
    Normalised complex-amplitude-profiled visibility residuals for two
    orbit/antenna-time conventions, plus the corresponding topocentric
    directions at the MS timestamp.

Logic:
    Generate the same numeric Julian dates used in the simulator. Compare
    (1) Skyfield UTC input via ts.from_datetimes and (2) the simulator's use of
    ts.ut1_jd on those numeric dates. Use Skyfield ITRSPosition for antennas
    in BOTH cases, matching the simulator's frame path, so the only
    experimental variable is the time-scale interpretation.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from datetime import timezone

import numpy as np
from astropy.time import Time
from casacore.tables import table
from scipy.constants import c
from skyfield.api import Distance, load, wgs84
from skyfield.toposlib import ITRSPosition

from tabsim_gnss_fit import read_snapshot
from tabsim_gnss_diagnostics import read_saved_omm
from tabsim.tle import earth_satellite

LOGGER = logging.getLogger(__name__)


def score(observed: np.ndarray, model: np.ndarray) -> dict:
    """Compute profiled complex-amplitude residual and phase coherence.

    Args:
        observed: Complex visibilities for unflagged first-integration baselines.
        model: Model complex phasors for corresponding baselines.

    Returns:
        Dictionary with normalised squared residual and coherence.
    """
    vv = float(np.vdot(observed, observed).real)
    mm = float(np.vdot(model, model).real)
    amp = np.vdot(model, observed) / mm
    residual = observed - amp * model
    return {
        "J": float(np.vdot(residual, residual).real / vv),
        "coherence": float(abs(np.vdot(model, observed)) / np.sqrt(mm * vv)),
    }


def model_for_scale(d: dict, omm: dict, jd_utc: np.ndarray,
                    convention: str) -> tuple[np.ndarray, tuple[float, float]]:
    """Predict satellite phasors using a selected numerical JD interpretation.

    Args:
        d: MS snapshot with ITRF antenna positions and baseline data.
        omm: Matching saved orbit record.
        jd_utc: Numeric UTC Julian dates for integration fine samples.
        convention: 'utc' (correct interpretation) or 'simulator_ut1'.

    Returns:
        Tuple of model averaged baseline phasors, and Az/El at t0.
    """
    ts = load.timescale()
    times = (ts.from_datetimes(Time(jd_utc, format="jd", scale="utc").to_datetime(timezone=timezone.utc))
             if convention == "utc" else ts.ut1_jd(jd_utc))
    sat = earth_satellite(omm, ts)
    sat_xyz = sat.at(times).position.m.T
    antennas = np.stack(
        [
            ITRSPosition(Distance(m=position)).at(times).position.m.T
            for position in d["antenna_ecef"]
        ],
        axis=1
    )
    path_lengths = np.linalg.norm(
        sat_xyz[:, None, :] - antennas, axis=-1
    )
    delta = path_lengths[:, d["a1"]] - path_lengths[:, d["a2"]]
    model = np.exp(
        (-2j * np.pi * d["freq"] / c) * (
            delta + d["uvw"][None, :, 2]
        )
    ).mean(axis=0)
    # Evaluate topocentric direction at exact MS timestamp, not first fine sample.
    jd0 = float(Time(d["ms_time"] / 86400.0, format="mjd", scale="utc").jd)
    t0 = (ts.from_datetime(Time(jd0, format="jd", scale="utc").to_datetime(timezone=timezone.utc))
          if convention == "utc" else ts.ut1_jd(jd0))
    centre = np.mean(d["antenna_ecef"], axis=0)
    from astropy.coordinates import EarthLocation
    from astropy import units as u
    loc = EarthLocation.from_geocentric(*(centre * u.m))
    observer = wgs84.latlon(
        float(loc.lat.deg),
        float(loc.lon.deg),
        elevation_m=float(loc.height.to_value(u.m)),
    )
    alt, az, _ = (sat - observer).at(t0).altaz()
    return model, (float(az.degrees), float(alt.degrees))


def run(ms_path: str, orbit_path: str, samples: int, snapshot: int,
        output_path: str) -> None:
    """Compare simulation-compatible and physical-UTC time interpretations.

    Args:
        ms_path: Simulated Ghana Measurement Set.
        orbit_path: Orbit record saved alongside that exact MS.
        samples: Integration samples, normally 33.
        snapshot: Zero-indexed integration.
        output_path: JSON result destination.

    Returns:
        None. Prints and stores comparison results.
    """
    if samples < 1:
        raise ValueError("samples must be positive")
    d = read_snapshot(ms_path, snapshot, "RFI_MODEL_DATA")
    numeric_jd0 = float(Time(
        d["ms_time"] / 86400.0, format="mjd", scale="utc"
    ).jd)
    offsets = ((np.arange(samples) + 0.5) / samples - 0.5) * d["dt"]
    jd = numeric_jd0 + offsets / 86400.0
    omm = read_saved_omm(orbit_path)
    records = {}
    for name in ("utc", "simulator_ut1"):
        model, direction = model_for_scale(d, omm, jd, name)
        results = score(np.asarray(d["vis"]), model)
        results["az_deg"], results["el_deg"] = direction
        records[name] = results

    az_utc, el_utc = records["utc"]["az_deg"], records["utc"]["el_deg"]
    az_sim, el_sim = records["simulator_ut1"]["az_deg"], records["simulator_ut1"]["el_deg"]
    daz = ((az_sim - az_utc + 180) % 360) - 180
    east = daz * np.cos(np.deg2rad(el_utc)) * 3600
    north = (el_sim - el_utc) * 3600
    records["simulator_minus_utc_arcsec"] = {
        "east": float(east), "north": float(north),
        "radial": float(np.hypot(east, north))
    }
    for name in ("utc", "simulator_ut1"):
        item = records[name]
        print(
            f"{name:17s} J={item['J']:.12g} "
            f"C={item['coherence']:.12f} "
            f"Az={item['az_deg']:.9f} El={item['el_deg']:.9f}"
        )
    print("simulator - UTC arcsec:", records["simulator_minus_utc_arcsec"])
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    LOGGER.info("Saved %s", out)


def main() -> None:
    """Parse command line and run the frame-consistent time-scale test."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ms", required=True)
    parser.add_argument("--orbits", required=True)
    parser.add_argument("--samples", type=int, default=33)
    parser.add_argument("--snapshot", type=int, default=0)
    parser.add_argument("--output", default="results/tabsim_validation/ghana_audit/time_scale_check.json")
    args = parser.parse_args()
    run(args.ms, args.orbits, args.samples, args.snapshot, args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
