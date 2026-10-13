"""Shared helpers for the Stage-5 tests."""

from __future__ import annotations

import json

import numpy as np

from tart_gnss.stage5 import synthetic as syn


def times(n=40, step=10.0, t0=syn.T0_UNIX_S, jitter=0.0, rng=None):
    """Observation epochs (Unix s), optionally jittered."""
    t = t0 + step * np.arange(n)
    if jitter and rng is not None:
        t = t + rng.uniform(-jitter, jitter, n)
    return t


def truth(t, extra_accel=None):
    """Exact two-body GCRS truth from the example initial state at T0."""
    return syn.truth_states(syn.example_initial_state(), syn.T0_UNIX_S, t, extra_accel)


def frame(rows):
    """Stage-4 rows as a DataFrame with parsed times (the loader's view)."""
    import pandas as pd

    df = pd.DataFrame(rows)
    df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    df["quality_flags"] = df["quality_flags"].map(list)
    return df


def run_track(rows, dynamics="two-body", sigma_a=None, policy="strict", smoother=True, **kw):
    """Build measurements from rows and run filter + smoother; return the TrackingResult."""
    from tart_gnss.stage5.cli import track
    from tart_gnss.stage5.filter import FilterConfig
    from tart_gnss.stage5.loader import build_measurements

    mset = build_measurements(frame(rows), policy)
    cfg = FilterConfig(dynamics=dynamics, process_accel_sigma_m_s2=sigma_a, **kw)
    return track(mset, cfg, smoother=smoother, policy=policy)


def run_cli(tmp_path, rows, *extra, **kw):
    """Write a synthetic Stage-4 run, run ``tart-track``; return (code, trajectory, innovations, manifest)."""
    from tart_gnss.stage5 import cli
    from tart_gnss.stage5.output import read_table

    m = syn.write_stage4_run(tmp_path / "stage4", rows, **kw)
    out = tmp_path / "stage5"
    code = cli.main(["--stage4-manifest", str(m), "--output-dir", str(out), *extra])
    found = list(out.glob("*/NORAD_*/stage5_manifest.json"))
    if not found:
        return code, None, None, None
    d = found[0].parent
    return (code, read_table(d / "trajectory.csv"), read_table(d / "innovations.csv"),
            json.loads(found[0].read_text()))


def smoothed_states(traj):
    """``(N, 6)`` smoothed GCRS states from a trajectory table."""
    cols = [f"{c}_gcrs_smoothed_m" for c in "xyz"] + [f"v{c}_gcrs_smoothed_m_s" for c in "xyz"]
    return traj[cols].to_numpy(dtype=float)
