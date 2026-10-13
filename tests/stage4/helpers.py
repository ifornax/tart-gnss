"""Shared helpers for the Stage-4 tests."""

from __future__ import annotations

import numpy as np

from tart_gnss.stage4.geometry import geodetic_to_ecef, lm_covariance_to_ecef, tangent_covariance
from tart_gnss.stage4.synthetic import EXAMPLE_POSITION_M, EXAMPLE_STATIONS, exact_lm, isotropic_lm_covariance

#: A fourth (real) station used where four lines are needed.
NAMIBIA = (-22.612053, 17.856784, 1655.0)


def lines(truth=EXAMPLE_POSITION_M, stations=EXAMPLE_STATIONS, sigma_deg=0.1, sigmas=None):
    """Station positions, exact sight vectors and isotropic angular covariances (ECEF)."""
    r, s, c = [], [], []
    for i, (sid, (lat, lon, alt)) in enumerate(stations.items()):
        ri = geodetic_to_ecef(lat, lon, alt)
        d = np.asarray(truth) - ri
        sig = np.radians(sigmas[i] if sigmas is not None else sigma_deg)
        l, m = exact_lm(truth, lat, lon, alt)
        r.append(ri)
        s.append(d / np.linalg.norm(d))
        c.append(lm_covariance_to_ecef(l, m, isotropic_lm_covariance(l, m, sig), lat, lon))
    return np.array(r), np.array(s), np.array(c)


def perturb(s, c, rng):
    """Redraw sight vectors from their tangent covariances."""
    out = []
    for si, ci in zip(s, c):
        e, c2 = tangent_covariance(si, ci)
        v = si + e @ rng.multivariate_normal(np.zeros(2), c2)
        out.append(v / np.linalg.norm(v))
    return np.array(out)


def run_cli(tmp_path, rows, *extra, stations=EXAMPLE_STATIONS, **kw):
    """Write a synthetic Stage-3 run, run ``tart-triangulate``; return (code, positions, contributions, manifest)."""
    import json

    from tart_gnss.stage4 import cli
    from tart_gnss.stage4.output import read_table
    from tart_gnss.stage4.synthetic import write_stage3_run

    m = write_stage3_run(tmp_path / "stage3", rows, stations=stations, **kw)
    out = tmp_path / "stage4"
    code = cli.main(["--stage3-manifest", str(m), "--output-dir", str(out), "--no-plots", *extra])
    found = list(out.glob("*/NORAD_*/stage4_manifest.json"))
    if not found:
        return code, None, None, None
    d = found[0].parent
    pos = read_table(d / "positions.csv")
    con = read_table(d / "station_contributions.csv")
    return code, pos, con, json.loads(found[0].read_text())


def truth_at(pos, truth):
    """Truth positions at the reference epochs of a positions table."""
    import pandas as pd

    t = (pos["time_utc"] - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds()
    return np.array([truth(x) for x in t])
