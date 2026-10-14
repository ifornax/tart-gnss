"""Shared helpers for the Stage-6 tests."""

from __future__ import annotations

import json

import numpy as np

from tart_gnss.stage6 import synthetic as syn
from tart_gnss.stage6.rtn import rtn_bases
from tart_gnss.stage6.tle_reference import TLEReference
from tart_gnss.tle import tle_epoch_utc

L1, L2 = syn.EXAMPLE_TLE
TLE_EPOCH_S = tle_epoch_utc(L1, L2).timestamp()
#: Shortly after the TLE epoch (fresh reference).
T_FRESH = TLE_EPOCH_S + 600.0


def reference(line1=L1, line2=L2) -> TLEReference:
    """A TLEReference without file provenance."""
    return TLEReference(norad_id=int(line1[2:7]), name="TEST", line1=line1, line2=line2,
                        epoch_utc=tle_epoch_utc(line1, line2), source_kind="tle-file")


def times(n=30, step=20.0, t0=T_FRESH):
    """Epochs (Unix s)."""
    return t0 + step * np.arange(n)


def truth(t):
    """TLE GCRS states by the independent (non-Skyfield) route."""
    return syn.independent_gcrs_states(L1, L2, t)


def inject_rtn(states, dr_rtn=(0.0, 0.0, 0.0), dv_rtn=(0.0, 0.0, 0.0)):
    """Add RTN offsets (defined on the reference state's axes) to GCRS states."""
    q = rtn_bases(states[:, :3], states[:, 3:])
    out = states.copy()
    out[:, :3] += np.einsum("nji,j->ni", q, np.asarray(dr_rtn, float))
    out[:, 3:] += np.einsum("nji,j->ni", q, np.asarray(dv_rtn, float))
    return out


def stage5_run(tmp_path, t, states, sigma_pos=1e3, sigma_vel=1.0, **kw):
    """Write a synthetic Stage-5 run; return its manifest path."""
    cov = np.repeat(syn.isotropic_state_covariance(sigma_pos, sigma_vel)[None], len(t), axis=0)
    return syn.write_stage5_run(tmp_path / "stage5", t, states, cov, **kw)


def tle_file(tmp_path, sets=None, name="ref.tle"):
    """Write a TLE file (default: the example set)."""
    return syn.write_tle_file(tmp_path / name, sets or [(L1, L2)])


def run_cli(tmp_path, manifest, *extra):
    """Run tart-validate; return (code, residuals, summary, manifest) (Nones if nothing was written)."""
    from tart_gnss.stage6 import cli
    from tart_gnss.stage6.output import read_table

    out = tmp_path / "stage6"
    code = cli.main(["--stage5-manifest", str(manifest), "--output-dir", str(out), *extra])
    found = list(out.glob("*/NORAD_*/stage6_manifest.json"))
    if not found:
        return code, None, None, None
    d = found[0].parent
    return (code, read_table(d / "residuals.csv"), json.loads((d / "summary.json").read_text()),
            json.loads(found[0].read_text()))
