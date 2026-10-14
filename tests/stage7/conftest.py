"""Module-level synthetic diagnosis chains (built once, reused by the Stage-7 tests)."""

import json

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")

from tart_gnss.stage4.synthetic import EXAMPLE_STATIONS  # noqa: E402
from tart_gnss.stage7.synthetic import Scenario, build_chain  # noqa: E402

NAMIBIA = (-22.612053, 17.856784, 1655.0)
_BLOCKS = (np.arange(60) // 5) % 2 == 0


def scenarios():
    """The scenarios of specification section 35 (and controls)."""
    rng = np.random.default_rng(3)
    return {
        "clean": Scenario("clean"),
        "bias": Scenario("bias", bias_deg={"madagascar": (0.03, 0.02)}),
        "calibration": Scenario("calibration", calibration_bias=("madagascar", (0.03, 0.02))),
        "timing": Scenario("timing", time_offsets_s={"mauritius": rng.uniform(1.0, 5.0, 36)},
                           time_match_tolerance_s=6.0),
        "geometry": Scenario("geometry", n_epochs=60, seed=2, subsets={"south_africa": _BLOCKS, "namibia": _BLOCKS},
                             min_stations=2, stations={**EXAMPLE_STATIONS, "namibia": NAMIBIA}),
        "dynamics": Scenario("dynamics", accel_m_s2=(0.35, 0.35, 0.0)),
        "stale": Scenario("stale", tau_s=10.0, tle_epoch_shift_d=20.0),
        "network": Scenario("network", tau_s=10.0),
        "nostage2": Scenario("nostage2", tau_s=10.0, with_stage2=False),
        "random": Scenario("random", random_metadata=True, seed=5),
    }


class Run:
    """A diagnosed scenario: chain paths, exit code, summary and tables."""

    def __init__(self, name, sc, base, *extra):
        from tart_gnss.stage7 import cli
        from tart_gnss.stage7.output import read_table

        self.chain = build_chain(sc, base / name)
        out = base / name / "stage7"
        self.code = cli.main(["--stage6-manifest", str(self.chain["stage6"]), "--output-dir", str(out), "--no-plots",
                              *extra])
        d = next(out.glob("*/NORAD_*"))
        self.dir = d
        self.summary = json.loads((d / "summary.json").read_text())
        self.manifest = json.loads((d / "stage7_manifest.json").read_text())
        self.diagnosis = read_table(d / "diagnosis.csv")
        self.influence = read_table(d / "station_influence.csv")
        self.pairs = read_table(d / "pairwise.csv")
        import pandas as pd

        self.correlations = pd.read_csv(d / "correlations.csv")

    @property
    def tags(self):
        return set(self.summary["overall_diagnostic_tags"])

    @property
    def level(self):
        return self.summary["overall_evidence_level"]


@pytest.fixture(scope="session")
def runs(tmp_path_factory):
    """Lazily built scenario runs, keyed by name."""
    base = tmp_path_factory.mktemp("stage7")
    cache = {}
    sc = scenarios()

    def get(name):
        if name not in cache:
            cache[name] = Run(name, sc[name], base)
        return cache[name]

    return get
