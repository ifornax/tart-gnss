"""Transparent rule-based diagnostic confidence (no probabilities).

Six boolean components are evaluated per epoch:

``sample_support``          at least ``min_samples`` diagnosed epochs in the run
``effect_size``             anomalous: Stage-6 D2 >= 2 x threshold; normal: D2 <= threshold / 2
``persistence``             anomalous: part of a persistent run (see attribution); normal: true
``cross_method_agreement``  the station test's pairwise and reference methods agree, or (network
                            consistency) both the leave-one-out and the pairwise tests are clean;
                            normal: true
``reference_quality``       Stage-6 ``precision_validation`` is true and the reference uncertainty is
                            known (empirical TLE ensemble)
``data_completeness``       Stages 2-5 linked, calibration metadata present, all links verified

``HIGH`` needs all six; ``MODERATE`` needs at least four with ``reference_quality``
and ``data_completeness`` not both false; otherwise ``LOW``. A missing
reference or incomplete data therefore never yields ``HIGH``.
"""

from __future__ import annotations

from typing import Dict

LEVELS = ("LOW", "MODERATE", "HIGH")
COMPONENTS = ("sample_support", "effect_size", "persistence", "cross_method_agreement", "reference_quality",
              "data_completeness")


def confidence(components: Dict[str, bool]) -> str:
    """Map the six components to LOW / MODERATE / HIGH."""
    n = sum(bool(components.get(c)) for c in COMPONENTS)
    if n == len(COMPONENTS):
        return "HIGH"
    if n >= 4 and (components.get("reference_quality") or components.get("data_completeness")):
        return "MODERATE"
    return "LOW"
