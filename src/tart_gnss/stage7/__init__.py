"""Stage 7: sensor-network diagnosis and systematic-effect evidence (``tart-diagnose``).

Stage 7 reads a Stage-6 run, follows the recorded provenance back through
Stages 5, 4, 3 and 2 (verifying the recorded SHA-256 at every link), aligns
the outputs by time and station, and examines each layer separately:
leave-one-station-out and pairwise re-triangulation (reusing the Stage-4
solver), station direction residuals against the network-derived Stage-5
trajectory, calibration / timing / geometry / trajectory-model correlations,
and TLE-reference quality. It assigns non-exclusive evidence tags and a
rule-based evidence level (0-5) per epoch and overall. These are evidence
patterns, never physical-cause labels, and no upstream output is modified.
"""

from tart_gnss.stage7.errors import (  # noqa: F401
    AlignmentError,
    CorrelationAnalysisError,
    DiagnosisError,
    InfluenceAnalysisError,
    InvalidStage6ManifestError,
    ProvenanceChainError,
    Stage7Error,
)
