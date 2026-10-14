"""Stage 6: TLE trajectory validation and residual analysis (``tart-validate``).

Stage 6 is the first stage in which a TLE enters the scientific comparison.
It propagates a selected element set with SGP4 to the exact epochs of the
TLE-independent Stage-5 trajectory, transforms it into the same GCRS frame,
and reports ``observed - TLE`` residuals in Cartesian and radial/along-track/
cross-track (RTN) components, together with the Stage-5 measurement
uncertainty rotated into RTN, measurement-conditional Mahalanobis distances,
optional empirical reference dispersion from an ensemble of nearby TLEs, and
TLE-age based precision-validation status. The TLE is never treated as exact
truth, never fitted to the observations, and the Stage-5 trajectory is never
modified.
"""

from tart_gnss.stage6.errors import (  # noqa: F401
    FrameTransformError,
    InvalidStage5ManifestError,
    InvalidStateCovarianceError,
    MissingTrajectoryError,
    ReferenceCovarianceError,
    RTNGeometryError,
    Stage6Error,
    TLENoradMismatchError,
    TLEPropagationError,
    TLEReferenceError,
)
