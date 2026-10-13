"""Stage 3: per-station GNSS satellite direction estimation from TART visibilities.

The primary estimator fits a point-source model directly to the complex
visibilities (:mod:`tart_gnss.stage3.estimator`); a dirty-image peak finder
(:mod:`tart_gnss.stage3.imaging`) is an independent validation path. Each
station is processed independently; nothing here combines stations.
"""

from tart_gnss.stage3.errors import (  # noqa: F401
    CoordinateTransformError,
    DirectionFitError,
    InsufficientBaselinesError,
    InvalidManifestError,
    InvalidVisibilityDataError,
    MissingMeasurementSetError,
    Stage3Error,
)
