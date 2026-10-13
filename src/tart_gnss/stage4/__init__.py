"""Stage 4: multi-station 3-D localisation from measured TART lines of sight.

Stage 4 reads the per-station directions measured by Stage 3
(``tart-localise``), matches them in time, converts each to an Earth-fixed
(ITRS/ECEF) sight line from the station's frozen coordinates, and intersects
the lines in a least-squares sense (:mod:`tart_gnss.stage4.triangulation`).
The primary estimator weights each line by its Stage-3 angular covariance
scaled to the station range. TLE predictions never enter the fit.
"""

from tart_gnss.stage4.errors import (  # noqa: F401
    CovarianceError,
    InsufficientStationsError,
    InvalidStage3ManifestError,
    MissingDirectionDataError,
    SingularGeometryError,
    Stage4Error,
    TimeMatchError,
    TriangulationError,
)
