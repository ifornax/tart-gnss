"""Stage 5: TLE-independent statistical trajectory tracking (``tart-track``).

Stage 5 reads the epoch-by-epoch 3-D positions and full covariances produced
by Stage 4 (``tart-triangulate``), rotates them from ITRS/ECEF into the
inertial GCRS frame (:mod:`tart_gnss.stage5.frames`), and estimates a
continuous position/velocity state with an extended Kalman filter
(:mod:`tart_gnss.stage5.filter`) and a Rauch-Tung-Striebel smoother
(:mod:`tart_gnss.stage5.smoother`) under constant-velocity, two-body or J2
dynamics (:mod:`tart_gnss.stage5.dynamics`). TLE predictions never enter the
estimate: not as initial state, prior, process model or regularisation target.
"""

from tart_gnss.stage5.errors import (  # noqa: F401
    DynamicsError,
    FilterDivergenceError,
    FrameTransformError,
    InitialisationError,
    InvalidCovarianceError,
    InvalidStage4ManifestError,
    MissingPositionDataError,
    SmootherError,
    Stage5Error,
)
