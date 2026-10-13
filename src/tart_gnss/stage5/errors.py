"""Exceptions raised by Stage 5. Each one names a specific failure."""


class Stage5Error(Exception):
    """Base class for all Stage-5 errors."""


class InvalidStage4ManifestError(Stage5Error):
    """The Stage-4 manifest is unreadable, inconsistent, or not usable for tracking."""


class MissingPositionDataError(Stage5Error):
    """The Stage-4 position table is missing, unreadable, or lacks required columns."""


class InvalidCovarianceError(Stage5Error):
    """A measurement or state covariance is non-finite, asymmetric or not positive definite."""


class FrameTransformError(Stage5Error):
    """An ITRS/ECEF <-> GCRS transformation could not be computed."""


class InitialisationError(Stage5Error):
    """No TLE-independent initial state can be formed from the Stage-4 measurements."""


class DynamicsError(Stage5Error):
    """State propagation (ODE integration) failed or produced non-finite values."""


class FilterDivergenceError(Stage5Error):
    """The forward filter produced a non-finite or non-positive-semidefinite covariance."""


class SmootherError(Stage5Error):
    """The backward RTS smoother failed."""
