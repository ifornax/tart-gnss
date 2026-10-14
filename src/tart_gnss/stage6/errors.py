"""Exceptions raised by Stage 6. Each one names a specific failure."""


class Stage6Error(Exception):
    """Base class for all Stage-6 errors."""


class InvalidStage5ManifestError(Stage6Error):
    """The Stage-5 manifest is unreadable, inconsistent, or not usable for validation."""


class MissingTrajectoryError(Stage6Error):
    """The Stage-5 trajectory table is missing, unreadable, or lacks required columns."""


class TLEReferenceError(Stage6Error):
    """The reference TLE cannot be found, read or verified."""


class TLENoradMismatchError(TLEReferenceError):
    """The supplied TLE belongs to a different NORAD object than the Stage-5 trajectory."""


class TLEPropagationError(Stage6Error):
    """SGP4 propagation of the reference TLE failed."""


class FrameTransformError(Stage6Error):
    """A TEME/GCRS/ITRS transformation failed or is inconsistent."""


class InvalidStateCovarianceError(Stage6Error):
    """A Stage-5 state covariance is non-finite or not positive definite."""


class RTNGeometryError(Stage6Error):
    """The RTN basis is undefined (zero position, or position parallel to velocity)."""


class ReferenceCovarianceError(Stage6Error):
    """The empirical TLE reference covariance cannot be formed."""
