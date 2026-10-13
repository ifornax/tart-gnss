"""Exceptions raised by Stage 4. Each one names a specific failure."""


class Stage4Error(Exception):
    """Base class for all Stage-4 errors."""


class InvalidStage3ManifestError(Stage4Error):
    """The Stage-3 manifest is unreadable, inconsistent, or not usable for triangulation."""


class MissingDirectionDataError(Stage4Error):
    """The Stage-3 direction table is missing, unreadable, or lacks required columns."""


class TimeMatchError(Stage4Error):
    """Station direction estimates cannot be matched in time as requested."""


class InsufficientStationsError(Stage4Error):
    """Fewer stations than required contribute to a solution."""


class TriangulationError(Stage4Error):
    """The sight-line intersection could not be computed (invalid inputs)."""


class SingularGeometryError(TriangulationError):
    """The sight-line geometry is rank deficient (e.g. parallel lines)."""


class CovarianceError(Stage4Error):
    """A directional or position covariance is singular, non-finite or not positive definite."""
