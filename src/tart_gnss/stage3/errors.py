"""Exceptions raised by Stage 3. Each one names a specific failure."""


class Stage3Error(Exception):
    """Base class for all Stage-3 errors."""


class InvalidManifestError(Stage3Error):
    """The Stage-2 manifest is unreadable, inconsistent, or not valid for the request."""


class MissingMeasurementSetError(Stage3Error):
    """A Measurement Set referenced by the Stage-2 manifest cannot be found."""


class InvalidVisibilityDataError(Stage3Error):
    """A Measurement Set is unreadable, or its contents violate the adopted conventions."""


class InsufficientBaselinesError(Stage3Error):
    """Too few valid baselines remain to estimate a direction."""


class DirectionFitError(Stage3Error):
    """The direction optimisation could not be carried out."""


class CoordinateTransformError(Stage3Error):
    """A coordinate conversion was requested for an invalid direction."""
