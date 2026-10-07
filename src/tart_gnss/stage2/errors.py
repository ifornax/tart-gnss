"""Exceptions raised by Stage 2. Each one names a specific failure."""


class Stage2Error(Exception):
    """Base class for all Stage-2 errors."""


class Stage1SelectionError(Stage2Error):
    """The Stage-1 results cannot be read, or the selection is invalid or ambiguous."""


class ArchiveUnavailableError(Stage2Error):
    """The TART archive could not be reached or returned an error."""


class ObservationNotFoundError(Stage2Error):
    """No archived observation intersects the requested interval."""


class DownloadError(Stage2Error):
    """A file could not be retrieved, or failed an integrity or provenance check."""


class ObservationValidationError(Stage2Error):
    """A native observation file is unreadable or internally inconsistent."""


class NoCommonTimeError(Stage2Error):
    """The station data intervals have no common (non-zero) overlap."""


class MeasurementSetConversionError(Stage2Error):
    """tart2ms failed, or the Measurement Set it produced failed validation."""
