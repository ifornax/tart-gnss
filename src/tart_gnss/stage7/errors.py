"""Exceptions raised by Stage 7. Each one names a specific failure."""


class Stage7Error(Exception):
    """Base class for all Stage-7 errors."""


class InvalidStage6ManifestError(Stage7Error):
    """The Stage-6 manifest is unreadable or not a Stage-6 manifest."""


class ProvenanceChainError(Stage7Error):
    """An upstream manifest or table cannot be found or fails its recorded hash."""


class AlignmentError(Stage7Error):
    """Stage outputs cannot be aligned in time."""


class InfluenceAnalysisError(Stage7Error):
    """Leave-one-out or pairwise analysis failed."""


class CorrelationAnalysisError(Stage7Error):
    """A correlation analysis could not be computed."""


class DiagnosisError(Stage7Error):
    """Evidence assignment failed."""
