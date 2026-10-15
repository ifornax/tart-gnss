"""Exceptions raised by the campaign layer. Each one names a specific failure."""


class CampaignError(Exception):
    """Base class for all campaign errors."""


class ManifestDiscoveryError(CampaignError):
    """No Stage-7 manifest could be found or an input path is invalid."""


class InvalidStage7ManifestError(CampaignError):
    """A Stage-7 manifest is unreadable or not a Stage-7 manifest."""


class CampaignSchemaError(CampaignError):
    """A Stage-7 run uses an unsupported (e.g. future) schema version."""


class InsufficientCampaignDataError(CampaignError):
    """Too few usable runs for the requested analysis."""


class AggregationError(CampaignError):
    """Campaign tables could not be built."""


class BootstrapError(CampaignError):
    """A bootstrap confidence interval could not be computed."""
