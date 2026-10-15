"""Campaign analysis across completed Stage-7 runs (``tart-campaign``).

The campaign layer turns many completed single-pass pipeline runs (one
satellite x one common observing window each) into population statistics:
localisation and trajectory performance, TLE-relative residuals, station
stability, constellation, satellite and station-set comparisons, geometry,
elevation, SNR and TLE-age dependence, evidence-level and diagnostic-tag
rates and pass-to-pass repeatability. It consumes completed outputs only:
it never runs Stages 1-7. Runs, not epochs, are the independent units of
every inference (epochs within a pass are strongly correlated).
"""

from tart_gnss.campaign.errors import (  # noqa: F401
    AggregationError,
    BootstrapError,
    CampaignError,
    CampaignSchemaError,
    InsufficientCampaignDataError,
    InvalidStage7ManifestError,
    ManifestDiscoveryError,
)
