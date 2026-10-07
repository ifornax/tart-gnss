"""tart_gnss: Stage 1 of the TART GNSS satellite-localisation pipeline.

Stage 1 predicts GNSS satellite positions from TLEs (SGP4 via Skyfield) and
identifies the time windows during which a satellite is simultaneously above a
minimum elevation at several TART radio-telescope stations.

Later stages (not implemented here) will retrieve the corresponding TART
visibilities, image each station independently, measure the satellite
direction and infer the satellite position with Bayesian methods.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
