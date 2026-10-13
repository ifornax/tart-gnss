"""tart_gnss: the TART GNSS satellite-localisation pipeline (Stages 1-3).

Stage 1 predicts GNSS satellite positions from TLEs (SGP4 via Skyfield) and
identifies the time windows during which a satellite is simultaneously above a
minimum elevation at several TART radio-telescope stations.

Stage 2 (:mod:`tart_gnss.stage2`) retrieves the corresponding TART
visibilities and prepares one Measurement Set per station. Stage 3
(:mod:`tart_gnss.stage3`) measures the satellite direction at each station
independently from the complex visibilities. Multi-station 3-D localisation
(Stage 4 onwards) is not implemented here.
"""

__version__ = "0.3.0"

__all__ = ["__version__"]
