"""tart_gnss: the TART GNSS satellite-localisation pipeline (Stages 1-4).

Stage 1 predicts GNSS satellite positions from TLEs (SGP4 via Skyfield) and
identifies the time windows during which a satellite is simultaneously above a
minimum elevation at several TART radio-telescope stations.

Stage 2 (:mod:`tart_gnss.stage2`) retrieves the corresponding TART
visibilities and prepares one Measurement Set per station. Stage 3
(:mod:`tart_gnss.stage3`) measures the satellite direction at each station
independently from the complex visibilities. Stage 4
(:mod:`tart_gnss.stage4`) intersects the measured lines of sight to estimate
the 3-D (ITRS/ECEF) position with its covariance. Trajectory inference and
TLE comparison (Stages 5-6) are not implemented here.

Every public command is named ``tart-*`` (see ``tests/test_cli_names.py``).
"""

__version__ = "0.4.0"

__all__ = ["__version__"]
