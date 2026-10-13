"""tart_gnss: the TART GNSS satellite-localisation pipeline (Stages 1-5).

Stage 1 predicts GNSS satellite positions from TLEs (SGP4 via Skyfield) and
identifies the time windows during which a satellite is simultaneously above a
minimum elevation at several TART radio-telescope stations.

Stage 2 (:mod:`tart_gnss.stage2`) retrieves the corresponding TART
visibilities and prepares one Measurement Set per station. Stage 3
(:mod:`tart_gnss.stage3`) measures the satellite direction at each station
independently from the complex visibilities. Stage 4
(:mod:`tart_gnss.stage4`) intersects the measured lines of sight to estimate
the 3-D (ITRS/ECEF) position with its covariance. Stage 5
(:mod:`tart_gnss.stage5`) turns those positions into a continuous
position/velocity trajectory (GCRS, full 6x6 covariance) with an extended
Kalman filter and RTS smoother, independently of any TLE. TLE comparison
(Stage 6) is not implemented here.

Every public command is named ``tart-*`` (see ``tests/test_cli_names.py``).
"""

__version__ = "0.5.0"

__all__ = ["__version__"]
