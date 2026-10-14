"""tart_gnss: the TART GNSS satellite-localisation pipeline (Stages 1-7).

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
Kalman filter and RTS smoother, independently of any TLE. Stage 6
(:mod:`tart_gnss.stage6`) is the first stage where a TLE enters: it
propagates a selected element set with SGP4 to the Stage-5 epochs and reports
``observed - TLE`` residuals (Cartesian, RTN) against the Stage-5 uncertainty.
Stage 7 (:mod:`tart_gnss.stage7`) traces the run back through Stages 2-6 and
reports evidence about where a discrepancy comes from (station, instrument,
timing, geometry, trajectory model, reference) with transparent rules, never
asserting a physical cause.

Every public command is named ``tart-*`` (see ``tests/test_cli_names.py``).
"""

__version__ = "0.7.0"

__all__ = ["__version__"]
