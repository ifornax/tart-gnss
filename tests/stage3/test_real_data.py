"""Optional real-data checks (skipped unless the 2026-06-25 Stage-2 Measurement Sets are present).

These are the non-circular checks of the visibility sign convention: synthetic
round trips use the same forward model as the estimator and cannot detect a
convention error, but a mirrored convention would put real satellites at
azimuth + 180 deg. Set ``TART_GNSS_STAGE2_DIR`` to the directory holding
``20260625/NORAD_40534/manifest.json`` (default: ``../sat_position/data/observations``).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from .helpers import HAS_CASACORE

ROOT = Path(os.environ.get("TART_GNSS_STAGE2_DIR",
                           Path(__file__).resolve().parents[3] / "sat_position" / "data" / "observations"))
MANIFEST = ROOT / "20260625" / "NORAD_40534" / "manifest.json"

pytestmark = pytest.mark.skipif(not (HAS_CASACORE and MANIFEST.is_file()),
                                reason="real 2026-06-25 TART Measurement Sets not available")


def _brightest(station, t_offset_s):
    from tart_gnss.stage3 import loader
    from tart_gnss.stage3.imaging import image_direction

    doc, _ = loader.load_stage2_manifest(MANIFEST)
    si = loader.station_input(doc, station, MANIFEST)
    ds = loader.load_station_ms(si.ms_path, station, si.latitude_deg, si.longitude_deg)
    ut = np.unique(ds.time_mjd_seconds)
    t = ut[np.argmin(np.abs(ut - ut[0] - t_offset_s))]
    sel = loader.select_samples(ds, np.flatnonzero(ds.time_mjd_seconds == t), si.zero_gain_antennas)
    return image_direction(sel.samples)  # blind: whole hemisphere


def _sep(az1, el1, az2, el2):
    from tart_gnss.stage3.geometry import angular_separation_deg, azel_to_unit

    return float(angular_separation_deg(azel_to_unit(az1, el1), azel_to_unit(az2, el2)))


def test_south_africa_brightest_source_is_target_not_mirror():
    img = _brightest("south_africa", 0.0)
    # TLE prediction for GPS PRN 26 (NORAD 40534) at 09:13:08.7 UTC: Az 176.5, El 71.0
    assert _sep(img.az_deg, img.el_deg, 176.5, 71.0) < 3.0
    assert _sep(img.az_deg, img.el_deg, 356.5, 71.0) > 20.0


def test_mauritius_brightest_source_is_prn28_not_mirror():
    img = _brightest("mauritius", 338.0)
    # GPS PRN 28 predicted at Az 106.2, El 66.6 (09:18:59 UTC); the mirror would be Az 286.
    assert _sep(img.az_deg, img.el_deg, 106.2, 66.6) < 3.0
    assert _sep(img.az_deg, img.el_deg, 286.2, 66.6) > 20.0
