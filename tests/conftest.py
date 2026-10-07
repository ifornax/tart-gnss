from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture
def sample_tle_path() -> Path:
    """Four real GNSS TLEs (GPS, Galileo, BeiDou GEO, GLONASS) from CelesTrak, 2026-10-07."""
    return DATA / "sample_gnss.tle"


@pytest.fixture
def station_yaml(tmp_path) -> Path:
    """Station file with three synthetic test stations (NOT real TART sites)."""
    path = tmp_path / "stations.yaml"
    path.write_text(
        """
stations:
  alpha:
    name: "Test station alpha"
    latitude_deg: -20.0
    longitude_deg: 57.0
    altitude_m: 0.0
    source: "synthetic test value"
  beta:
    name: "Test station beta"
    latitude_deg: -19.0
    longitude_deg: 47.0
    altitude_m: 1000.0
    source: "synthetic test value"
  gamma:
    name: "Test station gamma"
    tart_id: "zz-gamma"
    latitude_deg: -33.0
    longitude_deg: 26.0
    altitude_m: 500.0
    source: "synthetic test value"
""",
        encoding="utf-8",
    )
    return path
