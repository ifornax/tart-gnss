import pytest

from tart_gnss.stage3 import synthetic as syn


@pytest.fixture
def layout():
    """24-antenna synthetic layout, baselines up to ~4 m (not a real TART array)."""
    return syn.tart_like_layout()
