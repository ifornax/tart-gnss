"""Elevation-threshold visibility logic, generic in the number of stations.

A satellite is visible from station i at time t when E_i(t) >= E_min, and is
*commonly* visible when this holds for every selected station.
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
from skyfield.api import EarthSatellite
from skyfield.toposlib import GeographicPosition

from tart_gnss import propagation

DEFAULT_MIN_ELEVATION_DEG = 20.0


def station_visibility(elevations_deg, min_elevation_deg: float = DEFAULT_MIN_ELEVATION_DEG) -> np.ndarray:
    """Per-station visibility flags ``E_i(t) >= E_min``.

    Parameters
    ----------
    elevations_deg
        Array-like of shape ``(n_stations, n_times)``, degrees.
    min_elevation_deg
        Elevation threshold, degrees.

    Returns
    -------
    numpy.ndarray of bool, same shape as the input. NaN elevations (failed
    propagation) are treated as not visible.
    """
    el = np.asarray(elevations_deg, dtype=float)
    with np.errstate(invalid="ignore"):
        return np.where(np.isnan(el), False, el >= min_elevation_deg)


def common_visibility(elevations_deg, min_elevation_deg: float = DEFAULT_MIN_ELEVATION_DEG) -> np.ndarray:
    """Simultaneous visibility: logical AND over stations of ``E_i(t) >= E_min``.

    Parameters
    ----------
    elevations_deg
        Array-like of shape ``(n_stations, n_times)`` (any number of stations).

    Returns
    -------
    numpy.ndarray of bool, shape ``(n_times,)``.
    """
    flags = station_visibility(elevations_deg, min_elevation_deg)
    if flags.ndim == 1:
        flags = flags[:, np.newaxis]
    return np.all(flags, axis=0)


def make_common_visibility_predicate(
    satellite: EarthSatellite,
    stations: Sequence[GeographicPosition],
    start,
    min_elevation_deg: float = DEFAULT_MIN_ELEVATION_DEG,
) -> Callable[[float], bool]:
    """Return ``f(offset_s) -> bool`` evaluating common visibility at one epoch.

    Used to refine window boundaries between coarse samples.

    Parameters
    ----------
    satellite
        Skyfield satellite.
    stations
        Station geographic positions.
    start
        Aware UTC datetime defining offset zero.
    min_elevation_deg
        Elevation threshold, degrees.
    """
    ts = propagation.get_timescale()

    def predicate(offset_s: float) -> bool:
        t = propagation.times_from_offsets(start, offset_s, ts)
        el = propagation.elevations(satellite, stations, t)
        return bool(common_visibility(el.reshape(-1, 1), min_elevation_deg)[0])

    return predicate
