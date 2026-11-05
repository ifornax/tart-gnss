"""Archive-aware candidate discovery (``tart-archive-scan``).

Stage 1 predicts *geometric* common visibility; Stage 2 checks the archive for
one selected window. This package lists the TART visibility archive once per
``(station, UTC day)``, turns the filenames into conservative per-station
*candidate* coverage, intersects the stations, and matches the shared
candidate segments with every Stage-1 window.

Evidence classes are kept separate:

1. ``GEOMETRICALLY_VISIBLE``: Stage-1 SGP4 prediction only;
2. ``ARCHIVE_COVERAGE_CANDIDATE``: filename/object-metadata evidence of
   coverage, not verified from sample timestamps (all this scanner claims);
3. ``VALIDATED_THREE_STATION_DATA``: Stage-2 native validation;
4. ``RADIO_DETECTED``: Stage-3 detection, never inferred here.

The scanner never downloads native files and never runs ``tart-retrieve``.

Modules: :mod:`.inventory` (listing), :mod:`.cache` (station/day snapshots),
:mod:`.intervals` (candidate segments and intersections), :mod:`.matching`
(Stage-1 windows), :mod:`.ranking`, :mod:`.output` and :mod:`.cli`; network-wide
discovery (``--all-stations``) in :mod:`.network` (eligibility, physical sites,
geodesic baselines, sweep, triads) and :mod:`.network_run`.
"""

#: Evidence class of every scanner result (see the module docstring).
EVIDENCE_CLASS = "ARCHIVE_COVERAGE_CANDIDATE"

#: Coverage confidence of a filename-derived segment.
COVERAGE_CONFIDENCE = "FILENAME_CANDIDATE"
