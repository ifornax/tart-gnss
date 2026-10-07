"""Stage 2: TART observation retrieval, synchronisation and Measurement-Set preparation.

Stage 2 takes one Stage-1 common-visibility window and then:

1. finds the matching native TART visibility files in each station's archive;
2. downloads them with SHA-256 provenance;
3. checks their actual timestamps and metadata;
4. measures the real three-station temporal overlap and coverage;
5. converts each station independently to a CASA Measurement Set with ``tart2ms``;
6. writes ``manifest.json``, the interface to Stage 3.

Stage 2 does not image, detect sources or infer positions.
"""

from tart_gnss.stage2.errors import (  # noqa: F401
    ArchiveUnavailableError,
    DownloadError,
    MeasurementSetConversionError,
    NoCommonTimeError,
    ObservationNotFoundError,
    ObservationValidationError,
    Stage1SelectionError,
    Stage2Error,
)
