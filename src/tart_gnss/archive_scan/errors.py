"""Exceptions raised by the archive scanner. Each one names a specific failure."""


class ArchiveScanError(Exception):
    """Base class for all archive-scan errors."""


class ArchiveScanInputError(ArchiveScanError):
    """Invalid arguments, station selection, interval or Stage-1 input (exit 2)."""


class InventoryCacheError(ArchiveScanError):
    """A cached station/day listing is missing, corrupt or has an incompatible schema."""
