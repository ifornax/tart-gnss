"""Discovery of Stage-7 runs: explicit manifests and/or recursive search for ``stage7_manifest.json``.

Only files named exactly ``stage7_manifest.json`` are discovered; runs are
never inferred from CSV or Parquet files. Duplicate paths (after resolving)
are listed once.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

from tart_gnss.campaign.errors import ManifestDiscoveryError

MANIFEST_NAME = "stage7_manifest.json"


def discover_stage7_runs(input_root: Path) -> List[Path]:
    """Discover Stage-7 manifests recursively (sorted)."""
    root = Path(input_root)
    if not root.is_dir():
        raise ManifestDiscoveryError(f"input root is not a directory: {root}")
    return sorted(p for p in root.rglob(MANIFEST_NAME) if p.is_file())


def collect_manifests(input_roots: Optional[Sequence[Path]] = None,
                      manifests: Optional[Sequence[Path]] = None) -> List[Path]:
    """Union of discovered and explicit manifests, without duplicates."""
    out: List[Path] = []
    seen = set()
    for r in input_roots or []:
        for p in discover_stage7_runs(r):
            key = p.resolve()
            if key not in seen:
                seen.add(key)
                out.append(p)
    for m in manifests or []:
        p = Path(m)
        if not p.is_file():
            raise ManifestDiscoveryError(f"Stage-7 manifest not found: {p}")
        key = p.resolve()
        if key not in seen:
            seen.add(key)
            out.append(p)
    if not out:
        raise ManifestDiscoveryError("no stage7_manifest.json found in the given inputs")
    return out
