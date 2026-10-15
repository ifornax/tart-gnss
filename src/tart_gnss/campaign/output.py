"""Campaign tables (Parquet + CSV), reporting tables, ``summary.json`` and ``campaign_manifest.json``."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional

from tart_gnss.campaign.schema import LIST_COLUMNS
from tart_gnss.stage4.output import _clean

CAMPAIGN_SCHEMA_VERSION = 1


def _csv_ready(df):
    c = df.copy()
    for col in c.columns:
        if col in LIST_COLUMNS or c[col].map(lambda v: isinstance(v, (list, dict))).any():
            c[col] = c[col].map(lambda v: ";".join(map(str, v)) if isinstance(v, list) else
                                (json.dumps(_clean(v)) if isinstance(v, dict) else v))
    return c


def write_table(df, out_dir: Path, stem: str, parquet: bool = True) -> Dict[str, Optional[str]]:
    """Write ``<stem>.csv`` and (if possible) ``<stem>.parquet``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{stem}.csv"
    _csv_ready(df).to_csv(p, index=False, float_format="%.10g")
    paths: Dict[str, Optional[str]] = {"csv": str(p), "parquet": None}
    if parquet:
        try:
            import pyarrow  # noqa: F401

            q = out_dir / f"{stem}.parquet"
            dfp = df.copy()
            for col in dfp.columns:
                if dfp[col].map(lambda v: isinstance(v, dict)).any():
                    dfp[col] = dfp[col].map(lambda v: json.dumps(_clean(v)) if isinstance(v, dict) else v)
            dfp.to_parquet(q, index=False)
            paths["parquet"] = str(q)
        except ImportError:
            pass
    return paths


def read_table(path):
    """Read a campaign table; ``;``-joined list columns become lists in CSV."""
    import pandas as pd

    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
        for col in LIST_COLUMNS:
            if col in df.columns:
                df[col] = df[col].map(lambda f: list(f) if f is not None else [])
        return df
    df = pd.read_csv(path)
    for col in LIST_COLUMNS:
        if col in df.columns:
            df[col] = df[col].map(lambda f: [x for x in str(f).split(";") if x] if isinstance(f, str) else [])
    return df


def write_json(doc, path: Path) -> str:
    """Indented JSON with NaN/Inf as null."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(doc), indent=2, default=str) + "\n", encoding="utf-8")
    return str(path)
