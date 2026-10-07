#!/usr/bin/env python3
"""Convenience wrapper: ``python common_gnss.py ...`` == ``common-gnss ...``.

Requires the package to be installed (``pip install .``) or ``src/`` on the path.
"""

import sys
from pathlib import Path

try:
    from tart_gnss.cli import main
except ImportError:  # running from a source checkout without installation
    sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
    from tart_gnss.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
