"""Public command names follow the permanent ``tart-*`` convention (a stable public API)."""

import importlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = {
    "tart-common-gnss": "tart_gnss.cli:main",
    "tart-retrieve": "tart_gnss.stage2.cli:main",
    "tart-stations-update": "tart_gnss.registry_cli:main_update",
    "tart-stations-list": "tart_gnss.registry_cli:main_list",
    "tart-localise": "tart_gnss.stage3.cli:main",
    "tart-triangulate": "tart_gnss.stage4.cli:main",
}
LEGACY = ("common-gnss", "retrieve-tart", "update-tart-stations", "list-tart-stations")


def _entry_points():
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    block = re.search(r"^\[project\.scripts\]\n(.*?)(?=^\[|\Z)", text, re.S | re.M).group(1)
    return dict(re.findall(r'^([\w.-]+)\s*=\s*"([^"]+)"', block, re.M))


def test_entry_points_are_canonical():
    eps = _entry_points()
    assert eps == CANONICAL
    assert all(name.startswith("tart-") for name in eps)


def test_prog_matches_entry_point_name():
    for name, target in CANONICAL.items():
        module, func = target.split(":")
        src = Path(importlib.util.find_spec(module).origin).read_text(encoding="utf-8")
        assert f'prog="{name}"' in src, (name, module)


def test_no_legacy_names_in_active_text():
    # Historical provenance (results/, docs/examples/, the report's note on renamed commands) is exempt.
    allowed = re.compile(r"old names|former names|provenance|[Ll]egacy")
    pattern = re.compile(r"(?<![\w-])(" + "|".join(map(re.escape, LEGACY)) + r")(?![\w-])")
    files = [ROOT / "README.md", ROOT / "pyproject.toml", *ROOT.glob("*.md"), *(ROOT / "docs").glob("*.md"),
             *(ROOT / "src").rglob("*.py"), *(ROOT / "tests").rglob("*.py"), *(ROOT / "scripts").glob("*.py")]
    bad = []
    for f in set(files):
        if f.name == "test_cli_names.py":
            continue
        lines = f.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines):
            context = " ".join(lines[max(0, i - 6): i + 6])
            if pattern.search(line) and not allowed.search(context):
                bad.append(f"{f.relative_to(ROOT)}:{i + 1}: {line.strip()}")
    assert not bad, "\n".join(bad)
