"""Stored JSON files to check against a description of ``gbdraw.contract_shapes``."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def stored_json_files() -> list[Path]:
    """The JSON files of the folders that hold Sessions and requests: fixtures,
    test inputs, documentation Sessions and the Gallery Sessions."""

    return sorted(
        path
        for folder in (
            "tests/fixtures", "tests/web/fixtures", "tests/test_inputs", "docs/images",
            "gbdraw/web/gallery/sessions",
        )
        for path in (REPO_ROOT / folder).rglob("*.json*")
        if path.suffix in {".json", ".gz"}
    )


def load_json(path: Path) -> object:
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)
