"""Python's side of the shared Legend row edit vectors.

``tests/fixtures/legend_row_edit_vectors.json`` (``tools/generate_legend_row_edit_vectors.py``)
holds, per case, the Legend row edits as config leaves, Python's Legend without
and with them, and the feature Legend rows Python draws with them. Python
redraws every case here and must draw the stored rows; the Web's executor must
make the same rows from the un-edited Legend
(``tests/web/legend-row-edit-vectors.test.mjs``).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = REPO_ROOT / "tests" / "fixtures" / "legend_row_edit_vectors.json"
GENERATOR = REPO_ROOT / "tools" / "generate_legend_row_edit_vectors.py"


def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_legend_row_edit_vectors", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


GENERATOR_MODULE = _generator()
STORED = {case["id"]: case for case in json.loads(VECTORS.read_text(encoding="utf-8"))["cases"]}


def test_the_vectors_hold_every_case_the_generator_draws() -> None:
    assert list(STORED) == [case.id for case in GENERATOR_MODULE.CASES]


@pytest.mark.parametrize("case", GENERATOR_MODULE.CASES, ids=lambda case: case.id)
def test_python_draws_the_stored_legend_rows(case, tmp_path: Path) -> None:
    drawn = GENERATOR_MODULE.vector(case, tmp_path)
    stored = STORED[case.id]
    assert drawn["expected"] == stored["expected"]
    assert drawn["unedited"] == stored["unedited"]
    assert drawn["edited"] == stored["edited"]
    assert drawn == stored

