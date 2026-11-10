"""Python's side of the shared feature density row vectors.

``tests/fixtures/feature_density_slot_vectors.json``
(``tools/generate_feature_density_slot_vectors.py``) holds the rows Python gives
a stack with feature density tracks. The stored file must equal what the
Python owners give now, so a change to ``with_feature_density_slots`` or to a
default stack fails here until the vectors and the Web follow it
(``tests/web/feature-density-slot-vectors.test.mjs``).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = REPO_ROOT / "tests" / "fixtures" / "feature_density_slot_vectors.json"
GENERATOR = REPO_ROOT / "tools" / "generate_feature_density_slot_vectors.py"


def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_feature_density_slot_vectors", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_stored_vectors_are_what_python_gives_now() -> None:
    assert VECTORS.read_text(encoding="utf-8") == _generator().render()


def test_the_vectors_cover_both_modes_and_every_track_count() -> None:
    cases = json.loads(VECTORS.read_text(encoding="utf-8"))["cases"]
    assert {(case["mode"], len(case["trackIds"])) for case in cases} == {
        (mode, count) for mode in ("circular", "linear") for count in (0, 1, 2)
    }
    assert any(case["kind"] == "stack" and case["axisIndex"] is not None for case in cases)
