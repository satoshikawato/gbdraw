"""Python's side of the shared ring drawn-position vectors.

``tests/fixtures/ring_drawn_position_vectors.json``
(``tools/generate_ring_drawn_position_vectors.py``) holds the series places of
the comparison rings Python draws, outermost first, for a few stacks. The stored
file must equal what ``_drawn_ring_colors`` gives now, so a change to Python's
drawn order fails here until the vectors and the Web follow it
(``tests/web/ring-drawn-position-vectors.test.mjs``).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS = REPO_ROOT / "tests" / "fixtures" / "ring_drawn_position_vectors.json"
GENERATOR = REPO_ROOT / "tools" / "generate_ring_drawn_position_vectors.py"


def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_ring_drawn_position_vectors", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_stored_vectors_are_what_python_draws_now() -> None:
    assert VECTORS.read_text(encoding="utf-8") == _generator().render()
