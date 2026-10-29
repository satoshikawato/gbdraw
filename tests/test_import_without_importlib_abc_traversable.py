"""gbdraw must import where ``importlib.abc.Traversable`` no longer exists (Python 3.14+)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

_PROBE = """
import importlib.abc as abc_module

# Python 3.12 and 3.13 still serve the name through a module __getattr__; 3.14 has neither.
abc_module.__dict__.pop("Traversable", None)
abc_module.__getattr__ = lambda name: (_ for _ in ()).throw(AttributeError(name))
import gbdraw, gbdraw.cli
"""


@pytest.mark.skipif(sys.version_info < (3, 11), reason="3.10 reads Traversable from importlib.abc")
def test_gbdraw_imports_without_importlib_abc_traversable() -> None:
    result = subprocess.run(
        [sys.executable, "-c", _PROBE],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
