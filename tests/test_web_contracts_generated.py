from __future__ import annotations

import functools
import importlib.util
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "tools" / "generate_web_contracts.py"


@functools.cache
def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_web_contracts", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The `browser` marker puts this test in the jobs a Web-only pull request runs.
@pytest.mark.browser
def test_generated_web_contracts_match_python_source() -> None:
    result = subprocess.run(
        [sys.executable, "tools/generate_web_contracts.py", "--check"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_a_changed_python_constant_makes_the_generated_module_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_io as session_io

    generator = _generator()
    current = generator.render_web_contracts_module()
    assert f"export const SESSION_VERSION = {session_io.CURRENT_SESSION_VERSION};\n" in current

    monkeypatch.setattr(session_io, "CURRENT_SESSION_VERSION", session_io.CURRENT_SESSION_VERSION + 1)

    assert generator.render_web_contracts_module() != current


def test_generated_web_contracts_module_imports_nothing() -> None:
    # `contracts/` is a leaf directory of the Web layers (R13).
    source = _generator().render_web_contracts_module()

    assert re.search(r"\bimport\b|\bfrom\s+['\"]", source) is None
