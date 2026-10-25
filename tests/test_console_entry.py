"""The installed ``gbdraw`` command sets ``OPENBLAS_NUM_THREADS`` before NumPy loads; the Python API does not."""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parents[1]

# A meta-path spy that records the variable at the moment NumPy is first searched for.
_SPY = textwrap.dedent(
    """
    import os, sys

    class _Spy:
        seen = []
        def find_spec(self, name, path=None, target=None):
            if name == "numpy" and not self.seen:
                self.seen.append(os.environ.get("OPENBLAS_NUM_THREADS"))
            return None

    sys.meta_path.insert(0, _Spy)
    """
)


def _run(body: str, *, openblas: str | None) -> str:
    env = {key: value for key, value in os.environ.items() if key != "OPENBLAS_NUM_THREADS"}
    env["PYTHONPATH"] = str(REPO_ROOT)
    if openblas is not None:
        env["OPENBLAS_NUM_THREADS"] = openblas
    script = _SPY.replace("sys.meta_path.insert(0, _Spy)", "sys.meta_path.insert(0, _Spy())") + body
    done = subprocess.run(
        [sys.executable, "-c", script], env=env, cwd=REPO_ROOT, capture_output=True, text=True, check=False
    )
    assert done.returncode == 0, done.stderr
    return done.stdout.strip().splitlines()[-1]


_CONSOLE = textwrap.dedent(
    """
    import sys
    sys.argv = ["gbdraw", "--version"]
    import _gbdraw_console
    try:
        _gbdraw_console.main()
    except SystemExit as exc:
        assert exc.code in (0, None), exc.code
    print(_Spy.seen)
    """
)


@pytest.mark.parametrize(("preset", "expected"), [(None, "['1']"), ("3", "['3']")])
def test_console_entry_sets_openblas_threads_only_when_unset(preset: str | None, expected: str) -> None:
    assert _run(_CONSOLE, openblas=preset) == expected


def test_python_api_import_leaves_openblas_threads_alone() -> None:
    body = "import gbdraw\ngbdraw.Diagram\nprint(_Spy.seen)\n"
    assert _run(body, openblas=None) == "[None]"


def test_installed_command_and_packaging_name_the_console_module() -> None:
    config = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert config["project"]["scripts"]["gbdraw"] == "_gbdraw_console:main"
    assert "_gbdraw_console" in config["tool"]["setuptools"]["py-modules"]
