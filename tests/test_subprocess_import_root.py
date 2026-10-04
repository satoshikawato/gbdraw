"""Subprocesses started by the tests must import the checkout under test."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_subprocess_imports_gbdraw_from_checkout_under_test(tmp_path: Path) -> None:
    """A subprocess outside the repository root must not import another checkout."""

    result = subprocess.run(
        [sys.executable, "-c", "import gbdraw; print(gbdraw.__file__)"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    imported = Path(result.stdout.strip()).resolve()
    assert imported.is_relative_to(PROJECT_ROOT), (
        f"subprocess imported {imported}, outside {PROJECT_ROOT}"
    )
