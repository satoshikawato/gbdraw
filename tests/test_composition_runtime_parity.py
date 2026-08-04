from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
NODE_SPECS = (
    "tests/web/composition-layout.test.mjs",
    "tests/web/composition-runtime-parity.test.mjs",
    "tests/web/legend-layout-actions.test.mjs",
)


def test_web_composition_contracts_and_python_runtime_parity() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    environment = os.environ.copy()
    environment["PYTHON"] = sys.executable
    result = subprocess.run(
        [node, "--test", *NODE_SPECS],
        cwd=REPO_ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
