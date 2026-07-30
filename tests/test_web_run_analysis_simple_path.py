from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_web_run_analysis_simple_path_integration() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run(
        [node, "tests/web/run-analysis-simple-path.test.mjs"],
        check=True,
        cwd=REPO_ROOT,
    )
