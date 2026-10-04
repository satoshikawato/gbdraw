"""Keep complete-record browser regressions in the required PR contract suite."""

import subprocess
from pathlib import Path

import pytest


@pytest.mark.browser
# One pytest test runs the whole Playwright suite (3-5 minutes on CI); each
# Playwright test keeps its own timeout.
@pytest.mark.timeout(900)
def test_linear_comparison_browser_contracts():
    result = subprocess.run(
        ["npm", "run", "test:web:comparison-contracts", "--", "--workers=2"],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
