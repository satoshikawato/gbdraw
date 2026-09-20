"""Keep complete-record browser regressions in the required PR contract suite."""

import subprocess
from pathlib import Path

import pytest


@pytest.mark.browser
@pytest.mark.parametrize("shard", ("1/2", "2/2"))
def test_linear_comparison_browser_contracts(shard):
    result = subprocess.run(
        ["npm", "run", "test:web:comparison-contracts", "--", f"--shard={shard}"],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
