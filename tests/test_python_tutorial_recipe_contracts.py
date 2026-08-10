from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from docs.recipes.run_python_scenarios import SCENARIO_IDS


pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = "docs/recipes/run_python_scenarios.py"
TUTORIAL_SCENARIO_IDS = tuple(
    scenario_id
    for scenario_id in SCENARIO_IDS
    # Full-data LOSATP recipes remain available through the manual scenario runner.
    if scenario_id.startswith("T-PY-")
    and scenario_id not in {"T-PY-01", "T-PY-05", "T-PY-07"}
)


@pytest.mark.parametrize(
    "scenario_id",
    TUTORIAL_SCENARIO_IDS,
)
def test_python_tutorial_recipe_regenerates_from_a_clean_external_context(
    scenario_id: str,
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not existing_pythonpath
        else os.pathsep.join((str(REPO_ROOT), existing_pythonpath))
    )

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / RUNNER),
            "--scenario",
            scenario_id,
            "--check",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert f"{scenario_id}: verified" in result.stdout
    assert list(tmp_path.iterdir()) == []
