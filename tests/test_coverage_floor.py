"""Line-coverage floor of the ``gbdraw`` package.

``Core (Python 3.14)`` runs its tests under pytest-cov and fails when the total
line coverage of ``gbdraw/`` (without ``gbdraw/web/``) falls below
``fail_under`` in ``[tool.coverage.report]`` of ``pyproject.toml``
(``docs/internal/SELECTIVE_CI.md``, "Line coverage"). ``MAX_UNCOVERED_PERCENT``
holds that floor as the share of lines that may stay uncovered. It is a
registered design-rule allowlist (``tools/web-design-rule-guards.json``,
COVERAGE-FLOOR, kind ``py-count``), which reads a raised count as a looser
rule: raising the floor lowers it and is free, and lowering the floor raises
it, an authority-only change that needs the Owner's approval.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


PROJECT_ROOT = Path(__file__).resolve().parent.parent
TESTS_WORKFLOW = PROJECT_ROOT / ".github" / "workflows" / "test.yml"
MAX_UNCOVERED_PERCENT = 7


def test_the_coverage_floor_is_the_registered_one() -> None:
    pyproject = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    floor = pyproject["tool"]["coverage"]["report"]["fail_under"]
    assert floor == 100 - MAX_UNCOVERED_PERCENT, (
        f"fail_under is {floor}; MAX_UNCOVERED_PERCENT in this file allows {100 - MAX_UNCOVERED_PERCENT}. "
        "Raising the floor lowers MAX_UNCOVERED_PERCENT; lowering it needs the Owner's approval."
    )


def test_the_core_job_checks_the_floor() -> None:
    workflows = sorted((PROJECT_ROOT / ".github" / "workflows").glob("*.y*ml"))
    overriding = [path.name for path in workflows if "--cov-fail-under" in path.read_text(encoding="utf-8")]
    assert not overriding, f"{overriding} pass --cov-fail-under, which replaces fail_under in pyproject.toml"
    core = re.search(r"\n  core:\n[\s\S]*?(?=\n  [A-Za-z0-9_-]+:\n)", TESTS_WORKFLOW.read_text(encoding="utf-8"))
    step = core and re.search(r"- name: Run core tests\n(?P<run>[\s\S]*?)\n\n", core[0])
    assert step, "test.yml has no `Run core tests` step in its `core` job"
    assert "--cov " in step["run"], "the `core` job runs its tests without --cov"
