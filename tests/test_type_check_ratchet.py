"""Type-check ratchet for the ``gbdraw`` package.

mypy, pinned by the ``typecheck`` extra, checks ``gbdraw/`` except
``gbdraw/web/`` with the ``[tool.mypy]`` settings in ``pyproject.toml``. A
file's type debt is its number of mypy errors plus its number of
``# type: ignore`` comments. ``TYPE_DEBT_BASELINE`` records the debt of every
file that has any. An entry may only go down; a file that is not listed,
including a new file, has no debt. Suppressions that would hide errors from
this count are rejected. Raising an entry or relaxing the configuration needs
the Owner's approval.
Plan: ``docs/internal/PYTHON_TYPE_CHECK_RATCHET_PLAN_2026-10-06.md``.
"""

from __future__ import annotations

import io
import json
import os
import re
import subprocess
import sys
import tokenize
from importlib import metadata
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = PROJECT_ROOT / "pyproject.toml"
TYPE_IGNORE = re.compile(r"#\s*type:\s*ignore\b")
MYPY_INLINE_CONFIG = re.compile(r"#\s*mypy:")

EXPECTED_MYPY_CONFIG = {
    "python_version": "3.10",
    "platform": "linux",
    "files": ["gbdraw"],
    "exclude": ["^gbdraw/web/"],
    "no_site_packages": True,
    "check_untyped_defs": True,
    "overrides": [
        {
            "module": [
                "BCBio.*",
                "Bio.*",
                "fontTools.*",
                "numpy.*",
                "pandas.*",
                "svgwrite.*",
                "tomli",
            ],
            "ignore_missing_imports": True,
        }
    ],
}

TYPE_DEBT_BASELINE: dict[str, int] = {
    "gbdraw/analysis/collinearity.py": 10,
    "gbdraw/analysis/collinearity_units.py": 2,
    "gbdraw/analysis/conservation.py": 9,
    "gbdraw/analysis/depth.py": 3,
    "gbdraw/analysis/depth_tracks.py": 7,
    "gbdraw/analysis/gc.py": 4,
    "gbdraw/analysis/ortholog_paths.py": 14,
    "gbdraw/analysis/protein_artifacts.py": 8,
    "gbdraw/analysis/protein_colinearity.py": 60,
    "gbdraw/annotations/feature_underlays.py": 1,
    "gbdraw/annotations/io.py": 6,
    "gbdraw/annotations/models.py": 8,
    "gbdraw/annotations/planning.py": 2,
    "gbdraw/annotations/resolve.py": 5,
    "gbdraw/api/__init__.py": 3,
    "gbdraw/api/config.py": 3,
    "gbdraw/api/diagram.py": 48,
    "gbdraw/api/io.py": 9,
    "gbdraw/api/options.py": 17,
    "gbdraw/api/prepared.py": 3,
    "gbdraw/api/record_planning.py": 10,
    "gbdraw/api/render.py": 3,
    "gbdraw/api/request_render.py": 89,
    "gbdraw/api/requests.py": 3,
    "gbdraw/api/session_compat.py": 17,
    "gbdraw/api/tracks.py": 1,
    "gbdraw/canvas/circular.py": 5,
    "gbdraw/canvas/linear.py": 3,
    "gbdraw/circular.py": 16,
    "gbdraw/cli.py": 1,
    "gbdraw/cli_utils/common.py": 1,
    "gbdraw/cli_utils/session.py": 7,
    "gbdraw/comparisons/circular_losat.py": 1,
    "gbdraw/comparisons/linear_losat.py": 1,
    "gbdraw/comparisons/losat_jobs.py": 8,
    "gbdraw/config/models/canvas.py": 1,
    "gbdraw/config/models/labels.py": 4,
    "gbdraw/config/models/objects.py": 3,
    "gbdraw/config/models/root.py": 1,
    "gbdraw/config/modify.py": 5,
    "gbdraw/config/toml.py": 1,
    "gbdraw/configurators/blast.py": 2,
    "gbdraw/configurators/depth.py": 1,
    "gbdraw/configurators/features.py": 2,
    "gbdraw/configurators/gc.py": 1,
    "gbdraw/configurators/legend.py": 19,
    "gbdraw/core/record_metadata.py": 4,
    "gbdraw/definition_line_styles.py": 1,
    "gbdraw/diagrams/circular/assemble.py": 39,
    "gbdraw/diagrams/circular/builders.py": 8,
    "gbdraw/diagrams/circular/positioning.py": 2,
    "gbdraw/diagrams/circular/presets.py": 6,
    "gbdraw/diagrams/circular/radial_layout.py": 34,
    "gbdraw/diagrams/linear/assemble.py": 30,
    "gbdraw/diagrams/linear/builders.py": 12,
    "gbdraw/diagrams/linear/orthogroup_alignment.py": 7,
    "gbdraw/diagrams/linear/positioning.py": 2,
    "gbdraw/diagrams/linear/precalc.py": 13,
    "gbdraw/diagrams/linear/track_slots.py": 3,
    "gbdraw/features/coordinates.py": 11,
    "gbdraw/features/factory.py": 3,
    "gbdraw/features/objects.py": 1,
    "gbdraw/features/overrides.py": 10,
    "gbdraw/features/placement.py": 18,
    "gbdraw/features/source.py": 1,
    "gbdraw/features/tracks.py": 2,
    "gbdraw/interface.py": 19,
    "gbdraw/io/cli_tables.py": 14,
    "gbdraw/io/genome.py": 3,
    "gbdraw/io/record_select.py": 1,
    "gbdraw/io/regions.py": 2,
    "gbdraw/labels/circular.py": 6,
    "gbdraw/labels/circular_candidates.py": 1,
    "gbdraw/labels/circular_radial.py": 2,
    "gbdraw/labels/circular_types.py": 1,
    "gbdraw/labels/coordinates.py": 1,
    "gbdraw/labels/linear.py": 14,
    "gbdraw/labels/policy.py": 1,
    "gbdraw/layout/composition.py": 7,
    "gbdraw/layout/linear.py": 4,
    "gbdraw/layout/linear_multi_record.py": 1,
    "gbdraw/layout/record_coordinates.py": 5,
    "gbdraw/layout/record_placement.py": 1,
    "gbdraw/layout/scalar_axis.py": 1,
    "gbdraw/layout/similarity_alignment.py": 3,
    "gbdraw/layout/spatial.py": 1,
    "gbdraw/layout/text_geometry.py": 3,
    "gbdraw/legend/circular_layout.py": 1,
    "gbdraw/legend/table.py": 2,
    "gbdraw/linear.py": 15,
    "gbdraw/linear_comparison.py": 18,
    "gbdraw/losat_setup.py": 9,
    "gbdraw/render/drawers/circular/annotations.py": 4,
    "gbdraw/render/drawers/circular/conservation.py": 5,
    "gbdraw/render/drawers/circular/definition.py": 2,
    "gbdraw/render/drawers/circular/features.py": 1,
    "gbdraw/render/drawers/circular/gc_content.py": 1,
    "gbdraw/render/drawers/circular/labels.py": 9,
    "gbdraw/render/drawers/linear/annotations.py": 6,
    "gbdraw/render/drawers/linear/features.py": 1,
    "gbdraw/render/groups/circular/axis.py": 2,
    "gbdraw/render/groups/circular/conservation.py": 6,
    "gbdraw/render/groups/circular/definition.py": 4,
    "gbdraw/render/groups/circular/depth.py": 6,
    "gbdraw/render/groups/circular/gc_content.py": 8,
    "gbdraw/render/groups/circular/gc_skew.py": 8,
    "gbdraw/render/groups/circular/labels.py": 8,
    "gbdraw/render/groups/circular/legend.py": 20,
    "gbdraw/render/groups/circular/seq_record.py": 12,
    "gbdraw/render/groups/circular/ticks.py": 7,
    "gbdraw/render/groups/linear/definition.py": 5,
    "gbdraw/render/groups/linear/depth.py": 6,
    "gbdraw/render/groups/linear/gc_content.py": 10,
    "gbdraw/render/groups/linear/gc_skew.py": 4,
    "gbdraw/render/groups/linear/legend.py": 2,
    "gbdraw/render/groups/linear/length_bar.py": 2,
    "gbdraw/render/groups/linear/pairwise_match.py": 4,
    "gbdraw/render/groups/linear/plot_title.py": 1,
    "gbdraw/render/groups/linear/seq_record.py": 2,
    "gbdraw/render/interactive_context.py": 9,
    "gbdraw/render/interactive_svg.py": 15,
    "gbdraw/render/label_binding.py": 1,
    "gbdraw/render/patterns.py": 2,
    "gbdraw/session.py": 2,
    "gbdraw/session_io.py": 28,
    "gbdraw/session_request_codec.py": 56,
    "gbdraw/svg/circular_ticks.py": 1,
    "gbdraw/svg/circular_tracks.py": 1,
    "gbdraw/svg/linear_tracks.py": 3,
    "gbdraw/tracks/__init__.py": 4,
    "gbdraw/web_support/config_overrides.py": 1,
    "gbdraw/web_support/error_adapter.py": 8,
    "gbdraw/web_support/feature_catalog.py": 16,
    "gbdraw/web_support/feature_metadata.py": 2,
    "gbdraw/web_support/orthogroup_metadata.py": 5,
    "gbdraw/web_support/request_render.py": 4,
    "gbdraw/web_support/rule_matching.py": 4,
    "gbdraw/web_support/similarity_alignment.py": 3,
}


def _pyproject() -> dict:
    with PYPROJECT.open("rb") as handle:
        return tomllib.load(handle)


def _checked_files() -> list[str]:
    package = PROJECT_ROOT / "gbdraw"
    return sorted(
        path.relative_to(PROJECT_ROOT).as_posix()
        for path in package.rglob("*.py")
        if not path.relative_to(package).as_posix().startswith("web/")
    )


def _tokens(path: str) -> list[tokenize.TokenInfo]:
    source = (PROJECT_ROOT / path).read_text(encoding="utf-8")
    return list(tokenize.generate_tokens(io.StringIO(source).readline))


def _installed_mypy_matches_the_pin() -> None:
    pins = _pyproject()["project"]["optional-dependencies"]["typecheck"]
    assert len(pins) == 1 and pins[0].startswith("mypy=="), pins
    pinned = pins[0].removeprefix("mypy==")
    try:
        installed = metadata.version("mypy")
    except metadata.PackageNotFoundError:
        installed = None
    assert installed == pinned, (
        f"mypy {pinned} is required, found {installed}: "
        "run `pip install -e '.[typecheck]'`"
    )


def _mypy_errors() -> dict[str, list[str]]:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy",
            "--config-file",
            str(PYPROJECT),
            "--cache-dir",
            os.devnull,
            "--output",
            "json",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    # Exit code 1 means "errors found", which always come with JSON lines.
    assert result.returncode == 0 or (
        result.returncode == 1 and result.stdout.strip()
    ), f"mypy failed with exit code {result.returncode}:\n{result.stderr}{result.stdout}"
    errors: dict[str, list[str]] = {}
    for line in result.stdout.splitlines():
        diagnostic = json.loads(line)
        if diagnostic["severity"] != "error":
            continue
        path = diagnostic["file"].replace("\\", "/")
        errors.setdefault(path, []).append(
            f"{path}:{diagnostic['line']}: {diagnostic['message']}"
            f"  [{diagnostic['code']}]"
        )
    return errors


def test_mypy_is_installed_at_the_pinned_version() -> None:
    _installed_mypy_matches_the_pin()


def test_mypy_configuration_is_the_registered_one() -> None:
    assert _pyproject()["tool"]["mypy"] == EXPECTED_MYPY_CONFIG, (
        "[tool.mypy] in pyproject.toml differs from EXPECTED_MYPY_CONFIG. "
        "A configuration change is a rule change: update both, and get the "
        "Owner's approval when it relaxes a check."
    )


def test_no_suppression_outside_the_count() -> None:
    """Reject suppressions that hide more than the one line they are counted on.

    A ``# mypy:`` comment changes the configuration of a file, a
    ``# type: ignore`` on a line of its own at the top of a module silences the
    whole module, and ``@no_type_check`` silences a whole function.
    """
    found = []
    for path in _checked_files():
        for token in _tokens(path):
            location = f"{path}:{token.start[0]}: {token.line.strip()}"
            if token.type == tokenize.COMMENT and (
                MYPY_INLINE_CONFIG.match(token.string)
                or (
                    TYPE_IGNORE.search(token.string)
                    and token.line.lstrip().startswith("#")
                )
            ):
                found.append(location)
            elif token.type == tokenize.NAME and token.string == "no_type_check":
                found.append(location)
    assert not found, (
        "`# mypy:` comments, `# type: ignore` on a line of its own, and "
        "`no_type_check` are not used:\n" + "\n".join(found)
    )


def test_type_debt_matches_the_baseline() -> None:
    _installed_mypy_matches_the_pin()
    files = _checked_files()
    errors = _mypy_errors()
    outside = sorted(set(errors) - set(files))
    assert not outside, f"mypy reported files outside the checked set: {outside}"
    ignores: dict[str, list[str]] = {}
    for path in files:
        for token in _tokens(path):
            if token.type == tokenize.COMMENT and TYPE_IGNORE.search(token.string):
                ignores.setdefault(path, []).append(
                    f"{path}:{token.start[0]}: {token.string}"
                )
    assert list(TYPE_DEBT_BASELINE) == sorted(TYPE_DEBT_BASELINE), (
        "keep TYPE_DEBT_BASELINE sorted by path"
    )

    checked = set(files)
    problems: list[str] = []
    for path in sorted(checked | set(TYPE_DEBT_BASELINE)):
        baseline = TYPE_DEBT_BASELINE.get(path, 0)
        if path not in checked:
            problems.append(f"{path} no longer exists: remove its entry")
            continue
        found = errors.get(path, []) + ignores.get(path, [])
        debt = len(found)
        detail = (
            f"{len(errors.get(path, []))} mypy errors + "
            f"{len(ignores.get(path, []))} `# type: ignore`"
        )
        if debt > baseline:
            problems.append(
                f"{path}: type debt {debt} ({detail}) is above its baseline "
                f"{baseline}. Fix the new errors; do not add `# type: ignore`. "
                "Raising the baseline needs the Owner's approval.\n    "
                + "\n    ".join(found)
            )
        elif debt < baseline:
            target = (
                "remove its entry" if debt == 0 else f"lower its entry to {debt}"
            )
            problems.append(
                f"{path}: type debt {debt} ({detail}) is below its baseline "
                f"{baseline}: {target} in TYPE_DEBT_BASELINE"
            )
    assert not problems, "\n".join(problems)
