from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from docs.recipes._scenario_support import (
    PUBLISHED_IMAGE_ROOT,
    extract_executable_block,
    load_chapter,
)


pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_MANIFEST = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
SCENARIO_IDS = ("T-CLI-01", "T-CLI-02", "T-PY-01")
RUNNER_BY_KIND = {
    "cli-recipe": "docs/recipes/run_cli_scenarios.py",
    "python-recipe": "docs/recipes/run_python_scenarios.py",
}


def _chapters() -> dict[str, dict[str, object]]:
    manifest = json.loads(SCENARIO_MANIFEST.read_text(encoding="utf-8"))
    return {
        chapter["id"]: chapter
        for chapter in manifest["scenarios"]
        if chapter["id"] in SCENARIO_IDS
    }


def test_onboarding_pages_match_their_approved_manifest_entries() -> None:
    chapters = _chapters()

    assert set(chapters) == set(SCENARIO_IDS)
    for scenario_id, chapter in chapters.items():
        destination = REPO_ROOT / str(chapter["destination"])
        source = destination.read_text(encoding="utf-8")
        execution = chapter["execution"]
        output_name = execution["expected_outputs"][0]

        assert source.startswith("[Home]")
        assert f"# {chapter['title']}" in source
        assert "## What you'll need" in source
        assert "## Step 2:" in source
        assert "## What you built" not in source
        assert output_name in _section(source, "Step 2")
        assert "tests/test_inputs" not in source
        assert "http://" not in source
        assert (PUBLISHED_IMAGE_ROOT / scenario_id.lower() / output_name).is_file()


def test_tutorials_omit_replay_verification_and_summary_sections() -> None:
    tutorial_root = REPO_ROOT / "docs" / "TUTORIALS"
    for path in sorted(tutorial_root.rglob("*.md")):
        source = path.read_text(encoding="utf-8")
        assert "## Step 3: Verify the replay" not in source, path
        assert "## What you built" not in source, path


def test_each_documented_recipe_has_one_marked_owner() -> None:
    public_sources = [
        path.read_text(encoding="utf-8")
        for path in (REPO_ROOT / "docs").rglob("*.md")
        if "internal" not in path.relative_to(REPO_ROOT / "docs").parts
    ]

    for scenario_id in SCENARIO_IDS:
        start = f"<!-- executable:{scenario_id}:start -->"
        end = f"<!-- executable:{scenario_id}:end -->"
        assert sum(source.count(start) for source in public_sources) == 1
        assert sum(source.count(end) for source in public_sources) == 1


def test_recipe_runners_extract_the_literal_documented_blocks() -> None:
    for scenario_id in ("T-CLI-01", "T-CLI-02"):
        chapter = load_chapter(
            scenario_id,
            expected_kind="cli-recipe",
            runner_path=RUNNER_BY_KIND["cli-recipe"],
        )
        recipe = extract_executable_block(chapter, language="bash")
        assert recipe.startswith(f"gbdraw {chapter['settings']['mode']} \\\n")
        assert chapter["execution"]["expected_outputs"][0].removesuffix(".svg") in recipe

    circular_chapter = load_chapter(
        "T-CLI-01",
        expected_kind="cli-recipe",
        runner_path=RUNNER_BY_KIND["cli-recipe"],
    )
    circular_recipe = extract_executable_block(circular_chapter, language="bash")
    assert "--qualifier_priority cds_gene_qualifier_priority.tsv" in circular_recipe
    assert "--labels out" in circular_recipe
    assert "--track_type middle" in circular_recipe
    assert '--species "<i>Homo sapiens</i>"' in circular_recipe

    python_chapter = load_chapter(
        "T-PY-01",
        expected_kind="python-recipe",
        runner_path=RUNNER_BY_KIND["python-recipe"],
    )
    python_recipe = extract_executable_block(python_chapter, language="python")
    compile(python_recipe, python_chapter["destination"], "exec")
    assert "diagram = draw_circular(record, options=options)" in python_recipe
    assert "saved_path = diagram.save(output_path)" in python_recipe


def test_onboarding_recipes_regenerate_from_an_external_clean_context(
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not existing_pythonpath
        else os.pathsep.join((str(REPO_ROOT), existing_pythonpath))
    )

    chapters = _chapters()
    for scenario_id in SCENARIO_IDS:
        execution = chapters[scenario_id]["execution"]
        assert isinstance(execution, dict)
        runner = RUNNER_BY_KIND[str(execution["kind"])]
        result = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / runner),
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
        assert "verified" in result.stdout

    assert list(tmp_path.iterdir()) == []


def _section(source: str, step: str) -> str:
    match = re.search(
        rf"^## {re.escape(step)}:[^\n]*\n(?P<body>.*?)(?=^## |\Z)",
        source,
        re.MULTILINE | re.DOTALL,
    )
    assert match is not None
    return match.group("body")
