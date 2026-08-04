from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPO_ROOT / "docs/scenarios/manifest.json"
SCENARIO_ID_RE = re.compile(r"^(?:T|H|R|E|A)-[A-Z]+-\d{2}$")
ROLE_COUNTS = {
    "tutorial": 14,
    "how-to": 33,
    "reference": 10,
    "explanation": 6,
    "auxiliary": 2,
}
ALLOWED_SURFACES = {"gui", "cli", "python", "cross-surface", "gallery"}
ALLOWED_PRIORITIES = {"P0", "P1", "P2"}
ALLOWED_TIERS = {"core", "extended", "nightly"}
ALLOWED_EXECUTION_KINDS = {
    "playwright",
    "cli-recipe",
    "python-recipe",
    "contract",
}
SCREENSHOT_BUDGETS = {
    "T-GUI-01": 6,
    "T-GUI-02": 4,
    "T-GUI-03": 5,
    "T-GUI-04": 5,
    "T-GUI-05": 5,
    "T-GUI-06": 5,
    "T-GUI-08": 5,
    "T-GUI-09": 6,
    "H-GUI-01": 3,
    "H-GUI-02": 3,
    "H-GUI-03": 4,
    "H-GUI-04": 3,
    "H-GUI-05": 3,
    "H-GUI-06": 4,
    "H-GUI-07": 4,
    "H-GUI-08": 4,
    "H-GUI-09": 4,
    "H-GUI-10": 4,
    "H-GUI-11": 4,
    "H-GUI-12": 3,
    "H-GUI-13": 5,
    "H-GUI-14": 3,
    "H-GUI-15": 2,
}


def _manifest() -> dict[str, Any]:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def _repo_path(relative_path: str) -> Path:
    path = (REPO_ROOT / relative_path).resolve()
    assert path.is_relative_to(REPO_ROOT), f"path escapes repository: {relative_path}"
    return path


def test_chapter_plan_gate_is_approved_with_a_bounded_fixture_budget() -> None:
    manifest = _manifest()

    assert manifest["schema_version"] == 1
    assert manifest["approval"]["status"] == "approved"
    assert _repo_path(manifest["plan_source"]).is_file()

    budgets = manifest["package_budgets"]
    assert budgets == {
        "core_uncompressed_bytes": 1_048_576,
        "extended_incremental_uncompressed_bytes": 10_485_760,
        "gallery_nightly_bundled": False,
    }


def test_chapter_census_matches_the_reviewed_plan() -> None:
    chapters = _manifest()["chapters"]

    assert len(chapters) == 65
    assert Counter(chapter["role"] for chapter in chapters) == ROLE_COUNTS

    ids = [chapter["id"] for chapter in chapters]
    destinations = [chapter["destination"] for chapter in chapters]
    assert len(ids) == len(set(ids))
    assert len(destinations) == len(set(destinations))

    assert [chapter["id"] for chapter in chapters[:14]] == [
        "T-GUI-01",
        "T-GUI-02",
        "T-GUI-03",
        "T-GUI-04",
        "T-GUI-05",
        "T-GUI-06",
        "T-GUI-08",
        "T-GUI-09",
        "T-CLI-01",
        "T-CLI-02",
        "T-CLI-03",
        "T-CLI-04",
        "T-CLI-05",
        "T-PY-01",
    ]


def test_every_chapter_records_its_destination_inputs_and_executable_proof() -> None:
    manifest = _manifest()
    known_capabilities = {
        capability
        for capabilities in manifest["capability_groups"].values()
        for capability in capabilities
    }
    known_fixtures = set(manifest["fixtures"])

    for chapter in manifest["chapters"]:
        chapter_id = chapter["id"]
        assert SCENARIO_ID_RE.fullmatch(chapter_id), chapter_id
        assert chapter["title"].strip()
        assert chapter["surface"] in ALLOWED_SURFACES
        assert chapter["priority"] in ALLOWED_PRIORITIES
        assert chapter["capabilities"], chapter_id
        assert set(chapter["capabilities"]) <= known_capabilities
        assert set(chapter["canonical_for"]) <= set(chapter["capabilities"])
        assert set(chapter["fixtures"]) <= known_fixtures
        assert isinstance(chapter["settings"], dict)

        if chapter["role"] == "how-to":
            assert chapter["title"].startswith("How to "), chapter_id

        destination = _repo_path(chapter["destination"])
        assert "internal" not in destination.relative_to(REPO_ROOT).parts
        for source in chapter["sources"]:
            assert _repo_path(source).exists(), f"{chapter_id}: missing source {source}"

        execution = chapter["execution"]
        assert execution["kind"] in ALLOWED_EXECUTION_KINDS
        assert execution["tier"] in ALLOWED_TIERS
        assert execution["path"].strip()
        assert isinstance(execution["expected_outputs"], list)
        assert execution["assertions"], chapter_id

        status = chapter["status"]
        assert status["review"] == "approved"
        assert status["implementation"] in {"planned", "implemented", "verified"}
        if status["implementation"] != "planned":
            assert destination.is_file(), f"{chapter_id}: destination is not implemented"
            assert _repo_path(execution["path"]).is_file(), (
                f"{chapter_id}: executable proof is not implemented"
            )


def test_each_major_capability_has_exactly_one_canonical_owner() -> None:
    manifest = _manifest()
    groups = manifest["capability_groups"]
    expected_groups = {
        "entry-points",
        "inputs",
        "layout",
        "comparisons",
        "quantitative-and-annotation-tracks",
        "feature-presentation",
        "interactive-work",
        "outputs",
        "programmatic-use",
    }
    assert set(groups) == expected_groups

    vocabulary = [capability for capabilities in groups.values() for capability in capabilities]
    owners = [
        capability
        for chapter in manifest["chapters"]
        for capability in chapter["canonical_for"]
    ]
    assert len(vocabulary) == len(set(vocabulary))
    assert Counter(owners) == Counter({capability: 1 for capability in vocabulary})


def test_screenshots_are_owned_purposeful_and_within_the_reviewed_budget() -> None:
    chapters = _manifest()["chapters"]
    owned_paths: list[str] = []

    for chapter in chapters:
        chapter_id = chapter["id"]
        screenshots = chapter["screenshots"]
        if chapter_id in SCREENSHOT_BUDGETS:
            assert 0 < len(screenshots) <= SCREENSHOT_BUDGETS[chapter_id]
        else:
            assert screenshots == [], f"{chapter_id}: screenshots are not planned"

        for screenshot in screenshots:
            path = screenshot["path"]
            assert chapter["surface"] == "gui"
            assert path.startswith(f"docs/images/{chapter_id.lower()}/")
            assert path.endswith(".png")
            assert screenshot["reason"].strip()
            assert screenshot["alt"].strip()
            owned_paths.append(path)

            if chapter["status"]["implementation"] == "verified":
                assert _repo_path(path).is_file(), f"missing verified screenshot: {path}"

    assert len(owned_paths) == len(set(owned_paths))

    image_root = REPO_ROOT / "docs/images"
    actual_paths = {
        path.relative_to(REPO_ROOT).as_posix()
        for path in image_root.rglob("*.png")
    } if image_root.exists() else set()
    cli_png_outputs = {
        f"docs/images/{chapter['id'].lower()}/{output}"
        for chapter in chapters
        if chapter["surface"] == "cli"
        for output in chapter["execution"]["expected_outputs"]
        if output.endswith(".png")
    }
    owned_png_paths = set(owned_paths) | cli_png_outputs
    assert actual_paths <= owned_png_paths, (
        f"unowned documentation PNGs: {sorted(actual_paths - owned_png_paths)}"
    )


def test_tutorials_declare_an_early_visible_result_on_their_real_surface() -> None:
    tutorials = [
        chapter for chapter in _manifest()["chapters"] if chapter["role"] == "tutorial"
    ]
    expected_kind = {
        "gui": "playwright",
        "cli": "cli-recipe",
        "python": "python-recipe",
    }

    for chapter in tutorials:
        assertion = next(
            value
            for value in chapter["execution"]["assertions"]
            if value.startswith("first_result_step=")
        )
        step = int(assertion.rsplit("=", 1)[1])
        assert step <= 3
        if chapter["id"] in {
            "T-GUI-05",
            "T-GUI-06",
            "T-GUI-08",
            "T-GUI-09",
            "T-CLI-03",
            "T-CLI-04",
            "T-CLI-05",
        }:
            assert step == 2
        assert chapter["execution"]["kind"] == expected_kind[chapter["surface"]]


def test_first_gui_tutorial_keeps_the_fixed_accepted_path() -> None:
    chapter = next(
        chapter
        for chapter in _manifest()["chapters"]
        if chapter["id"] == "T-GUI-01"
    )

    assert chapter["fixtures"] == ["human-mitochondrion"]
    assert chapter["settings"] == {
        "mode": "circular",
        "output_prefix": "human_mitochondrion",
        "species": "<i>Homo sapiens</i>",
        "track_preset": "middle",
        "separate_strands": True,
        "hide_gc_content": False,
        "hide_gc_skew": False,
        "label_mode": "out",
        "cds_label_qualifier": "gene",
        "legend_position": "right",
    }
    assert chapter["execution"]["expected_outputs"] == ["human_mitochondrion.svg"]
    assert len(chapter["screenshots"]) == 6


def test_public_scenario_contracts_do_not_depend_on_internal_or_test_inputs() -> None:
    manifest_text = MANIFEST_PATH.read_text(encoding="utf-8")

    assert "tests/test_inputs" not in manifest_text
    for chapter in _manifest()["chapters"]:
        assert not chapter["destination"].startswith("docs/internal/")
        assert not chapter["execution"]["path"].startswith("docs/internal/")


def test_explanation_chapters_are_present_and_keep_surface_boundaries_clear() -> None:
    explanations = [
        chapter
        for chapter in _manifest()["chapters"]
        if chapter["role"] == "explanation"
    ]

    assert len(explanations) == 6
    for chapter in explanations:
        path = _repo_path(chapter["destination"])
        assert path.is_file(), chapter["id"]
        source = path.read_text(encoding="utf-8")
        assert source.startswith("[Documentation home]")
        assert f"# {chapter['title']}" in source
        assert "Web UX profile" not in source
        assert "tests/test_inputs" not in source
        assert "lambda_two_contigs" not in source

    comparison = _repo_path(
        next(
            chapter["destination"]
            for chapter in explanations
            if chapter["id"] == "E-COMPARISON-01"
        )
    ).read_text(encoding="utf-8")
    assert "LOSATN and TLOSATX run in the web app" in comparison
    assert "not phylogenetic orthogroups" in comparison


def test_explanation_pages_do_not_duplicate_runnable_recipe_authority() -> None:
    for chapter in _manifest()["chapters"]:
        if chapter["role"] != "explanation":
            continue
        source = _repo_path(chapter["destination"]).read_text(encoding="utf-8")
        assert "```bash" not in source
        assert "```python" not in source
