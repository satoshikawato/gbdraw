from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any

from docs.recipes._scenario_support import parse_translate_chain


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPO_ROOT / "docs/scenarios/manifest.json"
SCENARIO_ID_RE = re.compile(r"^(?:T|H|R|A)-[A-Z]+-\d{2}$")
ALLOWED_ROLES = {"tutorial", "evidence", "reference", "auxiliary"}
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
    "T-GUI-04": 6,
    "T-GUI-05": 5,
    "T-GUI-06": 5,
    "T-GUI-08": 6,
    "T-GUI-09": 6,
    "T-GUI-10": 2,
    "T-GUI-12": 2,
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


def _svg_tree(path: Path) -> bytes:
    return ET.tostring(ET.parse(path).getroot())


def test_recipe_transform_parser_accepts_composed_layout_translations() -> None:
    assert parse_translate_chain("translate(16.0,-31.5)") == (16.0, -31.5)
    assert parse_translate_chain(
        "translate(16.0,-31.5) translate(75.0,61.0)"
    ) == (91.0, 29.5)
    assert parse_translate_chain("translate(4 5) translate(-1)") == (3.0, 5.0)
    assert parse_translate_chain("") is None
    assert parse_translate_chain("translate(1,2) scale(2)") is None


def test_scenario_plan_gate_is_approved_with_a_bounded_fixture_budget() -> None:
    manifest = _manifest()

    assert manifest["schema_version"] == 2
    assert manifest["approval"]["status"] == "approved"
    assert _repo_path(manifest["plan_source"]).is_file()

    budgets = manifest["package_budgets"]
    assert budgets == {
        "core_uncompressed_bytes": 1_048_576,
        "extended_incremental_uncompressed_bytes": 10_485_760,
        "gallery_nightly_bundled": False,
    }


def test_scenario_registry_has_unique_ids_and_public_destinations() -> None:
    scenarios = _manifest()["scenarios"]

    assert set(Counter(scenario["role"] for scenario in scenarios)) <= ALLOWED_ROLES
    ids = [scenario["id"] for scenario in scenarios]
    destinations = [
        scenario["destination"]
        for scenario in scenarios
        if "destination" in scenario
    ]
    assert len(ids) == len(set(ids))
    assert len(destinations) == len(set(destinations))


def test_every_scenario_records_inputs_and_executable_proof() -> None:
    manifest = _manifest()
    known_capabilities = {
        capability
        for capabilities in manifest["capability_groups"].values()
        for capability in capabilities
    }
    known_fixtures = set(manifest["fixtures"])

    for chapter in manifest["scenarios"]:
        chapter_id = chapter["id"]
        assert SCENARIO_ID_RE.fullmatch(chapter_id), chapter_id
        assert chapter["title"].strip()
        assert chapter["role"] in ALLOWED_ROLES
        assert chapter["surface"] in ALLOWED_SURFACES
        assert chapter["priority"] in ALLOWED_PRIORITIES
        assert chapter["capabilities"], chapter_id
        assert set(chapter["capabilities"]) <= known_capabilities
        assert ("canonical" + "_for") not in chapter
        assert set(chapter["fixtures"]) <= known_fixtures
        assert isinstance(chapter["settings"], dict)

        for source in chapter["sources"]:
            assert _repo_path(source).exists(), f"{chapter_id}: missing source {source}"

        execution = chapter["execution"]
        assert execution["kind"] in ALLOWED_EXECUTION_KINDS
        assert execution["tier"] in ALLOWED_TIERS
        assert execution["path"].strip()
        assert isinstance(execution["expected_outputs"], list)
        assert execution["assertions"], chapter_id

        destination = None
        if chapter["role"] == "evidence":
            assert "destination" not in chapter
            if execution["kind"] in {"cli-recipe", "python-recipe"}:
                assert execution["source"] == "docs/internal/SCENARIO_EVIDENCE.md"
                assert _repo_path(execution["source"]).is_file()
        else:
            destination = _repo_path(chapter["destination"])
            assert "internal" not in destination.relative_to(REPO_ROOT).parts

        status = chapter["status"]
        assert status["review"] == "approved"
        assert status["implementation"] in {"planned", "implemented", "verified"}
        if status["implementation"] != "planned":
            if destination is not None:
                assert destination.is_file(), (
                    f"{chapter_id}: destination is not implemented"
                )
            assert _repo_path(execution["path"]).is_file(), (
                f"{chapter_id}: executable proof is not implemented"
            )


def _markdown_heading_anchors(source: str) -> set[str]:
    anchors: set[str] = set()
    counts: Counter[str] = Counter()
    for line in source.splitlines():
        match = re.match(r"^#{1,6}\s+(.+?)\s*$", line)
        if match is None:
            continue
        heading = re.sub(r"[<>`*_]", "", match.group(1)).lower()
        anchor = re.sub(r"[^a-z0-9 _-]", "", heading)
        anchor = re.sub(r"[ -]+", "-", anchor).strip("-")
        suffix = counts[anchor]
        counts[anchor] += 1
        anchors.add(f"{anchor}-{suffix}" if suffix else anchor)
    return anchors


def test_each_major_capability_has_one_public_owner_and_execution_evidence() -> None:
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
    owners = [capability for values in manifest["public_owners"].values() for capability in values]
    exercised = {
        capability
        for scenario in manifest["scenarios"]
        for capability in scenario["capabilities"]
    }
    assert len(vocabulary) == len(set(vocabulary))
    assert Counter(owners) == Counter({capability: 1 for capability in vocabulary})
    assert exercised == set(vocabulary)

    technical_index = (REPO_ROOT / "docs/REFERENCE/README.md").read_text(
        encoding="utf-8"
    )
    for target in manifest["public_owners"]:
        relative_path, _, anchor = target.partition("#")
        path = _repo_path(relative_path)
        assert path.is_file(), target
        assert "internal" not in path.relative_to(REPO_ROOT).parts
        assert f"({path.name})" in technical_index, target
        if anchor:
            source = path.read_text(encoding="utf-8")
            assert anchor in _markdown_heading_anchors(source), target


def test_screenshots_are_owned_purposeful_and_within_the_reviewed_budget() -> None:
    chapters = _manifest()["scenarios"]
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
        chapter for chapter in _manifest()["scenarios"] if chapter["role"] == "tutorial"
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
            "T-CLI-05",
            "T-CLI-06",
            "T-PY-02",
        }:
            assert step == 2
        assert chapter["execution"]["kind"] == expected_kind[chapter["surface"]]


def test_tutorial_projects_define_one_figure_with_intentional_surface_variants() -> None:
    manifest = _manifest()
    policy = manifest["tutorial_project_policy"]
    projects = manifest["tutorial_projects"]
    tutorials = {
        chapter["id"]: chapter
        for chapter in manifest["scenarios"]
        if chapter["role"] == "tutorial"
    }

    allowed_surfaces = set(policy["allowed_surfaces"])
    assert allowed_surfaces == {"gui", "cli", "python"}
    assert "materially different reader journey" in policy["rule"]
    assert "does not require a separate public page" in policy["evidence"]

    owned_variants: dict[str, str] = {}
    for project_id, project in projects.items():
        assert project["figure"].strip(), project_id
        assert project["parity"] == "verified"
        variants = project["variants"]
        assert variants
        assert set(variants) <= allowed_surfaces
        implemented = {
            surface: scenario_id
            for surface, scenario_id in variants.items()
            if scenario_id is not None
        }
        assert implemented, project_id
        for surface, scenario_id in implemented.items():
            chapter = tutorials[scenario_id]
            assert chapter["surface"] == surface
            assert chapter["project_id"] == project_id
            assert scenario_id not in owned_variants
            owned_variants[scenario_id] = project_id
            source = _repo_path(chapter["destination"]).read_text(encoding="utf-8")
            assert "## Choose how to build this figure" in source

    assert set(owned_variants) == set(tutorials)

    chloroplast = projects["gallery-chloroplast-map"]
    assert chloroplast["parity"] == "verified"
    shared_keys = {
        "mode",
        "record_id",
        "feature_types",
        "annotation_set",
        "track_order",
        "labels",
        "label_offsets",
        "legend",
    }
    shared_settings = [
        {
            key: tutorials[scenario_id]["settings"][key]
            for key in shared_keys
        }
        for scenario_id in chloroplast["variants"].values()
    ]
    assert shared_settings[0] == shared_settings[1] == shared_settings[2]


def test_verified_tutorial_non_browser_renderers_publish_the_same_svg_tree() -> None:
    first_circular = [
        REPO_ROOT / "docs/images/t-cli-01/human_mitochondrion.svg",
        REPO_ROOT / "docs/images/t-py-01/python_human_mitochondrion.svg",
    ]
    chloroplast = [
        REPO_ROOT / "gbdraw/web/gallery/sources/tobacco-chloroplast.svg",
        REPO_ROOT / "docs/images/t-cli-06/cli_annotated_chloroplast.svg",
        REPO_ROOT / "docs/images/t-py-02/python_annotated_chloroplast.svg",
    ]

    for project in (first_circular, chloroplast):
        trees = [_svg_tree(path) for path in project]
        assert all(tree == trees[0] for tree in trees[1:])


def test_first_gui_tutorial_keeps_the_fixed_accepted_path() -> None:
    chapter = next(
        chapter
        for chapter in _manifest()["scenarios"]
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


def test_public_destinations_and_execution_paths_keep_their_boundaries() -> None:
    manifest_text = MANIFEST_PATH.read_text(encoding="utf-8")

    assert "tests/test_inputs" not in manifest_text
    for chapter in _manifest()["scenarios"]:
        if "destination" in chapter:
            assert not chapter["destination"].startswith("docs/internal/")
        assert not chapter["execution"]["path"].startswith("docs/internal/")
        source = chapter["execution"].get("source")
        if source is not None:
            assert source == "docs/internal/SCENARIO_EVIDENCE.md"


def test_retired_public_categories_are_not_scenarios() -> None:
    roles = {scenario["role"] for scenario in _manifest()["scenarios"]}
    assert ("how" + "-to") not in roles
    assert ("expla" + "nation") not in roles
