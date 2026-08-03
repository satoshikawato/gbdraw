from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from docs.recipes._scenario_support import (
    PUBLISHED_IMAGE_ROOT,
    extract_executable_block,
    load_chapter,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = "docs/recipes/run_cli_scenarios.py"
SCENARIOS = {
    "H-CLI-02": (
        "use-input-tables.md",
        [
            "record_table.svg",
            "comparison_table.svg",
            "conservation_table.svg",
            "annotation_table.svg",
            "track_table.svg",
        ],
    ),
    "H-CLI-09": (
        "add-depth-gc-and-skew-tracks.md",
        ["cli_quantitative_tracks.svg"],
    ),
    "H-CLI-10": (
        "add-region-annotations-and-track-slots.md",
        ["cli_annotations_slots.svg"],
    ),
    "H-CLI-11": (
        "set-colors-labels-visibility-shapes-and-strokes.md",
        ["cli_feature_presentation.svg"],
    ),
    "H-CLI-12": (
        "save-and-regenerate-sessions.md",
        ["cli_session.json", "cli_session.json.gz", "cli_session_roundtrip.svg"],
    ),
    "H-CLI-13": (
        "export-static-and-interactive-outputs.md",
        [
            "cli_export.svg",
            "cli_export.interactive.svg",
            "cli_export.png",
            "cli_export.pdf",
            "cli_export.eps",
            "cli_export.ps",
        ],
    ),
}


def _fixture_manifest() -> dict[str, object]:
    return json.loads(
        (REPO_ROOT / "gbdraw/web/tutorial-data/manifest.json").read_text(
            encoding="utf-8"
        )
    )


def _recipe(scenario_id: str) -> tuple[dict[str, object], str, str]:
    chapter = load_chapter(
        scenario_id,
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    source = (REPO_ROOT / chapter["destination"]).read_text(encoding="utf-8")
    return chapter, source, extract_executable_block(chapter, language="bash")


def test_how_to_manifest_and_navigation_are_complete() -> None:
    index = (REPO_ROOT / "docs/HOW_TO/CLI/README.md").read_text(encoding="utf-8")

    for scenario_id, (filename, expected_outputs) in SCENARIOS.items():
        chapter, source, recipe = _recipe(scenario_id)
        assert chapter["execution"]["expected_outputs"] == expected_outputs
        assert chapter["status"] == {
            "implementation": "verified",
            "review": "approved",
        }
        assert "## Prerequisites" in source
        assert "## Verification" in source or "## Verify the export set" in source
        assert "## Troubleshooting" in source
        assert recipe
        assert f"]({filename})" in index


def test_input_tables_use_whole_sources_and_exact_table_kinds() -> None:
    _, source, recipe = _recipe("H-CLI-02")
    manifest = _fixture_manifest()
    losatn = manifest["files"]["lambda-de3-losatn"]["expectedSemantics"]
    mt_comparison = manifest["fixtures"]["metazoan-mitochondria-comparison"]

    assert recipe.count("gbdraw ") == 5
    assert "--records_table tables/records.tsv" in recipe
    assert "--comparisons_table tables/comparisons.tsv" in recipe
    assert "--conservation_table tables/conservation.tsv" in recipe
    assert "--annotation_table nicotiana-tabacum-regions.tsv" in recipe
    assert "--circular_track_table tables/tracks.tsv" in recipe
    assert "--gbk NC_001416.gb NC_042057.1.gb" in recipe
    assert "--gbk HmmtDNA.gbk" in recipe
    assert "--conservation_reference subject" in recipe
    assert "--qualifier_priority cds_gene_qualifier_priority.tsv" in recipe
    assert "--labels out" in recipe
    assert "--label_font_size 10" in recipe
    assert "four independent, naturally\ncircular RefSeq records" in source
    assert "none is cropped, split, or synthesized" in source
    assert (
        "Every CDS label in `record_table.svg` and `conservation_table.svg` "
        "comes from\nits `gene` qualifier"
    ) in source
    assert "No product description is used as a CDS label" in source
    assert manifest["fixtures"]["metazoan-mitochondria-four"][
        "expectedSemantics"
    ]["recordsAreWholeCanonicalSources"]
    assert losatn["queryRecordId"] == "NC_001416.1"
    assert losatn["subjectRecordId"] == "NC_042057.1"
    assert losatn["queryCoordinateRange"] == [1, 48_502]
    assert losatn["subjectCoordinateRange"] == [1, 42_925]
    assert losatn["rawRowCount"] == 6
    assert mt_comparison["expectedSemantics"]["sourceTopologies"] == ["circular"] * 4
    assert mt_comparison["expectedSemantics"]["referenceSide"] == "subject"
    assert mt_comparison["expectedSemantics"]["retainedRows"] == 106


def test_quantitative_and_annotation_recipes_pin_track_semantics() -> None:
    _, depth_source, depth_recipe = _recipe("H-CLI-09")
    _, annotation_source, annotation_recipe = _recipe("H-CLI-10")
    manifest = _fixture_manifest()

    assert "complete 606,194 bp circular `AP027133.1`" in depth_source
    assert "--depth_window 1" in depth_recipe
    assert "--depth_step 1000" in depth_recipe
    assert "--depth_min 0" in depth_recipe
    assert "--depth_max 80" in depth_recipe
    assert depth_recipe.count("--circular_track_slot") == 6
    assert manifest["fixtures"]["depth-1kb"]["expectedSemantics"][
        "binCount"
    ] == 607

    assert "complete 155,943 bp circular `NC_001879.2`" in annotation_source
    assert "LSC, IRb, SSC,\nand IRa" in annotation_source
    assert annotation_recipe.count("--circular_track_slot") == 5
    assert "set_id=plastome_regions" in annotation_recipe
    assert "--linear_track_slot" in annotation_source
    assert manifest["fixtures"]["tobacco-plastome-regions"][
        "expectedSemantics"
    ]["annotationIds"] == ["lsc", "irb", "ssc", "ira"]


def test_presentation_recipe_uses_one_complete_mitochondrion() -> None:
    _, source, recipe = _recipe("H-CLI-11")

    assert "complete circular 16,569 bp `NC_012920.1`" in source
    assert "--palette colorblind" in recipe
    assert "--table tables/presentation_colors.tsv" in recipe
    assert "--label_whitelist tables/presentation_labels.tsv" in recipe
    assert "--label_table tables/presentation_label_overrides.tsv" in recipe
    assert "--feature_visibility_table HmmtDNA_feature_visibility.tsv" in recipe
    assert recipe.count("--feature_shape") == 4
    assert "--resolve_overlaps" in recipe
    assert "--block_stroke_width 1.5" in recipe
    assert "--axis_stroke_width 4" in recipe
    assert "BGC000" not in source


def test_session_and_export_recipes_make_overwrite_and_formats_explicit() -> None:
    _, session_source, session_recipe = _recipe("H-CLI-12")
    _, export_source, export_recipe = _recipe("H-CLI-13")

    assert session_recipe.count("gbdraw circular") == 2
    assert "--session_output cli_session.json" in session_recipe
    assert "--session cli_session.json" in session_recipe
    assert "--session_output cli_session.json.gz" in session_recipe
    assert session_recipe.count("--overwrite") == 1
    assert "version 40" in session_source
    assert "schema 5" in session_source
    assert "byte-identical" in session_source

    assert export_recipe.count("gbdraw circular") == 1
    assert "-f svg,interactive_svg,png,pdf,eps,ps" in export_recipe
    for suffix in ("svg", "interactive.svg", "png", "pdf", "eps", "ps"):
        assert f"`cli_export.{suffix}`" in export_source


@pytest.mark.parametrize("scenario_id", SCENARIOS)
def test_recipe_regenerates_from_a_clean_external_directory(
    scenario_id: str,
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        value
        for value in (str(REPO_ROOT), environment.get("PYTHONPATH"))
        if value
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
        timeout=240,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert f"{scenario_id}: verified" in result.stdout
    assert list(tmp_path.iterdir()) == []
    output_dir = PUBLISHED_IMAGE_ROOT / scenario_id.lower()
    assert all((output_dir / name).is_file() for name in SCENARIOS[scenario_id][1])
