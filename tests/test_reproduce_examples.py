from __future__ import annotations

import json
from functools import cache
import re
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree

import pytest

from tools.reproduce_examples import PROJECT_ROOT, Reproducer
from tools.reproduce_examples_manifest import (
    MANUALLY_MANAGED_FIGURES,
    UNREFERENCED_FIGURE_RETENTION,
    CliRecipe,
    CompositeRecipe,
    FigureSpec,
    build_figure_specs,
    load_palette_names,
)
from tools.generate_palette_page import OUTPUT_PATH as PALETTE_PAGE, render_palette_page
from tools.generate_palette_explorer_assets import (
    JSON_OUTPUT_PATH as PALETTE_EXPLORER_JSON,
    SVG_OUTPUT_PATH as PALETTE_EXPLORER_SVG,
    palette_payload,
    validate_svg,
)


pytestmark = pytest.mark.gallery


PUBLIC_MARKDOWN = (
    PROJECT_ROOT / "README.md",
    *sorted(
        path
        for path in (PROJECT_ROOT / "docs").rglob("*.md")
        if "internal" not in path.relative_to(PROJECT_ROOT / "docs").parts
    ),
    PROJECT_ROOT / "examples" / "color_palette_examples.md",
)
MARKDOWN_TARGET_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
HTML_IMAGE_RE = re.compile(r'<img\b[^>]*\bsrc=["\']([^"\']+)', re.IGNORECASE)
SELF_GITHUB_PREFIX = "https://github.com/satoshikawato/gbdraw/blob/main/"


@cache
def _figure_specs() -> dict[str, FigureSpec]:
    return build_figure_specs()


def _local_target(markdown_path: Path, raw_target: str) -> Path | None:
    target = raw_target.strip().split()[0].strip("<>")
    if target.startswith(SELF_GITHUB_PREFIX):
        return PROJECT_ROOT / target.removeprefix(SELF_GITHUB_PREFIX).split("#", 1)[0]
    if target.startswith("#") or "://" in target or target.startswith(("mailto:", "data:")):
        return None
    return (markdown_path.parent / target.split("#", 1)[0].split("?", 1)[0]).resolve()


@cache
def _markdown_targets(markdown_path: Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    source = markdown_path.read_text(encoding="utf-8")
    image_targets = (
        *re.findall(r"!\[[^\]]*\]\(([^)]+)\)", source),
        *HTML_IMAGE_RE.findall(source),
    )
    return (
        (*MARKDOWN_TARGET_RE.findall(source), *HTML_IMAGE_RE.findall(source)),
        image_targets,
    )


def _local_image_references() -> set[str]:
    references: set[str] = set()
    for markdown_path in PUBLIC_MARKDOWN:
        _all_targets, image_targets = _markdown_targets(markdown_path)
        for raw_target in image_targets:
            target = _local_target(markdown_path, raw_target)
            if target is not None:
                references.add(target.relative_to(PROJECT_ROOT).as_posix())
    return references


def _documentation_scenario_artifacts() -> set[str]:
    manifest = json.loads(
        (PROJECT_ROOT / "docs" / "scenarios" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    artifacts: set[str] = set()
    for chapter in manifest["scenarios"]:
        artifacts.update(item["path"] for item in chapter["screenshots"])
        if chapter["execution"]["kind"] not in {"cli-recipe", "python-recipe"}:
            continue
        for output_name in chapter["execution"]["expected_outputs"]:
            if Path(output_name).suffix.lower() in {".png", ".svg"}:
                artifacts.add(f"docs/images/{chapter['id'].lower()}/{output_name}")
    return artifacts


def test_manifest_counts_and_unique_paths() -> None:
    figures = _figure_specs()

    docs_and_readme = [spec for spec in figures.values() if "palettes" not in spec.groups]
    palette_circular = [figure_id for figure_id in figures if figure_id.startswith("palette_circular_")]
    palette_linear = [figure_id for figure_id in figures if figure_id.startswith("palette_linear_")]

    assert len(docs_and_readme) == 53
    assert palette_circular == [
        "palette_circular_default",
        "palette_circular_ajisai",
        "palette_circular_soft_pastels",
    ]
    assert palette_linear == [
        "palette_linear_default",
        "palette_linear_ajisai",
        "palette_linear_soft_pastels",
    ]
    assert len(figures) == 53 + 6
    assert "gbdraw_social_preview" not in figures

    output_paths = [spec.output_path for spec in figures.values()]
    assert len(output_paths) == len(set(output_paths))


def test_showcase_comparisons_keep_worked_example_context() -> None:
    figures = _figure_specs()

    track_recipe = figures["track_layout_separate_strands"].recipe
    assert isinstance(track_recipe, CompositeRecipe)
    assert len(track_recipe.panels) == 6
    assert [
        panel.recipe.extra_args[panel.recipe.extra_args.index("--track_type") + 1]
        for panel in track_recipe.panels
        if panel.recipe
    ] == ["tuckin", "middle", "spreadout"] * 2

    label_recipe = figures["label_font_size_comparison"].recipe
    assert isinstance(label_recipe, CompositeRecipe)
    assert figures["label_font_size_comparison"].required_inputs == ("HmmtDNA.gbk",)
    assert label_recipe.tile_size == (3400, 2200)
    assert [
        panel.recipe.extra_args[panel.recipe.extra_args.index("--label_font_size") + 1]
        for panel in label_recipe.panels
        if panel.recipe
    ] == ["12", "18", "24"]

    definition_recipe = figures["definition_font_size_comparison"].recipe
    assert isinstance(definition_recipe, CompositeRecipe)
    assert [
        panel.recipe.extra_args[panel.recipe.extra_args.index("--plot_title") + 1]
        for panel in definition_recipe.panels
        if panel.recipe
    ] == [
        "--definition_font_size 20",
        "--definition_font_size 28",
        "--definition_font_size 36",
    ]

    offset_recipe = figures["outer_label_offset_comparison"].recipe
    assert isinstance(offset_recipe, CompositeRecipe)
    assert len(offset_recipe.panels) == 9
    assert all(
        "--plot_title" in panel.recipe.extra_args
        for panel in offset_recipe.panels
        if panel.recipe
    )

    skew_recipe = figures["skew_comparison"].recipe
    assert isinstance(skew_recipe, CompositeRecipe)
    assert len(skew_recipe.panels) == 12
    assert skew_recipe.columns == 4
    assert [
        panel.recipe.extra_args[panel.recipe.extra_args.index("--nt") + 1]
        for panel in skew_recipe.panels
        if panel.recipe
    ] == [
        "GC",
        "CG",
        "AG",
        "GA",
        "CT",
        "TC",
        "TG",
        "GT",
        "CA",
        "AC",
        "AT",
        "TA",
    ]

    window_recipe = figures["window_step_comparison"].recipe
    assert isinstance(window_recipe, CompositeRecipe)
    assert [
        (
            panel.recipe.extra_args[panel.recipe.extra_args.index("--window") + 1],
            "--step",
            panel.recipe.extra_args[panel.recipe.extra_args.index("--step") + 1],
        )
        for panel in window_recipe.panels
        if panel.recipe
    ] == [
        ("100000", "--step", "10000"),
        ("10000", "--step", "1000"),
        ("1000", "--step", "100"),
    ]

    for recipe in (track_recipe, label_recipe, offset_recipe, window_recipe, skew_recipe):
        for panel in recipe.panels:
            assert panel.recipe is not None
            assert "--no-gc" not in panel.recipe.extra_args
            assert "--no-skew" not in panel.recipe.extra_args
            assert "none" not in panel.recipe.extra_args


def test_non_layout_linear_showcases_use_on_axis_features_and_line_styling() -> None:
    figures = _figure_specs()
    on_axis_figure_ids = (
        "tutorial_2_pairwise_blast",
        "linear_multi_record",
        "tutorial_protein_pairwise",
        "tutorial_5_records_table",
        "tutorial_7_linear_layout",
        "tutorial_7_definition_lines",
    )

    for figure_id in on_axis_figure_ids:
        recipe = figures[figure_id].recipe
        assert isinstance(recipe, CliRecipe)
        assert "--track_layout" not in recipe.extra_args
        assert recipe.extra_args[recipe.extra_args.index("--block_stroke_color") + 1] == "gray"
        assert recipe.extra_args[recipe.extra_args.index("--block_stroke_width") + 1] == "1"
        assert recipe.extra_args[recipe.extra_args.index("--line_stroke_color") + 1] == "lightgray"
        assert recipe.extra_args[recipe.extra_args.index("--line_stroke_width") + 1] == "2"

    below_recipe = figures["tutorial_7_track_layout_below"].recipe
    assert isinstance(below_recipe, CliRecipe)
    assert below_recipe.extra_args[
        below_recipe.extra_args.index("--track_layout") + 1
    ] == "below"


def test_palette_manifest_stays_in_sync_with_palette_file() -> None:
    figures = _figure_specs()
    palette_names = load_palette_names()

    for palette_name in ("default", "ajisai", "soft_pastels"):
        assert f"palette_circular_{palette_name}" in figures
        assert f"palette_linear_{palette_name}" in figures
    assert len(palette_names) == 55
    assert "palettes_combined_image_1" not in figures
    assert "palettes_combined_image_2" not in figures


def test_palette_page_is_generated_from_palette_file() -> None:
    rendered = render_palette_page()
    assert PALETTE_PAGE.read_text(encoding="utf-8") == rendered
    assert '<span style="color:' not in rendered
    assert r'$\textcolor{#54bcf8}{\blacksquare}$' in rendered


def test_palette_explorer_uses_one_semantic_circular_svg() -> None:
    payload = json.loads(PALETTE_EXPLORER_JSON.read_text(encoding="utf-8"))
    assert payload == palette_payload()
    assert len(payload["palettes"]) == 55
    counts = validate_svg(PALETTE_EXPLORER_SVG)
    assert counts["CDS"] > 0
    assert counts["skew_high"] > 0
    assert counts["skew_low"] > 0
    assert counts["gc_content"] > 0


def test_public_markdown_local_targets_exist() -> None:
    missing: list[str] = []
    for markdown_path in PUBLIC_MARKDOWN:
        targets, _image_targets = _markdown_targets(markdown_path)
        for raw_target in targets:
            target = _local_target(markdown_path, raw_target)
            if target is not None and not target.exists():
                missing.append(f"{markdown_path.relative_to(PROJECT_ROOT)} -> {raw_target}")
    assert missing == []


def test_public_figures_have_reproduction_inventory_coverage() -> None:
    references = _local_image_references()
    manifest_paths = {spec.output_path for spec in _figure_specs().values()}
    scenario_paths = _documentation_scenario_artifacts()
    manual_paths = set(MANUALLY_MANAGED_FIGURES)
    retained_unreferenced = set(UNREFERENCED_FIGURE_RETENTION)

    assert references - manifest_paths - scenario_paths - manual_paths == set()
    assert retained_unreferenced <= manifest_paths - references
    assert all(reason.strip() for reason in MANUALLY_MANAGED_FIGURES.values())
    assert all(reason.strip() for reason in UNREFERENCED_FIGURE_RETENTION.values())
    assert all((PROJECT_ROOT / path).exists() for path in manual_paths)
    assert "examples/gbdraw_social_preview.png" in manual_paths


def test_alias_resolution_and_support_asset_materialization(tmp_path: Path) -> None:
    reproducer = Reproducer(
        project_root=PROJECT_ROOT,
        output_root=tmp_path / "out",
        figures=_figure_specs(),
    )
    try:
        aliased = reproducer.resolve_input("NC_000913.gbk", "ecoli_k12_plot", {}, dry_run=False)
        support = reproducer.resolve_input(
            "feature_specific_colors.tsv",
            "MjeNMV_feature_specifc_colors_with_labels",
            {},
            dry_run=False,
        )
        assert aliased is not None
        assert aliased.name == "MG1655.gbk"
        assert support is not None
        assert support.exists()
        assert "wsv.*-like protein" in support.read_text()

        label_override = reproducer.resolve_input(
            "label_override.tsv",
            "tutorial_3_label_override",
            {},
            dry_run=False,
        )
        assert label_override is not None
        assert label_override.parent == tmp_path / "out" / "_support_assets"
        assert "LC738868.1\tCDS\tlabel" in label_override.read_text()

        payload = reproducer.report_payload()
        assert payload["aliases_used"]
        assert payload["aliases_used"][0]["requested"] == "NC_000913.gbk"
    finally:
        reproducer.close()


def test_label_override_showcase_uses_its_manifest_owned_table(tmp_path: Path) -> None:
    reproducer = Reproducer(
        project_root=PROJECT_ROOT,
        output_root=tmp_path / "out",
        figures=_figure_specs(),
    )
    try:
        assert reproducer.render_figure("tutorial_3_label_override") is True
        svg = reproducer.output_path_for("tutorial_3_label_override").read_text()
        assert ">gustavus-like protein<" in svg
        assert ">protein gustavus-like protein<" not in svg
    finally:
        reproducer.close()


def test_missing_report_structure_for_missing_figure(tmp_path: Path) -> None:
    figure_id = "missing_figure"
    missing_filename = "definitely_missing.gb"
    figures = dict(_figure_specs())
    figures[figure_id] = FigureSpec(
        figure_id=figure_id,
        output_path="examples/missing_figure.svg",
        groups=("docs",),
        required_inputs=(missing_filename,),
        recipe=CliRecipe(subcommand="linear", gbk_files=(missing_filename,)),
        description="Synthetic missing-input figure for report validation.",
    )
    reproducer = Reproducer(
        project_root=PROJECT_ROOT,
        output_root=tmp_path / "out",
        figures=figures,
    )
    try:
        assert reproducer.render_figure(figure_id) is False
        report_path = tmp_path / "missing_inputs.json"
        reproducer.write_report(report_path)

        payload = json.loads(report_path.read_text())
        assert set(payload) == {
            "generated",
            "skipped_missing_inputs",
            "failed",
            "aliases_used",
            "missing_inputs",
        }
        assert payload["generated"] == []
        assert payload["skipped_missing_inputs"] == [figure_id]
        missing_entry = next(entry for entry in payload["missing_inputs"] if entry["filename"] == missing_filename)
        assert missing_entry["figures"] == [figure_id]
        assert missing_entry["could_derive_if_base_inputs_present"] is False
    finally:
        reproducer.close()


def test_smoke_render_subset_via_cli(tmp_path: Path) -> None:
    output_root = tmp_path / "reproduced"
    report_path = tmp_path / "report.json"
    cmd = [
        sys.executable,
        str(PROJECT_ROOT / "tools" / "reproduce_examples.py"),
        "--output-root",
        str(output_root),
        "--missing-report",
        str(report_path),
        "--figure",
        "ecoli_k12_plot",
        "--figure",
        "majani",
        "--figure",
        "palette_circular_default",
        "--figure",
        "palette_linear_default",
        "--figure",
        "window_step_comparison",
    ]
    result = subprocess.run(
        cmd,
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        timeout=3600,
    )

    assert result.returncode == 0, result.stderr or result.stdout
    assert (output_root / "examples" / "ecoli_k12_plot.svg").exists()
    assert (output_root / "examples" / "majani.svg").exists()
    assert (output_root / "examples" / "AP027078_tuckin_separate_strands_default.svg").exists()
    assert (output_root / "examples" / "hepatoplasmataceae_default.svg").exists()
    assert (output_root / "examples" / "window_step_comparison.png").exists()

    payload = json.loads(report_path.read_text())
    assert payload["failed"] == []
    assert payload["skipped_missing_inputs"] == []
    assert "ecoli_k12_plot" in payload["generated"]
    assert "majani" in payload["generated"]


@pytest.fixture(scope="module")
def reproduced_arrow_geometry_variants(
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[str, Path]:
    figure_ids = (
        "tutorial_9_arrow_geometry_circular",
        "tutorial_9_arrow_geometry_linear",
    )
    reproducer = Reproducer(
        project_root=PROJECT_ROOT,
        output_root=tmp_path_factory.mktemp("arrow-geometry") / "out",
        figures=_figure_specs(),
    )
    try:
        for figure_id in figure_ids:
            assert reproducer.render_figure(figure_id) is True
        return {
            figure_id: reproducer.output_path_for(figure_id)
            for figure_id in figure_ids
        }
    finally:
        reproducer.close()


def test_gallery_session_arrow_geometry_variants_reproduce_tracked_svgs(
    reproduced_arrow_geometry_variants: dict[str, Path],
) -> None:
    for figure_id, generated in reproduced_arrow_geometry_variants.items():
        tracked = PROJECT_ROOT / _figure_specs()[figure_id].output_path
        assert generated.read_text(encoding="utf-8") == tracked.read_text(
            encoding="utf-8"
        )


def test_gallery_session_arrow_geometry_variants_only_change_feature_paths(
    reproduced_arrow_geometry_variants: dict[str, Path],
) -> None:
    pairs = (
        (
            PROJECT_ROOT / "gbdraw/web/gallery/sources/HmmtDNA_ATskew.svg",
            reproduced_arrow_geometry_variants[
                "tutorial_9_arrow_geometry_circular"
            ],
            15,
        ),
        (
            PROJECT_ROOT
            / "gbdraw/web/gallery/sources/BGC0000708-BGC0000713.svg",
            reproduced_arrow_geometry_variants[
                "tutorial_9_arrow_geometry_linear"
            ],
            152,
        ),
    )
    for source_path, variant_path, expected_changed_paths in pairs:
        source_elements = list(ElementTree.parse(source_path).getroot().iter())
        variant_elements = list(ElementTree.parse(variant_path).getroot().iter())
        assert len(source_elements) == len(variant_elements)

        changed_paths = 0
        for source, variant in zip(source_elements, variant_elements, strict=True):
            assert source.tag == variant.tag
            assert source.text == variant.text
            assert source.tail == variant.tail
            source_attributes = dict(source.attrib)
            variant_attributes = dict(variant.attrib)
            source_path_data = source_attributes.pop("d", None)
            variant_path_data = variant_attributes.pop("d", None)
            assert source_attributes == variant_attributes
            if source_path_data != variant_path_data:
                assert source.tag.endswith("path")
                changed_paths += 1

        assert changed_paths == expected_changed_paths
