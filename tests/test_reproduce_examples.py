from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

from tools.reproduce_examples import PROJECT_ROOT, Reproducer
from tools.reproduce_examples_manifest import (
    MANUALLY_MANAGED_FIGURES,
    UNREFERENCED_FIGURE_RETENTION,
    CliRecipe,
    FigureSpec,
    build_figure_specs,
    load_palette_names,
)
from tools.generate_palette_page import OUTPUT_PATH as PALETTE_PAGE, render_palette_page


PUBLIC_MARKDOWN = (
    PROJECT_ROOT / "README.md",
    *sorted((PROJECT_ROOT / "docs").rglob("*.md")),
    PROJECT_ROOT / "examples" / "color_palette_examples.md",
)
MARKDOWN_TARGET_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
HTML_IMAGE_RE = re.compile(r'<img\b[^>]*\bsrc=["\']([^"\']+)', re.IGNORECASE)
SELF_GITHUB_PREFIX = "https://github.com/satoshikawato/gbdraw/blob/main/"


def _local_target(markdown_path: Path, raw_target: str) -> Path | None:
    target = raw_target.strip().split()[0].strip("<>")
    if target.startswith(SELF_GITHUB_PREFIX):
        return PROJECT_ROOT / target.removeprefix(SELF_GITHUB_PREFIX).split("#", 1)[0]
    if target.startswith("#") or "://" in target or target.startswith(("mailto:", "data:")):
        return None
    return (markdown_path.parent / target.split("#", 1)[0].split("?", 1)[0]).resolve()


def _local_image_references() -> set[str]:
    references: set[str] = set()
    for markdown_path in PUBLIC_MARKDOWN:
        source = markdown_path.read_text(encoding="utf-8")
        targets = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", source)
        targets.extend(HTML_IMAGE_RE.findall(source))
        for raw_target in targets:
            target = _local_target(markdown_path, raw_target)
            if target is not None:
                references.add(target.relative_to(PROJECT_ROOT).as_posix())
    return references


def test_manifest_counts_and_unique_paths() -> None:
    figures = build_figure_specs()

    docs_and_readme = [spec for spec in figures.values() if "palettes" not in spec.groups]
    palette_circular = [figure_id for figure_id in figures if figure_id.startswith("palette_circular_")]
    palette_linear = [figure_id for figure_id in figures if figure_id.startswith("palette_linear_")]

    assert len(docs_and_readme) == 57
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
    assert len(figures) == 57 + 6 + 2

    output_paths = [spec.output_path for spec in figures.values()]
    assert len(output_paths) == len(set(output_paths))


def test_palette_manifest_stays_in_sync_with_palette_file() -> None:
    figures = build_figure_specs()
    palette_names = load_palette_names()

    for palette_name in ("default", "ajisai", "soft_pastels"):
        assert f"palette_circular_{palette_name}" in figures
        assert f"palette_linear_{palette_name}" in figures
    assert len(palette_names) == 55
    assert "palettes_combined_image_1" in figures
    assert "palettes_combined_image_2" in figures


def test_palette_page_is_generated_from_palette_file() -> None:
    assert PALETTE_PAGE.read_text(encoding="utf-8") == render_palette_page()


def test_public_markdown_local_targets_exist() -> None:
    missing: list[str] = []
    for markdown_path in PUBLIC_MARKDOWN:
        source = markdown_path.read_text(encoding="utf-8")
        targets = MARKDOWN_TARGET_RE.findall(source) + HTML_IMAGE_RE.findall(source)
        for raw_target in targets:
            target = _local_target(markdown_path, raw_target)
            if target is not None and not target.exists():
                missing.append(f"{markdown_path.relative_to(PROJECT_ROOT)} -> {raw_target}")
    assert missing == []


def test_public_figures_have_reproduction_inventory_coverage() -> None:
    references = _local_image_references()
    manifest_paths = {spec.output_path for spec in build_figure_specs().values()}
    manual_paths = set(MANUALLY_MANAGED_FIGURES)
    retained_unreferenced = set(UNREFERENCED_FIGURE_RETENTION)

    assert references - manifest_paths - manual_paths == set()
    assert manifest_paths - references == retained_unreferenced
    assert all(reason.strip() for reason in MANUALLY_MANAGED_FIGURES.values())
    assert all(reason.strip() for reason in UNREFERENCED_FIGURE_RETENTION.values())
    assert all((PROJECT_ROOT / path).exists() for path in manual_paths)


def test_alias_resolution_and_support_asset_materialization(tmp_path: Path) -> None:
    reproducer = Reproducer(
        project_root=PROJECT_ROOT,
        output_root=tmp_path / "out",
        figures=build_figure_specs(),
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

        payload = reproducer.report_payload()
        assert payload["aliases_used"]
        assert payload["aliases_used"][0]["requested"] == "NC_000913.gbk"
    finally:
        reproducer.close()


def test_missing_report_structure_for_missing_figure(tmp_path: Path) -> None:
    figure_id = "missing_figure"
    missing_filename = "definitely_missing.gb"
    figures = dict(build_figure_specs())
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
