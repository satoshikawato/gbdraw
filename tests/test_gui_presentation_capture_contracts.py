from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest
from Bio import SeqIO
from PIL import Image

from docs.capture.config import chapter_for

pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
FLOW_PATH = (
    REPO_ROOT / "docs" / "capture" / "flows" / "how_to" / "presentation.py"
)
INDEX_PATH = REPO_ROOT / "gbdraw" / "web" / "index.html"
FIXTURE_ROOT = (
    REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion"
)
GENBANK_PATH = FIXTURE_ROOT / "HmmtDNA.gbk"
VISIBILITY_PATH = FIXTURE_ROOT / "HmmtDNA_feature_visibility.tsv"
COLOR_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "tobacco-plastome-regions"
    / "modified_default_colors.tsv"
)
PRIORITY_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "shared"
    / "cds_gene_qualifier_priority.tsv"
)
SCREENSHOTS = {
    "H-GUI-11": (
        REPO_ROOT / "docs" / "images" / "h-gui-11" / "style-settings.png",
        REPO_ROOT / "docs" / "images" / "h-gui-11" / "style-result.png",
    ),
    "H-GUI-12": (
        REPO_ROOT
        / "docs"
        / "images"
        / "h-gui-12"
        / "presentation-settings.png",
        REPO_ROOT
        / "docs"
        / "images"
        / "h-gui-12"
        / "presentation-result.png",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_presentation_flows_use_one_complete_natural_mitochondrial_record() -> None:
    assert GENBANK_PATH.stat().st_size == 64_640
    assert _sha256(GENBANK_PATH) == (
        "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f"
    )
    records = list(SeqIO.parse(GENBANK_PATH, "genbank"))
    assert len(records) == 1
    record = records[0]
    assert record.id == "NC_012920.1"
    assert len(record) == 16_569
    assert record.annotations["topology"] == "circular"
    assert "complete" in record.description.lower()


def test_presentation_rule_tables_are_frozen_and_gene_first() -> None:
    assert COLOR_PATH.stat().st_size == 95
    assert _sha256(COLOR_PATH) == (
        "e48654dfc5225c8c1eec251f773fc07892228dee906cb1e105e4d24cb5ae8bc1"
    )
    assert PRIORITY_PATH.read_bytes() == b"CDS\tgene\n"
    assert _sha256(PRIORITY_PATH) == (
        "1d3c787c5768b52191cc7e8a325cd00fdc4b1e9b495d37e580c3a69dec21b87a"
    )

    assert VISIBILITY_PATH.stat().st_size == 212
    assert _sha256(VISIBILITY_PATH) == (
        "aca42671c08de16dfbb2317d8b8f53801d1601e956c04aab6de5deaed2a5f348"
    )
    with VISIBILITY_PATH.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "record_id": "NC_012920.1",
            "feature_type": "D-loop",
            "qualifier": "location",
            "value": r"^0\.\.16569$",
            "action": "show",
        },
        {
            "record_id": "NC_012920.1",
            "feature_type": "CDS",
            "qualifier": "product",
            "value": r"^cytochrome c oxidase subunit I$",
            "action": "off",
        },
        {
            "record_id": "*",
            "feature_type": "CDS",
            "qualifier": "product",
            "value": r"^ATP synthase F0 subunit 6$",
            "action": "exclude_matching",
        },
    ]


def test_presentation_flows_use_visible_controls_without_state_injection() -> None:
    source = FLOW_PATH.read_text(encoding="utf-8")
    required = (
        'get_by_role("button", name="Circular", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'get_by_label("GenBank/DDBJ File", exact=True).set_input_files',
        'get_by_label("Palette", exact=True)',
        'get_by_label("Override File (-d)", exact=True).set_input_files',
        'get_by_role("button", name="Whitelist", exact=True)',
        'get_by_role("textbox", name="Qual", exact=True).fill("gene")',
        'get_by_label("Priority File (TSV)", exact=True).set_input_files',
        'get_by_label("Plot Title", exact=True)',
        'get_by_label("Legend Position", exact=True)',
        'get_by_label("Rendering for CDS", exact=True)',
        'get_by_label("Rendering for rRNA", exact=True)',
        'get_by_label("Rendering for tRNA", exact=True)',
        'name="Add feature visibility rule", exact=True',
        'f"Visibility rule {rule_number} record ID"',
        'get_by_label("Block Stroke Width", exact=True)',
        'get_by_label("Line Stroke Width", exact=True)',
        'get_by_label("Resolve Overlaps", exact=True)',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
        "download.save_as",
    )
    for fragment in required:
        assert fragment in source

    for forbidden in (
        "page.locator(",
        "get_by_text(",
        "localStorage",
        "sessionStorage",
        "add_init_script",
        ".runAnalysis(",
        ".downloadSVG(",
        "set_content(",
        "region_start",
        "region_end",
        "part_",
    ):
        assert forbidden not in source


def test_presentation_capture_has_semantic_svg_guards() -> None:
    source = FLOW_PATH.read_text(encoding="utf-8")
    for fragment in (
        'if report.get("featureElementCount") != 37:',
        "missing_labels = STYLE_LABELS - texts",
        'raise AssertionError("CDS labels came from product instead of gene")',
        "legend_order != EXPECTED_LEGEND_ORDER",
        "_normalized_visibility_rules(report) != EXPECTED_VISIBILITY_RULES",
        "if len(underlays) != 22:",
        'raise AssertionError("The off rule did not remove COX1")',
        'raise AssertionError("Exclude from matching incorrectly hid ATP6")',
        'raise AssertionError("D-loop or rRNA did not render as a rectangle")',
        'raise AssertionError("A CDS did not render with arrow-head geometry")',
        "if len(radial_bands) < 2:",
        'raise AssertionError("The legend is not positioned to the right of the record")',
        "assert_static_svg_safety(report)",
    ):
        assert fragment in source


def test_presentation_accessibility_labels_are_public_ui_contracts() -> None:
    source = INDEX_PATH.read_text(encoding="utf-8")
    for label in (
        'summary aria-label="Colors"',
        'aria-label="Palette"',
        'summary aria-label="Features"',
        'aria-label="Arrow head length ratio"',
        'aria-label="Arrow shaft width ratio"',
        'aria-label="Add feature visibility rule"',
        "`Visibility rule ${i + 1} record ID`",
        "`Visibility rule ${i + 1} feature type`",
        "`Visibility rule ${i + 1} qualifier`",
        "`Visibility rule ${i + 1} action`",
        "`Visibility rule ${i + 1} value regex`",
        'aria-label="Block Stroke Width"',
        'aria-label="Line Stroke Width"',
    ):
        assert label in source


def test_presentation_evidence_matches_scenarios_and_screenshots() -> None:
    for scenario_id in ("H-GUI-11", "H-GUI-12"):
        chapter = chapter_for(scenario_id)
        assert chapter["role"] == "evidence"
        assert "destination" not in chapter
        assert tuple(
            REPO_ROOT / screenshot["path"] for screenshot in chapter["screenshots"]
        ) == SCREENSHOTS[scenario_id]


def test_committed_presentation_screenshots_are_full_pinned_viewports() -> None:
    for paths in SCREENSHOTS.values():
        for path in paths:
            assert path.is_file(), path
            assert path.stat().st_size > 150_000
            with Image.open(path) as image:
                assert image.size == (1440, 900)
                assert image.mode == "RGB"
