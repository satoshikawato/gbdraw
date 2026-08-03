from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
DATA_ROOT = WEB_ROOT / "tutorial-data"
INDEX_PATH = WEB_ROOT / "index.html"
LAYOUT_FLOW_PATH = (
    REPO_ROOT / "docs" / "capture" / "flows" / "how_to" / "layouts.py"
)
COMPARISON_FLOW_PATH = (
    REPO_ROOT
    / "docs"
    / "capture"
    / "flows"
    / "how_to"
    / "nucleotide_comparisons.py"
)
MANIFEST_PATH = REPO_ROOT / "docs" / "scenarios" / "manifest.json"

PAGES = {
    "H-GUI-03": REPO_ROOT
    / "docs"
    / "HOW_TO"
    / "GUI"
    / "arrange-linear-records-regions-and-orientation.md",
    "H-GUI-04": REPO_ROOT
    / "docs"
    / "HOW_TO"
    / "GUI"
    / "use-uploaded-blast-results.md",
    "H-GUI-05": REPO_ROOT / "docs" / "HOW_TO" / "GUI" / "use-tlosatx.md",
    "H-GUI-06": REPO_ROOT
    / "docs"
    / "HOW_TO"
    / "GUI"
    / "add-circular-similarity-rings.md",
}
SCREENSHOT_NAMES = {
    "H-GUI-03": ("record-layout.png", "orientation-result.png"),
    "H-GUI-04": ("comparison-plan.png", "comparison-result.png"),
    "H-GUI-05": ("tlosatx-settings.png", "tlosatx-result.png"),
    "H-GUI-06": ("ring-settings.png", "ring-result.png", "hsp-popup.png"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _chapter(scenario_id: str) -> dict[str, object]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    return next(
        chapter for chapter in manifest["chapters"] if chapter["id"] == scenario_id
    )


def _blast_rows(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            row
            for row in csv.reader(handle, delimiter="\t")
            if row and not row[0].startswith("#")
        ]


def test_linear_how_tos_use_whole_lambda_and_de3_records_only() -> None:
    expected = (
        (
            DATA_ROOT / "lambda" / "NC_001416.gb",
            176_723,
            "4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7",
            "NC_001416.1",
            48_502,
        ),
        (
            DATA_ROOT / "de3" / "NC_042057.1.gb",
            111_686,
            "288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09",
            "NC_042057.1",
            42_925,
        ),
    )
    for path, size, digest, record_id, length in expected:
        assert path.stat().st_size == size
        assert _sha256(path) == digest
        records = list(SeqIO.parse(path, "genbank"))
        assert len(records) == 1
        assert records[0].id == record_id
        assert len(records[0]) == length
        assert records[0].annotations["topology"] == "linear"
        assert "complete" in records[0].description.lower()

    flow_source = (
        LAYOUT_FLOW_PATH.read_text(encoding="utf-8")
        + COMPARISON_FLOW_PATH.read_text(encoding="utf-8")
    )
    assert "NC_000866" not in flow_source
    assert "lambda_two_contigs" not in flow_source
    assert "NC_001416_part" not in flow_source


def test_frozen_tlosatx_evidence_has_397_rows_and_seven_qualified_links() -> None:
    path = DATA_ROOT / "lambda-de3-comparison" / "lambda-de3.tlosatx.tsv"
    assert path.stat().st_size == 28_303
    assert _sha256(path) == (
        "483e98f8b3dce172523cf00c82eb9d47e4faee13e4781b06f3d58ea4fb63532d"
    )
    rows = _blast_rows(path)
    assert len(rows) == 397
    assert {len(row) for row in rows} == {12}
    assert {(row[0], row[1]) for row in rows} == {
        ("NC_001416.1", "NC_042057.1")
    }
    qualified = [row for row in rows if int(row[3]) >= 1_000]
    assert len(qualified) == 7
    assert min(int(row[3]) for row in qualified) == 1_016
    assert max(int(row[3]) for row in qualified) == 1_647


def test_circular_ring_flow_uses_four_complete_natural_mitochondria() -> None:
    fixtures = (
        (
            DATA_ROOT / "human-mitochondrion" / "HmmtDNA.gbk",
            "NC_012920.1",
            16_569,
            "Homo sapiens",
            13,
        ),
        (
            DATA_ROOT / "metazoan-mitochondria-four" / "NC_002333.2.gb",
            "NC_002333.2",
            16_596,
            "Danio rerio",
            13,
        ),
        (
            DATA_ROOT / "metazoan-mitochondria-four" / "NC_024511.2.gb",
            "NC_024511.2",
            19_524,
            "Drosophila melanogaster",
            13,
        ),
        (
            DATA_ROOT / "metazoan-mitochondria-four" / "NC_001328.1.gb",
            "NC_001328.1",
            13_794,
            "Caenorhabditis elegans",
            12,
        ),
    )
    for path, record_id, length, organism, cds_gene_count in fixtures:
        records = list(SeqIO.parse(path, "genbank"))
        assert len(records) == 1
        record = records[0]
        assert record.id == record_id
        assert len(record) == length
        assert record.annotations["topology"] == "circular"
        assert record.annotations["organism"] == organism
        assert "complete" in record.description.lower()
        cds_features = [feature for feature in record.features if feature.type == "CDS"]
        assert len(cds_features) == cds_gene_count
        assert all(feature.qualifiers.get("gene") for feature in cds_features)

    priority = DATA_ROOT / "shared" / "cds_gene_qualifier_priority.tsv"
    assert priority.read_text(encoding="utf-8") == "CDS\tgene\n"
    human = SeqIO.read(fixtures[0][0], "genbank")
    assert sum(
        feature.type == "CDS" and bool(feature.qualifiers.get("gene"))
        for feature in human.features
    ) == 13


def test_h_gui_03_to_06_flows_use_visible_controls_without_state_injection() -> None:
    layout = LAYOUT_FLOW_PATH.read_text(encoding="utf-8")
    comparisons = COMPARISON_FLOW_PATH.read_text(encoding="utf-8")
    for fragment in (
        'get_by_role("button", name="Linear", exact=True)',
        'get_by_test_id("linear-genbank-1").set_input_files',
        'f"Definition for sequence {index}"',
        '"Arrange linear records in rows"',
        '"Linear scale style"',
        '"Ruler on Axis"',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
    ):
        assert fragment in layout
    for fragment in (
        'name="Upload BLAST TSV", exact=True',
        '"BLAST TSV", exact=True',
        'name="TLOSATX", exact=True',
        'name="LOSAT execution", exact=True',
        '"Pairwise Comparisons", exact=True',
        '"Separate Strands", exact=True',
        "expect(separate_strands).not_to_be_checked()",
        'name=re.compile(r"Add Seq")',
        '"Priority File (TSV)", exact=True',
        'name="Pairwise match details", exact=True',
        "page.expect_file_chooser",
        "page.expect_download",
    ):
        assert fragment in comparisons

    for source in (layout, comparisons):
        for forbidden in (
            "page.locator(",
            "get_by_text(",
            "localStorage",
            "sessionStorage",
            "add_init_script",
            ".runAnalysis(",
            ".downloadSVG(",
            "set_content(",
        ):
            assert forbidden not in source


def test_nucleotide_capture_accessibility_labels_are_in_the_public_ui() -> None:
    source = INDEX_PATH.read_text(encoding="utf-8")
    for label in (
        "Definition for sequence ${idx + 1}",
        "Region start for sequence ${idx + 1}",
        "Region end for sequence ${idx + 1}",
        "Reverse complement for sequence ${idx + 1}",
        'aria-label="Arrange linear records in rows"',
        "Linear record row for sequence ${idx + 1}",
        'aria-label="Linear record gap"',
        "TLOSATX gencode for sequence ${idx + 1}",
        'aria-label="Linear comparison minimum alignment length"',
        'aria-label="Pairwise Comparisons"',
        'aria-label="Circular reference gencode"',
        "Comparison ring label ${row.index + 1}",
        "Comparison subject gencode ${row.index + 1}",
        'aria-label="Circular comparison minimum alignment length"',
        'aria-label="Circular comparison ring width"',
        'aria-label="Circular comparison ring gap"',
    ):
        assert label in source


def test_nucleotide_pages_and_manifest_use_the_executable_scenarios() -> None:
    expected_values = {
        "H-GUI-03": (
            "NC_001416.1",
            "48,502 bp",
            "NC_042057.1",
            "42,925 bp",
            "linear_regions_orientation.svg",
        ),
        "H-GUI-04": (
            "397",
            "7",
            "uploaded_comparison.svg",
            "lambda-de3.tlosatx.tsv",
        ),
        "H-GUI-05": (
            "Serial",
            "397",
            "lambda-de3.tlosatx.tsv",
            "lambda-de3-tlosatx.svg",
        ),
        "H-GUI-06": (
            "NC_012920.1",
            "NC_002333.2",
            "NC_024511.2",
            "NC_001328.1",
            "cds_gene_qualifier_priority.tsv",
            "circular_similarity_rings.svg",
            "comparison_spans.fasta",
        ),
    }
    for scenario_id, page_path in PAGES.items():
        page = page_path.read_text(encoding="utf-8")
        chapter = _chapter(scenario_id)
        assert chapter["destination"] == str(page_path.relative_to(REPO_ROOT))
        assert tuple(
            Path(screenshot["path"]).name for screenshot in chapter["screenshots"]
        ) == SCREENSHOT_NAMES[scenario_id]
        for value in expected_values[scenario_id]:
            assert value in page

    h3 = _chapter("H-GUI-03")
    assert set(h3["fixtures"]) == {"lambda", "de3"}
    h6 = _chapter("H-GUI-06")
    assert set(h6["fixtures"]) == {
        "human-mitochondrion",
        "metazoan-mitochondria-four",
    }
    assert h6["settings"]["mode"] == "circular"
    assert h6["settings"]["separate_strands"] is False
    assert "lambda" not in h6["fixtures"]
    assert "de3" not in h6["fixtures"]


def test_nucleotide_screenshots_are_full_pinned_viewports() -> None:
    for scenario_id, names in SCREENSHOT_NAMES.items():
        root = REPO_ROOT / "docs" / "images" / scenario_id.lower()
        for name in names:
            path = root / name
            assert path.is_file(), path
            with Image.open(path) as image:
                assert image.size == (1440, 900)
                assert image.mode == "RGB"
