from __future__ import annotations

import hashlib
import json
import runpy
from pathlib import Path

from Bio import SeqIO
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
CAPTURE_ROOT = REPO_ROOT / "docs" / "capture"
WEB_INDEX_PATH = REPO_ROOT / "gbdraw" / "web" / "index.html"
CONFIG_PATH = CAPTURE_ROOT / "config.py"
WEB_CAPTURE_PATH = CAPTURE_ROOT / "flows" / "web_capture.py"
CIRCULAR_FLOW_PATH = (
    CAPTURE_ROOT / "flows" / "tutorials" / "gui_first_circular.py"
)
LINEAR_FLOW_PATH = CAPTURE_ROOT / "flows" / "tutorials" / "gui_first_linear.py"
GUI_INPUTS_FLOW_PATH = CAPTURE_ROOT / "flows" / "how_to" / "inputs.py"
GUI_CIRCULAR_LAYOUT_FLOW_PATH = CAPTURE_ROOT / "flows" / "how_to" / "layouts.py"
GUI_LOSATN_FLOW_PATH = CAPTURE_ROOT / "flows" / "tutorials" / "gui_losatn.py"
SEMANTICS_PATH = CAPTURE_ROOT / "assertions" / "svg_semantics.py"
DOWNLOADS_PATH = CAPTURE_ROOT / "assertions" / "downloads.py"
RUNNER_PATH = CAPTURE_ROOT / "run_all.py"
README_PATH = CAPTURE_ROOT / "README.md"
MANIFEST_PATH = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
CIRCULAR_TUTORIAL_PATH = (
    REPO_ROOT / "docs" / "TUTORIALS" / "GUI" / "first-circular-genome-diagram.md"
)
LINEAR_TUTORIAL_PATH = (
    REPO_ROOT / "docs" / "TUTORIALS" / "GUI" / "first-linear-genome-diagram.md"
)
GUI_INPUTS_HOW_TO_PATH = (
    REPO_ROOT
    / "docs"
    / "HOW_TO"
    / "GUI"
    / "use-genbank-and-gff3-fasta-inputs.md"
)
GUI_CIRCULAR_LAYOUT_HOW_TO_PATH = (
    REPO_ROOT
    / "docs"
    / "HOW_TO"
    / "GUI"
    / "arrange-multiple-circular-records.md"
)
GUI_LOSATN_TUTORIAL_PATH = (
    REPO_ROOT / "docs" / "TUTORIALS" / "GUI" / "compare-genomes-losatn.md"
)
GUI_HOW_TO_INDEX_PATH = REPO_ROOT / "docs" / "HOW_TO" / "GUI" / "README.md"
CIRCULAR_FIXTURE_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "human-mitochondrion"
    / "HmmtDNA.gbk"
)
LINEAR_FIXTURE_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "lambda"
    / "NC_001416.gb"
)
LINEAR_LABEL_RULE_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "shared"
    / "cds_gene_qualifier_priority.tsv"
)
GUI_INPUTS_GFF3_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "lambda-gff3"
    / "NC_001416.gff3"
)
GUI_INPUTS_FASTA_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "lambda-gff3"
    / "NC_001416.fna"
)
GUI_LOSATN_DE3_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "de3"
    / "NC_042057.1.gb"
)
GUI_LOSATN_REFERENCE_TSV_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "lambda-de3-comparison"
    / "lambda-de3.losatn.tsv"
)
GUI_CIRCULAR_LAYOUT_FIXTURE_PATHS = (
    CIRCULAR_FIXTURE_PATH,
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "metazoan-mitochondria-four"
    / "NC_002333.2.gb",
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "metazoan-mitochondria-four"
    / "NC_024511.2.gb",
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "metazoan-mitochondria-four"
    / "NC_001328.1.gb",
)
CIRCULAR_SCREENSHOT_ROOT = REPO_ROOT / "docs" / "images" / "t-gui-01"
LINEAR_SCREENSHOT_ROOT = REPO_ROOT / "docs" / "images" / "t-gui-02"
GUI_INPUTS_SCREENSHOT_ROOT = REPO_ROOT / "docs" / "images" / "h-gui-01"
GUI_LOSATN_SCREENSHOT_ROOT = REPO_ROOT / "docs" / "images" / "t-gui-03"
GUI_CIRCULAR_LAYOUT_SCREENSHOT_ROOT = REPO_ROOT / "docs" / "images" / "h-gui-02"
CIRCULAR_SCREENSHOTS = {
    "01-input-ready.png": "Circular GenBank input showing HmmtDNA.gbk selected",
    "02-first-diagram.png": "First circular human mitochondrial genome diagram",
    "03-publication-label.png": "Circular preview labeled Homo sapiens",
    "04-layout-settings.png": (
        "Circular layout controls set to Middle, Labels Out, and Legend Right"
    ),
    "04-finished-diagram.png": (
        "Finished circular human mitochondrial genome diagram with external labels "
        "and a right legend"
    ),
    "05-export-svg.png": "SVG download button below the finished result preview",
}
LINEAR_SCREENSHOTS = {
    "01-input-ready.png": "Linear GenBank input showing NC_001416.gb selected",
    "02-first-diagram.png": "First linear Lambda genome diagram",
    "03-layout-settings.png": "Linear layout controls configured for labels and a ruler",
    "04-finished-diagram.png": (
        "Finished linear Lambda genome diagram with labels and ruler"
    ),
}
GUI_INPUTS_SCREENSHOTS = {
    "genbank-input.png": "Lambda GenBank file selected in the web app",
    "gff3-fasta-input.png": "Matched Lambda GFF3 and FASTA files selected",
    "id-error.png": "GFF3 and FASTA record-ID validation error",
}
GUI_LOSATN_SCREENSHOTS = {
    "01-input-ready.png": "Two complete GenBank records selected for a Linear diagram",
    "02-first-diagram.png": "Two complete records in a plain linear diagram",
    "03-losatn-settings.png": (
        "LOSATN selected with deterministic single-thread settings"
    ),
    "04-comparison-result.png": (
        "Linear genome comparison with nucleotide similarity links"
    ),
    "05-match-popup.png": "LOSATN match details popup in the result preview",
}
GUI_CIRCULAR_LAYOUT_SCREENSHOTS = {
    "grid-settings.png": "Circular multi-record grid settings",
    "grid-result.png": (
        "Four complete mitochondrial genomes with gene labels arranged in a grid"
    ),
}


def _chapter(scenario_id: str) -> dict[str, object]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    return next(
        chapter for chapter in manifest["chapters"] if chapter["id"] == scenario_id
    )


def test_capture_environment_is_pinned_and_loopback_only() -> None:
    config = runpy.run_path(str(CONFIG_PATH))

    assert config["PYTHON_PLAYWRIGHT_VERSION"] == "1.61.0"
    assert config["NODE_PLAYWRIGHT_VERSION"] == "1.61.1"
    assert config["CHROMIUM_VERSION"] == "149.0.7827.55"
    assert config["CHROMIUM_REVISION"] == "v1228"
    assert (config["VIEWPORT_WIDTH"], config["VIEWPORT_HEIGHT"]) == (1440, 900)
    assert config["DEVICE_SCALE_FACTOR"] == 1
    assert config["LOCALE"] == "en-US"
    assert config["TIMEZONE_ID"] == "UTC"
    assert config["COLOR_SCHEME"] == "light"
    assert config["REDUCED_MOTION"] == "reduce"
    assert config["LOCAL_HOST"] == "127.0.0.1"
    assert config["ISOLATION_HEADERS"]["Cross-Origin-Embedder-Policy"] == (
        "require-corp"
    )
    assert config["IMPLEMENTED_SCENARIO_IDS"] == (
        "T-GUI-01",
        "T-GUI-02",
        "H-GUI-01",
        "T-GUI-03",
        "H-GUI-02",
        "H-GUI-03",
        "H-GUI-04",
        "H-GUI-05",
        "H-GUI-06",
        "T-GUI-04",
        "T-GUI-05",
        "T-GUI-06",
        "T-GUI-08",
        "T-GUI-09",
        "H-GUI-07",
        "H-GUI-08",
        "H-GUI-09",
        "H-GUI-10",
        "H-GUI-11",
        "H-GUI-12",
        "H-GUI-13",
        "H-GUI-14",
        "H-GUI-15",
    )
    assert config["FIRST_CIRCULAR_SCREENSHOT_NAMES"] == tuple(CIRCULAR_SCREENSHOTS)
    assert config["FIRST_LINEAR_SCREENSHOT_NAMES"] == tuple(LINEAR_SCREENSHOTS)
    assert config["GUI_INPUTS_SCREENSHOT_NAMES"] == tuple(GUI_INPUTS_SCREENSHOTS)
    assert config["GUI_LOSATN_SCREENSHOT_NAMES"] == tuple(GUI_LOSATN_SCREENSHOTS)
    assert config["GUI_CIRCULAR_LAYOUT_SCREENSHOT_NAMES"] == tuple(
        GUI_CIRCULAR_LAYOUT_SCREENSHOTS
    )
    assert config["GUI_LINEAR_LAYOUT_SCREENSHOT_NAMES"] == (
        "record-layout.png",
        "orientation-result.png",
    )
    assert config["GUI_UPLOADED_COMPARISON_SCREENSHOT_NAMES"] == (
        "comparison-plan.png",
        "comparison-result.png",
    )
    assert config["GUI_TLOSATX_SCREENSHOT_NAMES"] == (
        "tlosatx-settings.png",
        "tlosatx-result.png",
    )
    assert config["GUI_CIRCULAR_RINGS_SCREENSHOT_NAMES"] == (
        "ring-settings.png",
        "ring-result.png",
        "hsp-popup.png",
    )


def test_first_circular_capture_uses_the_frozen_bundled_fixture() -> None:
    fixture_bytes = CIRCULAR_FIXTURE_PATH.read_bytes()

    assert len(fixture_bytes) == 64_640
    assert hashlib.sha256(fixture_bytes).hexdigest() == (
        "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f"
    )


def test_first_linear_capture_uses_the_complete_frozen_lambda_fixture() -> None:
    fixture_bytes = LINEAR_FIXTURE_PATH.read_bytes()
    label_rule_bytes = LINEAR_LABEL_RULE_PATH.read_bytes()

    assert len(fixture_bytes) == 176_723
    assert hashlib.sha256(fixture_bytes).hexdigest() == (
        "4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7"
    )
    assert len(label_rule_bytes) == 9
    assert hashlib.sha256(label_rule_bytes).hexdigest() == (
        "1d3c787c5768b52191cc7e8a325cd00fdc4b1e9b495d37e580c3a69dec21b87a"
    )
    assert label_rule_bytes == b"CDS\tgene\n"


def test_gui_inputs_capture_uses_the_complete_frozen_gff3_fasta_pair() -> None:
    gff3_bytes = GUI_INPUTS_GFF3_PATH.read_bytes()
    fasta_bytes = GUI_INPUTS_FASTA_PATH.read_bytes()

    assert len(gff3_bytes) == 36_794
    assert hashlib.sha256(gff3_bytes).hexdigest() == (
        "d53e05de87933104cd26111bca42006cce9b5e903fb5b187740f963b3a2098cb"
    )
    assert len(fasta_bytes) == 49_253
    assert hashlib.sha256(fasta_bytes).hexdigest() == (
        "80897a7ee6b8aaffbab5442e0daad292592ac74701dbdf35af4b400ae0770ef3"
    )
    assert fasta_bytes.startswith(b">NC_001416.1 ")


def test_gui_losatn_capture_uses_only_the_qualified_complete_record_pair() -> None:
    lambda_bytes = LINEAR_FIXTURE_PATH.read_bytes()
    de3_bytes = GUI_LOSATN_DE3_PATH.read_bytes()
    reference_bytes = GUI_LOSATN_REFERENCE_TSV_PATH.read_bytes()

    assert len(lambda_bytes) == 176_723
    assert hashlib.sha256(lambda_bytes).hexdigest() == (
        "4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7"
    )
    assert len(de3_bytes) == 111_686
    assert hashlib.sha256(de3_bytes).hexdigest() == (
        "288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09"
    )
    assert len(reference_bytes) == 436
    assert hashlib.sha256(reference_bytes).hexdigest() == (
        "703b0ac749669152a8e6d5fa6fb246cf5973cb1a3e0ca9db2f346ee33628317c"
    )
    rows = [line.split(b"\t") for line in reference_bytes.splitlines()]
    assert len(rows) == 6
    assert all(len(row) == 12 for row in rows)
    assert {(row[0], row[1]) for row in rows} == {
        (b"NC_001416.1", b"NC_042057.1")
    }


def test_gui_circular_layout_uses_four_complete_natural_mitochondrial_records() -> None:
    expected = (
        (
            "NC_012920.1",
            16_569,
            "Homo sapiens",
            64_640,
            "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f",
        ),
        (
            "NC_002333.2",
            16_596,
            "Danio rerio",
            55_541,
            "94ab35da6f81abc2595fcd425c23585ed78d9396b5143918d9d1025d8a4d2140",
        ),
        (
            "NC_024511.2",
            19_524,
            "Drosophila melanogaster",
            72_340,
            "79fa36199682961919c4a11f3a8fc50c9e598e68b867eb25e847bce1aa1c4229",
        ),
        (
            "NC_001328.1",
            13_794,
            "Caenorhabditis elegans",
            39_227,
            "8de5f7cf3686f493ee5b8068dba39b31c5d02e70d997063f65ed19d0fa859a59",
        ),
    )

    for path, (record_id, length, organism, size, digest) in zip(
        GUI_CIRCULAR_LAYOUT_FIXTURE_PATHS,
        expected,
        strict=True,
    ):
        contents = path.read_bytes()
        assert len(contents) == size
        assert hashlib.sha256(contents).hexdigest() == digest
        records = list(SeqIO.parse(path, "genbank"))
        assert len(records) == 1
        record = records[0]
        assert record.id == record_id
        assert len(record) == length
        assert record.annotations["organism"] == organism
        assert record.annotations["topology"] == "circular"
        assert "complete" in record.description.lower()


def test_first_circular_flow_uses_accessible_real_actions_without_state_shortcuts() -> None:
    source = CIRCULAR_FLOW_PATH.read_text(encoding="utf-8")
    shared_source = WEB_CAPTURE_PATH.read_text(encoding="utf-8")

    for required in (
        'get_by_role("button", name="Circular", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'get_by_label("GenBank/DDBJ File", exact=True).set_input_files',
        'get_by_label("Output Prefix", exact=True)',
        'get_by_label("Species", exact=True)',
        'get_by_label("Track Preset", exact=True)',
        'get_by_label("Separate Strands", exact=True)',
        'get_by_label("Hide GC Content", exact=True)',
        'get_by_label("Hide GC Skew", exact=True)',
        'get_by_label("Legend Position", exact=True)',
        'get_by_label("Label Mode", exact=True)',
        'get_by_label("Priority File (TSV)", exact=True).set_input_files',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
        "download.save_as",
    ):
        assert required in source

    for required_value in (
        'fill("human_mitochondrion")',
        'fill("<i>Homo sapiens</i>")',
        'select_option("middle")',
        'select_option("out")',
        'select_option("right")',
    ):
        assert required_value in source

    assert 'get_by_role("button", name="Generate Diagram", exact=True)' in shared_source
    assert 'get_by_role("region", name="Result Preview", exact=True)' in shared_source
    assert 'context.route("**/*", enforce_local_requests)' in shared_source
    assert "external_requests.append(request_url)" in shared_source

    for forbidden in (
        "page.locator(",
        "get_by_text(",
        "localStorage",
        "sessionStorage",
        "add_init_script",
        ".runAnalysis(",
        ".downloadSVG(",
        ".downloadPNG(",
        "set_content(",
    ):
        assert forbidden not in source
        assert forbidden not in shared_source


def test_first_linear_flow_uses_accessible_real_actions_without_state_shortcuts() -> None:
    source = LINEAR_FLOW_PATH.read_text(encoding="utf-8")

    for required in (
        'get_by_role("button", name="Linear", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'name="No comparison", exact=True',
        'get_by_label("GenBank File", exact=True).set_input_files',
        'get_by_label("Output Prefix", exact=True)',
        'get_by_label("Track Layout", exact=True)',
        'get_by_label("Separate Strands", exact=True)',
        'get_by_label("Labels", exact=True)',
        'get_by_label("Show Labels", exact=True)',
        'get_by_label("Priority File (TSV)", exact=True).set_input_files',
        'get_by_label("Legend Position", exact=True)',
        'get_by_label("Axis & Scale", exact=True)',
        'get_by_label("Show Coordinate Scale (Linear)", exact=True)',
        'get_by_label("Linear scale style", exact=True)',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
        "download.save_as",
    ):
        assert required in source

    for required_value in (
        'fill("lambda_linear")',
        'select_option("middle")',
        'select_option("all")',
        'select_option("left")',
        'select_option("ruler")',
    ):
        assert required_value in source

    for forbidden in (
        "page.locator(",
        "get_by_text(",
        "localStorage",
        "sessionStorage",
        "add_init_script",
        ".runAnalysis(",
        ".downloadSVG(",
        ".downloadPNG(",
        "set_content(",
        "NC_001416_part",
        "region_start",
        "region_end",
    ):
        assert forbidden not in source


def test_gui_inputs_flow_uses_whole_files_real_actions_and_actual_error() -> None:
    source = GUI_INPUTS_FLOW_PATH.read_text(encoding="utf-8")

    for required in (
        'get_by_role("button", name="Linear", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'name="GFF3 + FASTA", exact=True',
        'name="No comparison", exact=True',
        'get_by_label("GenBank File", exact=True).set_input_files',
        'get_by_label("GFF3", exact=True).set_input_files',
        'get_by_label("FASTA", exact=True).set_input_files',
        'get_by_role("button", name="SVG", exact=True)',
        'get_by_role("alert", name="Generation Error", exact=True)',
        "page.expect_download",
        "download.save_as",
        "assert_matching_lambda_input_semantics",
        'b">MISMATCHED_ID"',
        "written_sequence != sequence",
    ):
        assert required in source

    for message_part in (
        "No matching FASTA record found for GFF record NC_001416.1.",
        "Please ensure that all GFF records have corresponding FASTA entries.",
    ):
        assert message_part in source

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
        "NC_001416_part",
    ):
        assert forbidden not in source


def test_gui_circular_layout_flow_uses_complete_records_and_public_controls() -> None:
    source = GUI_CIRCULAR_LAYOUT_FLOW_PATH.read_text(encoding="utf-8")
    circular_source = source[
        source.index("def capture_gui_circular_layout(") : source.index(
            "\ndef capture_gui_linear_layout("
        )
    ]

    for required in (
        'get_by_role("button", name="Circular", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'get_by_label("GenBank/DDBJ File", exact=True).set_input_files',
        'get_by_label("Output Prefix", exact=True)',
        'get_by_label("Multi-Record Canvas", exact=True)',
        'get_by_label("Record Size Mode", exact=True)',
        'get_by_label("Min Radius Ratio", exact=True)',
        'get_by_label("Column Gap Ratio", exact=True)',
        'get_by_label("Row Gap Ratio", exact=True)',
        'get_by_label("Title & Legend", exact=True)',
        'get_by_label("Plot Title", exact=True)',
        'get_by_label("Plot Title Position", exact=True)',
        'get_by_label("Definition Font Size", exact=True)',
        '"Keep Full Definition with Plot Title"',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
        "download.save_as",
        "SeqIO.parse",
        'str(record.annotations.get("topology", "")).lower() != "circular"',
        '"complete" not in record.description.lower()',
        'combined_bytes = b"".join(source_bytes)',
        'str(combined.seq) != str(source.seq)',
        "assert_gui_circular_layout_svg",
        "assert_gui_circular_layout_svg_download",
        "window.getSelection()?.removeAllRanges()",
        "selection_range_count != 0",
    ):
        assert required in source

    for required_value in (
        'fill("multi_record_circular")',
        'select_option("equal")',
        'fill("0.75")',
        'fill("0.40")',
        'fill("0.08")',
        'fill("Complete metazoan mitochondrial genomes")',
        'select_option("top")',
        'fill("20")',
        'select_option("out")',
        "FIRST_LINEAR_LABEL_RULE_PATH",
        'fill("16")',
    ):
        assert required_value in source

    for record_id in (
        "NC_012920.1",
        "NC_002333.2",
        "NC_024511.2",
        "NC_001328.1",
    ):
        assert record_id in source

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
        "BGC0000708",
        "BGC0000713",
    ):
        assert forbidden not in circular_source


def test_gui_losatn_flow_runs_the_real_serial_one_thread_journey() -> None:
    source = GUI_LOSATN_FLOW_PATH.read_text(encoding="utf-8")
    index_source = WEB_INDEX_PATH.read_text(encoding="utf-8")

    assert 'aria-label="Pairwise Match Height"' in index_source

    for required in (
        'get_by_role("button", name="Linear", exact=True)',
        'name="No comparison", exact=True',
        'get_by_role("button", name="Add sequence", exact=True)',
        'get_by_test_id("linear-genbank-1").set_input_files',
        'get_by_test_id("linear-genbank-2").set_input_files',
        'name="Run LOSAT", exact=True',
        'get_by_role("radio", name="LOSATN", exact=True)',
        'name="LOSAT execution", exact=True',
        'name="LOSAT total threads", exact=True',
        'name="LOSAT parallel runs", exact=True',
        'name="LOSAT threads per run", exact=True',
        'name="LOSATN task", exact=True',
        'name="Raw LOSAT filename", exact=True',
        'get_by_label("Pairwise Match", exact=True)',
        '"Pairwise Match Height", exact=True',
        'pairwise_match_height.fill("120")',
        'select_option("serial")',
        'select_option("1")',
        'select_option("megablast")',
        'name="Save Raw LOSAT TSV", exact=True',
        'get_by_role("button", name="SVG", exact=True)',
        'name="Pairwise match 1", exact=True',
        'first_match.press("Enter")',
        'get_by_role("dialog", name="Pairwise match details", exact=True)',
        "page.expect_download",
        "assert_gui_losatn_tsv_download",
        "assert_gui_losatn_svg_download",
    ):
        assert required in source

    first_result = source.index("first_report = generate_and_inspect")
    configure_losat = source.index('name="Run LOSAT", exact=True')
    final_result = source.index("final_report = generate_and_inspect")
    tsv_download = source.index('name="Save Raw LOSAT TSV", exact=True')
    svg_download = source.index('name="SVG", exact=True')
    popup = source.index('name="Pairwise match 1", exact=True')
    assert first_result < configure_losat < final_result
    assert final_result < tsv_download < svg_download < popup

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
        "NC_001416_part",
    ):
        assert forbidden not in source
    assert "tlosatx" not in source.lower()


def test_download_contract_parses_and_validates_all_static_svgs() -> None:
    downloads_source = DOWNLOADS_PATH.read_text(encoding="utf-8")
    semantics_source = SEMANTICS_PATH.read_text(encoding="utf-8")

    assert 'EXPECTED_FIRST_CIRCULAR_SVG = "human_mitochondrion.svg"' in (
        downloads_source
    )
    assert 'EXPECTED_FIRST_LINEAR_SVG = "lambda_linear.svg"' in downloads_source
    assert 'EXPECTED_GUI_INPUTS_SVG = "lambda_gff3.svg"' in downloads_source
    assert 'EXPECTED_GUI_LOSATN_TSV = "lambda-de3.losatn.tsv"' in downloads_source
    assert 'EXPECTED_GUI_LOSATN_SVG = "lambda-de3-losatn.svg"' in downloads_source
    assert (
        'EXPECTED_GUI_CIRCULAR_LAYOUT_SVG = "multi_record_circular.svg"'
        in downloads_source
    )
    assert downloads_source.count("inspect_svg_file(path)") == 5
    assert "assert_finished_circular_svg(report)" in downloads_source
    assert "assert_finished_linear_svg(report)" in downloads_source
    for required in (
        '"NC_001416.1"',
        '"48,502 bp"',
        "feature_count != 73",
        '"5 kbp"',
        '"45 kbp"',
        '{"A", "B", "J", "int"}',
        '{"+": 47, "-": 26, "undefined": 0}',
        "Rendered Lambda feature IDs do not match generated CDS metadata",
        "EXPECTED_GUI_LOSATN_MATCHES",
        "featureElementCount\") != 130",
        "Rendered LOSATN link metadata does not match the qualified six-row TSV",
        "assert_gui_circular_layout_svg",
        '"NC_002333.2"',
        '"NC_024511.2"',
        '"NC_001328.1"',
        'featureElementCount") != 147',
        "Expected a 2 by 2 Circular grid",
        "Circular record bounds overlap on the shared grid",
        "assert_static_svg_safety(report)",
    ):
        assert required in semantics_source


def test_runner_regenerates_and_checks_every_implemented_scenario() -> None:
    source = RUNNER_PATH.read_text(encoding="utf-8")

    assert '"--scenario"' in source
    assert '"--tier"' in source
    assert '"--check"' in source
    assert '"T-GUI-01": ScenarioCapture' in source
    assert '"T-GUI-02": ScenarioCapture' in source
    assert '"H-GUI-01": ScenarioCapture' in source
    assert '"T-GUI-03": ScenarioCapture' in source
    assert '"H-GUI-02": ScenarioCapture' in source
    for scenario_id in ("H-GUI-03", "H-GUI-04", "H-GUI-05", "H-GUI-06"):
        assert f'"{scenario_id}": ScenarioCapture' in source
    assert "first_circular_screenshot_paths" in source
    assert "first_linear_screenshot_paths" in source
    assert "gui_inputs_screenshot_paths" in source
    assert "gui_losatn_screenshot_paths" in source
    assert "gui_circular_layout_screenshot_paths" in source
    assert "gui_linear_layout_screenshot_paths" in source
    assert "gui_uploaded_comparison_screenshot_paths" in source
    assert "gui_tlosatx_screenshot_paths" in source
    assert "gui_circular_rings_screenshot_paths" in source
    assert "capture_first_circular" in source
    assert "capture_first_linear" in source
    assert "capture_gui_inputs" in source
    assert "capture_gui_losatn" in source
    assert "capture_gui_circular_layout" in source
    assert "capture_gui_linear_layout" in source
    assert "capture_gui_uploaded_comparison" in source
    assert "capture_gui_tlosatx" in source
    assert "capture_gui_circular_rings" in source
    assert 'tier="extended"' in source
    assert "TIER_RANK" in source
    assert "ImageChops.difference" in source
    assert "for name, committed_path in committed_paths.items()" in source
    assert source.count("run_scenario as run_cli_scenario") == 1
    assert source.count("run_scenario as run_python_scenario") == 1
    assert "_run_recipes(cli_ids, run_cli_scenario, check=args.check)" in source
    assert "_run_recipes(python_ids, run_python_scenario, check=args.check)" in source


def test_runner_all_orchestrates_gui_cli_and_python_with_check(monkeypatch) -> None:
    from docs.capture import run_all

    gui_calls: list[tuple[str, ...]] = []
    cli_calls: list[tuple[str, bool]] = []
    python_calls: list[tuple[str, bool]] = []

    monkeypatch.setattr(
        run_all,
        "_check",
        lambda scenario_ids, _paths: gui_calls.append(scenario_ids),
    )
    monkeypatch.setattr(
        run_all,
        "run_cli_scenario",
        lambda scenario_id, *, check: cli_calls.append((scenario_id, check)) or (),
    )
    monkeypatch.setattr(
        run_all,
        "run_python_scenario",
        lambda scenario_id, *, check: python_calls.append((scenario_id, check)) or (),
    )

    assert run_all.main(["--tier", "core", "--check"]) == 0
    assert gui_calls == [
        ("T-GUI-01", "T-GUI-02", "H-GUI-01", "H-GUI-15")
    ]
    assert cli_calls == [(scenario_id, True) for scenario_id in run_all.CLI_SCENARIO_IDS]
    assert python_calls == [
        (scenario_id, True) for scenario_id in run_all.PYTHON_SCENARIO_IDS
    ]


def test_runner_focused_scenarios_execute_only_the_selected_surface(monkeypatch) -> None:
    from docs.capture import run_all

    calls: list[tuple[str, str, bool]] = []
    monkeypatch.setattr(
        run_all,
        "_capture_many",
        lambda scenario_ids, _paths: calls.append(("gui", scenario_ids[0], False)),
    )
    monkeypatch.setattr(
        run_all,
        "_check",
        lambda scenario_ids, _paths: calls.append(("gui", scenario_ids[0], True)),
    )
    monkeypatch.setattr(
        run_all,
        "run_cli_scenario",
        lambda scenario_id, *, check: calls.append(("cli", scenario_id, check)) or (),
    )
    monkeypatch.setattr(
        run_all,
        "run_python_scenario",
        lambda scenario_id, *, check: calls.append(("python", scenario_id, check))
        or (),
    )

    assert run_all.main(["--scenario", "T-GUI-01"]) == 0
    assert run_all.main(["--scenario", "H-CLI-03", "--check"]) == 0
    assert run_all.main(["--scenario", "H-PY-02"]) == 0
    assert calls == [
        ("gui", "T-GUI-01", False),
        ("cli", "H-CLI-03", True),
        ("python", "H-PY-02", False),
    ]


def test_t_gui_01_manifest_owns_the_complete_verified_journey() -> None:
    chapter = _chapter("T-GUI-01")
    manifest_screenshots = {
        Path(item["path"]).name: item["alt"] for item in chapter["screenshots"]
    }

    assert chapter["execution"]["path"] == (
        "docs/capture/flows/tutorials/gui_first_circular.py"
    )
    assert chapter["execution"]["expected_outputs"] == ["human_mitochondrion.svg"]
    assert chapter["settings"]["cds_label_qualifier"] == "gene"
    assert {
        "cds_labels_use_gene=true",
        "cds_product_labels_absent=true",
    } <= set(chapter["execution"]["assertions"])
    assert chapter["destination"] == (
        "docs/TUTORIALS/GUI/first-circular-genome-diagram.md"
    )
    assert manifest_screenshots == CIRCULAR_SCREENSHOTS
    assert chapter["status"] == {"implementation": "verified", "review": "approved"}


def test_t_gui_02_manifest_owns_the_complete_verified_journey() -> None:
    chapter = _chapter("T-GUI-02")
    manifest_screenshots = {
        Path(item["path"]).name: item["alt"] for item in chapter["screenshots"]
    }

    assert chapter["fixtures"] == ["lambda"]
    assert chapter["settings"] == {
        "mode": "linear",
        "comparison_mode": "none",
        "ruler": True,
    }
    assert chapter["execution"]["path"] == (
        "docs/capture/flows/tutorials/gui_first_linear.py"
    )
    assert chapter["execution"]["expected_outputs"] == ["lambda_linear.svg"]
    assert chapter["destination"] == (
        "docs/TUTORIALS/GUI/first-linear-genome-diagram.md"
    )
    assert manifest_screenshots == LINEAR_SCREENSHOTS
    assert chapter["status"] == {"implementation": "verified", "review": "approved"}


def test_h_gui_01_manifest_owns_the_complete_verified_journey() -> None:
    chapter = _chapter("H-GUI-01")
    manifest_screenshots = {
        Path(item["path"]).name: item["alt"] for item in chapter["screenshots"]
    }

    assert chapter["fixtures"] == ["lambda", "lambda-gff3"]
    assert chapter["settings"] == {"mode": "linear", "external_network": False}
    assert chapter["execution"]["path"] == "docs/capture/flows/how_to/inputs.py"
    assert chapter["execution"]["expected_outputs"] == ["lambda_gff3.svg"]
    assert chapter["destination"] == (
        "docs/HOW_TO/GUI/use-genbank-and-gff3-fasta-inputs.md"
    )
    assert manifest_screenshots == GUI_INPUTS_SCREENSHOTS
    assert chapter["status"] == {"implementation": "verified", "review": "approved"}


def test_h_gui_02_manifest_owns_the_complete_verified_journey() -> None:
    chapter = _chapter("H-GUI-02")
    manifest_screenshots = {
        Path(item["path"]).name: item["alt"] for item in chapter["screenshots"]
    }

    assert chapter["fixtures"] == [
        "human-mitochondrion",
        "metazoan-mitochondria-four",
    ]
    assert chapter["settings"] == {"mode": "circular", "grouping": "grid"}
    assert chapter["execution"]["path"] == "docs/capture/flows/how_to/layouts.py"
    assert chapter["execution"]["tier"] == "extended"
    assert chapter["execution"]["expected_outputs"] == [
        "multi_record_circular.svg"
    ]
    assert chapter["execution"]["assertions"] == [
        "record_count=4",
        "record_ids_match=true",
        "all_source_topologies=circular",
        "all_records_are_complete=true",
        "organism_labels_match=true",
        "feature_count=147",
        "labels_all_records=true",
        "cds_labels_use_gene=true",
        "cds_product_labels_absent=true",
        "grid_topology_match=true",
        "records_do_not_overlap=true",
        "preview_text_selection=none",
    ]
    assert chapter["destination"] == (
        "docs/HOW_TO/GUI/arrange-multiple-circular-records.md"
    )
    assert manifest_screenshots == GUI_CIRCULAR_LAYOUT_SCREENSHOTS
    assert chapter["status"] == {"implementation": "verified", "review": "approved"}


def test_t_gui_03_manifest_owns_the_complete_verified_journey() -> None:
    chapter = _chapter("T-GUI-03")
    manifest_screenshots = {
        Path(item["path"]).name: item["alt"] for item in chapter["screenshots"]
    }

    assert chapter["fixtures"] == ["lambda", "de3", "lambda-de3-comparison"]
    assert chapter["settings"] == {
        "mode": "linear",
        "program": "losatn",
        "task": "megablast",
        "scheduling": "serial",
        "threads": 1,
        "pairwise_match_height": 120,
    }
    assert chapter["execution"]["path"] == (
        "docs/capture/flows/tutorials/gui_losatn.py"
    )
    assert chapter["execution"]["tier"] == "extended"
    assert chapter["execution"]["expected_outputs"] == [
        "lambda-de3.losatn.tsv",
        "lambda-de3-losatn.svg",
    ]
    assert chapter["execution"]["assertions"] == [
        "first_result_step=2",
        "comparison.rows=6",
        "pairwise_match_height=120",
        "svg.links=6",
        "download.endpoints_match=true",
    ]
    assert chapter["destination"] == (
        "docs/TUTORIALS/GUI/compare-genomes-losatn.md"
    )
    assert manifest_screenshots == GUI_LOSATN_SCREENSHOTS
    assert chapter["status"] == {"implementation": "verified", "review": "approved"}


def test_circular_tutorial_follows_steps_and_defers_related_links() -> None:
    tutorial = CIRCULAR_TUTORIAL_PATH.read_text(encoding="utf-8")
    headings = [
        "## What you'll need",
        "## Step 1: Load the bundled mitochondrial genome",
        "## Step 2: Generate the first diagram",
        "## Step 3: Add a publication label",
        "## Step 4: Make the feature map easier to read",
        "## Step 5: Export the SVG",
        "## What you built",
        "## Next steps",
    ]

    positions = [tutorial.index(heading) for heading in headings]
    assert positions == sorted(positions)
    assert tutorial.index("04-finished-diagram.png") < positions[0]
    for name, alt in CIRCULAR_SCREENSHOTS.items():
        image = f"![{alt}](../../images/t-gui-01/{name})"
        assert image in tutorial

    for value in (
        "`human_mitochondrion`",
        "`<i>Homo sapiens</i>`",
        "| Track Preset | Middle |",
        "| Separate Strands | On |",
        "| Hide GC Content | Off |",
        "| Hide GC Skew | Off |",
        "| Label Mode | Out |",
        "| Priority File (TSV) | `cds_gene_qualifier_priority.tsv` |",
        "| Legend Position | Right |",
        "`human_mitochondrion.svg`",
        "all 13 CDS labels from `gene`",
    ):
        assert value in tutorial

    what_you_built = tutorial.index("## What you built")
    for related_target in (
        "first-linear-genome-diagram.md",
        "../../HOW_TO/GUI/style-features-labels-titles-and-legends.md",
        "compare-genomes-losatn.md",
        "../../HOW_TO/GUI/save-restore-undo-and-reproduce-work.md",
    ):
        assert tutorial.index(related_target) > what_you_built


def test_linear_tutorial_shows_the_step_two_result_and_defers_related_links() -> None:
    tutorial = LINEAR_TUTORIAL_PATH.read_text(encoding="utf-8")
    headings = [
        "## What you'll need",
        "## Step 1: Load the bundled Lambda genome",
        "## Step 2: Generate the first diagram",
        "## Step 3: Add concise labels and a ruler",
        "## Step 4: Regenerate and export the SVG",
        "## What you built",
        "## Next steps",
    ]

    positions = [tutorial.index(heading) for heading in headings]
    assert positions == sorted(positions)
    assert tutorial.index("04-finished-diagram.png") < positions[0]
    for name, alt in LINEAR_SCREENSHOTS.items():
        image = f"![{alt}](../../images/t-gui-02/{name})"
        assert image in tutorial

    for value in (
        "`NC_001416.1`",
        "`48,502 bp`",
        "without cropping or splitting it",
        "| Output Prefix | `lambda_linear` |",
        "| Track Layout | Features on axis |",
        "| Separate Strands | On |",
        "| Show Labels | All Records |",
        "| Priority File (TSV) | `cds_gene_qualifier_priority.tsv` |",
        "| Show Coordinate Scale | On |",
        "| Scale Style | Ruler (Ticks) |",
        "| Legend Position | Left |",
        "`lambda_linear.svg`",
        "all 73 displayed CDS features",
    ):
        assert value in tutorial

    what_you_built = tutorial.index("## What you built")
    for related_target in (
        "../../HOW_TO/GUI/arrange-linear-records-regions-and-orientation.md",
        "../../HOW_TO/GUI/style-features-labels-titles-and-legends.md",
        "compare-genomes-losatn.md",
        "../../REFERENCE/output-formats-and-export.md",
    ):
        assert tutorial.index(related_target) > what_you_built
        assert (LINEAR_TUTORIAL_PATH.parent / related_target).resolve().is_file()


def test_gui_inputs_how_to_documents_both_paths_and_the_actual_error() -> None:
    how_to = GUI_INPUTS_HOW_TO_PATH.read_text(encoding="utf-8")
    headings = [
        "## Before you start",
        "## Use the GenBank input",
        "## Use a matched GFF3 and FASTA pair",
        "## Diagnose a record-ID mismatch",
        "## Verify the result",
        "## Related guides",
    ]
    positions = [how_to.index(heading) for heading in headings]
    assert positions == sorted(positions)

    for name, alt in GUI_INPUTS_SCREENSHOTS.items():
        assert f"![{alt}](../../images/h-gui-01/{name})" in how_to

    for value in (
        "`NC_001416.gb`",
        "`NC_001416.gff3`",
        "`NC_001416.fna`",
        "complete `NC_001416.1` record",
        "do not crop or split the 48,502 bp sequence",
        "| Rendered CDS features | 73 |",
        "| Positive-strand CDS | 47 |",
        "| Negative-strand CDS | 26 |",
        "`lambda_gff3.svg`",
        "No matching FASTA record found for GFF record NC_001416.1.",
    ):
        assert value in how_to

    for related_target in (
        "../../GFF3_FASTA.md",
        "../../REFERENCE/input-formats-and-tsv-schemas.md",
        "../../TUTORIALS/GUI/first-linear-genome-diagram.md",
    ):
        assert (GUI_INPUTS_HOW_TO_PATH.parent / related_target).resolve().is_file()

    assert "./use-genbank-and-gff3-fasta-inputs.md" in (
        GUI_HOW_TO_INDEX_PATH.read_text(encoding="utf-8")
    )


def test_gui_circular_layout_how_to_uses_only_complete_circular_records() -> None:
    how_to = GUI_CIRCULAR_LAYOUT_HOW_TO_PATH.read_text(encoding="utf-8")
    headings = [
        "## Before you start",
        "## Load all four complete records",
        "## Set the 2 by 2 grid",
        "## Generate and verify the grid",
        "## Troubleshooting",
        "## Related guides",
    ]
    positions = [how_to.index(heading) for heading in headings]
    assert positions == sorted(positions)

    for name, alt in GUI_CIRCULAR_LAYOUT_SCREENSHOTS.items():
        assert f"![{alt}](../../images/h-gui-02/{name})" in how_to

    for value in (
        "complete, naturally",
        "does not crop a record, split a sequence into artificial contigs",
        "*Homo sapiens*",
        "*Danio rerio*",
        "*Drosophila melanogaster*",
        "*Caenorhabditis elegans*",
        "`NC_012920.1`",
        "`NC_002333.2`",
        "`NC_024511.2`",
        "`NC_001328.1`",
        "16,569 bp",
        "16,596 bp",
        "19,524 bp",
        "13,794 bp",
        "Put the four unchanged",
        "GenBank entries in one file",
        "do not concatenate their nucleotide sequences",
        "| Record Size Mode | Equal |",
        "| Min Radius Ratio | `0.75` |",
        "| Column Gap Ratio | `0.40` |",
        "| Row Gap Ratio | `0.08` |",
        "| Label Mode | Out |",
        "`cds_gene_qualifier_priority.tsv`",
        "CDS labels use concise gene",
        "rather than the longer `product`",
        "descriptions",
        "Every circle has feature labels",
        "**Keep Full Definition with Plot Title**",
        "147 displayed mitochondrial features",
        "Circular mode",
        "must not be used to disguise a partial genomic region",
        "`multi_record_circular.svg`",
    ):
        assert value in how_to

    assert "BGC0000708" not in how_to
    assert "BGC0000713" not in how_to
    assert "./arrange-multiple-circular-records.md" in (
        GUI_HOW_TO_INDEX_PATH.read_text(encoding="utf-8")
    )


def test_gui_losatn_tutorial_preserves_the_approved_five_step_journey() -> None:
    tutorial = GUI_LOSATN_TUTORIAL_PATH.read_text(encoding="utf-8")
    headings = [
        "## What you'll need",
        "## Step 1: Load both complete genomes",
        "## Step 2: Generate the map without comparison links",
        "## Step 3: Configure LOSATN",
        "## Step 4: Run LOSATN and download the evidence",
        "## Step 5: Inspect one nucleotide match",
        "## What you built",
        "## Next steps",
    ]
    positions = [tutorial.index(heading) for heading in headings]
    assert positions == sorted(positions)
    assert tutorial.index("04-comparison-result.png") < positions[0]

    for name, alt in GUI_LOSATN_SCREENSHOTS.items():
        assert f"![{alt}](../../images/t-gui-03/{name})" in tutorial

    for value in (
        "`NC_001416.1` (48,502 bp)",
        "`NC_042057.1` (42,925 bp)",
        "| Execution | Serial |",
        "| Total threads | 1 |",
        "| Parallel runs | 1 run |",
        "| Threads per run | Fixed at 1 |",
        "| Task | `megablast` |",
        "| Pairwise Match Height | `120` |",
        "`lambda-de3.losatn.tsv`",
        "`lambda-de3-losatn.svg`",
        "six qualified `megablast` rows",
        "Lambda 1..21231 to DE3 20081..41311",
    ):
        assert value in tutorial

    what_you_built = tutorial.index("## What you built")
    for related_target in (
        "../../HOW_TO/GUI/use-uploaded-blast-results.md",
        "../../HOW_TO/GUI/arrange-linear-records-regions-and-orientation.md",
        "../../REFERENCE/output-formats-and-export.md",
        "../../REFERENCE/input-formats-and-tsv-schemas.md",
    ):
        assert tutorial.index(related_target) > what_you_built
        assert (GUI_LOSATN_TUTORIAL_PATH.parent / related_target).resolve().is_file()

    assert "tlosatx" not in tutorial.lower()


def test_all_committed_gui_screenshots_match_the_display_contract() -> None:
    for root, screenshots in (
        (CIRCULAR_SCREENSHOT_ROOT, CIRCULAR_SCREENSHOTS),
        (LINEAR_SCREENSHOT_ROOT, LINEAR_SCREENSHOTS),
        (GUI_INPUTS_SCREENSHOT_ROOT, GUI_INPUTS_SCREENSHOTS),
        (GUI_LOSATN_SCREENSHOT_ROOT, GUI_LOSATN_SCREENSHOTS),
        (GUI_CIRCULAR_LAYOUT_SCREENSHOT_ROOT, GUI_CIRCULAR_LAYOUT_SCREENSHOTS),
    ):
        for name in screenshots:
            path = root / name
            assert path.is_file(), name
            with Image.open(path) as screenshot:
                assert screenshot.format == "PNG", name
                assert screenshot.mode == "RGB", name
                assert screenshot.size == (1440, 900), name


def test_capture_readme_records_versions_regeneration_and_download_checks() -> None:
    readme = README_PATH.read_text(encoding="utf-8")

    assert "Python Playwright 1.61.0" in readme
    assert "Node Playwright 1.61.1" in readme
    assert "Chrome for Testing 149.0.7827.55" in readme
    assert "Playwright Chromium revision v1228" in readme
    assert "writes every manifest-owned screenshot" in readme
    assert "each file is checked for its name, structure, biological endpoints" in readme
    assert "--scenario T-GUI-01" in readme
    assert "--scenario T-GUI-02" in readme
    assert "--scenario H-GUI-01" in readme
    assert "--scenario T-GUI-03" in readme
    assert "--scenario H-GUI-02" in readme
    for scenario_id in ("H-GUI-03", "H-GUI-04", "H-GUI-05", "H-GUI-06"):
        assert f"--scenario {scenario_id}" in readme
    assert "--tier extended" in readme
    assert "--check" in readme
