from __future__ import annotations

import ast
import hashlib
import json
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
BGC_ROOT = WEB_ROOT / "tutorial-data" / "aminoglycoside-bgc-five"
MANIFEST_PATH = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
FIXTURE_MANIFEST_PATH = WEB_ROOT / "tutorial-data" / "manifest.json"
COMMON_FLOW_PATH = REPO_ROOT / "docs" / "capture" / "flows" / "bgc_losatp.py"
TUTORIAL_FLOW_PATH = (
    REPO_ROOT
    / "docs"
    / "capture"
    / "flows"
    / "tutorials"
    / "gui_losatp_groups.py"
)
HOW_TO_FLOW_PATH = (
    REPO_ROOT / "docs" / "capture" / "flows" / "how_to" / "comparisons.py"
)
HEPATOPLASMA_FLOW_PATH = (
    REPO_ROOT
    / "docs"
    / "capture"
    / "flows"
    / "hepatoplasmataceae_losatp.py"
)
GALLERY_SVG_PATH = (
    WEB_ROOT / "gallery" / "examples" / "BGC0000708-BGC0000713.svg"
)
APP_SETUP_PATH = WEB_ROOT / "js" / "app" / "app-setup.js"
RUN_ANALYSIS_PATH = WEB_ROOT / "js" / "app" / "run-analysis.js"
SVG_ACTIONS_PATH = WEB_ROOT / "js" / "app" / "feature-editor" / "svg-actions.js"

FIXTURES = (
    (
        "BGC0000708.gbk",
        105_185,
        "9a5f971c5ed8c406b20574fb50aac567609deb787eb1e8d4635050aa264a04b0",
        "BGC0000708",
        40_579,
        30,
    ),
    (
        "BGC0000709.gbk",
        121_511,
        "4b66b7e4b78d429d12176e1e36d0e48178c562a9d128d4308b38753af9995255",
        "BGC0000709",
        50_466,
        38,
    ),
    (
        "BGC0000711.gbk",
        72_291,
        "32393648f6a91166444331b83687f1b9b7b24c60553a7ddcb677dfe207736789",
        "BGC0000711",
        30_837,
        21,
    ),
    (
        "BGC0000712.gbk",
        134_734,
        "705104a0daa5c44981b0a1e5352d3e56f012dd2e3ae94c98c85cd0ee9198bf94",
        "BGC0000712",
        48_169,
        40,
    ),
    (
        "BGC0000713.gbk",
        79_197,
        "bf182663de453f4a3fc30ed0aa8f040a164eeab1c98e604983844994996e58fb",
        "BGC0000713",
        31_892,
        26,
    ),
)
EXPECTED_RECORD_IDS = tuple(contract[3] for contract in FIXTURES)
EXPECTED_CDS_COUNTS = {
    contract[3]: contract[5]
    for contract in FIXTURES
}
EXPECTED_ADJACENT_PAIRS = {
    ("BGC0000708", "BGC0000709"),
    ("BGC0000709", "BGC0000711"),
    ("BGC0000711", "BGC0000712"),
    ("BGC0000712", "BGC0000713"),
}
EXPECTED_FIRST_RECORD_LABELS = {
    "livA",
    "livB",
    "livC",
    "livD",
    "livE",
    "livF",
    "livG",
    "livH",
    "livI",
    "livK",
    "livL",
    "livM",
    "livN",
    "livO",
    "livP",
    "livQ",
    "livS",
    "livT",
    "livU",
    "livV",
    "livW",
    "livX",
    "livY",
    "livZ",
}
HEPATOPLASMATACEAE_RECORD_IDS = {
    "AP027078.1",
    "AP027131.1",
    "AP027133.1",
    "AP027132.1",
    "NZ_CP006932.1",
}
PAGES = {
    "T-GUI-04": (
        REPO_ROOT / "docs" / "TUTORIALS" / "GUI" / "compare-proteins-losatp.md"
    ),
    "H-GUI-07": (
        REPO_ROOT
        / "docs"
        / "HOW_TO"
        / "GUI"
        / "create-protein-similarity-groups.md"
    ),
    "H-GUI-08": (
        REPO_ROOT
        / "docs"
        / "HOW_TO"
        / "GUI"
        / "draw-collinear-protein-blocks.md"
    ),
    "T-GUI-08": (
        REPO_ROOT
        / "docs"
        / "TUTORIALS"
        / "GUI"
        / "compare-proteins-losatp-collinear.md"
    ),
}
SCREENSHOT_NAMES = {
    "T-GUI-04": (
        "01-input-ready.png",
        "02-first-diagram.png",
        "03-losatp-settings.png",
        "04-align-og1.png",
        "05-comparison-result.png",
        "06-match-popup.png",
    ),
    "H-GUI-07": (
        "group-settings.png",
        "group-result.png",
        "group-popup.png",
    ),
    "H-GUI-08": (
        "collinear-settings.png",
        "collinear-result.png",
        "block-popup.png",
    ),
    "T-GUI-08": (
        "01-input-ready.png",
        "02-first-diagram.png",
        "03-collinear-settings.png",
        "04-collinear-result.png",
        "05-block-popup.png",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _chapter(scenario_id: str) -> dict[str, object]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    return next(
        chapter for chapter in manifest["chapters"] if chapter["id"] == scenario_id
    )


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _normalized_text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def _module_literal(path: Path, variable: str) -> object:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == variable for target in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"Missing literal assignment {variable} in {path}")


def test_bgc_capture_uses_five_frozen_whole_native_linear_records() -> None:
    total_cds = 0
    for filename, size, digest, record_id, length, cds_count in FIXTURES:
        path = BGC_ROOT / filename
        assert path.stat().st_size == size
        assert _sha256(path) == digest
        records = list(SeqIO.parse(path, "genbank"))
        assert len(records) == 1
        record = records[0]
        assert record.id == record_id
        assert len(record) == length
        assert str(record.annotations.get("topology", "")).lower() == "linear"
        assert sum(feature.type == "CDS" for feature in record.features) == cds_count
        total_cds += cds_count

    assert total_cds == 155


def test_bgc_fixture_manifest_keeps_native_sources_and_scenario_ownership() -> None:
    manifest = json.loads(FIXTURE_MANIFEST_PATH.read_text(encoding="utf-8"))
    fixture = manifest["fixtures"]["aminoglycoside-bgc-five"]
    semantics = fixture["expectedSemantics"]
    assert semantics["recordCount"] == 5
    assert semantics["recordIds"] == list(EXPECTED_RECORD_IDS)
    assert semantics["displayedFeatureCount"] == 155
    assert semantics["displayedFeatureTypes"] == ["CDS"]
    assert semantics["recordsAreWholeCanonicalSources"] is True
    assert semantics["tutorialScenarioId"] == "T-GUI-04"
    assert semantics["similarityGroupScenarioId"] == "H-GUI-07"
    assert semantics["collinearScenarioId"] == "H-CLI-08"
    assert {"T-GUI-04", "H-GUI-07"} <= set(
        fixture["scenarioIds"]
    )
    assert "H-GUI-08" not in fixture["scenarioIds"]

    files = manifest["files"]
    for filename, size, digest, record_id, length, cds_count in FIXTURES:
        relative_path = f"aminoglycoside-bgc-five/{filename}"
        entry = next(
            value for value in files.values() if value["relativePath"] == relative_path
        )
        assert entry["role"] == "raw"
        assert entry["derivation"] is None
        assert entry["sizeBytes"] == size
        assert entry["sha256"] == digest
        assert {"T-GUI-04", "H-GUI-07"} <= set(
            entry["scenarioIds"]
        )
        assert "H-GUI-08" not in entry["scenarioIds"]
        assert len(entry["records"]) == 1
        record = entry["records"][0]
        assert record["id"] == record_id
        assert record["topology"] == "linear"
        assert record["length"] == length
        assert record["displayedFeatureCount"] == cds_count
        assert record["cdsCount"] == cds_count


def test_hepatoplasmataceae_fixture_owns_all_collinear_tutorial_variants() -> None:
    manifest = json.loads(FIXTURE_MANIFEST_PATH.read_text(encoding="utf-8"))
    fixture = manifest["fixtures"]["hepatoplasmataceae-five"]
    assert set(fixture["scenarioIds"]) == {
        "T-GUI-08",
        "T-CLI-10",
        "T-PY-07",
        "H-GUI-08",
    }
    assert set(fixture["expectedSemantics"]["recordIds"]) == (
        HEPATOPLASMATACEAE_RECORD_IDS
    )
    assert fixture["expectedSemantics"]["recordsAreWholeCanonicalSources"] is True
    semantics = fixture["expectedSemantics"]
    assert semantics["galleryCollinearSearchScope"] == "adjacent"
    assert semantics["gallerySearchedRecordPairCount"] == 4
    assert semantics["allRecordHowToSearchScope"] == "all"
    assert semantics["allRecordHowToPairCount"] == 10

    for file_id in (*fixture["fileIds"], *fixture["fileReferences"]):
        assert {"T-GUI-08", "H-GUI-08"} <= set(
            manifest["files"][file_id]["scenarioIds"]
        )


def test_protein_flows_separate_bgc_groups_from_hepatoplasma_collinear() -> None:
    tutorial = TUTORIAL_FLOW_PATH.read_text(encoding="utf-8")
    how_to = HOW_TO_FLOW_PATH.read_text(encoding="utf-8")
    common = COMMON_FLOW_PATH.read_text(encoding="utf-8")

    assert tutorial.count('mode="orthogroup"') == 1
    assert 'mode="collinear"' not in tutorial
    assert tutorial.count("include_plain_result=True") == 1
    assert tutorial.count("download_member_fasta=False") == 1
    assert tutorial.count("separate_strands=False") == 1
    assert how_to.count('mode="orthogroup"') == 1
    assert 'mode="collinear"' not in how_to
    assert how_to.count("include_plain_result=False") == 1
    assert how_to.count("capture_hepatoplasmataceae_collinear(") == 1
    assert 'output_prefix="hepatoplasmataceae_collinear"' in how_to

    for fragment in (
        'get_by_role("button", name="Linear", exact=True)',
        'page.locator(\'[data-capture="linear-blast-source"]\')',
        'fieldset[data-edge-key=',
        'get_attribute("data-linear-record-uid")',
        'get_by_role("radio", name="GenBank", exact=True)',
        '"radio", name="No comparison", exact=True',
        'get_by_role("button", name="Add sequence", exact=True)',
        'get_by_test_id(f"linear-genbank-{index}").set_input_files',
        'f"Definition for sequence {index}"',
        'f"Subtitle / title for sequence {index}"',
        'f"Reverse complement for sequence {index}"',
        'name="Linear sequence 5", exact=True',
        'to_contain_text("BGC0000713.gbk")',
        '"Reverse complement for sequence 5", exact=True',
        "to_be_in_viewport()",
        '"radio", name="Run LOSAT", exact=True',
        'get_by_role("radio", name="LOSATP", exact=True)',
        'get_by_label("LOSATP blastp mode", exact=True)',
        'name="Raw LOSAT filename for #1 to #2", exact=True',
        'name="Save Raw LOSAT TSV for #1 to #2", exact=True',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
        "download.save_as",
    ):
        assert fragment in common

    presentation = _module_literal(COMMON_FLOW_PATH, "BGC_RECORD_PRESENTATION")
    assert [entry[2] for entry in presentation] == [False, False, False, False, True]

    combined_source = tutorial + how_to + common
    for forbidden in (
        "Region start for sequence",
        "Region end for sequence",
        "region_start",
        "region_end",
        "localStorage",
        "sessionStorage",
        "add_init_script",
        ".runAnalysis(",
        ".downloadSVG(",
        "set_content(",
        "window.__GBDRAW_APP__ =",
    ):
        assert forbidden not in combined_source


def test_bgc_flow_pins_gallery_presentation_and_losatp_values() -> None:
    source = COMMON_FLOW_PATH.read_text(encoding="utf-8")
    for fragment in (
        'show_labels.select_option("first")',
        'page.get_by_label("Priority File (TSV)", exact=True).set_input_files',
        'feature_height.fill("75")',
        'block_stroke_width.fill("2")',
        'line_stroke_width.fill("2")',
        'scale_style.select_option("ruler")',
        'axis_stroke_width.fill("5")',
        'fit_complete_linear_preview(page, target_zoom="40%")',
        'if scenario_id == "H-GUI-07"',
        'expect(result_pair).to_contain_text("Raw result ready")',
        'execution.select_option("serial")',
        'total_threads.select_option("1")',
        'parallel_runs.select_option("1")',
        'page.get_by_label("LOSATP minimum bitscore", exact=True).fill("50")',
        'page.get_by_label("LOSATP maximum e-value", exact=True).fill("0.01")',
        'page.get_by_label("LOSATP minimum identity", exact=True).fill("30")',
        'page.get_by_label("LOSATP minimum alignment length", exact=True).fill("0")',
        'page.get_by_label("Collinear max unit gap", exact=True).fill("1")',
        'page.get_by_label("Collinear minimum block genes", exact=True).fill("2")',
        'page.get_by_label("Collinear evidence scope", exact=True).select_option',
        '"orientation_identity"',
        'first_pair = _linear_pair(page, 1, 2)',
        'name="Raw LOSAT filename for #1 to #2", exact=True',
        'name="Save Raw LOSAT TSV for #1 to #2", exact=True',
    ):
        assert fragment in source

    for forbidden in (
        "raw_names",
        "raw_buttons",
        "enabled_indexes",
        "raw_index",
        "to_have_count(4)",
        'name="Raw LOSAT filename", exact=True',
        'name="Save Raw LOSAT TSV", exact=True',
    ):
        assert forbidden not in source

    priority = WEB_ROOT / "tutorial-data" / "shared" / "cds_gene_qualifier_priority.tsv"
    assert priority.read_bytes() == b"CDS\tgene\n"
    first_record = SeqIO.read(BGC_ROOT / "BGC0000708.gbk", "genbank")
    first_record_gene_labels = {
        feature.qualifiers["gene"][0]
        for feature in first_record.features
        if feature.type == "CDS" and feature.qualifiers.get("gene")
    }
    assert first_record_gene_labels == EXPECTED_FIRST_RECORD_LABELS
    assert len(first_record_gene_labels) == 24


def test_gallery_svg_is_the_frozen_similarity_group_reference() -> None:
    root = ET.parse(GALLERY_SVG_PATH).getroot()
    elements = list(root.iter())
    feature_paths = [
        element
        for element in elements
        if element.get("data-gbdraw-feature-id")
    ]
    assert len(feature_paths) == 155
    assert {_local_name(element.tag) for element in feature_paths} == {"path"}
    assert Counter(
        element.get("data-gbdraw-record-id") for element in feature_paths
    ) == Counter(EXPECTED_CDS_COUNTS)
    assert {element.get("stroke-width") for element in feature_paths} == {"2.0"}
    assert all(
        "-18.75" in element.get("d", "")
        and "18.75" in element.get("d", "")
        for element in feature_paths
    )

    labels = {
        text
        for element in elements
        if _local_name(element.tag) == "text"
        and (text := _normalized_text(element)).startswith("liv")
    }
    assert labels == EXPECTED_FIRST_RECORD_LABELS
    assert len(labels) == 24

    orthogroup_links = [
        element
        for element in elements
        if element.get("data-gbdraw-pairwise-match-id")
        and element.get("data-match-kind") == "orthogroup"
    ]
    assert len(orthogroup_links) == 77
    endpoints = {
        (
            element.get("data-query-record-id"),
            element.get("data-subject-record-id"),
        )
        for element in orthogroup_links
    }
    assert endpoints == EXPECTED_ADJACENT_PAIRS
    assert ("BGC0000708", "BGC0000713") not in endpoints

    metadata = next(
        element
        for element in elements
        if _local_name(element.tag) == "metadata"
        and element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    payload = json.loads(metadata.text or "")
    orthogroups = payload["items"][0]["orthogroups"]
    assert len(orthogroups) == 23
    assert {group["id"] for group in orthogroups} == {
        f"og_{index}" for index in range(1, 24)
    }


def test_hepatoplasmataceae_collinear_guards_pin_evidence_and_span_fasta() -> None:
    source = HEPATOPLASMA_FLOW_PATH.read_text(encoding="utf-8")
    wrapper = HOW_TO_FLOW_PATH.read_text(encoding="utf-8")
    for fragment in (
        'page.locator(\'[data-capture="linear-blast-source"]\')',
        'page.get_by_label("Collinear evidence scope", exact=True).select_option(',
        'evidence_scope = "all" if starting_project == "all-vs-all" else "adjacent"',
        "_assert_gallery_adjacent_search(page, session_keys=session_cache_keys)",
        'page.get_by_label("Track Layout", exact=True).select_option("middle")',
        '"Separate Strands"',
        'page.get_by_label("Collinear color mode", exact=True).select_option',
        '"orientation_identity"',
        'telemetry.get("totalPairs") != 25',
        'telemetry.get("cacheHits") != 25',
        'telemetry.get("cacheMisses") != 0',
        'telemetry.get("uniqueJobs") != 0',
        'if endpoints != ADJACENT_PAIRS',
        'download_dir, "collinear_members.fasta"',
        'if len(records) != 2',
        'sequence not in source_sequences[source_id]',
        'capture_hepatoplasmataceae_collinear(',
        'output_prefix="hepatoplasmataceae_collinear"',
    ):
        assert fragment in source + wrapper
    assert (
        'get_by_role("radio", name="Run LOSAT", exact=True).first'
        not in source
    )


def test_protein_comparison_pages_and_manifest_state_the_complete_recipe() -> None:
    group_assertions = {
        "T-GUI-04": {
            "whole_linear_record_count=5",
            "displayed_features=155",
            "feature_height=75",
            "feature_stroke_width=2",
            "labels=first_record_gene",
            "separate_strands=false",
            "group_count=23",
            "svg.group_links=77",
            "adjacent_endpoints_only=true",
            "direct_0708_0713_edge=false",
        },
        "H-GUI-07": {
            "whole_linear_record_count=5",
            "displayed_features=155",
            "feature_height=75",
            "feature_stroke_width=2",
            "labels=first_record_gene",
            "separate_strands=false",
            "group_count=23",
            "svg.group_links=77",
            "adjacent_endpoints_only=true",
            "direct_0708_0713_edge=false",
        },
    }
    for scenario_id in ("T-GUI-04", "H-GUI-07"):
        page_path = PAGES[scenario_id]
        chapter = _chapter(scenario_id)
        page = page_path.read_text(encoding="utf-8")
        assert chapter["destination"] == str(page_path.relative_to(REPO_ROOT))
        assert chapter["fixtures"] == ["aminoglycoside-bgc-five"]
        assert chapter["settings"]["mode"] == "linear"
        assert chapter["settings"]["program"] == "losatp"
        assert chapter["settings"]["scheduling"] == "serial"
        assert chapter["settings"]["threads"] == 1
        assert chapter["settings"]["separate_strands"] is False
        assert group_assertions[scenario_id] <= set(
            chapter["execution"]["assertions"]
        )
        assert tuple(
            Path(screenshot["path"]).name for screenshot in chapter["screenshots"]
        ) == SCREENSHOT_NAMES[scenario_id]
        assert chapter["status"] == {
            "implementation": "verified",
            "review": "approved",
        }
        for record_id in EXPECTED_RECORD_IDS:
            assert record_id in page
        for value in (
            "Feature Height",
            "`75`",
            "Block Stroke Width",
            "Line Stroke Width",
            "Show Coordinate Scale",
            "Ruler",
            "Axis Stroke Width",
            "`5`",
            "40%",
            "First record",
            "`gene`",
        ):
            assert value in page

    tutorial = PAGES["T-GUI-04"].read_text(encoding="utf-8")
    groups = PAGES["H-GUI-07"].read_text(encoding="utf-8")
    assert "leave all optional Region fields blank" in tutorial
    assert "Reverse complement** only for `BGC0000713`" in tutorial
    assert "| Separate Strands | Off |" in tutorial
    assert "23 stable\ngroups and 77 displayed group links" in tutorial
    assert "comparison between sequence 1 and sequence 2" in tutorial
    assert "Raw LOSAT results" not in tutorial
    assert "without\n   setting regions" in groups
    assert "Turn off **Separate Strands**" in groups
    assert _chapter("T-GUI-04")["settings"]["separate_strands"] is False
    assert _chapter("H-GUI-07")["settings"]["separate_strands"] is False
    assert "separate_strands=False" in TUTORIAL_FLOW_PATH.read_text(encoding="utf-8")
    assert "separate_strands=False" in HOW_TO_FLOW_PATH.read_text(encoding="utf-8")
    assert "23\ngroups, and 77 links" in groups

    gallery_collinear = PAGES["T-GUI-08"].read_text(encoding="utf-8")
    gallery_chapter = _chapter("T-GUI-08")
    assert gallery_chapter["settings"]["comparison_scope"] == "adjacent"
    assert tuple(
        Path(screenshot["path"]).name
        for screenshot in gallery_chapter["screenshots"]
    ) == SCREENSHOT_NAMES["T-GUI-08"]
    assert "no Similarity-groups result is mixed" in gallery_collinear
    assert "adjacent-pair Collinear" in gallery_collinear

    collinear_path = PAGES["H-GUI-08"]
    collinear = collinear_path.read_text(encoding="utf-8")
    chapter = _chapter("H-GUI-08")
    assert chapter["destination"] == str(collinear_path.relative_to(REPO_ROOT))
    assert chapter["fixtures"] == ["hepatoplasmataceae-five"]
    assert chapter["settings"] == {
        "mode": "linear",
        "program": "losatp",
        "protein_mode": "collinear",
        "scheduling": "threaded",
        "threads": 32,
        "evidence_scope": "all",
        "track_preset": "middle",
        "separate_strands": True,
    }
    assert {
        "whole_record_count=5",
        "displayed_features=2994",
        "all_source_topologies=circular",
        "all_records_are_complete=true",
        "evidence_scope=all",
        "cache_hits=25",
        "cache_misses=0",
        "worker_jobs=0",
        "adjacent_display_pairs=true",
        "orientation_match=true",
        "color_mode=orientation_identity",
        "span_fasta_records=2",
    } <= set(chapter["execution"]["assertions"])
    assert chapter["execution"]["expected_outputs"] == [
        "hepatoplasmataceae_collinear.svg",
        "collinear_members.fasta",
    ]
    assert tuple(
        Path(screenshot["path"]).name for screenshot in chapter["screenshots"]
    ) == SCREENSHOT_NAMES["H-GUI-08"]
    for record_id in HEPATOPLASMATACEAE_RECORD_IDS:
        assert record_id in collinear
    for value in (
        "Hepatoplasmataceae",
        "All records",
        "25 directional and self LOSATP results",
        "Middle",
        "Separate Strands",
        "Ajisai",
        "40%",
        "hepatoplasmataceae_collinear.svg",
        "contain two non-empty nucleotide envelope sequences",
    ):
        assert value in collinear
    for bgc_id in EXPECTED_RECORD_IDS:
        assert bgc_id not in collinear


def test_protein_comparison_screenshots_are_full_pinned_viewports() -> None:
    for scenario_id, names in SCREENSHOT_NAMES.items():
        root = REPO_ROOT / "docs" / "images" / scenario_id.lower()
        for name in names:
            path = root / name
            assert path.is_file(), path
            with Image.open(path) as image:
                assert image.size == (1440, 900)
                assert image.mode == "RGB"


def test_protein_popup_state_uses_catalog_commit_path() -> None:
    app_setup = APP_SETUP_PATH.read_text(encoding="utf-8")
    run_analysis = RUN_ANALYSIS_PATH.read_text(encoding="utf-8")
    svg_actions = SVG_ACTIONS_PATH.read_text(encoding="utf-8")

    for fragment in (
        "const setLinearComparisonLosatFilename = (id, value) => {",
        "const deactivateLinearComparisonLosatFilename = (id) => {",
        "const updateResolvedLosatFilenameDraft = (edgeKey, updater) => {",
    ):
        start = app_setup.index(fragment)
        body = app_setup[start : start + 1_600]
        assert "replaceLinearComparisonPlan(next, { invalidate: false });" in body

    for fragment in (
        "collinearGroups,",
        "collinearGroups.value = [];",
        "kind: 'collinearityResult'",
        "typedResource: convertedPayload.collinearityResult",
        "setOrthogroupMetadata(candidateCommit.featureState.orthogroups);",
        "collinearGroups.value = Array.isArray(candidateCommit.featureState.collinearGroups)",
    ):
        assert fragment in run_analysis
    assert "Array.isArray(convertedPayload.collinearGroups)" not in run_analysis
    assert "const committedOrthogroups" not in run_analysis

    for fragment in (
        "collinearGroups,",
        "...(Array.isArray(orthogroups?.value) ? orthogroups.value : []),",
        "...(Array.isArray(collinearGroups?.value) ? collinearGroups.value : [])",
        "const payload = buildMatchPayload(matchElement, featureLookup);",
    ):
        assert fragment in svg_actions
