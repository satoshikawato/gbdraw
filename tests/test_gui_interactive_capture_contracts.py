from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from Bio import SeqIO
from PIL import Image

from docs.capture.config import chapter_for
from gbdraw.session_io import CURRENT_SESSION_VERSION

pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
CAPTURE_ROOT = REPO_ROOT / "docs" / "capture"
INTERACTIVE_FLOW = (
    CAPTURE_ROOT / "flows" / "how_to" / "interactive_sessions.py"
)
EXPORT_FLOW = CAPTURE_ROOT / "flows" / "how_to" / "exports.py"
HUMAN_HELPER = CAPTURE_ROOT / "flows" / "human_circular.py"
INDEX_PATH = REPO_ROOT / "gbdraw" / "web" / "index.html"
HUMAN_FIXTURE = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "human-mitochondrion"
    / "HmmtDNA.gbk"
)
BGC_ROOT = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "aminoglycoside-bgc-five"
)
SCREENSHOT_NAMES = {
    "H-GUI-13": (
        "search-popup.png",
        "editor.png",
        "edited-result.png",
        "match-popup.png",
        "group-popup.png",
    ),
    "H-GUI-14": (
        "history-actions.png",
        "session-download.png",
        "reloaded-result.png",
    ),
    "H-GUI-15": (
        "export-actions.png",
        "exported-result.png",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_interactive_flows_use_whole_natural_source_records() -> None:
    assert HUMAN_FIXTURE.stat().st_size == 64_640
    assert _sha256(HUMAN_FIXTURE) == (
        "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f"
    )
    human_records = list(SeqIO.parse(HUMAN_FIXTURE, "genbank"))
    assert len(human_records) == 1
    human = human_records[0]
    assert human.id == "NC_012920.1"
    assert len(human) == 16_569
    assert human.annotations["topology"] == "circular"
    assert "complete" in human.description.lower()

    expected_bgcs = {
        "BGC0000708": 40_579,
        "BGC0000709": 50_466,
        "BGC0000711": 30_837,
        "BGC0000712": 48_169,
        "BGC0000713": 31_892,
    }
    for record_id, length in expected_bgcs.items():
        records = list(SeqIO.parse(BGC_ROOT / f"{record_id}.gbk", "genbank"))
        assert len(records) == 1
        assert records[0].id == record_id
        assert len(records[0]) == length
        assert records[0].annotations["topology"] == "linear"


def test_interactive_flows_use_raw_inputs_and_visible_actions() -> None:
    source = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (HUMAN_HELPER, INTERACTIVE_FLOW, EXPORT_FLOW)
    )
    for required in (
        'get_by_role("button", name="Circular", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'get_by_label("GenBank/DDBJ File", exact=True).set_input_files',
        'get_by_role("searchbox", name="Search features", exact=True)',
        'name="Open active feature", exact=True',
        'name="Apply color and caption to selected features"',
        'name="Apply visibility to selected features", exact=True',
        'name="Apply stroke to selected features", exact=True',
        'name="Toggle layout edit mode", exact=True',
        'name="Save Session", exact=True',
        'name="Load Session", exact=True',
        'name="Reset Settings", exact=True',
        'name="Undo", exact=True',
        'name="Redo", exact=True',
        'button_name="Interactive SVG"',
        'button_name="PNG"',
        'button_name="PDF"',
        "page.expect_download",
        "page.expect_file_chooser",
    ):
        assert required in source

    for forbidden in (
        "gallery/sessions",
        "window.__GBDRAW_APP__ =",
        ".runAnalysis(",
        ".downloadSVG(",
        ".downloadInteractiveSVG(",
        ".downloadPNG(",
        ".downloadPDF(",
        "add_init_script",
        "set_content(",
        "localStorage",
        "sessionStorage",
        "region_start",
        "region_end",
    ):
        assert forbidden not in source


def test_h_gui_13_persists_real_edits_and_keeps_group_semantics_distinct() -> None:
    source = INTERACTIVE_FLOW.read_text(encoding="utf-8")
    for value in (
        'EDIT_TARGET_QUERY = "COX"',
        'EDIT_FILL = "#d81b60"',
        'EDIT_STROKE = "#5e35b1"',
        "EDIT_STROKE_WIDTH = 2.5",
        'EDIT_LEGEND = "Oxidative phosphorylation"',
        "_select_all_search_matches(page, 3)",
        "_assert_target_edit(target_style)",
        "The legend drag did not reach the committed SVG",
        "_inspect_downloaded_edit(svg_path, feature_ids)",
        'mode="orthogroup"',
        "assert_gui_bgc_similarity_groups_svg",
        '[data-match-kind=\"orthogroup\"][data-orthogroup-id=\"og_1\"]',
        '"Similarity group ID", "Members", "og_"',
        'query=r"^og_1$"',
        'field="orthogroup"',
        "use_regex=True",
        "expected_matches=5",
        '"Similarity group", "og_1", "Members", "Record coverage"',
        '"23 similarity groups", "og_1", "5 members", "5 records"',
        '{"text": group_text, "groupId": "og_1", "members": 5}',
        "_assert_closed_right_drawer_contract(page)",
        "_assert_closed_drawer_does_not_reserve_canvas(page)",
        "for width in (1280, VIEWPORT_WIDTH)",
        'drawer_state.get("visibility") != "hidden"',
        'drawer_state.get("pointerEvents") != "none"',
        'toggle_state.get("ariaExpanded") != "false"',
        'final_report["closedDrawerGeometry"]',
        "_reset_finished_preview_viewport(page, target_zoom=60)",
    ):
        assert value in source
    assert 'mode="collinear"' not in source


def test_h_gui_14_writes_current_gzip_session_and_uses_a_fresh_context() -> None:
    source = INTERACTIVE_FLOW.read_text(encoding="utf-8")
    assert CURRENT_SESSION_VERSION == 41
    for value in (
        "CURRENT_SESSION_VERSION = 41",
        "CURRENT_RENDER_REQUEST_SCHEMA = 6",
        'SESSION_FILENAME = f"{SESSION_TITLE}.gbdraw-session.json.gz"',
        'contents[:2] != b"\\x1f\\x8b"',
        'with gzip.open(path, "rt", encoding="utf-8")',
        'payload.get("format") != "gbdraw-session"',
        "base64.b64decode",
        "FIRST_CIRCULAR_FIXTURE_SHA256",
        "first_capture.close()",
        "second_capture = open_browser_capture",
        'dialog.message != "Session loaded successfully!"',
        "_assert_semantically_equivalent(source_report, restored_report)",
        "_assert_semantically_equivalent(source_report, exported)",
        "_stabilize_static_capture_surface(page)",
        "backdrop-filter: none !important",
        "_reset_finished_preview_viewport(page, target_zoom=60)",
        'final_report["restoredPreviewFrame"]',
        "_frame_finished_preview_with_legend(page)",
        'error: \'record and legend do not fit\'',
    ):
        assert value in source


def test_h_gui_15_validates_every_actual_export() -> None:
    source = EXPORT_FLOW.read_text(encoding="utf-8")
    for value in (
        'STATIC_SVG_NAME = f"{OUTPUT_PREFIX}.svg"',
        'INTERACTIVE_SVG_NAME = f"{OUTPUT_PREFIX}.interactive.svg"',
        'PNG_NAME = f"{OUTPUT_PREFIX}.png"',
        'PDF_NAME = f"{OUTPUT_PREFIX}.pdf"',
        "PNG_MAGIC",
        "assert_finished_circular_svg(report)",
        'root.attrib.get("data-gbdraw-interactive-svg") != "true"',
        "INTERACTIVE_ASSET_IDS.issubset(ids)",
        'payload.get("schema") != 3',
        'page.goto(path.resolve().as_uri(), wait_until="load")',
        'name="Expand feature search", exact=True',
        'fill("COX1")',
        'not contents.startswith(b"%PDF-")',
        'b"%%EOF"',
        'rb"/MediaBox',
        "PNG_DPI / 96",
    ):
        assert value in source


def test_interactive_accessibility_labels_are_public_ui_contracts() -> None:
    source = INDEX_PATH.read_text(encoding="utf-8") + (
        REPO_ROOT / "gbdraw" / "web" / "js" / "components.js"
    ).read_text(encoding="utf-8")
    for label in (
        ':aria-label="`Feature details: ${clickedFeature.label}`"',
        'aria-label="Close feature popup"',
        'aria-label="Feature legend name"',
        'aria-label="Feature visibility"',
        'aria-label="Feature stroke"',
        ":aria-label=\"ariaLabel + ' width'\"",
        'aria-label="Feature search status"',
        'aria-label="Selected feature count"',
        'aria-label="PNG DPI"',
        '.right-drawer[aria-hidden="true"]',
        "transform: translateX(100%);",
        "visibility: hidden;",
        "pointer-events: none;",
        ':aria-hidden="!showRightDrawer"',
        '.drawer-toggle[aria-expanded="true"]',
        ':aria-expanded="showRightDrawer"',
    ):
        assert label in source


def test_manifest_runner_and_screenshot_evidence_contracts_match() -> None:
    from docs.capture import config as capture_config
    from docs.capture import run_all

    for scenario_id, names in SCREENSHOT_NAMES.items():
        chapter = chapter_for(scenario_id)
        assert chapter["status"] == {
            "implementation": "verified",
            "review": "approved",
        }
        assert chapter["role"] == "evidence"
        assert "destination" not in chapter
        assert tuple(
            Path(screenshot["path"]).name for screenshot in chapter["screenshots"]
        ) == names
        assert capture_config.screenshot_names_for(scenario_id) == names
        assert scenario_id in run_all.CAPTURE_FUNCTIONS



def test_interactive_screenshots_are_deterministic_pngs() -> None:
    for scenario_id, names in SCREENSHOT_NAMES.items():
        for name in names:
            path = REPO_ROOT / "docs" / "images" / scenario_id.lower() / name
            assert path.is_file(), path
            with Image.open(path) as image:
                assert image.format == "PNG"
                assert image.size == (1440, 900)
                assert image.mode == "RGB"
            assert path.stat().st_size <= 2_500_000
