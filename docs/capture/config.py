"""Authoritative environment contract for documentation screenshots."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit


CAPTURE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = CAPTURE_ROOT.parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
MANIFEST_PATH = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
VENDORED_FONT_ROOT = WEB_ROOT / "vendor" / "fonts"

LOCAL_SCHEME = "http"
LOCAL_HOST = "127.0.0.1"
LOCAL_BASE_URL_TEMPLATE = f"{LOCAL_SCHEME}://{LOCAL_HOST}:{{port}}/"
ALLOWED_CAPTURE_HOSTS = frozenset({"127.0.0.1", "localhost"})
ISOLATION_HEADERS = {
    "Cross-Origin-Opener-Policy": "same-origin",
    "Cross-Origin-Embedder-Policy": "require-corp",
    "Cross-Origin-Resource-Policy": "same-origin",
    "Cache-Control": "no-store",
}

PYTHON_PLAYWRIGHT_VERSION = "1.61.0"
NODE_PLAYWRIGHT_VERSION = "1.61.1"
CHROMIUM_VERSION = "149.0.7827.55"
CHROMIUM_REVISION = "v1228"

VIEWPORT_WIDTH = 1440
VIEWPORT_HEIGHT = 900
DEVICE_SCALE_FACTOR = 1
LOCALE = "en-US"
TIMEZONE_ID = "UTC"
COLOR_SCHEME = "light"
REDUCED_MOTION = "reduce"

ACTION_TIMEOUT_MS = 30_000
NAVIGATION_TIMEOUT_MS = 30_000
WORKER_READY_TIMEOUT_MS = 180_000
GENERATION_TIMEOUT_MS = 180_000
SCREENSHOT_MAX_BYTES = 2_500_000

SUPPORTED_TIERS = ("core", "extended", "nightly")
FIRST_CIRCULAR_SCENARIO_ID = "T-GUI-01"
FIRST_LINEAR_SCENARIO_ID = "T-GUI-02"
GUI_INPUTS_SCENARIO_ID = "H-GUI-01"
GUI_LOSATN_SCENARIO_ID = "T-GUI-03"
GUI_CIRCULAR_LAYOUT_SCENARIO_ID = "H-GUI-02"
GUI_LINEAR_LAYOUT_SCENARIO_ID = "H-GUI-03"
GUI_UPLOADED_COMPARISON_SCENARIO_ID = "H-GUI-04"
GUI_TLOSATX_SCENARIO_ID = "H-GUI-05"
GUI_CIRCULAR_RINGS_SCENARIO_ID = "H-GUI-06"
GUI_LOSATP_GROUPS_SCENARIO_ID = "T-GUI-04"
GUI_ANNOTATED_CHLOROPLAST_SCENARIO_ID = "T-GUI-05"
GUI_PRECOMPUTED_CIRCULAR_RINGS_SCENARIO_ID = "T-GUI-06"
GUI_LOSATP_TUTORIAL_COLLINEAR_SCENARIO_ID = "T-GUI-08"
GUI_INTERACTIVE_HANDOFF_SCENARIO_ID = "T-GUI-09"
GUI_FEATURE_HIGHLIGHT_SCENARIO_ID = "T-GUI-10"
GUI_QUANTITATIVE_MAP_SCENARIO_ID = "T-GUI-12"
GUI_LOSATP_GROUPS_HOW_TO_SCENARIO_ID = "H-GUI-07"
GUI_LOSATP_COLLINEAR_SCENARIO_ID = "H-GUI-08"
GUI_QUANTITATIVE_TRACKS_SCENARIO_ID = "H-GUI-09"
GUI_ANNOTATION_TRACKS_SCENARIO_ID = "H-GUI-10"
GUI_STYLING_SCENARIO_ID = "H-GUI-11"
GUI_FEATURE_PRESENTATION_SCENARIO_ID = "H-GUI-12"
GUI_INTERACTIVE_EDITING_SCENARIO_ID = "H-GUI-13"
GUI_SESSION_REPRODUCTION_SCENARIO_ID = "H-GUI-14"
GUI_EXPORTS_SCENARIO_ID = "H-GUI-15"
IMPLEMENTED_SCENARIO_IDS = (
    FIRST_CIRCULAR_SCENARIO_ID,
    FIRST_LINEAR_SCENARIO_ID,
    GUI_INPUTS_SCENARIO_ID,
    GUI_LOSATN_SCENARIO_ID,
    GUI_CIRCULAR_LAYOUT_SCENARIO_ID,
    GUI_LINEAR_LAYOUT_SCENARIO_ID,
    GUI_UPLOADED_COMPARISON_SCENARIO_ID,
    GUI_TLOSATX_SCENARIO_ID,
    GUI_CIRCULAR_RINGS_SCENARIO_ID,
    GUI_LOSATP_GROUPS_SCENARIO_ID,
    GUI_ANNOTATED_CHLOROPLAST_SCENARIO_ID,
    GUI_PRECOMPUTED_CIRCULAR_RINGS_SCENARIO_ID,
    GUI_LOSATP_TUTORIAL_COLLINEAR_SCENARIO_ID,
    GUI_INTERACTIVE_HANDOFF_SCENARIO_ID,
    GUI_FEATURE_HIGHLIGHT_SCENARIO_ID,
    GUI_QUANTITATIVE_MAP_SCENARIO_ID,
    GUI_LOSATP_GROUPS_HOW_TO_SCENARIO_ID,
    GUI_LOSATP_COLLINEAR_SCENARIO_ID,
    GUI_QUANTITATIVE_TRACKS_SCENARIO_ID,
    GUI_ANNOTATION_TRACKS_SCENARIO_ID,
    GUI_STYLING_SCENARIO_ID,
    GUI_FEATURE_PRESENTATION_SCENARIO_ID,
    GUI_INTERACTIVE_EDITING_SCENARIO_ID,
    GUI_SESSION_REPRODUCTION_SCENARIO_ID,
    GUI_EXPORTS_SCENARIO_ID,
)
FIRST_CIRCULAR_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-publication-label.png",
    "04-layout-settings.png",
    "04-finished-diagram.png",
    "05-export-svg.png",
)
FIRST_LINEAR_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-layout-settings.png",
    "04-finished-diagram.png",
)
GUI_INPUTS_SCREENSHOT_NAMES = (
    "genbank-input.png",
    "gff3-fasta-input.png",
    "id-error.png",
)
GUI_LOSATN_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-losatn-settings.png",
    "04-comparison-result.png",
    "05-match-popup.png",
)
GUI_CIRCULAR_LAYOUT_SCREENSHOT_NAMES = (
    "grid-settings.png",
    "grid-result.png",
)
GUI_LINEAR_LAYOUT_SCREENSHOT_NAMES = (
    "record-layout.png",
    "orientation-result.png",
)
GUI_UPLOADED_COMPARISON_SCREENSHOT_NAMES = (
    "comparison-plan.png",
    "comparison-result.png",
)
GUI_TLOSATX_SCREENSHOT_NAMES = (
    "tlosatx-settings.png",
    "tlosatx-result.png",
)
GUI_CIRCULAR_RINGS_SCREENSHOT_NAMES = (
    "ring-settings.png",
    "ring-result.png",
    "hsp-popup.png",
)
GUI_LOSATP_GROUPS_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-losatp-settings.png",
    "04-align-og1.png",
    "05-comparison-result.png",
    "06-match-popup.png",
)
GUI_ANNOTATED_CHLOROPLAST_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-annotation-table.png",
    "04-track-settings.png",
    "05-finished-diagram.png",
)
GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-ring-settings.png",
    "04-ring-result.png",
    "05-hsp-popup.png",
)
GUI_LOSATP_TUTORIAL_COLLINEAR_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-collinear-settings.png",
    "04-collinear-result.png",
    "04-collinear-detail.png",
    "05-block-popup.png",
)
GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-interactive-export.png",
    "04-feature-search.png",
    "05-session-download.png",
    "06-reloaded-result.png",
)
GUI_FEATURE_HIGHLIGHT_SCREENSHOT_NAMES = (
    "presentation-settings.png",
    "presentation-result.png",
)
GUI_QUANTITATIVE_MAP_SCREENSHOT_NAMES = (
    "track-settings.png",
    "track-result.png",
)
GUI_LOSATP_GROUPS_HOW_TO_SCREENSHOT_NAMES = (
    "group-settings.png",
    "group-result.png",
    "group-popup.png",
)
GUI_LOSATP_COLLINEAR_SCREENSHOT_NAMES = (
    "collinear-settings.png",
    "collinear-result.png",
    "collinear-detail.png",
    "block-popup.png",
)
GUI_QUANTITATIVE_TRACKS_SCREENSHOT_NAMES = (
    "track-settings.png",
    "track-result.png",
)
GUI_ANNOTATION_TRACKS_SCREENSHOT_NAMES = (
    "slot-settings.png",
    "annotation-result.png",
)
GUI_STYLING_SCREENSHOT_NAMES = (
    "style-settings.png",
    "style-result.png",
)
GUI_FEATURE_PRESENTATION_SCREENSHOT_NAMES = (
    "presentation-settings.png",
    "presentation-result.png",
)
GUI_INTERACTIVE_EDITING_SCREENSHOT_NAMES = (
    "search-popup.png",
    "editor.png",
    "edited-result.png",
    "match-popup.png",
    "group-popup.png",
)
GUI_SESSION_REPRODUCTION_SCREENSHOT_NAMES = (
    "history-actions.png",
    "session-download.png",
    "reloaded-result.png",
)
GUI_EXPORTS_SCREENSHOT_NAMES = (
    "export-actions.png",
    "exported-result.png",
)
SCREENSHOT_NAMES_BY_SCENARIO = {
    FIRST_CIRCULAR_SCENARIO_ID: FIRST_CIRCULAR_SCREENSHOT_NAMES,
    FIRST_LINEAR_SCENARIO_ID: FIRST_LINEAR_SCREENSHOT_NAMES,
    GUI_INPUTS_SCENARIO_ID: GUI_INPUTS_SCREENSHOT_NAMES,
    GUI_LOSATN_SCENARIO_ID: GUI_LOSATN_SCREENSHOT_NAMES,
    GUI_CIRCULAR_LAYOUT_SCENARIO_ID: GUI_CIRCULAR_LAYOUT_SCREENSHOT_NAMES,
    GUI_LINEAR_LAYOUT_SCENARIO_ID: GUI_LINEAR_LAYOUT_SCREENSHOT_NAMES,
    GUI_UPLOADED_COMPARISON_SCENARIO_ID: (
        GUI_UPLOADED_COMPARISON_SCREENSHOT_NAMES
    ),
    GUI_TLOSATX_SCENARIO_ID: GUI_TLOSATX_SCREENSHOT_NAMES,
    GUI_CIRCULAR_RINGS_SCENARIO_ID: GUI_CIRCULAR_RINGS_SCREENSHOT_NAMES,
    GUI_LOSATP_GROUPS_SCENARIO_ID: GUI_LOSATP_GROUPS_SCREENSHOT_NAMES,
    GUI_ANNOTATED_CHLOROPLAST_SCENARIO_ID: (
        GUI_ANNOTATED_CHLOROPLAST_SCREENSHOT_NAMES
    ),
    GUI_PRECOMPUTED_CIRCULAR_RINGS_SCENARIO_ID: (
        GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES
    ),
    GUI_LOSATP_TUTORIAL_COLLINEAR_SCENARIO_ID: (
        GUI_LOSATP_TUTORIAL_COLLINEAR_SCREENSHOT_NAMES
    ),
    GUI_INTERACTIVE_HANDOFF_SCENARIO_ID: GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES,
    GUI_FEATURE_HIGHLIGHT_SCENARIO_ID: GUI_FEATURE_HIGHLIGHT_SCREENSHOT_NAMES,
    GUI_QUANTITATIVE_MAP_SCENARIO_ID: GUI_QUANTITATIVE_MAP_SCREENSHOT_NAMES,
    GUI_LOSATP_GROUPS_HOW_TO_SCENARIO_ID: (
        GUI_LOSATP_GROUPS_HOW_TO_SCREENSHOT_NAMES
    ),
    GUI_LOSATP_COLLINEAR_SCENARIO_ID: GUI_LOSATP_COLLINEAR_SCREENSHOT_NAMES,
    GUI_QUANTITATIVE_TRACKS_SCENARIO_ID: (
        GUI_QUANTITATIVE_TRACKS_SCREENSHOT_NAMES
    ),
    GUI_ANNOTATION_TRACKS_SCENARIO_ID: GUI_ANNOTATION_TRACKS_SCREENSHOT_NAMES,
    GUI_STYLING_SCENARIO_ID: GUI_STYLING_SCREENSHOT_NAMES,
    GUI_FEATURE_PRESENTATION_SCENARIO_ID: (
        GUI_FEATURE_PRESENTATION_SCREENSHOT_NAMES
    ),
    GUI_INTERACTIVE_EDITING_SCENARIO_ID: (
        GUI_INTERACTIVE_EDITING_SCREENSHOT_NAMES
    ),
    GUI_SESSION_REPRODUCTION_SCENARIO_ID: (
        GUI_SESSION_REPRODUCTION_SCREENSHOT_NAMES
    ),
    GUI_EXPORTS_SCENARIO_ID: GUI_EXPORTS_SCREENSHOT_NAMES,
}

FIRST_CIRCULAR_FIXTURE_PATH = (
    WEB_ROOT / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
)
FIRST_CIRCULAR_FIXTURE_SIZE = 64_640
FIRST_CIRCULAR_FIXTURE_SHA256 = (
    "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f"
)
FIRST_LINEAR_FIXTURE_PATH = WEB_ROOT / "tutorial-data" / "lambda" / "NC_001416.gb"
FIRST_LINEAR_FIXTURE_SIZE = 176_723
FIRST_LINEAR_FIXTURE_SHA256 = (
    "4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7"
)
FIRST_LINEAR_LABEL_RULE_PATH = (
    WEB_ROOT / "tutorial-data" / "shared" / "cds_gene_qualifier_priority.tsv"
)
FIRST_LINEAR_LABEL_RULE_SIZE = 9
FIRST_LINEAR_LABEL_RULE_SHA256 = (
    "1d3c787c5768b52191cc7e8a325cd00fdc4b1e9b495d37e580c3a69dec21b87a"
)
GUI_INPUTS_GFF3_PATH = WEB_ROOT / "tutorial-data" / "lambda-gff3" / "NC_001416.gff3"
GUI_INPUTS_GFF3_SIZE = 36_794
GUI_INPUTS_GFF3_SHA256 = (
    "d53e05de87933104cd26111bca42006cce9b5e903fb5b187740f963b3a2098cb"
)
GUI_INPUTS_FASTA_PATH = WEB_ROOT / "tutorial-data" / "lambda-gff3" / "NC_001416.fna"
GUI_INPUTS_FASTA_SIZE = 49_253
GUI_INPUTS_FASTA_SHA256 = (
    "80897a7ee6b8aaffbab5442e0daad292592ac74701dbdf35af4b400ae0770ef3"
)
GUI_LOSATN_DE3_FIXTURE_PATH = (
    WEB_ROOT / "tutorial-data" / "de3" / "NC_042057.1.gb"
)
GUI_LOSATN_DE3_FIXTURE_SIZE = 111_686
GUI_LOSATN_DE3_FIXTURE_SHA256 = (
    "288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09"
)
GUI_LOSATN_REFERENCE_TSV_PATH = (
    WEB_ROOT
    / "tutorial-data"
    / "lambda-de3-comparison"
    / "lambda-de3.losatn.tsv"
)
GUI_LOSATN_REFERENCE_TSV_SIZE = 436
GUI_LOSATN_REFERENCE_TSV_SHA256 = (
    "703b0ac749669152a8e6d5fa6fb246cf5973cb1a3e0ca9db2f346ee33628317c"
)
GUI_CIRCULAR_LAYOUT_FIXTURES = (
    (
        FIRST_CIRCULAR_FIXTURE_PATH,
        FIRST_CIRCULAR_FIXTURE_SIZE,
        FIRST_CIRCULAR_FIXTURE_SHA256,
        "NC_012920.1",
        16_569,
        "Homo sapiens",
    ),
    (
        WEB_ROOT
        / "tutorial-data"
        / "metazoan-mitochondria-four"
        / "NC_002333.2.gb",
        55_541,
        "94ab35da6f81abc2595fcd425c23585ed78d9396b5143918d9d1025d8a4d2140",
        "NC_002333.2",
        16_596,
        "Danio rerio",
    ),
    (
        WEB_ROOT
        / "tutorial-data"
        / "metazoan-mitochondria-four"
        / "NC_024511.2.gb",
        72_340,
        "79fa36199682961919c4a11f3a8fc50c9e598e68b867eb25e847bce1aa1c4229",
        "NC_024511.2",
        19_524,
        "Drosophila melanogaster",
    ),
    (
        WEB_ROOT
        / "tutorial-data"
        / "metazoan-mitochondria-four"
        / "NC_001328.1.gb",
        39_227,
        "8de5f7cf3686f493ee5b8068dba39b31c5d02e70d997063f65ed19d0fa859a59",
        "NC_001328.1",
        13_794,
        "Caenorhabditis elegans",
    ),
)
GUI_CIRCULAR_LAYOUT_COMBINED_FILENAME = "complete_metazoan_mitochondria.gb"
GUI_BGC_FIXTURES = (
    (
        WEB_ROOT / "tutorial-data" / "aminoglycoside-bgc-five" / "BGC0000708.gbk",
        105_185,
        "9a5f971c5ed8c406b20574fb50aac567609deb787eb1e8d4635050aa264a04b0",
        "BGC0000708",
        40_579,
    ),
    (
        WEB_ROOT / "tutorial-data" / "aminoglycoside-bgc-five" / "BGC0000709.gbk",
        121_511,
        "4b66b7e4b78d429d12176e1e36d0e48178c562a9d128d4308b38753af9995255",
        "BGC0000709",
        50_466,
    ),
    (
        WEB_ROOT / "tutorial-data" / "aminoglycoside-bgc-five" / "BGC0000711.gbk",
        72_291,
        "32393648f6a91166444331b83687f1b9b7b24c60553a7ddcb677dfe207736789",
        "BGC0000711",
        30_837,
    ),
    (
        WEB_ROOT / "tutorial-data" / "aminoglycoside-bgc-five" / "BGC0000712.gbk",
        134_734,
        "705104a0daa5c44981b0a1e5352d3e56f012dd2e3ae94c98c85cd0ee9198bf94",
        "BGC0000712",
        48_169,
    ),
    (
        WEB_ROOT / "tutorial-data" / "aminoglycoside-bgc-five" / "BGC0000713.gbk",
        79_197,
        "bf182663de453f4a3fc30ed0aa8f040a164eeab1c98e604983844994996e58fb",
        "BGC0000713",
        31_892,
    ),
)
GUI_HEPATOPLASMATACEAE_FIXTURES = (
    (
        WEB_ROOT / "tutorial-data" / "hepatoplasmataceae-five" / "AP027078.gb",
        1_344_275,
        "2b7f3fe01757416ed09a34c95c5d326269519f7323ff272e2950e8d6617a87c7",
        "AP027078.1",
        615_622,
        "Candidatus Tyloplasma litorale",
    ),
    (
        WEB_ROOT / "tutorial-data" / "hepatoplasmataceae-five" / "AP027131.gb",
        1_455_315,
        "4281dfd61d963e1264823ef8edef97f6072f567e443f650a97629bd683fba6da",
        "AP027131.1",
        662_108,
        "Candidatus Hepatoplasma vulgare",
    ),
    (
        WEB_ROOT / "tutorial-data" / "depth-1kb" / "AP027133.gb",
        1_344_094,
        "913af50dd9d37cc2107be5e46484b885c5d586fb414b4b501380fc8f17a659d6",
        "AP027133.1",
        606_194,
        "Candidatus Hepatoplasma scabrum",
    ),
    (
        WEB_ROOT / "tutorial-data" / "hepatoplasmataceae-five" / "AP027132.gb",
        1_419_813,
        "b251675fb9dc1853851204da0dff8b3dc7b46110292798e9a5792582383e3903",
        "AP027132.1",
        643_039,
        "Candidatus Hepatoplasma crinochetorum",
    ),
    (
        WEB_ROOT
        / "tutorial-data"
        / "hepatoplasmataceae-five"
        / "NZ_CP006932.gb",
        1_684_495,
        "c92f2df6b5b1eb911d2569af758677cae91e64af4af67e47b5ec0d2c866361a8",
        "NZ_CP006932.1",
        657_101,
        "Candidatus Hepatoplasma crinochetorum Av",
    ),
)


def load_manifest() -> dict[str, Any]:
    """Load the approved chapter manifest."""

    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def chapter_for(scenario_id: str) -> dict[str, Any]:
    """Return one approved scenario entry."""

    for chapter in load_manifest()["chapters"]:
        if chapter["id"] == scenario_id:
            return chapter
    raise KeyError(f"Unknown documentation scenario: {scenario_id}")


def _repo_path(relative_path: str) -> Path:
    path = (REPO_ROOT / relative_path).resolve()
    if not path.is_relative_to(REPO_ROOT):
        raise ValueError(f"Documentation path escapes the repository: {relative_path}")
    return path


def screenshot_paths_for(scenario_id: str) -> dict[str, Path]:
    """Resolve all images for one implemented scenario from the manifest."""

    expected_names = SCREENSHOT_NAMES_BY_SCENARIO.get(scenario_id)
    if expected_names is None:
        raise KeyError(f"No implemented screenshot contract for {scenario_id}")
    chapter = chapter_for(scenario_id)
    paths = {
        Path(screenshot["path"]).name: _repo_path(screenshot["path"])
        for screenshot in chapter["screenshots"]
    }
    if tuple(paths) != expected_names:
        raise ValueError(
            f"{scenario_id} screenshot names or order differ from the approved capture contract"
        )
    return paths


def first_circular_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-01 images from the approved manifest."""

    return screenshot_paths_for(FIRST_CIRCULAR_SCENARIO_ID)


def first_circular_screenshot_path() -> Path:
    """Resolve the Step 2 image for callers that only need the early result."""

    return first_circular_screenshot_paths()["02-first-diagram.png"]


def first_circular_tutorial_path() -> Path:
    """Resolve the Tutorial destination from the approved manifest."""

    return _repo_path(chapter_for(FIRST_CIRCULAR_SCENARIO_ID)["destination"])


def first_linear_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-02 images from the approved manifest."""

    return screenshot_paths_for(FIRST_LINEAR_SCENARIO_ID)


def first_linear_tutorial_path() -> Path:
    """Resolve the T-GUI-02 Tutorial destination from the approved manifest."""

    return _repo_path(chapter_for(FIRST_LINEAR_SCENARIO_ID)["destination"])


def gui_inputs_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-01 images from the approved manifest."""

    return screenshot_paths_for(GUI_INPUTS_SCENARIO_ID)


def gui_inputs_how_to_path() -> Path:
    """Resolve the H-GUI-01 How-to destination from the approved manifest."""

    return _repo_path(chapter_for(GUI_INPUTS_SCENARIO_ID)["destination"])


def gui_losatn_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-03 images from the approved manifest."""

    return screenshot_paths_for(GUI_LOSATN_SCENARIO_ID)


def gui_losatn_tutorial_path() -> Path:
    """Resolve the T-GUI-03 Tutorial destination from the approved manifest."""

    return _repo_path(chapter_for(GUI_LOSATN_SCENARIO_ID)["destination"])


def gui_circular_layout_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-02 images from the approved manifest."""

    return screenshot_paths_for(GUI_CIRCULAR_LAYOUT_SCENARIO_ID)


def gui_circular_layout_how_to_path() -> Path:
    """Resolve the H-GUI-02 How-to destination from the approved manifest."""

    return _repo_path(chapter_for(GUI_CIRCULAR_LAYOUT_SCENARIO_ID)["destination"])


def gui_linear_layout_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-03 images from the approved manifest."""

    return screenshot_paths_for(GUI_LINEAR_LAYOUT_SCENARIO_ID)


def gui_uploaded_comparison_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-04 images from the approved manifest."""

    return screenshot_paths_for(GUI_UPLOADED_COMPARISON_SCENARIO_ID)


def gui_tlosatx_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-05 images from the approved manifest."""

    return screenshot_paths_for(GUI_TLOSATX_SCENARIO_ID)


def gui_circular_rings_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-06 images from the approved manifest."""

    return screenshot_paths_for(GUI_CIRCULAR_RINGS_SCENARIO_ID)


def gui_losatp_groups_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-04 images from the approved manifest."""

    return screenshot_paths_for(GUI_LOSATP_GROUPS_SCENARIO_ID)


def gui_annotated_chloroplast_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-05 images from the approved manifest."""

    return screenshot_paths_for(GUI_ANNOTATED_CHLOROPLAST_SCENARIO_ID)


def gui_precomputed_circular_rings_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-06 images from the approved manifest."""

    return screenshot_paths_for(GUI_PRECOMPUTED_CIRCULAR_RINGS_SCENARIO_ID)


def gui_losatp_tutorial_collinear_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-08 images from the approved manifest."""

    return screenshot_paths_for(GUI_LOSATP_TUTORIAL_COLLINEAR_SCENARIO_ID)


def gui_interactive_handoff_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-09 images from the approved manifest."""

    return screenshot_paths_for(GUI_INTERACTIVE_HANDOFF_SCENARIO_ID)


def gui_feature_highlight_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-10 images from the approved manifest."""

    return screenshot_paths_for(GUI_FEATURE_HIGHLIGHT_SCENARIO_ID)


def gui_quantitative_map_screenshot_paths() -> dict[str, Path]:
    """Resolve all T-GUI-12 images from the approved manifest."""

    return screenshot_paths_for(GUI_QUANTITATIVE_MAP_SCENARIO_ID)


def gui_losatp_groups_how_to_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-07 images from the approved manifest."""

    return screenshot_paths_for(GUI_LOSATP_GROUPS_HOW_TO_SCENARIO_ID)


def gui_losatp_collinear_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-08 images from the approved manifest."""

    return screenshot_paths_for(GUI_LOSATP_COLLINEAR_SCENARIO_ID)


def gui_quantitative_tracks_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-09 images from the approved manifest."""

    return screenshot_paths_for(GUI_QUANTITATIVE_TRACKS_SCENARIO_ID)


def gui_annotation_tracks_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-10 images from the approved manifest."""

    return screenshot_paths_for(GUI_ANNOTATION_TRACKS_SCENARIO_ID)


def gui_styling_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-11 images from the approved manifest."""

    return screenshot_paths_for(GUI_STYLING_SCENARIO_ID)


def gui_feature_presentation_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-12 images from the approved manifest."""

    return screenshot_paths_for(GUI_FEATURE_PRESENTATION_SCENARIO_ID)


def gui_interactive_editing_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-13 images from the approved manifest."""

    return screenshot_paths_for(GUI_INTERACTIVE_EDITING_SCENARIO_ID)


def gui_session_reproduction_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-14 images from the approved manifest."""

    return screenshot_paths_for(GUI_SESSION_REPRODUCTION_SCENARIO_ID)


def gui_exports_screenshot_paths() -> dict[str, Path]:
    """Resolve all H-GUI-15 images from the approved manifest."""

    return screenshot_paths_for(GUI_EXPORTS_SCENARIO_ID)


def validate_capture_base_url(base_url: str) -> None:
    """Reject production or other non-loopback capture targets."""

    parsed = urlsplit(base_url)
    if parsed.scheme != LOCAL_SCHEME or parsed.hostname not in ALLOWED_CAPTURE_HOSTS:
        raise ValueError(
            "Documentation capture is restricted to an HTTP loopback server; "
            f"received {base_url!r}"
        )
    if parsed.port is None:
        raise ValueError(f"Documentation capture requires an explicit local port: {base_url!r}")
