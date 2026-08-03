from __future__ import annotations

import inspect
import json
import re
import subprocess
import sys
from pathlib import Path

import gbdraw
from gbdraw.annotations.io import ANNOTATION_TABLE_COLUMNS
from gbdraw.io.cli_tables import (
    _CIRCULAR_TRACK_COLUMNS,
    _COMPARISON_COLUMNS,
    _CONSERVATION_COLUMNS,
    _RECORDS_COLUMNS,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION, SUPPORTED_SESSION_VERSIONS
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_ROOT = REPO_ROOT / "docs" / "REFERENCE"
REFERENCE_FILES = {
    "web-app.md",
    "command-line.md",
    "input-formats-and-tsv-schemas.md",
    "comparison-programs-thresholds-and-results.md",
    "palettes-feature-rules-labels-shapes-and-tracks.md",
    "python-api.md",
    "typed-requests.md",
    "session-and-request-compatibility.md",
    "output-formats-and-export.md",
    "interactive-svg-and-semantic-hooks.md",
    "tutorial-fixture-provenance.md",
}
DOCUMENTATION_PAGES = (
    *(REFERENCE_ROOT / name for name in REFERENCE_FILES),
    *(REPO_ROOT / "docs" / "EXPLANATION").glob("*.md"),
    REPO_ROOT / "docs" / "GALLERY.md",
    REPO_ROOT / "docs" / "PALETTE_EXPLORER.md",
)
MARKDOWN_LINK_RE = re.compile(r"(?<!!)\[[^\]\n]+\]\(([^)\n]+)\)")


def _read(name: str) -> str:
    return (REFERENCE_ROOT / name).read_text(encoding="utf-8")


def _documented_code_tokens(source: str) -> set[str]:
    return set(re.findall(r"`([^`]+)`", source))


def _local_markdown_target(source: Path, raw_target: str) -> Path | None:
    target = raw_target.strip()
    if target.startswith("<"):
        target = target[1:].split(">", 1)[0]
    else:
        target = target.split(maxsplit=1)[0]
    if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.IGNORECASE) or target.startswith("//"):
        return None
    path_part = target.split("#", 1)[0].split("?", 1)[0]
    return source if not path_part else (source.parent / path_part).resolve()


def test_reference_census_is_complete() -> None:
    assert {path.name for path in REFERENCE_ROOT.glob("*.md")} >= (
        REFERENCE_FILES | {"README.md"}
    )


def test_session_reference_matches_current_implementation_constants() -> None:
    source = _read("session-and-request-compatibility.md")

    assert f"session version {CURRENT_SESSION_VERSION}" in source
    assert f"`renderRequest` schema {CANONICAL_REQUEST_SCHEMA}" in source
    assert SUPPORTED_SESSION_VERSIONS == frozenset(
        {27, 28, 29, 30, 31, 32, 33, 39, 40}
    )
    assert SUPPORTED_CANONICAL_REQUEST_SCHEMAS == frozenset({1, 2, 5})
    assert "27–33 and 39–40" in source
    assert "1, 2, and 5" in source


def test_input_reference_lists_parser_owned_table_columns() -> None:
    tokens = _documented_code_tokens(_read("input-formats-and-tsv-schemas.md"))

    for columns in (
        _RECORDS_COLUMNS,
        _COMPARISON_COLUMNS,
        _CONSERVATION_COLUMNS,
        _CIRCULAR_TRACK_COLUMNS,
        ANNOTATION_TABLE_COLUMNS,
    ):
        assert set(columns) <= tokens
    assert {"reference_name", "position", "depth"} <= tokens


def test_command_reference_entry_points_exist_in_current_help() -> None:
    source = _read("command-line.md")
    assert {"gbdraw circular", "gbdraw linear", "gbdraw gui"} <= (
        _documented_code_tokens(source)
    )

    for mode, required_option in (("circular", "--gbk"), ("linear", "--gbk")):
        result = subprocess.run(
            [sys.executable, "-m", "gbdraw.cli", mode, "--help"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        assert required_option in result.stdout


def test_first_web_tutorial_controls_have_stable_accessible_names() -> None:
    index = (REPO_ROOT / "gbdraw" / "web" / "index.html").read_text(encoding="utf-8")
    reference = _read("web-app.md")
    accessible_names = {
        "Circular",
        "Linear",
        "GenBank/DDBJ File",
        "Output Prefix",
        "Species",
        "Track Preset",
        "Separate Strands",
        "Hide GC Content",
        "Hide GC Skew",
        "Label Mode",
        "Legend Position",
        "Generate Diagram",
        "Result Preview",
        "SVG",
    }

    assert accessible_names <= set(re.findall(r"\*\*([^*]+)\*\*", reference))
    for name in accessible_names:
        assert name in index, name


def test_web_reference_defaults_match_the_browser_profile() -> None:
    reference = _read("web-app.md")
    profile = (REPO_ROOT / "gbdraw" / "web" / "js" / "web-ux-profile.js").read_text(
        encoding="utf-8"
    )
    state = (REPO_ROOT / "gbdraw" / "web" / "js" / "state.js").read_text(
        encoding="utf-8"
    )

    assert "separateStrands: true" in profile
    assert "legend: 'left'" in profile
    assert "legend: 'bottom'" in profile
    assert "track_type: 'tuckin'" in state
    assert "linear_track_layout: 'middle'" in state
    assert "| Separate strands | On | On |" in reference
    assert "| Legend | Left | Bottom |" in reference
    assert "| Feature placement | Tuckin preset | Features on axis |" in reference


def test_output_reference_covers_supported_public_format_tokens() -> None:
    source = _read("output-formats-and-export.md")
    tokens = {"SVG", "Interactive SVG", "PNG", "PDF", "EPS", "PS"}

    for token in tokens:
        assert re.search(rf"^\| {re.escape(token)} \|", source, re.MULTILINE)
    assert "<prefix>.interactive.svg" in source


def test_python_reference_core_signatures_match_the_public_api() -> None:
    source = _read("python-api.md")

    for name in ("read_genbank", "read_gff", "draw_circular", "draw_linear"):
        assert name in gbdraw.__all__
        signature = str(inspect.signature(getattr(gbdraw, name))).replace("'", "")
        assert f"{name}{signature}" in source
    for name in gbdraw.__all__:
        assert f"`{name}`" in source or name in {
            "read_genbank",
            "read_gff",
            "draw_circular",
            "draw_linear",
        }


def test_public_reference_rejects_artificially_split_lambda_fixture() -> None:
    public_reference = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(REFERENCE_ROOT.glob("*.md"))
    )

    assert "lambda_two_contigs" not in public_reference
    assert "lambda_left" not in public_reference
    assert "lambda_right" not in public_reference
    assert "naturally single sequence is never divided" in public_reference


def test_reference_explanation_and_auxiliary_local_links_resolve() -> None:
    missing: list[str] = []
    for source in DOCUMENTATION_PAGES:
        for raw_target in MARKDOWN_LINK_RE.findall(source.read_text(encoding="utf-8")):
            target = _local_markdown_target(source, raw_target)
            if target is not None and not target.exists():
                missing.append(f"{source.relative_to(REPO_ROOT)} -> {raw_target}")
    assert missing == []


def test_comparison_and_gallery_examples_use_whole_lambda_de3_and_five_bgcs() -> None:
    comparison = _read("comparison-programs-thresholds-and-results.md")
    gallery = (REPO_ROOT / "docs" / "GALLERY.md").read_text(encoding="utf-8")
    gallery_examples = json.loads(
        (REPO_ROOT / "gbdraw" / "web" / "gallery" / "examples.json").read_text(
            encoding="utf-8"
        )
    )
    fixture_manifest = json.loads(
        (REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )

    assert "NC_001416.1" in comparison
    assert "NC_042057.1" in comparison
    assert re.search(r"six\s+raw rows", comparison)
    assert re.search(r"397\s+raw rows", comparison)
    assert "T4" not in comparison
    assert fixture_manifest["fixtures"]["lambda-de3-comparison"]["expectedSemantics"][
        "recordsAreWholeCanonicalSources"
    ] is True

    assert "BGC0000708, BGC0000709, BGC0000711, BGC0000712, and BGC0000713" in gallery
    assert "LOSATP <em>Similarity groups</em>" in gallery
    assert "phylogenetic orthology" in gallery
    bgc_example = next(
        example
        for example in gallery_examples
        if example["id"] == "BGC0000708-BGC0000713"
    )
    assert bgc_example["featureSources"] == [
        "BGC0000708.gbk",
        "BGC0000709.gbk",
        "BGC0000711.gbk",
        "BGC0000712.gbk",
        "BGC0000713.gbk",
    ]
    assert bgc_example["workflow"] == "LOSATP similarity groups and color rules"
