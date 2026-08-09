"""Authoritative environment contract for documentation screenshots."""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit


CAPTURE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = CAPTURE_ROOT.parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
MANIFEST_PATH = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
TUTORIAL_DATA_MANIFEST_PATH = WEB_ROOT / "tutorial-data" / "manifest.json"
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


@cache
def load_manifest() -> dict[str, Any]:
    """Load the approved chapter manifest."""

    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


@cache
def load_tutorial_data_manifest() -> dict[str, Any]:
    """Load the tutorial fixture manifest."""

    return json.loads(TUTORIAL_DATA_MANIFEST_PATH.read_text(encoding="utf-8"))


def chapter_for(scenario_id: str) -> dict[str, Any]:
    """Return one approved scenario entry."""

    for chapter in load_manifest()["chapters"]:
        if chapter["id"] == scenario_id:
            return chapter
    raise KeyError(f"Unknown documentation scenario: {scenario_id}")


def scenario_ids_for(execution_kind: str) -> tuple[str, ...]:
    """Return verified scenario IDs for one executable surface."""

    return tuple(
        chapter["id"]
        for chapter in load_manifest()["chapters"]
        if chapter.get("execution", {}).get("kind") == execution_kind
        and chapter.get("status", {}).get("implementation") == "verified"
    )


def scenario_tier(scenario_id: str) -> str:
    """Return the manifest-owned capture tier for one scenario."""

    return str(chapter_for(scenario_id)["execution"]["tier"])


def supported_tiers() -> tuple[str, ...]:
    """Return execution tiers in their manifest-defined order."""

    return tuple(
        dict.fromkeys(
            str(chapter["execution"]["tier"])
            for chapter in load_manifest()["chapters"]
            if chapter.get("execution", {}).get("tier")
        )
    )


def _repo_path(relative_path: str) -> Path:
    path = (REPO_ROOT / relative_path).resolve()
    if not path.is_relative_to(REPO_ROOT):
        raise ValueError(f"Documentation path escapes the repository: {relative_path}")
    return path


def screenshot_names_for(scenario_id: str) -> tuple[str, ...]:
    """Return screenshot filenames in manifest order."""

    return tuple(
        Path(screenshot["path"]).name
        for screenshot in chapter_for(scenario_id)["screenshots"]
    )


def screenshot_paths_for(scenario_id: str) -> dict[str, Path]:
    """Resolve all screenshot paths for one scenario from the manifest."""

    screenshots = chapter_for(scenario_id)["screenshots"]
    paths = {
        Path(screenshot["path"]).name: _repo_path(screenshot["path"])
        for screenshot in screenshots
    }
    if len(paths) != len(screenshots):
        raise ValueError(f"Duplicate screenshot filename for {scenario_id}")
    return paths


def tutorial_file_for(file_id: str) -> dict[str, Any]:
    """Return one tutorial fixture file entry."""

    try:
        return load_tutorial_data_manifest()["files"][file_id]
    except KeyError:
        raise KeyError(f"Unknown tutorial fixture file: {file_id}") from None


def tutorial_file_identity(file_id: str) -> tuple[Path, int, str]:
    """Resolve a tutorial fixture path, byte size, and digest."""

    manifest = load_tutorial_data_manifest()
    entry = tutorial_file_for(file_id)
    path = _repo_path(f"{manifest['canonicalRoot']}/{entry['relativePath']}")
    return path, int(entry["sizeBytes"]), str(entry["sha256"])


def genbank_fixtures_for(
    fixture_id: str,
) -> tuple[tuple[Path, int, str, str, int, str], ...]:
    """Resolve ordered single-record GenBank metadata for one fixture set."""

    manifest = load_tutorial_data_manifest()
    try:
        fixture = manifest["fixtures"][fixture_id]
    except KeyError:
        raise KeyError(f"Unknown tutorial fixture: {fixture_id}") from None

    entries = []
    for file_id in (*fixture.get("fileIds", ()), *fixture.get("fileReferences", ())):
        entry = tutorial_file_for(file_id)
        if entry.get("inputType") == "genbank":
            entries.append(entry)
    by_record_id: dict[str, dict[str, Any]] = {}
    for entry in entries:
        records = entry.get("records", ())
        if len(records) != 1:
            raise ValueError(
                f"{fixture_id} requires one record per GenBank fixture file"
            )
        by_record_id[str(records[0]["id"])] = entry

    record_ids = fixture.get("expectedSemantics", {}).get("recordIds")
    if record_ids:
        try:
            entries = [by_record_id[str(record_id)] for record_id in record_ids]
        except KeyError as exc:
            raise ValueError(
                f"{fixture_id} record order references an unknown GenBank record"
            ) from exc

    canonical_root = str(manifest["canonicalRoot"])
    resolved = []
    for entry in entries:
        record = entry["records"][0]
        resolved.append(
            (
                _repo_path(f"{canonical_root}/{entry['relativePath']}"),
                int(entry["sizeBytes"]),
                str(entry["sha256"]),
                str(record["id"]),
                int(record["length"]),
                str(record["organism"]),
            )
        )
    return tuple(resolved)


# Capture flows still consume these descriptive names; their values come from
# the two manifests above rather than restating screenshot or fixture metadata.
FIRST_CIRCULAR_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-01")
FIRST_LINEAR_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-02")
GUI_INPUTS_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-01")
GUI_LOSATN_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-03")
GUI_CIRCULAR_LAYOUT_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-02")
GUI_LINEAR_LAYOUT_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-03")
GUI_UPLOADED_COMPARISON_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-04")
GUI_TLOSATX_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-05")
GUI_CIRCULAR_RINGS_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-06")
GUI_LOSATP_GROUPS_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-04")
GUI_ANNOTATED_CHLOROPLAST_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-05")
GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-06")
GUI_LOSATP_TUTORIAL_COLLINEAR_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-08")
GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-09")
GUI_FEATURE_HIGHLIGHT_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-10")
GUI_QUANTITATIVE_MAP_SCREENSHOT_NAMES = screenshot_names_for("T-GUI-12")
GUI_LOSATP_GROUPS_HOW_TO_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-07")
GUI_LOSATP_COLLINEAR_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-08")
GUI_QUANTITATIVE_TRACKS_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-09")
GUI_ANNOTATION_TRACKS_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-10")
GUI_STYLING_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-11")
GUI_FEATURE_PRESENTATION_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-12")
GUI_INTERACTIVE_EDITING_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-13")
GUI_SESSION_REPRODUCTION_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-14")
GUI_EXPORTS_SCREENSHOT_NAMES = screenshot_names_for("H-GUI-15")

(
    FIRST_CIRCULAR_FIXTURE_PATH,
    FIRST_CIRCULAR_FIXTURE_SIZE,
    FIRST_CIRCULAR_FIXTURE_SHA256,
) = tutorial_file_identity("human-mitochondrion-genbank")
(
    FIRST_LINEAR_FIXTURE_PATH,
    FIRST_LINEAR_FIXTURE_SIZE,
    FIRST_LINEAR_FIXTURE_SHA256,
) = tutorial_file_identity("lambda-genbank")
(
    FIRST_LINEAR_LABEL_RULE_PATH,
    FIRST_LINEAR_LABEL_RULE_SIZE,
    FIRST_LINEAR_LABEL_RULE_SHA256,
) = tutorial_file_identity("shared-cds-gene-qualifier-priority")
(
    GUI_INPUTS_GFF3_PATH,
    GUI_INPUTS_GFF3_SIZE,
    GUI_INPUTS_GFF3_SHA256,
) = tutorial_file_identity("lambda-full-record-gff3")
(
    GUI_INPUTS_FASTA_PATH,
    GUI_INPUTS_FASTA_SIZE,
    GUI_INPUTS_FASTA_SHA256,
) = tutorial_file_identity("lambda-full-record-fasta")
(
    GUI_LOSATN_DE3_FIXTURE_PATH,
    GUI_LOSATN_DE3_FIXTURE_SIZE,
    GUI_LOSATN_DE3_FIXTURE_SHA256,
) = tutorial_file_identity("de3-genbank")
(
    GUI_LOSATN_REFERENCE_TSV_PATH,
    GUI_LOSATN_REFERENCE_TSV_SIZE,
    GUI_LOSATN_REFERENCE_TSV_SHA256,
) = tutorial_file_identity("lambda-de3-losatn")

GUI_CIRCULAR_LAYOUT_FIXTURES = genbank_fixtures_for("metazoan-mitochondria-four")
GUI_CIRCULAR_LAYOUT_COMBINED_FILENAME = "complete_metazoan_mitochondria.gb"
GUI_BGC_FIXTURES = genbank_fixtures_for("aminoglycoside-bgc-five")
GUI_HEPATOPLASMATACEAE_FIXTURES = genbank_fixtures_for("hepatoplasmataceae-five")


def validate_capture_base_url(base_url: str) -> None:
    """Reject production or other non-loopback capture targets."""

    parsed = urlsplit(base_url)
    if parsed.scheme != LOCAL_SCHEME or parsed.hostname not in ALLOWED_CAPTURE_HOSTS:
        raise ValueError(
            "Documentation capture is restricted to an HTTP loopback server; "
            f"received {base_url!r}"
        )
    if parsed.port is None:
        raise ValueError(
            f"Documentation capture requires an explicit local port: {base_url!r}"
        )
