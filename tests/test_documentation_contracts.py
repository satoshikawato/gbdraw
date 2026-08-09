from __future__ import annotations

import re
from pathlib import Path

from gbdraw.session_io import CURRENT_SESSION_VERSION, SUPPORTED_SESSION_VERSIONS
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
RELEASE_NOTES = REPO_ROOT / "docs" / "RELEASE_NOTES_0.14.0b0.md"
SESSION_COMPATIBILITY = (
    REPO_ROOT / "docs" / "REFERENCE" / "session-and-request-compatibility.md"
)
COMPATIBILITY_HISTORY = REPO_ROOT / "docs" / "SESSION_COMPATIBILITY.md"
BROWSER_ACCEPTANCE = REPO_ROOT / "tests" / "run_losat_cache_browser_acceptance.py"
CURRENT_TASK_DOCS = (
    "docs/CLI_Reference.md",
    "docs/FAQ.md",
    "docs/HOW_TO/GUI/save-restore-undo-and-reproduce-work.md",
    "docs/HOW_TO/CLI/save-and-regenerate-sessions.md",
    "docs/HOW_TO/PYTHON/build-typed-requests-and-round-trip-sessions.md",
)


def _parse_documented_versions(value: str) -> frozenset[int]:
    versions: set[int] = set()
    for token in value.replace("and", ",").split(","):
        token = token.strip()
        if not token:
            continue
        if "–" in token:
            start, end = (int(part) for part in token.split("–", maxsplit=1))
            versions.update(range(start, end + 1))
        else:
            versions.add(int(token))
    return frozenset(versions)


def test_session_compatibility_table_matches_implementation() -> None:
    text = SESSION_COMPATIBILITY.read_text(encoding="utf-8")
    session_row = re.search(
        r"^\| gbdraw session \| (\d+) \| ([^|]+) \|$", text, re.MULTILINE
    )
    request_row = re.search(
        r"^\| Canonical `renderRequest` \| (\d+) \| ([^|]+) \|$",
        text,
        re.MULTILINE | re.IGNORECASE,
    )

    assert session_row is not None
    assert request_row is not None
    assert int(session_row.group(1)) == CURRENT_SESSION_VERSION
    assert _parse_documented_versions(session_row.group(2)) == SUPPORTED_SESSION_VERSIONS
    assert int(request_row.group(1)) == CANONICAL_REQUEST_SCHEMA
    assert (
        _parse_documented_versions(request_row.group(2))
        == SUPPORTED_CANONICAL_REQUEST_SCHEMAS
    )


def test_release_session_history_and_acceptance_use_current_authority() -> None:
    release_notes = RELEASE_NOTES.read_text(encoding="utf-8")
    acceptance_source = BROWSER_ACCEPTANCE.read_text(encoding="utf-8")

    assert f"## Python/Web session version {CURRENT_SESSION_VERSION}" in release_notes
    assert re.search(
        rf"gbdraw 0\.14\.0b0 writes session version {CURRENT_SESSION_VERSION} and\s+"
        rf"canonical `renderRequest`\s+schema {CANONICAL_REQUEST_SCHEMA}",
        release_notes,
    )
    assert "Session version 39 introduced compact runtime handles" in release_notes
    assert "Current session version 39" not in release_notes
    assert "Current writers emit version 39" not in release_notes
    assert "from gbdraw.session_io import CURRENT_SESSION_VERSION" in acceptance_source
    assert re.search(
        r"^CURRENT_SESSION_VERSION\s*=\s*\d+", acceptance_source, re.MULTILINE
    ) is None


def test_current_task_docs_delegate_persisted_format_details_to_one_authority() -> None:
    forbidden_details = (
        "__gbdraw_legacy_spacing",
        "schema-2 protein",
        "versions 34–38",
        "request schemas 3–4",
        "ui.layoutPreferences",
        "h_[a-z2-7]{26}",
    )

    for relative_path in CURRENT_TASK_DOCS:
        text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
        assert (
            "SESSION_COMPATIBILITY.md" in text
            or "session-and-request-compatibility.md" in text
        )
        for detail in forbidden_details:
            assert detail not in text, f"{relative_path} repeats {detail}"

    authority = COMPATIBILITY_HISTORY.read_text(encoding="utf-8")
    for detail in forbidden_details:
        assert detail in authority


def test_current_compatibility_reference_and_history_have_distinct_roles() -> None:
    current = SESSION_COMPATIBILITY.read_text(encoding="utf-8")
    history = COMPATIBILITY_HISTORY.read_text(encoding="utf-8")

    assert "this page is the current support authority" in current
    assert "# Session and request compatibility history" in history
    assert "canonical owner of current support" in history
    assert "This page is the current compatibility reference" not in history


def test_generated_cli_inventory_delegates_current_semantics() -> None:
    inventory = (REPO_ROOT / "docs/CLI_Reference.md").read_text(encoding="utf-8")
    semantics = (REPO_ROOT / "docs/REFERENCE/command-line.md").read_text(
        encoding="utf-8"
    )

    assert "# Generated command-line option inventory" in inventory
    assert "[command-line reference](./REFERENCE/command-line.md)" in inventory
    assert "canonical owner of current CLI semantics" in inventory
    assert "This page is the canonical owner of current command semantics" in semantics


def test_python_api_describes_its_four_output_forms() -> None:
    text = (REPO_ROOT / "docs/REFERENCE/python-api.md").read_text(encoding="utf-8")

    for token in ("Diagram", "Diagram.save()", "Diagram.to_svg()", "Diagram.to_bytes()"):
        assert token in text
    assert "[Typed request reference](typed-requests.md)" in text


def test_web_maintenance_guide_describes_typed_worker_rendering() -> None:
    text = (REPO_ROOT / "gbdraw/web/CLAUDE.md").read_text(encoding="utf-8")

    assert "canonical schema-5 render request" in text
    assert "workers/diagram-generation-worker.js" in text
    assert "typed request decoding and render_request()" in text
    assert "setup() ~" not in text
    assert "args ~" not in text
