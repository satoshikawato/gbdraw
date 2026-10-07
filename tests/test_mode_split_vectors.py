"""The Python split of a Session 27-44 draft against the vectors shared with the Web.

``tests/fixtures/sessions/mode-split-vectors.json`` (Phase E owns its cases) is
read by this file and by a node test against ``splitDraftIntoModes``. A case
holds a ``fixture`` (a Session 27-44 in the repository) or an ``input`` (an
inline Session 27-44, or a flat draft after the OV-102, OV-114 and OV-132
normalizers), the split ``context``, and JSON pointers into the Session 46
result: ``expect`` (pointer -> value) and ``expectAbsent``. A fixture case may
``omit`` fixture fields (JSON pointers removed before the migration). A case
may also hold ``expectedModes``, the JavaScript split's two slices, which the
Python slices must equal.

The split writes only migrated values (partial slices); Web Load fills the
rest, so a value the split does not write is asserted absent.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.session_io import (
    expand_session_feature_catalog,
    migrate_session_flat_draft,
    mode_split_palette_colors,
    session_depth_source_widths,
    session_mode,
    split_draft_into_modes,
)
from gbdraw.web_support.mode_scoped_settings import DIAGRAM_MODES, MODE_SCOPED_SETTINGS_REVISION

REPO_ROOT = Path(__file__).resolve().parents[1]
VECTORS_PATH = REPO_ROOT / "tests" / "fixtures" / "sessions" / "mode-split-vectors.json"
VECTORS = json.loads(VECTORS_PATH.read_text(encoding="utf-8"))
CASES = VECTORS["cases"]
_ABSENT = object()


def test_vectors_follow_this_registry_revision() -> None:
    assert VECTORS["schemaVersion"] == 2
    assert VECTORS["registryRevision"] == MODE_SCOPED_SETTINGS_REVISION


def _resolve(document: Any, pointer: str) -> Any:
    """The value at an RFC 6901 JSON pointer, or ``_ABSENT``."""

    current = document
    for token in pointer.split("/")[1:]:
        token = token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, dict) and token in current:
            current = current[token]
        elif isinstance(current, list) and token.isdigit() and int(token) < len(current):
            current = current[int(token)]
        else:
            return _ABSENT
    return current


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def _read_fixture(path: str) -> dict[str, Any]:
    data = (REPO_ROOT / path).read_bytes()
    return json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)


def _remove(document: Any, pointer: str) -> None:
    *parents, last = [token.replace("~1", "/").replace("~0", "~") for token in pointer.split("/")[1:]]
    for token in parents:
        document = document[int(token)] if isinstance(document, list) else document[token]
    if isinstance(document, list):
        del document[int(last)]
    else:
        del document[last]


def _split_case(case: dict[str, Any]) -> dict[str, Any]:
    context = case["context"]
    source = _read_fixture(case["fixture"]) if "fixture" in case else json.loads(json.dumps(case["input"]))
    for pointer in case.get("omit", []):
        _remove(source, pointer)
    if source.get("format") == "gbdraw-session":
        # A whole Session: the CLI's own context and migrations, then the split.
        session = expand_session_feature_catalog(source)
        web_files = session.get("webFiles")
        bindings = web_files.get("bindings") if isinstance(web_files, dict) else session.get("files")
        config = session.get("config") if isinstance(session.get("config"), dict) else {}
        derived = {
            "committedMode": session_mode(session),
            "modeProfiles": config.get("modeProfiles"),
            "depthSources": session_depth_source_widths(bindings),
        }
        for name, value in derived.items():
            if name in context:
                assert value == context[name], f"{case['name']}: context {name}"
        migrated = migrate_session_flat_draft(session).session
        result = split_draft_into_modes(
            migrated,
            committed_mode=derived["committedMode"],
            palette_colors=mode_split_palette_colors(migrated.get("config")),
        )
        result["version"] = 46
        return result
    return split_draft_into_modes(
        source,
        committed_mode=context["committedMode"],
        mode_profiles=context.get("modeProfiles"),
        depth_sources=context.get("depthSources"),
        palette_colors=context.get("paletteColors"),
    )


@pytest.mark.parametrize("case", CASES, ids=[case["name"][:80] for case in CASES])
def test_python_split_matches_the_shared_vectors(case: dict[str, Any]) -> None:
    result = _split_case(case)

    mismatches = []
    for pointer, expected in case.get("expect", {}).items():
        actual = _resolve(result, pointer)
        if actual is _ABSENT or _canonical(actual) != _canonical(expected):
            shown = "<absent>" if actual is _ABSENT else _canonical(actual)[:200]
            mismatches.append(f"{pointer}: expected {_canonical(expected)[:200]}, got {shown}")
    for pointer in case.get("expectAbsent", []):
        if _resolve(result, pointer) is not _ABSENT:
            mismatches.append(f"{pointer}: expected absent, got {_canonical(_resolve(result, pointer))[:200]}")
    if "expectedModes" in case:
        for mode in DIAGRAM_MODES:
            if _canonical(result["modes"][mode]) != _canonical(case["expectedModes"][mode]):
                mismatches.append(f"/modes/{mode}: differs from the JavaScript split")
    assert not mismatches, "\n".join(mismatches)
