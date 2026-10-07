"""The Python split of a Session 27-44 draft against the vectors shared with the Web.

``tests/fixtures/sessions/mode-split-vectors.json`` (Phase E owns its cases) is
read by this file and by a node test against ``splitDraftIntoModes``. A case
holds a ``fixture`` (a Session 27-44 in the repository) or an ``input`` (an
inline Session 27-44, or a flat draft after the OV-102, OV-114 and OV-132
normalizers), the split ``context``, and JSON pointers into the Session 46
result: ``expect`` (pointer -> value) and ``expectAbsent``. A case may also
hold ``expectedModes``, the JavaScript split's two slices, which the Python
slices must equal.

The split writes only migrated values; a value missing from a slice is that
mode's default. The expected values assume complete slices, so a missing
slice value is compared with the mode's default from the Web's own default
creators.
"""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path
from typing import Any

import pytest

from gbdraw.session_io import (
    expand_session_feature_catalog,
    migrate_session_flat_draft,
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
    assert VECTORS["schemaVersion"] == 1
    assert VECTORS["registryRevision"] == MODE_SCOPED_SETTINGS_REVISION


@pytest.fixture(scope="module")
def mode_defaults() -> dict[str, Any]:
    """Each mode's default slice values, from the Web's default creators."""

    script = (
        "import { createDefaultForm, createDefaultAdv, createDefaultLosat } from "
        "'./gbdraw/web/js/services/session-active-config-contract.js';\n"
        "import { createDefaultLayoutPreferences } from './gbdraw/web/js/services/layout-preferences.js';\n"
        "const layout = createDefaultLayoutPreferences();\n"
        # Until the Web exports one creator of a mode's default slice, the
        # values outside the form, adv and LOSAT creators are the state.js
        # initial values.
        "const slice = (mode) => ({\n"
        "  config: { form: createDefaultForm(), adv: createDefaultAdv(mode), losat: createDefaultLosat(),\n"
        "    losatProgram: 'blastn', linearRecordLayout: { enabled: true, recordGap: 24, rows: [] },\n"
        "    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] } },\n"
        "  features: { featureOverrides: {}, featureColorOverrides: {} },\n"
        "  editorState: {\n"
        "    legend: { entries: [], deletedEntries: [], colorOverrides: {}, strokeOverrides: {}, addedCaptions: [] },\n"
        "    featureStrokes: { overrides: {} }\n"
        "  },\n"
        "  ui: { layoutPreferences: layout[mode], selectedFeatureRecordIdx: 0, linearTypographyLinked: true }\n"
        "});\n"
        "process.stdout.write(JSON.stringify({ circular: slice('circular'), linear: slice('linear') }));\n"
    )
    result = subprocess.run(
        ["node", "--input-type=module", "-e", script],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(result.stdout)


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


def _split_case(case: dict[str, Any]) -> dict[str, Any]:
    context = case["context"]
    source = _read_fixture(case["fixture"]) if "fixture" in case else json.loads(json.dumps(case["input"]))
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
        result = split_draft_into_modes(migrated, committed_mode=derived["committedMode"])
        result["version"] = 46
        return result
    return split_draft_into_modes(
        source,
        committed_mode=context["committedMode"],
        mode_profiles=context.get("modeProfiles"),
        depth_sources=context.get("depthSources"),
    )


@pytest.mark.parametrize("case", CASES, ids=[case["name"][:80] for case in CASES])
def test_python_split_matches_the_shared_vectors(case: dict[str, Any], mode_defaults: dict[str, Any]) -> None:
    result = _split_case(case)

    mismatches = []
    for pointer, expected in case.get("expect", {}).items():
        actual = _resolve(result, pointer)
        parts = pointer.split("/")
        if actual is _ABSENT and len(parts) > 2 and parts[1] == "modes" and parts[2] in DIAGRAM_MODES:
            # Missing means default.
            actual = _resolve(mode_defaults[parts[2]], "/" + "/".join(parts[3:]))
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
