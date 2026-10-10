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

The last test is the design's composition oracle (§3.2): every Session 27-44
fixture through the chain the CLI save path runs.
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


# Design §3.2 (composition oracle): a Session 27-44 brought to Session 46 by
# the chain the CLI save path runs (upgrade_session_document: the 27-44
# migrations, then this split) holds one drawing, of its committed mode A.
# Drawing B is the 46 -> 47 step's (§3.1 rule 3), which reads B's signals
# from the 46: ui.mode (b), a binding in B's input slots (c), and B rows that
# A does not hold (d). F1-F3 carry them; a single-mode fixture carries none of
# b, c, or the keyed and Result-mode rows of d. The value rows of d need the
# Web defaults and are pinned by the expectations of the vectors above.
_SESSIONS = REPO_ROOT / "tests" / "fixtures" / "sessions"
_LEGACY_FIXTURES = sorted(
    path.relative_to(REPO_ROOT).as_posix()
    for path in _SESSIONS.iterdir()
    if path.name.endswith((".json", ".json.gz"))
    and _read_fixture(path.relative_to(REPO_ROOT).as_posix()).get("format") == "gbdraw-session"
    and _read_fixture(path.relative_to(REPO_ROOT).as_posix()).get("version", 46) <= 44
)
# The drawing B of each two-mode fixture (§3.3 E0) and the signals the 46 holds.
_TWO_MODE_FIXTURES = {
    "two-mode-project.v44.gbdraw-session.json.gz": ("circular", {"ui.mode", "binding", "rows"}),
    "inactive-class-m.v44.gbdraw-session.json.gz": ("linear", {"rows"}),
    "two-mode-thresholds.v42.gbdraw-session.json.gz": ("circular", {"rows"}),
}
_FIXTURE_VECTORS = {case["fixture"]: case for case in CASES if "fixture" in case and "omit" not in case}


def _bound(value: Any) -> bool:
    if isinstance(value, list):
        return any(_bound(item) for item in value)
    return value is not None


def _mode_bindings(bindings: Any, mode: str) -> bool:
    """A non-null binding in the mode's input slots; the empty Linear row is no binding."""

    bindings = bindings if isinstance(bindings, dict) else {}
    if mode == "circular":
        return any(_bound(value) for key, value in bindings.items() if key.startswith("c_"))
    rows = bindings.get("linearSeqs") if isinstance(bindings.get("linearSeqs"), list) else []
    return _bound(bindings.get("linearComparisons") or []) or any(
        isinstance(row, dict) and any(_bound(row.get(key)) for key in ("gb", "gff", "fasta", "depth", "blast"))
        for row in rows
    )


def _rows_only_in(document: dict[str, Any], mode: str, other: str) -> list[str]:
    """The keyed and Result-mode registry rows of ``mode``'s slice that ``other`` lacks."""

    from gbdraw.web_support.mode_scoped_settings import MODE_SCOPED_SETTINGS

    rows = []
    for row in MODE_SCOPED_SETTINGS:
        if row.key is None and row.migrate != "result-mode":
            continue
        pointer = "/" + "/".join([*row.domain.split("."), row.path])
        # A Session written without a draft (CLI, Python API) has no slices.
        slices = document.get("modes", {})
        value = _resolve(slices.get(mode, {}), pointer)
        if value is _ABSENT or value in ({}, [], None):
            continue
        kept = _resolve(slices.get(other, {}), pointer)
        if isinstance(value, dict) and isinstance(kept, dict):
            if any(key not in kept or _canonical(kept[key]) != _canonical(item) for key, item in value.items()):
                rows.append(pointer)
        elif isinstance(value, list) and isinstance(kept, list):
            # An annotation set without annotations is the default shell that
            # by-binding writes to both slices; it carries no use of B.
            used = {
                _canonical(item)
                for item in value
                if not (row.path == "annotationSets" and isinstance(item, dict) and not item.get("annotations"))
            }
            if not used <= {_canonical(item) for item in kept}:
                rows.append(pointer)
        elif kept is _ABSENT or _canonical(kept) != _canonical(value):
            rows.append(pointer)
    return rows


@pytest.mark.parametrize("fixture", _LEGACY_FIXTURES, ids=[Path(path).name for path in _LEGACY_FIXTURES])
def test_the_27_44_chain_composes_to_one_drawing_and_carries_drawing_b(fixture: str, tmp_path: Path) -> None:
    from gbdraw.session import upgrade_session_document
    from gbdraw.session_drawings import DRAWING_NAMES, read_session_drawings

    source = _read_fixture(fixture)
    name = Path(fixture).name
    request = source.get("renderRequest")
    # §3.2 rule 2: A is the request's mode; settings-only, ui.mode; 27-30, session_mode().
    if source["version"] <= 30:
        mode_a = session_mode(source)
    elif isinstance(request, dict):
        mode_a = request["mode"]
    else:
        mode_a = source["ui"]["mode"]
    if source["version"] <= 30:
        # No canonical request: the CLI replays it (build_session_json), and its
        # Session view is its one drawing of its declared mode.
        drawings = read_session_drawings(source)
        assert [(drawing.id, drawing.name, drawing.request) for drawing in drawings] == [
            (mode_a, DRAWING_NAMES[mode_a], None)
        ]
        return
    upgraded = upgrade_session_document(REPO_ROOT / fixture, temporary_directory=tmp_path).document.to_dict()

    assert upgraded["version"] == 46
    assert [(drawing.id, drawing.name, drawing.request is not None) for drawing in read_session_drawings(upgraded)] == [
        (mode_a, DRAWING_NAMES[mode_a], isinstance(request, dict))
    ]
    case = _FIXTURE_VECTORS.get(fixture)
    if case is not None:
        # The chain's slices are the split's (the vectors' Phase E oracle).
        for pointer, expected in case.get("expect", {}).items():
            if not pointer.startswith(("/renderRequest", "/version")):
                assert _canonical(_resolve(upgraded, pointer)) == _canonical(expected), pointer
        for pointer in case.get("expectAbsent", []):
            if not pointer.startswith(("/renderRequest", "/version")):
                assert _resolve(upgraded, pointer) is _ABSENT, pointer
        for mode in case.get("expectedModes", {}):
            assert _canonical(upgraded["modes"][mode]) == _canonical(case["expectedModes"][mode]), mode
    mode_b = "linear" if mode_a == "circular" else "circular"
    # §3.2 rule 4: the active values (the Result-mode registry rows: rendered-ID
    # per-feature edits, Legend editor state, feature strokes, the selected
    # record) go to A only, even where B's would equal A's.
    from gbdraw.web_support.mode_scoped_settings import MODE_SCOPED_SETTINGS

    slice_b = upgraded.get("modes", {}).get(mode_b, {})
    assert [
        pointer
        for pointer in (
            "/" + "/".join([*row.domain.split("."), row.path])
            for row in MODE_SCOPED_SETTINGS
            if row.migrate == "result-mode"
        )
        if _resolve(slice_b, pointer) is not _ABSENT
    ] == []
    signals = {
        signal
        for signal, present in (
            ("ui.mode", upgraded.get("ui", {}).get("mode") == mode_b),
            ("binding", _mode_bindings(upgraded.get("webFiles", {}).get("bindings"), mode_b)),
            ("rows", bool(_rows_only_in(upgraded, mode_b, mode_a))),
        )
        if present
    }
    if name in _TWO_MODE_FIXTURES:
        expected_b, expected_signals = _TWO_MODE_FIXTURES[name]
        assert mode_b == expected_b
        if "rows" in expected_signals and "rows" not in signals:
            # The value rows of rule d: B's slice holds a vector-pinned value
            # that A's does not.
            assert case is not None
            assert any(
                pointer.startswith(f"/modes/{mode_b}/")
                and _canonical(_resolve(upgraded, f"/modes/{mode_a}/" + pointer.split("/", 3)[3]))
                != _canonical(expected)
                for pointer, expected in case.get("expect", {}).items()
            )
            signals.add("rows")
        assert signals == expected_signals
    else:
        # Negative control: a single-mode fixture gives drawing B no signal.
        assert signals == set(), _rows_only_in(upgraded, mode_b, mode_a)
