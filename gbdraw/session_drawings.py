"""Where each drawing of a Session lives: the one owner of the document layout.

A drawing is one diagram of a Session: its mode, its committed request and
Results, and the editor state that belongs to it. Selection, rendering, the
CLI, output naming and the Session writers reach a drawing's parts only
through this module, so a new document layout changes this module alone.

Sessions 46-49 hold at most one drawing of each mode:

- the top-level committed set: ``renderRequest``, ``results``, their catalog,
  ``runMetadata`` and ``cliInvocation``;
- the other mode's committed set in ``otherModeResult``;
- each mode's draft (settings, per-feature edits, Legend edits and the
  mode's ``ui`` keys) in its slice ``modes[mode]``. A mode may have a slice
  without a drawing.

The drawings share every other field: the other ``ui`` keys, ``webFiles``,
``cliOptions`` and the LOSAT artifacts. A settings-only Session
(``renderRequest: null``) is one drawing of ``ui.mode`` without a request; a
Session 27-30 is one drawing of its declared mode without a canonical
request. Sessions 31-44 keep one flat draft (``config``, ``features``) for
their one drawing. A drawing's ID is its mode and its name is the mode's name.

Session 50 lists its drawings in ``drawings`` (tab order), each with its ID,
name and mode; the project keeps only the fields in :data:`PROJECT_FIELDS`,
the ``ui`` keys in :data:`PROJECT_UI_FIELDS` with ``ui.activeDrawingId``, and
``webFiles.resourceOriginalNames``. Readers see a Session 50 drawing as the
same one-drawing view that a Session 46-49 drawing gives
(:func:`_drawing_view`); the writer is its inverse (:func:`_drawing_of_view`).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from importlib import resources
from types import MappingProxyType
import json
from typing import Any, Callable, Literal, Mapping, Sequence, TypeAlias, TypedDict, cast, get_args

from gbdraw.exceptions import ValidationError
from gbdraw.session_contract import (
    DRAWING_ID_MAX_LENGTH,
    SESSION_ID_RE,
    GbdrawSession50,
    SessionDrawing,
    SessionProjectUi,
    SessionProjectWebFiles,
    keys_of,
)
from gbdraw.session_io import (
    CANONICAL_SESSION_MIN_VERSION,
    CURRENT_SESSION_VERSION,
    DRAWINGS_SESSION_MIN_VERSION,
    FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION,
    MODE_SCOPED_SESSION_MAX_VERSION,
    MODE_SCOPED_SESSION_MIN_VERSION,
    OTHER_MODE_RESULT_EDITOR_FIELDS,
    OTHER_MODE_RESULT_FIELDS,
    OTHER_MODE_RESULT_LEGEND_FIELDS,
    OTHER_MODE_RESULT_UI_FIELDS,
    WEB_BINDING_CARD_SLOTS,
    ROW_KEYED_LEGEND_SESSION_MIN_VERSION,
    WEB_BINDING_MODE_SLOTS,
    session_legend_edits_by_row_key,
    session_mode,
)
from gbdraw.web_support.mode_scoped_settings import (
    MODE_SCOPED_SETTINGS,
    SLICE_CONTAINERS,
    ModeScopedSetting,
)

DrawingMode = Literal["circular", "linear"]
# A committed request as a document stores it, of any supported schema.
StoredRequest: TypeAlias = Mapping[str, object]
# A Web draft config as a document stores it: an open object.
DraftConfig: TypeAlias = Mapping[str, object]
DRAWING_NAMES: Mapping[str, str] = MappingProxyType(
    {"circular": "Circular", "linear": "Linear"}
)

# The Web reader names a document it cannot lay out as a Session field error.
_LAYOUT_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}
# The fields of one committed Result set, named as at the top level: the set's
# own fields, and its part of the two objects it shares with the other fields.
_SET_FIELDS = OTHER_MODE_RESULT_FIELDS - {"editorState", "ui"}
_SET_EDITOR_FIELDS = OTHER_MODE_RESULT_EDITOR_FIELDS - {"legend"}

# Session 50 (declared in gbdraw.session_contract): the project fields beside
# ``ui``, ``webFiles`` and ``drawings``; the ``ui`` keys of the whole project
# (a drawing's ``ui`` holds every other view key and its mode slice's ``ui``
# keys); resource metadata, the only ``webFiles`` key of the project; and the
# fields a drawing keeps under the names a one-drawing view gives them.
PROJECT_FIELDS = keys_of(GbdrawSession50) - {"ui", "webFiles", "drawings"}
PROJECT_UI_FIELDS = keys_of(SessionProjectUi) - {"activeDrawingId"}
PROJECT_WEB_FILES_FIELDS = keys_of(SessionProjectWebFiles)
_DRAWING_FIELDS = keys_of(SessionDrawing) - {
    "id", "name", "mode", "config", "features", "editorState", "ui", "webFiles",
}
# The ``webFiles`` request metadata each mode's drawing keeps (D-15-79).
_DRAWING_WEB_FILES_FIELDS: Mapping[str, frozenset[str]] = MappingProxyType({
    "circular": frozenset({
        "circularOutputPrefixExplicit", "circularInputOriginalName",
        # D-15-83: the conservation inputs are the Circular drawing's.
        "conservationLosatFastaSources", "conservationSequenceSources", "conservationBlastSource",
    }),
    "linear": frozenset({"linearRecordMetadata"}),
})
# The view's input-type key of each mode; a drawing's ``ui`` names it ``inputType``.
_INPUT_TYPE_FIELDS: Mapping[str, str] = MappingProxyType({"circular": "cInputType", "linear": "lInputType"})
_SLICE_UI_FIELDS = SLICE_CONTAINERS["ui"]
_SLICE_EDITOR_FIELDS = SLICE_CONTAINERS["editorState"] - {"legend"}
_SLICE_LEGEND_FIELDS = SLICE_CONTAINERS["editorState.legend"]


@dataclass(frozen=True)
class SessionDrawingArtifacts:
    """One drawing of a validated Session, seen as a one-drawing Session.

    ``fields`` names the drawing's request, Results, editor state, draft and
    LOSAT artifacts as the top level of a one-drawing Session names them. Only
    :class:`gbdraw.session.SessionDocument` builds it, from a validated
    document; readers do not change it.
    """

    version: int
    fields: Mapping[str, Any]


@dataclass(frozen=True)
class SessionDrawingParts:
    """One drawing as the document holds it: identity, request and Results."""

    id: str
    name: str
    mode: DrawingMode
    request: Mapping[str, Any] | None
    results: tuple[Mapping[str, Any], ...]
    resource_ids: frozenset[str]
    artifacts: SessionDrawingArtifacts


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _only(source: Mapping[str, Any], fields: frozenset[str]) -> dict[str, Any]:
    return {key: value for key, value in source.items() if key in fields}


def _without(source: Mapping[str, Any], fields: frozenset[str]) -> dict[str, Any]:
    return {key: value for key, value in source.items() if key not in fields}


def _shared(view: Mapping[str, Any]) -> dict[str, Any]:
    """Every field of a one-drawing view outside its committed Result set."""

    editor = _mapping(view.get("editorState"))
    shared = _without(view, _SET_FIELDS | {"otherModeResult"})
    shared["editorState"] = {
        **_without(editor, OTHER_MODE_RESULT_EDITOR_FIELDS),
        "legend": _without(_mapping(editor.get("legend")), OTHER_MODE_RESULT_LEGEND_FIELDS),
    }
    shared["ui"] = _without(_mapping(view.get("ui")), OTHER_MODE_RESULT_UI_FIELDS)
    return shared


def _view(shared: Mapping[str, Any], result_set: Mapping[str, Any]) -> dict[str, Any]:
    """The shared fields with one committed Result set at the top level.

    A per-set field that the set lacks is absent, so a reader takes its
    default; it is never another set's value.
    """

    shared_editor = _mapping(shared.get("editorState"))
    set_editor = _mapping(result_set.get("editorState"))
    view = {**_without(shared, frozenset({"editorState", "ui"})), **_only(result_set, _SET_FIELDS)}
    view["editorState"] = {
        **_without(shared_editor, frozenset({"legend"})),
        **_only(set_editor, _SET_EDITOR_FIELDS),
        "legend": {
            **_mapping(shared_editor.get("legend")),
            **_only(_mapping(set_editor.get("legend")), OTHER_MODE_RESULT_LEGEND_FIELDS),
        },
    }
    view["ui"] = {
        **_mapping(shared.get("ui")),
        **_only(_mapping(result_set.get("ui")), OTHER_MODE_RESULT_UI_FIELDS),
    }
    return view


def drawing_draft_config(fields: Mapping[str, Any], mode: DrawingMode) -> Mapping[str, Any]:
    """The Web draft config of ``mode`` in a Session or one-drawing view.

    Session 46 keeps it in the mode's slice (``modes[mode].config``; absent
    means that mode's defaults); Sessions 44 and older keep one flat
    ``config``.
    """

    modes = fields.get("modes")
    if isinstance(modes, Mapping):
        return _mapping(_mapping(modes.get(mode)).get("config"))
    return _mapping(fields.get("config"))


def session_with_stored_settings_replaced(
    data: Mapping[str, object],
    *,
    request: Callable[[StoredRequest, str], StoredRequest],
    draft_config: Callable[[DraftConfig, str], DraftConfig],
) -> dict[str, object]:
    """``data`` with each committed request and draft config as the callbacks
    return them. Each callback also takes the value's path in the document.

    The drawings' requests (``renderRequest``, ``otherModeResult``) and the
    drafts (each slice in ``modes``, also a slice without a drawing, or the flat
    ``config`` of a Session 44 or older). Every other field keeps its value, so
    this serves a document of any canonical version, where
    :func:`replace_session_drawings` writes the current layout.
    """

    changed = dict(data)
    stored = data.get("renderRequest")
    if isinstance(stored, Mapping):
        changed["renderRequest"] = request(stored, "renderRequest")
    other = data.get("otherModeResult")
    if isinstance(other, Mapping) and isinstance(other.get("renderRequest"), Mapping):
        changed["otherModeResult"] = {
            **other,
            "renderRequest": request(other["renderRequest"], "otherModeResult.renderRequest"),
        }
    modes = data.get("modes")
    if isinstance(modes, Mapping):
        slices = dict(modes)
        for mode in get_args(DrawingMode):
            config = drawing_draft_config(data, mode)
            if config:
                slices[mode] = {**modes[mode], "config": draft_config(config, f"modes.{mode}.config")}
        changed["modes"] = slices
    else:
        # Without slices both modes read the one flat draft.
        config = drawing_draft_config(data, "linear")
        if config:
            changed["config"] = draft_config(config, "config")
    return changed


def _request_mode(view: Mapping[str, Any]) -> str | None:
    request = view.get("renderRequest")
    mode = request.get("mode") if isinstance(request, Mapping) else None
    return mode if mode in DRAWING_NAMES else None


def _parts(
    mode: str,
    view: Mapping[str, Any],
    version: int,
    *,
    canonical: bool,
    identity: tuple[str, str] | None = None,
) -> SessionDrawingParts:
    from gbdraw.session_resources import canonical_resource_ids

    request = view.get("renderRequest") if canonical else None
    request = request if isinstance(request, Mapping) else None
    results = view.get("results")
    drawing_id, name = identity or (mode, DRAWING_NAMES[mode])
    return SessionDrawingParts(
        id=drawing_id,
        name=name,
        mode=cast(DrawingMode, mode),
        request=request,
        results=(
            tuple(item for item in results if isinstance(item, Mapping))
            if isinstance(results, list)
            else ()
        ),
        resource_ids=frozenset(canonical_resource_ids(request)),
        artifacts=SessionDrawingArtifacts(version, MappingProxyType(dict(view))),
    )


def _with_key(target: dict[str, Any], key: str, value: Mapping[str, Any], present: bool) -> None:
    if present or value:
        target[key] = dict(value)


def _drawing_view(project: Mapping[str, Any], drawing: Mapping[str, Any]) -> dict[str, Any]:
    """One Session 50 drawing as the one-drawing view that Sessions 46-49 give.

    The inverse of :func:`_drawing_of_view`.
    """

    mode = str(drawing["mode"])
    view = {**_only(project, PROJECT_FIELDS), **_only(drawing, _DRAWING_FIELDS)}
    editor = _mapping(drawing.get("editorState"))
    legend = _mapping(editor.get("legend"))
    if "editorState" in drawing:
        view["editorState"] = _without(editor, _SLICE_EDITOR_FIELDS | {"legend"})
        _with_key(view["editorState"], "legend", _without(legend, _SLICE_LEGEND_FIELDS), "legend" in editor)
    drawing_slice = _only(drawing, frozenset({"config", "features"}))
    slice_editor = _only(editor, _SLICE_EDITOR_FIELDS)
    _with_key(slice_editor, "legend", _only(legend, _SLICE_LEGEND_FIELDS), False)
    _with_key(drawing_slice, "editorState", slice_editor, False)
    drawing_ui = _mapping(drawing.get("ui"))
    _with_key(drawing_slice, "ui", _only(drawing_ui, _SLICE_UI_FIELDS), False)
    view["modes"] = {mode: drawing_slice}
    view["ui"] = {
        **_without(_mapping(project.get("ui")), frozenset({"activeDrawingId"})),
        **_without(drawing_ui, _SLICE_UI_FIELDS | {"inputType"}),
        "mode": mode,
    }
    if "inputType" in drawing_ui:
        view["ui"][_INPUT_TYPE_FIELDS[mode]] = drawing_ui["inputType"]
    web_files = {**_mapping(project.get("webFiles")), **_mapping(drawing.get("webFiles"))}
    _with_key(view, "webFiles", web_files, False)
    return view


def read_session_drawings(data: Mapping[str, Any]) -> tuple[SessionDrawingParts, ...]:
    """The drawings of a validated Session, in document order."""

    version = int(data["version"])
    if version >= DRAWINGS_SESSION_MIN_VERSION:
        return tuple(
            _parts(
                str(drawing["mode"]),
                _drawing_view(data, drawing),
                version,
                canonical=True,
                identity=(str(drawing["id"]), str(drawing["name"])),
            )
            for drawing in data["drawings"]
        )
    if version < CANONICAL_SESSION_MIN_VERSION:
        declared = session_mode(data)
        if declared not in DRAWING_NAMES:
            return ()
        assert declared is not None
        return (_parts(declared, data, version, canonical=False),)
    top_mode = _request_mode(data)
    if data.get("renderRequest") is None:
        shown = _mapping(data.get("ui")).get("mode")
        if not isinstance(shown, str) or shown not in DRAWING_NAMES:
            return ()
        return (_parts(shown, data, version, canonical=True),)
    if top_mode is None:
        return ()
    drawings = [_parts(top_mode, _without(data, frozenset({"otherModeResult"})), version, canonical=True)]
    other = data.get("otherModeResult")
    if isinstance(other, Mapping):
        view = _view(_shared(data), other)
        other_mode = _request_mode(view)
        if other_mode is not None:
            drawings.append(_parts(other_mode, view, version, canonical=True))
    return tuple(drawings)


def active_session_drawing(
    data: Mapping[str, Any],
    drawings: Sequence[SessionDrawingParts],
) -> SessionDrawingParts | None:
    """The drawing the Web app opens: ``ui.activeDrawingId``'s from Session 50,
    the shown mode's before it; else the first."""

    ui = _mapping(data.get("ui"))
    active_id = ui.get("activeDrawingId")
    shown = ui.get("mode")
    for drawing in drawings:
        by_id = int(data["version"]) >= DRAWINGS_SESSION_MIN_VERSION
        if (drawing.id == active_id) if by_id else (drawing.mode == shown):
            return drawing
    return drawings[0] if drawings else None


def allocate_drawing_id(mode: str, taken_ids: Sequence[str]) -> str:
    """The ID a new drawing of ``mode`` gets beside ``taken_ids``.

    The mode itself when free, else ``<mode>-<n>`` with the smallest free
    ``n >= 2`` (B2-D4). ``gbdraw/web/js/services/drawing-view.js``
    ``allocateDrawingId`` takes the same rule; both read the ``ids`` vectors
    of ``tests/fixtures/sessions/drawings-vectors.json``.
    """

    if mode not in DRAWING_NAMES:
        raise ValidationError(
            f"A drawing mode must be circular or linear, not {mode!r}.",
            diagnostic=_LAYOUT_INVALID,
        )
    taken = frozenset(taken_ids)
    if mode not in taken:
        return mode
    number = 2
    while f"{mode}-{number}" in taken:
        number += 1
    return f"{mode}-{number}"


def allocate_session_drawing(
    mode: str,
    *,
    requested_id: str | None,
    requested_name: str | None,
    taken_ids: Sequence[str],
) -> tuple[str, str]:
    """The ID and name a new drawing of ``mode`` gets beside ``taken_ids``.

    A drawing's ID comes from its kind (:func:`allocate_drawing_id`), so a
    requested ID must be that one. The name is the mode's name, followed by
    the ID's number for a second drawing of the mode (``Linear 2``).
    """

    drawing_id = allocate_drawing_id(mode, taken_ids)
    number = drawing_id[len(mode) + 1:]
    name = DRAWING_NAMES[mode] + (f" {number}" if number else "")
    if requested_id not in (None, drawing_id) or requested_name not in (None, name):
        raise ValidationError(
            f"Session {CURRENT_SESSION_VERSION} names each drawing by its mode: "
            f"this {mode} drawing has ID {drawing_id!r} and name {name!r}.",
            diagnostic=_LAYOUT_INVALID,
        )
    return drawing_id, name


def _pruned(value: Any) -> Any:
    if isinstance(value, Mapping):
        pruned = {key: _pruned(item) for key, item in value.items()}
        return {key: item for key, item in pruned.items() if item not in (None, {}, [])}
    return value


def _view_mode(view: Mapping[str, Any]) -> DrawingMode:
    """The mode of a one-drawing view: its request's, else its shown mode."""

    mode = _request_mode(view) or _mapping(view.get("ui")).get("mode")
    if mode not in DRAWING_NAMES:
        raise ValidationError("A drawing's view names no diagram mode.", diagnostic=_LAYOUT_INVALID)
    return cast(DrawingMode, mode)


def _project_of_view(view: Mapping[str, Any]) -> dict[str, Any]:
    """The Session 50 project fields of a one-drawing view, without ``drawings``."""

    project = _only(view, PROJECT_FIELDS)
    project["ui"] = _only(_mapping(view.get("ui")), PROJECT_UI_FIELDS)
    _with_key(project, "webFiles", _only(_mapping(view.get("webFiles")), PROJECT_WEB_FILES_FIELDS), False)
    return project


# A one-drawing view holds no ``otherModeResult``: the 46-49 readers split it
# off before a view is built, and Session 50 has no place for it.
_VIEW_FIELDS = PROJECT_FIELDS | _DRAWING_FIELDS | {"editorState", "ui", "webFiles", "modes"}
_DRAWING_WEB_FILE_KEYS = PROJECT_WEB_FILES_FIELDS | {"bindings"} | frozenset().union(*_DRAWING_WEB_FILES_FIELDS.values())


def other_mode_fields(state: Mapping[str, Any], mode: str) -> list[str]:
    """The fields of a new drawing's ``state`` that belong to the other mode:
    its ``modes`` slice, input type, Web file bindings and request metadata,
    which :func:`_drawing_of_view` leaves out (D-15-125: a new build refuses
    them)."""

    (other,) = set(DRAWING_NAMES) - {mode}
    web_files = _mapping(state.get("webFiles"))
    return [
        *(f"modes.{key}" for key in sorted(set(_mapping(state.get("modes"))) - {mode})),
        *(f"ui.{key}" for key in (_INPUT_TYPE_FIELDS[other],) if key in _mapping(state.get("ui"))),
        *(f"webFiles.{key}" for key in sorted(_DRAWING_WEB_FILES_FIELDS[other] & set(web_files))),
        *(
            f"webFiles.bindings.{key}"
            for key in sorted(WEB_BINDING_MODE_SLOTS[other] & set(_mapping(web_files.get("bindings"))))
        ),
    ]


def _drawing_of_view(view: Mapping[str, Any], drawing_id: str, name: str) -> dict[str, Any]:
    """One Session 50 drawing from its one-drawing view; the inverse of :func:`_drawing_view`.

    The drawing keeps its own mode's slice, input slots and request metadata;
    the other mode's are not this drawing's and are left out.
    """

    mode = _view_mode(view)
    web_files = _mapping(view.get("webFiles"))
    unplaced = sorted({*(set(view) - _VIEW_FIELDS), *(f"webFiles.{key}" for key in set(web_files) - _DRAWING_WEB_FILE_KEYS)})
    if unplaced:
        raise ValidationError(
            "Session 50 has no place for " + ", ".join(unplaced) + ".",
            diagnostic=_LAYOUT_INVALID,
        )
    drawing_slice = _mapping(_mapping(view.get("modes")).get(mode))
    drawing: dict[str, Any] = {
        "id": drawing_id,
        "name": name,
        "mode": mode,
        **_only(view, _DRAWING_FIELDS),
        **_only(drawing_slice, frozenset({"config", "features"})),
    }
    editor = _mapping(view.get("editorState"))
    slice_editor = _mapping(drawing_slice.get("editorState"))
    if "editorState" in view or "editorState" in drawing_slice:
        merged = {**_without(editor, frozenset({"legend"})), **_without(slice_editor, frozenset({"legend"}))}
        _with_key(
            merged,
            "legend",
            {**_mapping(editor.get("legend")), **_mapping(slice_editor.get("legend"))},
            "legend" in editor or "legend" in slice_editor,
        )
        drawing["editorState"] = merged
    view_ui = _mapping(view.get("ui"))
    drawing["ui"] = {
        **_without(view_ui, PROJECT_UI_FIELDS | {"mode", "activeDrawingId", *_INPUT_TYPE_FIELDS.values()}),
        **_mapping(drawing_slice.get("ui")),
    }
    if _INPUT_TYPE_FIELDS[mode] in view_ui:
        drawing["ui"]["inputType"] = view_ui[_INPUT_TYPE_FIELDS[mode]]
    drawing_web_files = _only(web_files, _DRAWING_WEB_FILES_FIELDS[mode])
    if "bindings" in web_files:
        drawing_web_files["bindings"] = _only(
            _mapping(web_files["bindings"]),
            frozenset({"schema"}) | WEB_BINDING_CARD_SLOTS | WEB_BINDING_MODE_SLOTS[mode],
        )
    _with_key(drawing, "webFiles", drawing_web_files, False)
    return drawing


def _identity_problem(ids: Sequence[str], names: Sequence[str]) -> str | None:
    """Why ``ids`` and ``names`` cannot name the drawings of a Session 50, if they cannot."""

    if len(set(ids)) != len(ids):
        return "two drawings share an ID"
    folded = [name.casefold() for name in names]
    if "" in (name.strip() for name in names) or len(set(folded)) != len(folded):
        return "each drawing needs a name of its own"
    return None



def _unreferenced_resources(project: Mapping[str, Any], drawings: Sequence[Mapping[str, Any]]) -> list[str]:
    """The project's resources that no drawing names."""

    from gbdraw.session_resources import drawing_resource_ids

    referenced = set().union(*(drawing_resource_ids(drawing) for drawing in drawings))
    return sorted(set(_mapping(project.get("resources"))) - referenced)


def _layout_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic=_LAYOUT_INVALID)


def _unknown_keys(value: Mapping[str, Any], allowed: frozenset[str], where: str) -> None:
    unknown = sorted(str(key) for key in value if key not in allowed)
    if unknown:
        raise _layout_error(f"Session 50 {where} contains unclassified field(s): " + ", ".join(unknown) + ".")


def _checked_drawing(drawing: object) -> SessionDrawing:
    """One entry of ``drawings``, with its own keys only (``SessionDrawing``)."""

    if not isinstance(drawing, Mapping):
        raise _layout_error("Each Session 50 drawing must be an object.")
    _unknown_keys(drawing, keys_of(SessionDrawing), "drawing")
    missing = sorted(SessionDrawing.__required_keys__ - set(drawing))
    if missing:
        raise _layout_error("A Session 50 drawing requires " + ", ".join(missing) + ".")
    drawing_id, mode = drawing["id"], drawing["mode"]
    if (not isinstance(drawing_id, str) or len(drawing_id) > DRAWING_ID_MAX_LENGTH
            or not SESSION_ID_RE.fullmatch(drawing_id)):
        raise _layout_error(f"Invalid Session 50 drawing ID: {drawing_id!r}.")
    if not isinstance(drawing["name"], str) or not isinstance(mode, str) or mode not in DRAWING_NAMES:
        raise _layout_error(f"Session 50 drawing {drawing_id} needs a name and a diagram mode.")
    if drawing["renderRequest"] is not None and _request_mode(drawing) != mode:
        raise _layout_error(f"Session 50 drawing {drawing_id} has a request of another mode than {mode}.")
    ui = drawing.get("ui", {})
    if not isinstance(ui, Mapping):
        raise _layout_error(f"Session 50 drawing {drawing_id} has a ui that is not an object.")
    _unknown_keys(ui, frozenset(ui) - PROJECT_UI_FIELDS - {"mode", "activeDrawingId", *_INPUT_TYPE_FIELDS.values()},
                  f"drawing {drawing_id} ui")
    web_files = drawing.get("webFiles", {})
    if not isinstance(web_files, Mapping):
        raise _layout_error(f"Session 50 drawing {drawing_id} has webFiles that are not an object.")
    _unknown_keys(web_files, frozenset({"bindings"}) | _DRAWING_WEB_FILES_FIELDS[mode], f"drawing {drawing_id} webFiles")
    if drawing["renderRequest"] is None and (
        drawing["results"] != [] or drawing.get("runMetadata")
        or _mapping(drawing.get("editorState")).get("featureCatalog") is not None
    ):
        raise _layout_error(f"Session 50 drawing {drawing_id} has Results without a request.")
    return cast(SessionDrawing, drawing)


def checked_drawing_views(data: Mapping[str, Any]) -> tuple[dict[str, Any], ...]:
    """The one-drawing view of each drawing of a Session 50, after checking
    the layout: the closed project (``GbdrawSession50``), its drawings
    (``SessionDrawing``) with unique IDs and names, ``ui.activeDrawingId``,
    and that each resource is some drawing's. ``validate_session`` checks
    each view as it checks a Session 46-49."""

    _unknown_keys(data, keys_of(GbdrawSession50), "project")
    missing = sorted(GbdrawSession50.__required_keys__ - set(data))
    if missing:
        raise _layout_error("Session 50 requires " + ", ".join(missing) + ".")
    ui, web_files, entries = data["ui"], data.get("webFiles", {}), data["drawings"]
    if not isinstance(ui, Mapping) or not isinstance(web_files, Mapping):
        raise _layout_error("Session 50 ui and webFiles must be objects.")
    _unknown_keys(ui, keys_of(SessionProjectUi), "ui")
    _unknown_keys(web_files, keys_of(SessionProjectWebFiles), "webFiles")
    if not isinstance(entries, list) or not entries:
        raise _layout_error("Session 50 requires a non-empty drawings list.")
    drawings = [_checked_drawing(entry) for entry in entries]
    problem = _identity_problem([drawing["id"] for drawing in drawings], [drawing["name"] for drawing in drawings])
    if problem:
        raise _layout_error(f"Session 50 cannot hold these drawings: {problem}.")
    if ui.get("activeDrawingId") not in {drawing["id"] for drawing in drawings}:
        raise _layout_error("Session 50 ui.activeDrawingId must name one of its drawings.")
    if not isinstance(data["resources"], Mapping):
        raise _layout_error("Session 50 requires a canonical resources object.")
    unreferenced = _unreferenced_resources(data, [cast(Mapping[str, Any], drawing) for drawing in drawings])
    if unreferenced:
        raise _layout_error("No drawing names the resources " + ", ".join(unreferenced) + ".")
    return tuple(_drawing_view(data, cast(Mapping[str, Any], drawing)) for drawing in drawings)

def _drawings_document(
    views: Sequence[Mapping[str, Any]],
    *,
    ids: Sequence[str],
    names: Sequence[str],
    active: int,
    project_views: Sequence[Mapping[str, Any]],
    drop_unreferenced: bool = False,
) -> dict[str, Any]:
    """The Session 50 document of ``views``, one one-drawing view per drawing.

    ``project_views`` give the project fields: the first gives them, and the
    others may repeat but not change them. Each drawing keeps the LOSAT
    entries its own view carries, so an entry two drawings use is in both.
    A resource that no drawing names is refused, or left out with
    ``drop_unreferenced``.
    """

    if not views or len(ids) != len(views) or len(names) != len(views):
        raise ValidationError("Each drawing of a Session needs one ID and one name.", diagnostic=_LAYOUT_INVALID)
    problem = _identity_problem(ids, names)
    if problem:
        raise ValidationError(f"Session 50 cannot hold these drawings: {problem}.", diagnostic=_LAYOUT_INVALID)
    project = _project_of_view(project_views[0])
    envelope = frozenset({"format", "version", "createdAt"})
    first = _pruned(_without(project, envelope))
    for view in project_views[1:]:
        changed = sorted(
            key
            for key, value in _pruned(_without(_project_of_view(view), envelope)).items()
            if first.get(key) != value
        )
        if changed:
            raise ValidationError(
                "Session 50 drawings share the project's " + ", ".join(changed) + "; give them one value.",
                diagnostic=_LAYOUT_INVALID,
            )
    drawings = [_drawing_of_view(view, drawing_id, name) for view, drawing_id, name in zip(views, ids, names)]
    unreferenced = _unreferenced_resources(project, drawings)
    if unreferenced and drop_unreferenced:
        project["resources"] = _without(_mapping(project["resources"]), frozenset(unreferenced))
    elif unreferenced:
        raise ValidationError(
            "No drawing names the resources " + ", ".join(unreferenced) + "; leave them out of the Session.",
            diagnostic=_LAYOUT_INVALID,
        )
    project["version"] = DRAWINGS_SESSION_MIN_VERSION
    project["ui"]["activeDrawingId"] = ids[active]
    project["drawings"] = drawings
    return project


def write_session_drawings(
    views: Sequence[Mapping[str, Any]],
    *,
    ids: Sequence[str] | None = None,
    names: Sequence[str] | None = None,
    active: int | None = None,
) -> dict[str, Any]:
    """The Session 50 document fields that hold ``views``, one one-drawing view per drawing.

    The drawings are listed in ``drawings`` in the order of ``views``, with
    ``ids`` and ``names`` (allocated from each mode when omitted,
    :func:`allocate_session_drawing`); ``active`` is the index of
    ``ui.activeDrawingId`` (the first drawing when omitted).
    """

    if not views:
        raise ValidationError("A Session holds at least one drawing.", diagnostic=_LAYOUT_INVALID)
    if ids is None and names is None:
        ids, names = [], []
        for view in views:
            drawing_id, name = allocate_session_drawing(
                _view_mode(view), requested_id=None, requested_name=None, taken_ids=ids
            )
            ids.append(drawing_id)
            names.append(name)
    if ids is None or names is None:
        raise ValidationError("Give both the IDs and the names of the drawings, or neither.", diagnostic=_LAYOUT_INVALID)
    return _drawings_document(views, ids=ids, names=names, active=active or 0, project_views=views)


def replace_session_drawings(
    data: Mapping[str, Any],
    views: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """The Session 50 ``data`` with the drawings named in ``views`` replaced, in place and order.

    Each view is one drawing seen as a one-drawing Session, of the same mode.
    The other drawings stay as they are, and the replaced views give the
    project fields. A Session 46-49 moves to drawings before it is saved
    again (:func:`drawings_document_from_mode_scoped`).
    """

    if int(data["version"]) < DRAWINGS_SESSION_MIN_VERSION:
        raise ValidationError(
            "Only a Session 50 has drawings to replace; move this one to drawings first.",
            diagnostic=_LAYOUT_INVALID,
        )
    drawings = read_session_drawings(data)
    unknown = sorted(set(views) - {drawing.id for drawing in drawings})
    if unknown:
        raise ValidationError("Session has no drawing " + ", ".join(unknown) + " to replace.", diagnostic=_LAYOUT_INVALID)
    replaced = [views[drawing.id] for drawing in drawings if drawing.id in views]
    for drawing in drawings:
        if drawing.id in views and _view_mode(views[drawing.id]) != drawing.mode:
            raise ValidationError(f"The {drawing.id} drawing keeps its {drawing.mode} mode.", diagnostic=_LAYOUT_INVALID)
    active = active_session_drawing(data, drawings)
    return _drawings_document(
        [views.get(drawing.id, drawing.artifacts.fields) for drawing in drawings],
        ids=[drawing.id for drawing in drawings],
        names=[drawing.name for drawing in drawings],
        active=drawings.index(active) if active is not None else 0,
        project_views=replaced or [drawings[0].artifacts.fields],
    )


class ModeScopedSourceFacts(TypedDict):
    """What a Session 46-49 held that its Session 50 cannot say.

    These are the only facts that carry a source version past
    :func:`drawings_document_from_mode_scoped`: ``cardRowsAreFileRows`` (a
    Session 48 without a Web draft: the file cards' rows are file rows), and
    ``proteinResultsInDrawnFrame`` (a Session 46 or 47: stored protein
    comparison results are in the drawn frame of their records).
    """

    cardRowsAreFileRows: bool
    proteinResultsInDrawnFrame: bool


class _ModeSliceDefaultRow(TypedDict):
    domain: str
    path: str


class ModeSliceDefault(_ModeSliceDefaultRow, total=False):
    """One registry row's default in one mode; no ``value`` when Save leaves the row out at its default."""

    value: object


class UnarrangedEntries(TypedDict):
    """A list that arranges nothing: each entry holds exactly ``entryKeys``, ``positionKey`` its 1-based position."""

    entryKeys: list[str]
    positionKey: str


class ModeSliceNormalizeRule(TypedDict):
    """A row field the Web fills with run-time identities: a reader empties it when it arranges nothing."""

    domain: str
    path: str
    field: str
    emptyWhenUnarranged: UnarrangedEntries
    reason: str


class ModeSliceDefaults(TypedDict):
    """``gbdraw/data/mode-slice-defaults.json``: each mode's default settings for rule d of the move to drawings.

    The Web owns these defaults (a new drawing of ``createDrawingState``);
    ``tools/generate_mode_slice_defaults.mjs`` writes them and checks them
    with ``--check``, so Python reads the Web's generated defaults.
    """

    schema: int
    description: str
    registryRevision: int
    normalize: list[ModeSliceNormalizeRule]
    earlierDefaults: list[ModeSliceDefault]
    defaults: dict[DrawingMode, list[ModeSliceDefault]]


_ABSENT = object()


@cache
def _mode_slice_defaults() -> ModeSliceDefaults:
    text = resources.files("gbdraw.data").joinpath("mode-slice-defaults.json").read_text(encoding="utf-8")
    return cast(ModeSliceDefaults, json.loads(text))


def _normalized(row: ModeScopedSetting, value: object) -> object:
    """``value`` without the run-time identities the packaged rules name for ``row``."""

    for rule in _mode_slice_defaults()["normalize"]:
        if (rule["domain"], rule["path"]) != (row.domain, row.path) or not isinstance(value, Mapping):
            continue
        entries = value.get(rule["field"])
        shape = rule["emptyWhenUnarranged"]
        keys = sorted(shape["entryKeys"])
        if isinstance(entries, list) and all(
            isinstance(entry, Mapping)
            and sorted(entry) == keys
            and not isinstance(entry[shape["positionKey"]], bool)
            and entry[shape["positionKey"]] == position
            for position, entry in enumerate(entries, start=1)
        ):
            value = {**value, rule["field"]: []}
    return value


def _default_value(mode: DrawingMode, row: ModeScopedSetting) -> object:
    for default in _mode_slice_defaults()["defaults"][mode]:
        if (default["domain"], default["path"]) == (row.domain, row.path):
            return default.get("value", _ABSENT)
    return _ABSENT


def _is_earlier_default(row: ModeScopedSetting, value: object) -> bool:
    """Whether ``value`` is a default of ``row`` that a released writer stored."""

    return any(
        (earlier["domain"], earlier["path"]) == (row.domain, row.path)
        and _canonical(earlier.get("value")) == _canonical(value)
        for earlier in _mode_slice_defaults()["earlierDefaults"]
    )


def _row_value(mode_slice: Mapping[str, object], row: ModeScopedSetting) -> object:
    holder: object = mode_slice
    for part in row.domain.split("."):
        holder = holder.get(part, _ABSENT) if isinstance(holder, Mapping) else _ABSENT
    value = holder.get(row.path, _ABSENT) if isinstance(holder, Mapping) else _ABSENT
    return value if value is _ABSENT else _normalized(row, value)


def _json_number(value: object) -> object:
    """``value`` with each integral float as an int: JSON numbers have one
    kind, so ``1.0`` a Python writer wrote is the Web's ``1`` (OV-1668)."""

    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, Mapping):
        return {key: _json_number(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_number(item) for item in value]
    return value


def _canonical(value: object) -> str:
    return json.dumps(_json_number(value), sort_keys=True, separators=(",", ":"))


def _keyed_entries(value: object, key: str) -> dict[str, str]:
    """A keyed-row value's entries, canonical, by entry key (registry ``key``)."""

    if isinstance(value, Mapping):
        return {str(entry_key): _canonical(entry) for entry_key, entry in value.items()}
    fields = key[len("JSON["):-1].split(",") if key.startswith("JSON[") else [key]
    entries: dict[str, str] = {}
    for entry in value if isinstance(value, list) else ():
        identity = [entry.get(field) for field in fields] if isinstance(entry, Mapping) else [entry]
        entries[_canonical(identity)] = _canonical(entry)
    return entries


def _used_by_settings(
    slices: Mapping[str, object],
    mode: DrawingMode,
    other: DrawingMode,
) -> bool:
    """Rule d: a registry row of ``mode``'s slice that differs from ``mode``'s
    default (the Web's, :func:`_mode_slice_defaults`) and, unless the row is
    ``mode``'s alone, from ``other``'s value (entry by entry for a keyed row).
    An absent row is the default. Both values are compared without the
    run-time identities of the packaged ``normalize`` rules."""

    own = _mapping(slices.get(mode))
    shown = _mapping(slices.get(other))
    for row in MODE_SCOPED_SETTINGS:
        if row.modes not in ("both", mode):
            continue
        value = _row_value(own, row)
        if value is _ABSENT:
            continue
        default = _default_value(mode, row)
        if default is not _ABSENT and _canonical(value) == _canonical(default):
            continue
        if _is_earlier_default(row, value):
            continue
        if row.modes == mode:
            return True
        shown_value = _row_value(shown, row)
        if shown_value is _ABSENT:
            return True
        if row.key is None:
            if _canonical(value) != _canonical(shown_value):
                return True
            continue
        shown_entries = _keyed_entries(shown_value, row.key)
        if any(shown_entries.get(entry_key) != entry for entry_key, entry in _keyed_entries(value, row.key).items()):
            return True
    return False


def _binds_input(bindings: Mapping[str, object], mode: DrawingMode) -> bool:
    """Rule c: a binding of ``mode``'s input slots names a file (a Depth file
    alone counts; an empty ``linearSeqs`` row does not)."""

    from gbdraw.session_resources import canonical_resource_ids

    return any(canonical_resource_ids(bindings.get(slot)) for slot in WEB_BINDING_MODE_SLOTS[mode])


def _undrawn_view(data: Mapping[str, Any], mode: DrawingMode) -> dict[str, Any]:
    """The one-drawing view of ``mode`` in a Session 46-49 that holds no Result of it."""

    view = _shared(data)
    view["renderRequest"] = None
    view["results"] = []
    view["editorState"] = {**view["editorState"], "featureCatalog": None}
    view["ui"] = {**view["ui"], "mode": mode}
    return view


def _migrated_drawing_modes(
    shown: DrawingMode,
    other: DrawingMode,
    *,
    other_used: bool,
    result_modes: frozenset[DrawingMode],
) -> list[DrawingMode]:
    """The modes of a migrated Session's drawings, the drawing that opens first.

    The mode shown at Save always gets a drawing; ``other`` gets one when it
    is used. The shown drawing opens unless it has no stored Result and the
    other drawing has one: then the drawing with the Result opens (the
    figure the released app showed), and the drawing that opens comes
    first. The Web's ``migratedDrawingModes`` is the twin.
    """

    if not other_used:
        return [shown]
    if other in result_modes and shown not in result_modes:
        return [other, shown]
    return [shown, other]


def drawings_document_from_mode_scoped(
    data: Mapping[str, object],
    *,
    stored_result_modes: frozenset[str] = frozenset(),
) -> tuple[GbdrawSession50, ModeScopedSourceFacts]:
    """The Session 50 of a validated Session 46-49, and what only the source says.

    The rules of the Product Contract (PD-OI-090, the move of an earlier
    Session) as the design lists them: each mode's committed set goes to its
    mode's drawing (:func:`read_session_drawings`); the drawing of ``ui.mode``
    always gets a drawing; the other mode's drawing is made only
    when it is used: it has a Result, a non-null input binding of its mode,
    or a registry row that differs from its default and from the shown
    mode's value (:func:`_used_by_settings`, with the Web's packaged
    defaults). ``stored_result_modes`` names the modes whose Results the
    caller set aside (Web Load of a Session 27-39): they count as Results.
    The order and the drawing that opens are :func:`_migrated_drawing_modes`'s.
    Each drawing keeps its own mode's slice, input
    bindings, request metadata and the file cards' bindings, and a copy of
    the LOSAT artifacts (:func:`_drawing_of_view`); the IDs and names are the
    modes'. A Session 46 gets its Legend edits keyed by row key first
    (:func:`gbdraw.session_io.session_legend_edits_by_row_key`). Resources no
    drawing names are left out. The Web's ``upgradeModeScopedSession`` is the
    twin; ``tests/fixtures/sessions/drawings-vectors.json`` pins both.
    """

    version = data.get("version")
    if not isinstance(version, int) or not MODE_SCOPED_SESSION_MIN_VERSION <= version < DRAWINGS_SESSION_MIN_VERSION:
        raise ValidationError(
            f"Only a Session {MODE_SCOPED_SESSION_MIN_VERSION}-{MODE_SCOPED_SESSION_MAX_VERSION} moves to drawings; "
            f"this one is {version!r}.",
            diagnostic=_LAYOUT_INVALID,
        )
    source: Mapping[str, Any] = data
    if version < ROW_KEYED_LEGEND_SESSION_MIN_VERSION:
        source, _outcomes = session_legend_edits_by_row_key(data)
    committed = {drawing.mode: dict(drawing.artifacts.fields) for drawing in read_session_drawings(source)}
    if not committed:
        raise ValidationError("The Session names no drawing to move.", diagnostic=_LAYOUT_INVALID)
    shown_mode = _mapping(source.get("ui")).get("mode")
    shown = cast(DrawingMode, shown_mode) if shown_mode in DRAWING_NAMES else next(iter(committed))
    other: DrawingMode = "linear" if shown == "circular" else "circular"
    slices = _mapping(source.get("modes"))
    result_modes = frozenset(
        mode
        for mode in (shown, other)
        if mode in stored_result_modes or _mapping(committed.get(mode)).get("renderRequest") is not None
    )
    other_used = (
        other in committed
        or other in result_modes
        or _binds_input(_mapping(_mapping(source.get("webFiles")).get("bindings")), other)
        or (other in slices and _used_by_settings(slices, other, shown))
    )
    modes = _migrated_drawing_modes(shown, other, other_used=other_used, result_modes=result_modes)
    views = [committed.get(mode) or _undrawn_view(source, mode) for mode in modes]
    document = _drawings_document(
        views,
        ids=modes,
        names=[DRAWING_NAMES[mode] for mode in modes],
        active=0,
        project_views=views,
        drop_unreferenced=True,
    )
    has_draft = any(isinstance(_mapping(value).get("config"), Mapping) for value in slices.values())
    facts: ModeScopedSourceFacts = {
        "cardRowsAreFileRows": version >= FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION and not has_draft,
        "proteinResultsInDrawnFrame": version < FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION,
    }
    return cast(GbdrawSession50, document), facts


class SessionDrawingsLayout(TypedDict):
    """:func:`describe_session_drawings_layout`: the generated ``SESSION_DRAWINGS_LAYOUT``."""

    version: int
    drawingNames: dict[str, str]
    setFields: list[str]
    setEditorFields: list[str]
    projectFields: list[str]
    projectUiFields: list[str]
    projectWebFilesFields: list[str]
    drawingFields: list[str]
    viewFields: list[str]
    drawingWebFileKeys: list[str]
    drawingWebFilesFields: dict[str, list[str]]
    inputTypeFields: dict[str, str]
    bindingCardSlots: list[str]
    bindingModeSlots: dict[str, list[str]]
    resourceReferenceFields: list[str]
    resourceIdListFields: list[str]
    drawingResourceFields: list[str]


def describe_session_drawings_layout() -> SessionDrawingsLayout:
    """The field sets of the move of a Session 46-49 to drawings, for the Web.

    The Web's twin of :func:`drawings_document_from_mode_scoped`
    (``upgradeModeScopedSession``, ``gbdraw/web/js/services/session-drawings-upgrade.js``)
    reads them from the generated ``SESSION_DRAWINGS_LAYOUT``; this module and
    :mod:`gbdraw.session_resources` stay their owners. Sets are sorted lists.
    """

    from gbdraw.session_resources import (
        _DRAWING_RESOURCE_FIELDS,
        _RESOURCE_ID_LIST_FIELDS,
        _RESOURCE_REFERENCE_FIELDS,
    )

    def per_mode(fields: Mapping[str, frozenset[str]]) -> dict[str, list[str]]:
        return {mode: sorted(fields[mode]) for mode in DRAWING_NAMES}

    return {
        "version": DRAWINGS_SESSION_MIN_VERSION,
        "drawingNames": dict(DRAWING_NAMES),
        "setFields": sorted(_SET_FIELDS),
        "setEditorFields": sorted(_SET_EDITOR_FIELDS),
        "projectFields": sorted(PROJECT_FIELDS),
        "projectUiFields": sorted(PROJECT_UI_FIELDS),
        "projectWebFilesFields": sorted(PROJECT_WEB_FILES_FIELDS),
        "drawingFields": sorted(_DRAWING_FIELDS),
        "viewFields": sorted(_VIEW_FIELDS),
        "drawingWebFileKeys": sorted(_DRAWING_WEB_FILE_KEYS),
        "drawingWebFilesFields": per_mode(_DRAWING_WEB_FILES_FIELDS),
        "inputTypeFields": dict(_INPUT_TYPE_FIELDS),
        "bindingCardSlots": sorted(WEB_BINDING_CARD_SLOTS),
        "bindingModeSlots": per_mode(WEB_BINDING_MODE_SLOTS),
        "resourceReferenceFields": sorted(_RESOURCE_REFERENCE_FIELDS),
        "resourceIdListFields": list(_RESOURCE_ID_LIST_FIELDS),
        "drawingResourceFields": list(_DRAWING_RESOURCE_FIELDS),
    }


__all__ = [
    "DRAWING_NAMES",
    "PROJECT_FIELDS",
    "PROJECT_UI_FIELDS",
    "PROJECT_WEB_FILES_FIELDS",
    "DraftConfig",
    "DrawingMode",
    "SessionDrawingArtifacts",
    "SessionDrawingParts",
    "StoredRequest",
    "active_session_drawing",
    "allocate_drawing_id",
    "allocate_session_drawing",
    "ModeScopedSourceFacts",
    "SessionDrawingsLayout",
    "checked_drawing_views",
    "describe_session_drawings_layout",
    "drawing_draft_config",
    "drawings_document_from_mode_scoped",
    "other_mode_fields",
    "read_session_drawings",
    "replace_session_drawings",
    "session_with_stored_settings_replaced",
    "write_session_drawings",
]
