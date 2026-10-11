"""Where each drawing of a Session lives: the one owner of the document layout.

A drawing is one diagram of a Session: its mode, its committed request and
Results, and the editor state that belongs to it. Selection, rendering, the
CLI, output naming and the Session writers reach a drawing's parts only
through this module, so a new document layout changes this module alone.

Session 46, the current writer, holds at most one drawing of each mode:

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
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Literal, Mapping, Sequence, cast

from gbdraw.exceptions import ValidationError
from gbdraw.session_io import (
    CANONICAL_SESSION_MIN_VERSION,
    CURRENT_SESSION_VERSION,
    OTHER_MODE_RESULT_EDITOR_FIELDS,
    OTHER_MODE_RESULT_LEGEND_FIELDS,
    OTHER_MODE_RESULT_UI_FIELDS,
    session_mode,
)

DrawingMode = Literal["circular", "linear"]
DRAWING_NAMES: Mapping[str, str] = MappingProxyType(
    {"circular": "Circular", "linear": "Linear"}
)

# The Web reader names a document it cannot lay out as a Session field error.
_LAYOUT_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}
# The fields of one committed Result set, named as at the top level.
_SET_FIELDS = frozenset({"renderRequest", "results", "runMetadata", "cliInvocation"})
_SET_EDITOR_FIELDS = OTHER_MODE_RESULT_EDITOR_FIELDS - {"legend"}


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
    return {key: source[key] for key in fields if key in source}


def _without(source: Mapping[str, Any], fields: frozenset[str]) -> dict[str, Any]:
    return {key: value for key, value in source.items() if key not in fields}


def _result_set(view: Mapping[str, Any]) -> dict[str, Any]:
    """The committed Result set of a one-drawing view, as ``otherModeResult`` holds it."""

    editor = _mapping(view.get("editorState"))
    result = _only(view, _SET_FIELDS)
    result["editorState"] = {
        **_only(editor, _SET_EDITOR_FIELDS),
        "legend": _only(_mapping(editor.get("legend")), OTHER_MODE_RESULT_LEGEND_FIELDS),
    }
    result["ui"] = _only(_mapping(view.get("ui")), OTHER_MODE_RESULT_UI_FIELDS)
    return result


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
) -> SessionDrawingParts:
    from gbdraw.session_resources import canonical_resource_ids

    request = view.get("renderRequest") if canonical else None
    request = request if isinstance(request, Mapping) else None
    results = view.get("results")
    return SessionDrawingParts(
        id=mode,
        name=DRAWING_NAMES[mode],
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


def read_session_drawings(data: Mapping[str, Any]) -> tuple[SessionDrawingParts, ...]:
    """The drawings of a validated Session, in document order."""

    version = int(data["version"])
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
    """The drawing the Web app opens: the shown mode's, else the first."""

    shown = _mapping(data.get("ui")).get("mode")
    for drawing in drawings:
        if drawing.mode == shown:
            return drawing
    return drawings[0] if drawings else None


def allocate_session_drawing(
    mode: str,
    *,
    requested_id: str | None,
    requested_name: str | None,
    taken_ids: Sequence[str],
) -> tuple[str, str]:
    """The ID and name a new drawing of ``mode`` gets in the current layout."""

    if mode not in DRAWING_NAMES:
        raise ValidationError(
            f"A drawing mode must be circular or linear, not {mode!r}.",
            diagnostic=_LAYOUT_INVALID,
        )
    name = DRAWING_NAMES[mode]
    if requested_id not in (None, mode) or requested_name not in (None, name):
        raise ValidationError(
            f"Session {CURRENT_SESSION_VERSION} names each drawing by its mode: "
            f"a {mode} drawing has ID {mode!r} and name {name!r}.",
            diagnostic=_LAYOUT_INVALID,
        )
    if mode in taken_ids:
        raise ValidationError(
            f"Session {CURRENT_SESSION_VERSION} holds at most one drawing of each mode; "
            f"it already has a {mode} drawing.",
            diagnostic=_LAYOUT_INVALID,
        )
    return mode, name


_ARTIFACT_ENTRY_FIELDS = ("losatCache", "losatDerivedCache")
_LEGACY_ARTIFACT_FIELDS = ("proteinRawCandidates", "proteinDerivedEvidence")


def _pruned(value: Any) -> Any:
    if isinstance(value, Mapping):
        pruned = {key: _pruned(item) for key, item in value.items()}
        return {key: item for key, item in pruned.items() if item not in (None, {}, [])}
    return value


def _merged_entries(groups: Sequence[Any]) -> list[Any]:
    merged: list[Any] = []
    for entries in groups:
        for entry in entries if isinstance(entries, list) else ():
            if entry not in merged:
                merged.append(entry)
    return merged


def _without_slices(shared: Mapping[str, Any], modes: frozenset[str]) -> dict[str, Any]:
    """``shared`` without the mode slices of ``modes``."""

    slices = shared.get("modes")
    if not isinstance(slices, Mapping):
        return dict(shared)
    return {**shared, "modes": _without(slices, modes)}


def _merged_slices(
    shared_slices: Any,
    owners: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """``shared_slices`` with each owner view's own mode slice, in place."""

    slices = dict(_mapping(shared_slices))
    for view in owners:
        mode = _request_mode(view)
        if mode is None:
            continue
        own = _mapping(view.get("modes")).get(mode)
        if own is None:
            slices.pop(mode, None)
        else:
            slices[mode] = own
    return slices


def _merged_shared(
    views: Sequence[Mapping[str, Any]],
    *,
    kept: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """The fields that the drawings of a Session 46 share, from their views.

    Each drawing owns its mode's slice in ``modes``: it comes from that
    drawing's view, or from ``kept`` for a drawing that is not replaced. The
    first view gives every other shared field, and the others may repeat but
    not change it, except the LOSAT artifacts: their entries are merged
    without repeats, and the protein identity manifest is the one non-empty
    manifest among the views.
    """

    owners = [*views, *kept]
    owned = frozenset(mode for view in owners if (mode := _request_mode(view)) is not None)
    shared = _shared(views[0])
    artifact_fields = {*_ARTIFACT_ENTRY_FIELDS, "legacyArtifacts", "proteinIdentityManifest"}
    first = _pruned(_without_slices(_without(shared, frozenset(artifact_fields)), owned))
    for view in views[1:]:
        changed = sorted(
            key
            for key, value in _pruned(
                _without_slices(_without(_shared(view), frozenset(artifact_fields)), owned)
            ).items()
            if first.get(key) != value
        )
        if changed:
            raise ValidationError(
                f"Session {CURRENT_SESSION_VERSION} drawings share "
                + ", ".join(changed)
                + "; give these fields to one drawing only.",
                diagnostic=_LAYOUT_INVALID,
            )
    if len(owners) > 1:
        slices = _merged_slices(shared.get("modes"), owners)
        if slices or "modes" in shared:
            shared["modes"] = slices
    if len(views) == 1:
        return shared
    for field_name in _ARTIFACT_ENTRY_FIELDS:
        if any(field_name in view for view in views):
            shared[field_name] = {
                "entries": _merged_entries(
                    [_mapping(view.get(field_name)).get("entries") for view in views]
                )
            }
    legacy: dict[str, Any] = {}
    for field_name in _LEGACY_ARTIFACT_FIELDS:
        envelopes = [
            _mapping(_mapping(view.get("legacyArtifacts")).get(field_name)) for view in views
        ]
        entries = _merged_entries([envelope.get("entries") for envelope in envelopes])
        if entries:
            legacy[field_name] = {"schema": 1, "entries": entries}
    if legacy:
        shared["legacyArtifacts"] = legacy
    else:
        shared.pop("legacyArtifacts", None)
    manifests = []
    for view in views:
        manifest = view.get("proteinIdentityManifest")
        if (
            isinstance(manifest, Mapping)
            and any(manifest.get(key) for key in ("proteinSets", "recordAnalyses", "recordInstances"))
            and manifest not in manifests
        ):
            manifests.append(manifest)
    if len(manifests) > 1:
        raise ValidationError(
            f"Session {CURRENT_SESSION_VERSION} drawings share one protein identity manifest; "
            "two drawings carry different ones.",
            diagnostic=_LAYOUT_INVALID,
        )
    if manifests:
        shared["proteinIdentityManifest"] = manifests[0]
    return shared


def _assemble(
    shared: Mapping[str, Any],
    result_sets: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    data = _view(shared, result_sets[0])
    if len(result_sets) > 1:
        if not result_sets[1].get("results"):
            mode = _request_mode(result_sets[1]) or "second"
            raise ValidationError(
                f"Session {CURRENT_SESSION_VERSION} keeps a second drawing only with its Result; "
                f"give the {mode} drawing its Results or save it alone.",
                diagnostic=_LAYOUT_INVALID,
            )
        data["otherModeResult"] = dict(result_sets[1])
    return data


def write_session_drawings(
    views: Sequence[Mapping[str, Any]],
    *,
    active: int | None = None,
) -> dict[str, Any]:
    """The document fields that hold ``views``, one one-drawing view per drawing.

    Session 46 keeps a second drawing, of the other mode, only with its
    Result (``otherModeResult``). Each drawing's view gives its mode's slice
    in ``modes``; the drawings share every other field outside a committed
    Result set (see :func:`_merged_shared`). ``active`` is the index
    of the drawing the Web app opens (``ui.mode``).
    """

    if not views:
        raise ValidationError("A Session holds at least one drawing.", diagnostic=_LAYOUT_INVALID)
    if len(views) > 2:
        raise ValidationError(
            f"Session {CURRENT_SESSION_VERSION} holds at most one Circular and one Linear drawing.",
            diagnostic=_LAYOUT_INVALID,
        )
    modes = [_request_mode(view) for view in views]
    if len(views) == 2:
        if None in modes:
            raise ValidationError(
                f"Session {CURRENT_SESSION_VERSION} holds a drawing without a request only as its only drawing.",
                diagnostic=_LAYOUT_INVALID,
            )
        if modes[0] == modes[1]:
            raise ValidationError(
                f"Session {CURRENT_SESSION_VERSION} holds at most one drawing of each mode; "
                f"both drawings are {modes[0]}.",
                diagnostic=_LAYOUT_INVALID,
            )
    if len(views) == 1:
        data = _without(views[0], frozenset({"otherModeResult"}))
    else:
        data = _assemble(_merged_shared(views), [_result_set(view) for view in views])
    if active is not None:
        shown = modes[active] or _mapping(views[active].get("ui")).get("mode")
        data["ui"] = {**_mapping(data.get("ui")), "mode": shown}
    return data


def replace_session_drawings(
    data: Mapping[str, Any],
    views: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """``data`` with the drawings named in ``views`` replaced, in place and order.

    Each view is one drawing seen as a one-drawing Session, of the same mode.
    The other drawings keep their committed Result sets and their mode
    slices. Session 46 drawings share every other field outside a committed
    Result set, so the replaced drawings' views give those fields (see
    :func:`_merged_shared`).
    """

    drawings = read_session_drawings(data)
    unknown = sorted(set(views) - {drawing.id for drawing in drawings})
    if unknown:
        raise ValidationError("Session has no drawing " + ", ".join(unknown) + " to replace.", diagnostic=_LAYOUT_INVALID)
    replaced = [views[drawing.id] for drawing in drawings if drawing.id in views]
    for drawing in drawings:
        if drawing.id in views and _request_mode(views[drawing.id]) != drawing.mode:
            raise ValidationError(f"The {drawing.id} drawing keeps its {drawing.mode} mode.", diagnostic=_LAYOUT_INVALID)
    if len(drawings) == 1:
        return _without(replaced[0], frozenset({"otherModeResult"}))
    return _assemble(
        _merged_shared(
            replaced,
            kept=[drawing.artifacts.fields for drawing in drawings if drawing.id not in views],
        ),
        [
            _result_set(views[drawing.id] if drawing.id in views else drawing.artifacts.fields)
            for drawing in drawings
        ],
    )


__all__ = [
    "DRAWING_NAMES",
    "DrawingMode",
    "SessionDrawingArtifacts",
    "SessionDrawingParts",
    "active_session_drawing",
    "allocate_session_drawing",
    "drawing_draft_config",
    "read_session_drawings",
    "replace_session_drawings",
    "write_session_drawings",
]
