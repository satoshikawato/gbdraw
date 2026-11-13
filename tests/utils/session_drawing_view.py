"""A written Session 49 document's drawing as the one-drawing view (the 46-48 layout).

Built on the production reader (``_drawing_view``); the 49 layout is derived
nowhere else.
"""

from __future__ import annotations

from typing import Any, Mapping

from gbdraw.session_drawings import _drawing_view
from gbdraw.session_io import FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION, mode_slice_of_version


def drawing_view(
    document: Mapping[str, Any],
    drawing_id: str | None = None,
    *,
    mode: str | None = None,
) -> dict[str, Any]:
    """The view of the drawing with ``drawing_id``, or the only drawing of ``mode``.

    With neither, the document must hold exactly one drawing.
    """

    drawings = [
        drawing
        for drawing in document["drawings"]
        if (drawing_id is None or drawing["id"] == drawing_id)
        and (mode is None or drawing["mode"] == mode)
    ]
    assert len(drawings) == 1, [(d["id"], d["mode"]) for d in document["drawings"]]
    return _drawing_view(document, drawings[0])


def view_at(view: Mapping[str, Any], version: int) -> dict[str, Any]:
    """``view`` as a Session ``version`` holds it: no slice field a later Session added (``since``)."""

    modes = {mode: mode_slice_of_version(slice_, version) for mode, slice_ in view.get("modes", {}).items()}
    return catalogs_at({**view, "version": version, **({"modes": modes} if "modes" in view else {})}, version)


def catalogs_at(session: Mapping[str, Any], version: int) -> dict[str, Any]:
    """``session`` with each feature catalog (its own, the other mode's Result set's, each mode slice's) labelled
    with the schema a Session ``version`` writes: a catalog with no sequence (schema 6, D5) is a valid schema 5
    catalog. ``session`` is not changed."""

    expected = FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION.get(version)

    def labelled(holder: Any) -> Any:
        editor = holder.get("editorState") if isinstance(holder, Mapping) else None
        catalog = editor.get("featureCatalog") if isinstance(editor, Mapping) else None
        if expected is None or not isinstance(catalog, Mapping):
            return holder
        return {**holder, "editorState": {**editor, "featureCatalog": {**catalog, "schema": expected}}}

    held = dict(labelled(session))
    if isinstance(held.get("otherModeResult"), Mapping):
        held["otherModeResult"] = labelled(held["otherModeResult"])
    if isinstance(held.get("modes"), Mapping):
        held["modes"] = {mode: labelled(slice_) for mode, slice_ in held["modes"].items()}
    return held


def build_mode_scoped_session_json(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """``build_session_json`` as its checked mode-scoped Session 48, before the move to drawings.

    The writer builds the 46-48 layout and moves it to drawings last; here the
    move is skipped.
    """

    from unittest import mock

    from gbdraw import session_drawings
    from tests.utils.session_files import build_session_payload

    def keep_mode_scoped(data: Mapping[str, Any], *_args: Any, **_kwargs: Any) -> tuple[Mapping[str, Any], None]:
        return catalogs_at(data, data["version"]), None

    with mock.patch.object(session_drawings, "drawings_document_from_mode_scoped", keep_mode_scoped):
        return build_session_payload(*args, **kwargs)


def build_mode_scoped_session_document(request: Any) -> dict[str, Any]:
    """The checked mode-scoped Session 48 of one request (``build_session_document`` before its move to drawings),
    its resources' bytes embedded."""

    from gbdraw.session import _mode_scoped_document_of_request
    from tests.utils.session_files import embedded_session

    document = _mode_scoped_document_of_request(request, resolved=False, created_at=None)
    return embedded_session(document.to_dict(), document.resources.sources)


def upgrade_to_mode_scoped(source: Any, temporary_directory: Any = None) -> dict[str, Any]:
    """An older Session upgraded to the newest mode-scoped Session (48), before the move to drawings,
    its resources' bytes embedded."""

    from gbdraw.session import _upgrade_to_current, load_session_document

    from tests.utils.session_files import embedded_session

    document = source if hasattr(source, "to_dict") else load_session_document(source)
    upgraded = _upgrade_to_current(document, temporary_directory).document
    return embedded_session(upgraded.to_dict(), upgraded.resources.sources)


def drawings_document_of_view(view: Mapping[str, Any]) -> dict[str, Any]:
    """A one-drawing view (the 46-48 layout) as the one-drawing Session 49 document, by the production writer."""

    from gbdraw.session_drawings import _drawings_document, _view_mode

    mode = _view_mode(view)
    return _drawings_document(
        [view], ids=[mode], names=[mode.capitalize()], active=0, project_views=[view]
    )


# fixround C2
def session_view(document: Mapping[str, Any]) -> dict[str, Any]:
    """The one-drawing view of a 49, or the 27-48 document itself (a fixture an older writer made)."""

    return drawing_view(document) if "drawings" in document else dict(document)
