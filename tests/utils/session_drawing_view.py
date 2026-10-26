"""A written Session 49 document's drawing as the one-drawing view (the 46-48 layout).

Built on the production reader (``_drawing_view``); the 49 layout is derived
nowhere else.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator, Mapping

from gbdraw.session_drawings import _drawing_view


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


@contextmanager
def without_drawings_move() -> Iterator[None]:
    """Make ``build_session_json`` return the checked mode-scoped Session 48 it builds.

    The writer builds the 46-48 layout and moves it to drawings last; inside
    this context the move is skipped. For tests that need the 46-48 layout of
    what a writer builds (a CLI save included).
    """

    from unittest import mock

    from gbdraw import session_drawings

    def keep_mode_scoped(data: Mapping[str, Any], *_args: Any, **_kwargs: Any) -> tuple[Mapping[str, Any], None]:
        return data, None

    with mock.patch.object(session_drawings, "drawings_document_from_mode_scoped", keep_mode_scoped):
        yield


def build_mode_scoped_session_json(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """``build_session_json`` as its checked mode-scoped Session 48, before the move to drawings."""

    from gbdraw.session_io import build_session_json

    with without_drawings_move():
        return dict(build_session_json(*args, **kwargs))


def build_mode_scoped_session_document(request: Any) -> dict[str, Any]:
    """The checked mode-scoped Session 48 of one request (``build_session_document`` before its move to drawings)."""

    from gbdraw.session import _mode_scoped_document_of_request

    return _mode_scoped_document_of_request(request, resolved=False, created_at=None).to_dict()


def upgrade_to_mode_scoped(source: Any, temporary_directory: Any = None) -> dict[str, Any]:
    """An older Session upgraded to the newest mode-scoped Session (48), before the move to drawings."""

    from gbdraw.session import _upgrade_to_current, load_session_document

    document = source if hasattr(source, "to_dict") else load_session_document(source)
    return _upgrade_to_current(document, temporary_directory).document.to_dict()


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
