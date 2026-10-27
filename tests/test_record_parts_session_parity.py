"""Warm equals cold over the Gallery Sessions' G2 changes (SW2 design 3.3 item 1).

Each change of a deciding Session runs as change, restore, change on one warm
cache, as the page sends them: the catalog revision the page holds is passed
back. Every reply is compared with a render of the same request on an empty
cache: the SVG results, and the catalog the page holds after the reply (the
one the reply sent, or the one it kept). A change must draw other bytes than
its restore, or the step proves nothing.

Slow: a Vnig render takes about 3.5 s; each Session runs alone.
"""

from __future__ import annotations

import copy
import gzip
import json
from collections.abc import Callable
from pathlib import Path

import pytest

from gbdraw.api.prepared import PreparedBiologicalInputCache
from tests.test_web_render_preparation_cache import _render

pytestmark = pytest.mark.slow

SESSIONS = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "gallery" / "sessions"

Request = dict[str, object]
Change = Callable[[Request], None]


def _options(request: Request) -> dict[str, object]:
    options = request["diagramOptions"]
    assert isinstance(options, dict)
    return options


def _override(**values: object) -> Change:
    """Config overrides by path (``__`` for ``.``), as the Web form writes them."""

    def apply(request: Request) -> None:
        options = _options(request)
        overrides = options.get("configOverrides")
        assert isinstance(overrides, dict)
        for name, value in values.items():
            overrides[name.replace("__", ".")] = value

    return apply


def _option(name: str, value: object) -> Change:
    def apply(request: Request) -> None:
        _options(request)[name] = value

    return apply


def _record_label(label: str) -> Change:
    def apply(request: Request) -> None:
        records = request["records"]
        assert isinstance(records, list)
        presentation = records[0].get("presentation") or {}
        records[0]["presentation"] = {**presentation, "label": label}

    return apply


def _layout(name: str, value: object) -> Change:
    def apply(request: Request) -> None:
        layout = request.get("layout") or {}
        assert isinstance(layout, dict)
        request["layout"] = {**layout, name: value}

    return apply


def _tracks(name: str, value: object) -> Change:
    def apply(request: Request) -> None:
        options = _options(request)
        tracks = options.get("tracks") or {}
        assert isinstance(tracks, dict)
        options["tracks"] = {**tracks, name: value}

    return apply


def _both(*changes: Change) -> Change:
    def apply(request: Request) -> None:
        for change in changes:
            change(request)

    return apply


def _gc_toggle(request: Request) -> None:
    overrides = _options(request)["configOverrides"]
    assert isinstance(overrides, dict)
    overrides["canvas.show_gc"] = not bool(overrides.get("canvas.show_gc", True))


_STROKE = _override(objects__features__block_stroke_color="#1a5fb4")
# (row, prepare or None, change): the S-run rows; the values that draw new bytes (SW2 design 7, OV-1314).
_ROWS: dict[str, tuple[tuple[str, Change | None, Change], ...]] = {
    "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz": (
        ("definition font", None, _override(objects__definition__circular__font_size=22)),
        ("scale interval", None, _override(objects__scale__interval=500000)),
        ("window", None, _option("window", 40000)),
        ("stroke", None, _STROKE),
        ("custom stack", None, _tracks("centerReservedRadius", 60)),
        ("GC toggle", None, _gc_toggle),
    ),
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz": (
        ("record label", None, _record_label("Test record X")),
        (
            "definition font",
            None,
            _override(objects__definition__linear__font_size__short=20, objects__definition__linear__font_size__long=20),
        ),
        ("scale interval", None, _override(objects__scale__interval=200000)),
        ("window", None, _option("window", 2000)),
        ("comparison height", None, _override(canvas__linear__comparison_height=80)),
        ("stroke", None, _STROKE),
        ("GC toggle", None, _gc_toggle),
    ),
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz": (
        ("record label", None, _record_label("Test record X")),
        ("titles", None, _override(objects__definition__linear__line_styles__name__font_size=22)),
        ("scale interval", None, _override(objects__scale__interval=500000)),
        ("comparison height", None, _override(canvas__linear__comparison_height=80)),
        ("record gap", None, _layout("recordGapPx", 96)),
        ("stroke", None, _STROKE),
        ("GC toggle", None, _gc_toggle),
        (
            "window with GC",
            _override(canvas__show_gc=True),
            _both(_override(canvas__show_gc=True), _option("window", 40000)),
        ),
    ),
}


@pytest.mark.parametrize("session", list(_ROWS))
def test_the_g2_changes_of_a_gallery_session_render_warm_as_cold(session: str, tmp_path: Path) -> None:
    path = SESSIONS / session
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        document = json.load(handle)
    base = document["renderRequest"]
    cache = PreparedBiologicalInputCache()
    request_ids = iter(range(1, 10_000))
    cold_replies: dict[str, dict[str, object]] = {}
    # What the page holds: the catalog revision it passes back, and that catalog.
    held_revision: str | None = None
    held_catalog: object = None
    failures: list[str] = []

    def render(request: Request, label: str) -> dict[str, object]:
        nonlocal held_revision, held_catalog
        step = {**document, "renderRequest": request}
        reply, _ = _render(tmp_path, cache, step, request_id=next(request_ids), held_revision=held_revision)
        revision = reply["featureCatalogRevision"]
        if reply["featureCatalog"] is not None:
            held_catalog = reply["featureCatalog"]
        elif revision != held_revision:
            failures.append(f"{label}: no catalog sent, yet the revision changed")
        held_revision = revision if isinstance(revision, str) else None
        key = json.dumps(request, sort_keys=True)
        if key not in cold_replies:
            cold_replies[key], _ = _render(
                tmp_path, PreparedBiologicalInputCache(), step, request_id=next(request_ids)
            )
        cold = cold_replies[key]
        if reply["results"] != cold["results"]:
            failures.append(f"{label}: warm SVG differs from cold")
        if held_catalog != cold["featureCatalog"]:
            failures.append(f"{label}: the catalog the page holds differs from cold")
        return reply

    render(base, "base")
    for row, prepare, change in _ROWS[session]:
        start = copy.deepcopy(base)
        if prepare is not None:
            prepare(start)
            render(start, f"{row}: prepare")
        changed = copy.deepcopy(start)
        change(changed)
        first = render(changed, f"{row}: change")
        restored = render(start, f"{row}: restore")
        render(changed, f"{row}: change again")
        if first["results"] == restored["results"]:
            failures.append(f"{row}: the change draws the bytes of its restore")
    assert not failures, "\n".join(failures)
