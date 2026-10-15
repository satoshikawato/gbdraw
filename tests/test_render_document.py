"""The document writer: one serialization per request, bytes as svgwrite writes them."""

from __future__ import annotations

import io
import re
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.render as api_render
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
)
from gbdraw.api import LinearDiagramOptions, LinearTrackOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_records,
    build_linear_diagram,
)
from gbdraw.circular import circular_main
from gbdraw.io.record_select import parse_record_selector
from gbdraw.render.document import SvgDocument, write_svg_document
from gbdraw.render.group_content import (
    ContentGroup,
    GroupContent,
    GroupKey,
    GroupKind,
    group_manifest,
    place_content,
)
from gbdraw.render.patterns import hatch_pattern
from gbdraw.svg.elements import ClipPath, Group, Rect, Text
from gbdraw.svg.elements import Path as SvgPath
from gbdraw.render.export import has_cairosvg

INPUTS = Path(__file__).parent / "test_inputs"
HATCHED = RegionAnnotationStyle(fill="#88aaff", hatch=HatchStyle(cross=True))


def _linear_with_clip_and_hatch() -> Drawing:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotations = (
        RegionAnnotation("one", CoordinateSpan(None, 100, 300), lane=0, mark="band", style=HATCHED),
        RegionAnnotation("two", CoordinateSpan(None, 400, 600), lane=1),
    )
    return build_linear_diagram(
        [record],
        options=LinearDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", annotations),)),
            tracks=LinearTrackOptions(
                linear_track_slots=("regions:annotations@set_id=regions,h=8px,overflow=clip",),
                linear_track_axis_index=1,
            ),
        ),
    )


def _multi_record_circular_with_clip_and_hatch() -> Drawing:
    records = [
        SeqRecord(Seq("A" * 1000), id=f"record-{name}", name=f"record-{name}")
        for name in ("a", "b")
    ]
    annotations = tuple(
        RegionAnnotation(
            f"record-{number}-lane-{lane}",
            CoordinateSpan(parse_record_selector(f"#{number}"), lane * 100, lane * 100 + 80),
            lane=lane - 1,
            mark="band",
            style=HATCHED,
        )
        for number in (1, 2)
        for lane in (1, 2)
    )
    return assemble_circular_diagram_from_records(
        records,
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", annotations),)),
        circular_track_slots=[
            "regions:annotations@set_id=regions,w=8px,overflow=clip",
            "features:features",
            "ticks:ticks",
        ],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    )


def _circular_with_empty_label_groups(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> Drawing:
    drawings: list[Drawing] = []
    original = api_render.write_svg_document

    def capture(drawing: Drawing) -> SvgDocument:
        drawings.append(drawing)
        return original(drawing)

    monkeypatch.setattr(api_render, "write_svg_document", capture)
    circular_main([
        "--gbk", str(INPUTS / "HmmtDNA.gbk"),
        "-o", str(tmp_path / "embedded-only"),
        "-f", "svg",
        "--labels", "out",
        "--label_rendering", "embedded_only",
    ])
    assert len(drawings) == 1
    return drawings[0]


def _top_level(drawing: Drawing) -> list[tuple[str, object, int]]:
    return [
        (element.elementname, element.attribs.get("id"), len(element.elements))
        for element in drawing.elements
    ]


def _defs_kinds(drawing: Drawing) -> set[str]:
    return {element.elementname for element in drawing.defs.elements}


def _assert_writer_matches_svgwrite(drawing: Drawing) -> SvgDocument:
    document = write_svg_document(drawing)
    assert document.body == drawing.tostring()
    with_declaration = io.StringIO()
    drawing.write(with_declaration)
    assert document.prolog + document.body == with_declaration.getvalue()
    assert document.prolog == '<?xml version="1.0" encoding="utf-8" ?>\n'
    written = io.StringIO()
    document.write(written)
    assert written.getvalue() == with_declaration.getvalue()
    return document


def test_writer_matches_svgwrite_for_linear_clip_paths_and_hatch_patterns() -> None:
    drawing = _linear_with_clip_and_hatch()
    assert {"clipPath", "pattern"} <= _defs_kinds(drawing)
    _assert_writer_matches_svgwrite(drawing)


def test_writer_matches_svgwrite_for_a_multi_record_canvas() -> None:
    drawing = _multi_record_circular_with_clip_and_hatch()
    assert {"clipPath", "pattern"} <= _defs_kinds(drawing)
    assert sum(1 for kind, _id, _count in _top_level(drawing) if kind == "g") >= 2
    _assert_writer_matches_svgwrite(drawing)


def test_writer_matches_svgwrite_with_empty_top_level_groups(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    drawing = _circular_with_empty_label_groups(tmp_path, monkeypatch)
    empty = {
        element_id
        for kind, element_id, count in _top_level(drawing)
        if kind == "g" and count == 0
    }
    assert {"label_leaders", "label_text"} <= empty
    document = _assert_writer_matches_svgwrite(drawing)
    assert re.search(r'<g [^>]*id="label_leaders"[^>]*/>', document.body)
    written = (tmp_path / "embedded-only.svg").read_text(encoding="utf-8")
    assert written == document.prolog + document.body


def _count_serializations(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """Record each serialization: a call of the owner, or an svgwrite call outside it."""

    calls: list[int] = []
    owner = api_render.write_svg_document
    original_tostring = Drawing.tostring
    original_write = Drawing.write
    inside_owner = False

    def counting_owner(drawing: Drawing) -> SvgDocument:
        nonlocal inside_owner
        calls.append(id(drawing))
        inside_owner = True
        try:
            return owner(drawing)
        finally:
            inside_owner = False

    def direct_tostring(self: Drawing) -> str:
        if not inside_owner:
            calls.append(id(self))
        return original_tostring(self)

    def direct_write(self: Drawing, *args: Any, **kwargs: Any) -> Any:
        if not inside_owner:
            calls.append(id(self))
        return original_write(self, *args, **kwargs)

    monkeypatch.setattr(api_render, "write_svg_document", counting_owner)
    monkeypatch.setattr(Drawing, "tostring", direct_tostring)
    monkeypatch.setattr(Drawing, "write", direct_write)
    return calls


def test_save_figure_to_serializes_once_and_keeps_the_declaration_in_the_file_only(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    drawing = _linear_with_clip_and_hatch()
    expected = drawing.tostring()
    formats = ["svg", "interactive_svg"] + (["png"] if has_cairosvg() else [])
    calls = _count_serializations(monkeypatch)
    paths = api_render.save_figure_to(
        drawing,
        formats,
        output_dir=str(tmp_path),
        output_prefix="once",
    )
    assert len(calls) == 1
    svg_path = tmp_path / "once.svg"
    assert str(svg_path) in paths
    assert svg_path.read_text(encoding="utf-8") == (
        '<?xml version="1.0" encoding="utf-8" ?>\n' + expected
    )
    interactive = next(Path(path) for path in paths if path.endswith(".interactive.svg"))
    assert not interactive.read_text(encoding="utf-8").startswith("<?xml")


@pytest.mark.parametrize("fmt", ["svg", "interactive_svg"])
def test_render_to_bytes_serializes_once_without_a_declaration(
    fmt: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    drawing = _linear_with_clip_and_hatch()
    expected = drawing.tostring()
    calls = _count_serializations(monkeypatch)
    rendered = api_render.render_to_bytes(drawing, fmt)
    assert len(calls) == 1
    assert rendered.startswith(b"<svg")
    if fmt == "svg":
        assert rendered == expected.encode("utf-8")


# Placed group contents (design 2.4, 2.5, F13).


class _CountingPath(SvgPath):
    serializations = 0

    def get_xml(self):
        type(self).serializations += 1
        return super().get_xml()


def _built_group(group_id: str, *, children: int = 2) -> Group:
    group = Group(id=group_id, transform="translate(3,4)")
    group.attribs["data-kind"] = "test & <quoted>"
    for index in range(children):
        group.add(_CountingPath(d=f"M{index},0 L10,{index}", fill="#ff0000"))
    if children:
        group.add(Text("label & <text>", insert=(1.5, 2.25)))
    return group


def _clip_path(clip_id: str) -> ClipPath:
    clip = ClipPath(id=clip_id)
    clip.add(Rect(insert=(0, 0), size=(10, 5)))
    return clip


def _drawing() -> tuple[Drawing, list[ContentGroup]]:
    """Contents with children, an empty content (F13), a plain empty group, defs."""

    hatch = hatch_pattern(HatchStyle(angle=45, spacing=4, color="#000", width=1))
    drawing = Drawing(debug=False, size=("100px", "50px"))
    drawing.add(Group(id="plain_empty"))
    contents = (
        GroupContent.of(_built_group("first"), definitions=(hatch, _clip_path("clip-1"))),
        # The same hatch again is not added twice; a clip path always is.
        GroupContent.of(_built_group("second"), definitions=(hatch, _clip_path("clip-1"))),
        GroupContent.of(_built_group("empty_content", children=0)),
    )
    placed = []
    for index, content in enumerate(contents):
        group = place_content(drawing, content, GroupKey("record-1", GroupKind.FEATURES, index))
        group.translate(5, 6)
        group.attribs["data-request"] = str(index)
        drawing.add(group)
        placed.append(group)
    drawing.add(_built_group("plain_full"))
    return drawing, placed


def _svgwrite_file_text(drawing: Drawing) -> str:
    handle = io.StringIO()
    drawing.write(handle)
    return handle.getvalue()


def test_the_writer_gives_svgwrites_bytes_for_contents_definitions_and_empty_groups() -> None:
    drawing, _placed = _drawing()
    document = write_svg_document(drawing)

    assert document.body == drawing.tostring()
    assert document.prolog + document.body == _svgwrite_file_text(drawing)
    assert 'data-request="2" id="empty_content" transform="translate(3,4) translate(5,6)" />' in (
        document.body
    )
    assert '<g id="plain_empty" />' in document.body
    definition_ids = [item.attribs["id"] for item in drawing.defs.elements]
    assert definition_ids.count("clip-1") == 2
    assert len([i for i in definition_ids if i.startswith("gbdraw-hatch-")]) == 1


def test_a_placed_content_is_not_serialized_again() -> None:
    drawing, placed = _drawing()
    for group in placed:
        _ = group.content.children_text
    _CountingPath.serializations = 0
    first = write_svg_document(drawing)
    # A second request places the same contents in new groups.
    again = Drawing(debug=False, size=("100px", "50px"))
    for group in placed:
        fresh = place_content(again, group.content, group.key)
        fresh.translate(5, 6)
        again.add(fresh)
    second = write_svg_document(again)

    # Only "plain_full" (two paths) is serialized, once per write of the first drawing.
    assert _CountingPath.serializations == 2
    assert second.body == again.tostring()
    assert first.body.count("label &amp; &lt;text&gt;") == 3


def test_a_group_whose_children_changed_is_written_by_svgwrite() -> None:
    drawing, placed = _drawing()
    _ = placed[0].content.children_text
    placed[0].add(Rect(insert=(1, 1), size=(2, 2), id="added"))

    assert not placed[0].holds_content_children()
    assert write_svg_document(drawing).body == drawing.tostring()
    assert 'id="added"' in write_svg_document(drawing).body


def test_the_manifest_lists_placed_groups_in_document_order() -> None:
    drawing, placed = _drawing()
    manifest = group_manifest(drawing)

    assert [entry.key for entry in manifest] == [group.key for group in placed]
    assert [entry.dom_id for entry in manifest] == ["first", "second", "empty_content"]
    assert manifest[0].position == "translate(3,4) translate(5,6)"
    assert manifest[0].definition_ids == ("gbdraw-hatch-" + manifest[0].definition_ids[0][13:], "clip-1")
    assert manifest[0].content.revision == manifest[1].content.revision
    assert manifest[0].content.revision != manifest[2].content.revision


def test_an_empty_drawing_is_written_as_svgwrite_writes_it() -> None:
    drawing = Drawing(debug=False)
    drawing.elements.clear()
    document = write_svg_document(drawing)

    assert document.prolog + document.body == _svgwrite_file_text(drawing)
