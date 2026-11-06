"""Sequence search hits as an annotation set, and the annotation gradient style they draw with (PR-B)."""

from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.annotations import (
    AnnotationGradient,
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
    ResolvedRegionAnnotation,
)
from gbdraw.annotations.search_hits import (
    QUERY_POSITION_LEGEND_LABEL,
    SEARCH_GRADIENT_END_COLOR,
    SEARCH_GRADIENT_START_COLOR,
    SEARCH_HIT_METADATA_KEYS,
    SEARCH_QUERY_COLORS,
    SearchQuery,
    SearchSubject,
    origin_window_span,
    search_hit_annotation_set,
)
from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.requests import InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.core.color import interpolate_color
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.io.record_select import parse_record_selector
from gbdraw.mode_profiles import ComparisonThresholds
from gbdraw.render.drawers.annotation_marks import gradient_steps, mark_parts
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request
from gbdraw.web_support.sequence_search import sequence_search_annotation_set_json

THRESHOLDS = ComparisonThresholds(evalue=1e-2, bitscore=50.0, identity=70.0, alignment_length=0)
# query subject identity length mismatch gap qstart qend sstart send evalue bitscore
ROWS = "\n".join(
    "\t".join(row)
    for row in (
        ("qA", "rec1:1-5000", "100.000", "600", "0", "0", "1", "600", "1201", "1800", "0.0", "1083"),
        ("qB", "rec1:1-5000", "84.783", "46", "7", "0", "320", "365", "2213", "2168", "2.71e-09", "52.7"),
        ("qB", "rec2:1-3000", "90.000", "200", "0", "0", "1", "200", "101", "300", "1e-30", "300"),
        # Removed by the thresholds: bitscore 41.9 < 50, identity 60 < 70.
        ("qB", "rec1:1-5000", "60.000", "30", "0", "0", "1", "30", "10", "39", "1e-3", "41.9"),
        # Names no searched sequence (PD-OI-074): left out and counted.
        ("qA", "other", "99.000", "100", "0", "0", "1", "100", "1", "100", "1e-40", "200"),
    )
)
# qA is 600 bp and its hit covers all of it; qB is 400 bp.
QUERIES = (SearchQuery("qA", 600), SearchQuery("qB", 400))
SUBJECTS = (
    SearchSubject("rec1:1-5000", parse_record_selector("rec1"), "record-1"),
    SearchSubject("rec2:1-3000", parse_record_selector("rec2"), "record-2", source_start=1001),
)


def _search(**kwargs: object):
    return search_hit_annotation_set(
        read_comparison_table(io.StringIO(ROWS)),
        set_id="Search 1",
        subjects=SUBJECTS,
        thresholds=THRESHOLDS,
        program="blastn",
        queries=QUERIES,
        **kwargs,  # type: ignore[arg-type]
    )


def test_hits_become_one_highlight_row_each_with_query_colours_and_metadata() -> None:
    result = _search()
    assert result.annotation_set is not None
    rows = result.annotation_set.annotations
    assert [row.id for row in rows] == ["qA#1", "qB#1", "qB#2"]
    assert {row.mark for row in rows} == {"highlight"} and {row.label for row in rows} == {""}
    assert [(row.target.record.record_id, row.target.start, row.target.end) for row in rows] == [
        ("rec1", 1201, 1800),
        ("rec1", 2168, 2213),  # a minus-strand hit spans its lower to higher coordinate
        ("rec2", 1101, 1300),  # the second record, its sequence starting at source base 1001
    ]
    assert [row.style.fill for row in rows] == [SEARCH_QUERY_COLORS[0], SEARCH_QUERY_COLORS[1], SEARCH_QUERY_COLORS[1]]
    assert [row.legend_label for row in rows] == ["qA", "qB", "qB"]
    assert all(row.style.gradient is None for row in rows)
    assert tuple(rows[1].metadata) == SEARCH_HIT_METADATA_KEYS
    assert dict(rows[1].metadata) == {
        "query": "qB", "program": "blastn", "identity": "84.783", "alignment_length": "46",
        "evalue": "2.71e-09", "bitscore": "52.7",
    }
    assert [hit.strand for hit in result.hits] == ["+", "-", "+"]
    assert [hit.record_key for hit in result.hits] == ["record-1", "record-1", "record-2"]


def test_thresholds_and_unnamed_subjects_are_counted() -> None:
    result = _search()
    assert (result.tally.rows_read, result.tally.rows_kept) == (4, 3)
    assert result.rows_removed_by_thresholds == 1
    assert {item.name for item in result.tally.thresholds} == {"bitscore", "identity"}
    assert (result.unmatched_rows, result.unmatched_subject_ids) == (1, ("other",))


def test_no_hit_gives_no_set() -> None:
    strict = ComparisonThresholds(evalue=1e-2, bitscore=5000.0, identity=70.0, alignment_length=0)
    result = search_hit_annotation_set(
        read_comparison_table(io.StringIO(ROWS)), set_id="Search 1", subjects=SUBJECTS, thresholds=strict,
        program="blastn", queries=QUERIES,
    )
    assert result.annotation_set is None and result.hits == ()
    assert result.rows_removed_by_thresholds == 4


def test_query_colours_replace_the_default() -> None:
    result = _search(query_colors={"qB": "#123456"})
    assert result.annotation_set is not None
    assert [row.style.fill for row in result.annotation_set.annotations] == [SEARCH_QUERY_COLORS[0], "#123456", "#123456"]


def test_gradient_hits_carry_their_direction_and_one_legend_caption() -> None:
    result = _search(color_mode="gradient")
    assert result.annotation_set is not None
    rows = result.annotation_set.annotations
    assert [row.style.gradient.direction for row in rows] == ["forward", "reverse", "forward"]
    assert {(row.style.gradient.start_color, row.style.gradient.end_color) for row in rows} == {
        (SEARCH_GRADIENT_START_COLOR, SEARCH_GRADIENT_END_COLOR)
    }
    # Each hit's place in its query (D-44): qA covers all of it, qB 320-365 and 1-200 of 400.
    assert [(row.style.gradient.query_start, row.style.gradient.query_end) for row in rows] == [
        (0.0, 1.0), (319 / 400, 365 / 400), (0.0, 0.5),
    ]
    # The end colour is the plain fill an annotation TSV export writes (receipt).
    assert {row.style.fill for row in rows} == {SEARCH_GRADIENT_END_COLOR}
    assert {row.legend_label for row in rows} == {QUERY_POSITION_LEGEND_LABEL}


def test_ambiguous_subject_id_is_refused() -> None:
    twice = (SUBJECTS[0], SearchSubject("rec1:1-5000", parse_record_selector("rec2"), "record-2"))
    with pytest.raises(ValidationError) as error:
        search_hit_annotation_set(
            read_comparison_table(io.StringIO(ROWS)), set_id="S", subjects=twice, thresholds=THRESHOLDS,
            program="blastn", queries=QUERIES[:1],
        )
    assert error.value.diagnostic == {"code": "COMPARISON_INPUT", "reason": "RECORD_ID", "column": 2}


def test_worker_helper_returns_the_set_as_the_request_writes_it() -> None:
    payload = {
        "table": ROWS, "setId": "Search 1",
        "subjects": [
            {"subjectId": "rec1:1-5000", "record": {"kind": "recordId", "value": "rec1"}, "recordKey": "record-1", "sourceStart": 1,
             "wrapsOrigin": False, "recordLength": 5000},
            {"subjectId": "rec2:1-3000", "record": {"kind": "recordId", "value": "rec2"}, "recordKey": "record-2", "sourceStart": 1001,
             "wrapsOrigin": False, "recordLength": 4000},
        ],
        "thresholds": {"evalue": 1e-2, "bitscore": 50, "identity": 70, "alignmentLength": 0},
        "program": "blastn", "queries": [{"name": "qA", "length": 600}, {"name": "qB", "length": 400}],
        "colorMode": "gradient", "queryColors": {},
    }
    reply = json.loads(sequence_search_annotation_set_json(json.dumps(payload)))
    assert [row["id"] for row in reply["annotationSet"]["annotations"]] == ["qA#1", "qB#1", "qB#2"]
    gradients = [row["style"]["gradient"] for row in reply["annotationSet"]["annotations"]]
    # The whole query writes no query positions; a part writes both.
    assert gradients[0] == {
        "startColor": SEARCH_GRADIENT_START_COLOR, "endColor": SEARCH_GRADIENT_END_COLOR, "direction": "forward",
    }
    assert gradients[1] == {
        "startColor": SEARCH_GRADIENT_START_COLOR, "endColor": SEARCH_GRADIENT_END_COLOR, "direction": "reverse",
        "queryStart": 319 / 400, "queryEnd": 365 / 400,
    }
    assert (reply["rowsRead"], reply["rowsKept"], reply["unmatchedRows"]) == (4, 3, 1)
    assert reply["hits"][1]["strand"] == "-"


def test_gradient_style_validation_has_typed_reasons() -> None:
    assert AnnotationGradient("pink", "#5B1A4A").start_color == "#ffc0cb"
    with pytest.raises(ValidationError) as color:
        AnnotationGradient("not-a-colour", "#5b1a4a")
    assert color.value.diagnostic == {"code": "INPUT_INVALID", "reason": "COLOR"}
    with pytest.raises(ValidationError) as direction:
        AnnotationGradient("#ffffff", "#000000", "sideways")  # type: ignore[arg-type]
    assert direction.value.diagnostic == {"code": "INPUT_INVALID", "reason": "FIELDS"}
    with pytest.raises(ValidationError) as both:
        RegionAnnotationStyle(hatch=HatchStyle(), gradient=AnnotationGradient("#ffffff", "#000000"))
    assert both.value.diagnostic == {"code": "INPUT_INVALID", "reason": "FIELDS"}


@pytest.mark.parametrize(
    ("query_start", "query_end", "reason"),
    [
        (-0.1, 1.0, "NONNEGATIVE"),
        (0.0, 1.5, "ORDER"),
        (0.8, 0.2, "ORDER"),
        (float("nan"), 1.0, "FINITE"),
        ("half", 1.0, "FINITE"),
        (True, 1.0, "FINITE"),
    ],
)
def test_gradient_query_positions_are_fractions_in_query_order(query_start: object, query_end: object, reason: str) -> None:
    assert AnnotationGradient("#ffffff", "#000000", query_start=0.25, query_end=0.75).query_end == 0.75
    with pytest.raises(ValidationError) as error:
        AnnotationGradient("#ffffff", "#000000", query_start=query_start, query_end=query_end)  # type: ignore[arg-type]
    assert error.value.diagnostic == {"code": "INPUT_INVALID", "reason": reason}


def test_gradient_style_round_trips_through_the_request(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ACGT" * 500), id="rec1", annotations={"molecule_type": "DNA"})
    result = _search(color_mode="gradient")
    assert result.annotation_set is not None
    plain = RegionAnnotation("plain", CoordinateSpan(None, 10, 50), mark="band", style=RegionAnnotationStyle(fill="#88aaff"))
    annotations = AnnotationOptions(sets=(
        AnnotationSet("Search 1", tuple(row for row in result.annotation_set.annotations if row.target.record.record_id == "rec1")),
        AnnotationSet("plain", (plain,)),
    ))
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(annotations=annotations),
    )
    encoded = encode_canonical_request(request)
    sets = encoded.payload["diagramOptions"]["annotations"]["sets"]
    # A style without a gradient is written as before: no key.
    assert "gradient" not in sets[1]["annotations"][0]["style"]
    # A hit over a part of its query writes where it lies in the query.
    assert [
        (row["style"]["gradient"].get("queryStart"), row["style"]["gradient"].get("queryEnd"))
        for row in sets[0]["annotations"]
    ] == [(None, None), (319 / 400, 365 / 400)]
    resource_paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        target.write_bytes(resource.content or b"")
        resource_paths[resource.resource_id] = target
    decoded = decode_canonical_request(encoded.payload, resource_paths=resource_paths, output_directory=tmp_path / "out")
    assert decoded.options.annotations == annotations
    assert encode_canonical_request(decoded).payload == encoded.payload


def _resolved(
    direction: str, source_step: int, query_start: float = 0.0, query_end: float = 1.0
) -> tuple[ResolvedRegionAnnotation, AnnotationGradient]:
    gradient = AnnotationGradient(
        "#ffffff", "#000000", direction, query_start=query_start, query_end=query_end  # type: ignore[arg-type]
    )
    style = RegionAnnotationStyle(fill="#000000", gradient=gradient)
    annotation = ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=((100, 200),), midpoint_bp=150.0, span_bp=100, label="",
        mark="highlight", lane=None, style=style, source_step=source_step,  # type: ignore[arg-type]
    )
    return annotation, gradient


@pytest.mark.parametrize(
    ("direction", "source_step", "pale_first"),
    [("forward", 1, True), ("reverse", 1, False), ("forward", -1, False), ("reverse", -1, True)],
)
def test_gradient_steps_follow_the_hit_strand_and_the_record_direction(
    direction: str, source_step: int, pale_first: bool
) -> None:
    annotation, gradient = _resolved(direction, source_step)
    (part,) = mark_parts(annotation, join_origin=False)
    steps = gradient_steps(annotation, gradient, part, 100.0, 200.0, record_length=1000, length_px=40.0)
    assert len(steps) == 20 and steps[0][0] == 100.0
    first, last = steps[0][1], steps[-1][1]
    assert (first > last) == pale_first  # "#f..." sorts above "#0..."
    # A part of the region takes the colours of its place in the whole region.
    half = gradient_steps(annotation, gradient, part, 150.0, 200.0, record_length=1000, length_px=20.0)
    assert half[-1][1] == last


@pytest.mark.parametrize("direction", ["forward", "reverse"])
def test_a_partial_hit_draws_only_its_part_of_the_query_scale(direction: str) -> None:
    # The hit covers query positions 301-600 of 600: from the mid colour to the dark end (D-44).
    annotation, gradient = _resolved(direction, 1, query_start=0.5, query_end=1.0)
    (part,) = mark_parts(annotation, join_origin=False)
    steps = gradient_steps(annotation, gradient, part, 100.0, 200.0, record_length=1000, length_px=40.0)
    # Each of the 20 steps takes the colour at its middle: 1/40 of the span in from each end.
    near_mid = interpolate_color("#ffffff", "#000000", 0.5 + 0.5 / 40)
    near_end = interpolate_color("#ffffff", "#000000", 1.0 - 0.5 / 40)
    low_end, high_end = steps[0][1], steps[-1][1]
    assert (low_end, high_end) == ((near_mid, near_end) if direction == "forward" else (near_end, near_mid))


def _gradient_svg(mode: str) -> str:
    record = SeqRecord(Seq("ACGT" * 2500), id="rec1", annotations={"molecule_type": "DNA", "topology": "circular"})
    forward = RegionAnnotation(
        "fwd", CoordinateSpan(None, 1000, 4000), mark="highlight", legend_label=QUERY_POSITION_LEGEND_LABEL,
        style=RegionAnnotationStyle(fill="#5b1a4a", gradient=AnnotationGradient("#f6c9d8", "#5b1a4a", "forward")),
    )
    reverse = RegionAnnotation(
        "rev", CoordinateSpan(None, 6000, 8000), mark="highlight", legend_label=QUERY_POSITION_LEGEND_LABEL,
        style=RegionAnnotationStyle(fill="#5b1a4a", gradient=AnnotationGradient("#f6c9d8", "#5b1a4a", "reverse")),
    )
    annotations = AnnotationOptions(sets=(AnnotationSet("Search 1", (forward, reverse)),))
    if mode == "circular":
        drawing = build_circular_diagram(record, options=CircularDiagramOptions(annotations=annotations))
    else:
        drawing = build_linear_diagram([record], options=LinearDiagramOptions(annotations=annotations))
    return drawing.tostring()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_gradient_row_is_drawn_from_pale_to_dark_along_the_hit(mode: str) -> None:
    svg = _gradient_svg(mode)
    for annotation_id, pale_first in (("fwd", True), ("rev", False)):
        item = re.search(rf'data-gbdraw-annotation-id="{annotation_id}".*?</g>\s*</g>', svg, re.S)
        assert item is not None
        steps = re.search(r'<g opacity="0.2">(.*?)</g>', item.group(0), re.S)
        assert steps is not None, "the steps sit in one group at the fill opacity"
        fills = re.findall(r'fill="(#[0-9a-f]{6})"', steps.group(1))
        assert len(fills) > 10 and fills == (sorted(fills, reverse=True) if pale_first else sorted(fills))
        # The mark itself keeps its outline and has no fill of its own.
        assert 'fill="none"' in item.group(0)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_gradient_hits_give_one_query_position_legend_row(mode: str) -> None:
    svg = _gradient_svg(mode)
    # One row for both hits; a Linear Legend is drawn in its horizontal and its vertical form.
    assert svg.count(f">{QUERY_POSITION_LEGEND_LABEL}<") == (1 if mode == "circular" else 2)
    assert ">start<" in svg and ">end<" in svg


FIXTURE = Path(__file__).resolve().parent / "fixtures" / "sessions" / "python-search-hits-gradient.v51.gbdraw-session.json.gz"


def test_session_fixture_keeps_the_gradient_through_load_and_save(tmp_path: Path) -> None:
    from gbdraw.session import load_session_document, materialize_session

    stored = json.loads(gzip.decompress(FIXTURE.read_bytes()))["drawings"][0]["renderRequest"]
    with materialize_session(load_session_document(FIXTURE), output_directory=tmp_path) as materialized:
        request = decode_canonical_request(
            stored, resource_paths=materialized.resource_paths, output_directory=tmp_path
        )
        encoded = encode_canonical_request(request).payload["diagramOptions"]["annotations"]["sets"]
    sets = request.options.annotations.sets
    assert [row.style.gradient for row in sets[0].annotations] == [None, None, None]
    assert [row.style.gradient.direction for row in sets[1].annotations] == ["forward", "forward", "reverse"]
    assert [(row.style.gradient.query_start, row.style.gradient.query_end) for row in sets[1].annotations] == [
        (0.0, 1.0), (16 / 650, 171 / 650), (319 / 650, 365 / 650),
    ]
    assert encoded == stored["diagramOptions"]["annotations"]["sets"]


def test_queries_are_read_as_nucleotide_fasta_records(tmp_path: Path) -> None:
    from gbdraw.web_support.sequence_search import read_search_queries_json

    path = tmp_path / "queries.fasta"
    path.write_text(">q1 first\nacgtac\ngtn\n>q2\nTTTT\n", encoding="utf-8")
    reply = json.loads(read_search_queries_json(str(path)))
    assert reply == {"queries": [
        {"name": "q1", "length": 9, "fasta": ">q1\nACGTACGTN\n"},
        {"name": "q2", "length": 4, "fasta": ">q2\nTTTT\n"},
    ]}


@pytest.mark.parametrize(
    ("text", "code", "reason"),
    [
        (">p\nMKVLHHW\n", "INPUT_INVALID", "NOT_NUCLEOTIDE"),
        ("not a sequence\n", "NO_RECORDS", "FASTA_NO_ENTRIES"),
        (">q\nACGT\n>q\nACGT\n", "INPUT_INVALID", "UNIQUE_IDS"),
    ],
)
def test_a_query_that_cannot_be_searched_is_refused_with_its_reason(
    tmp_path: Path, text: str, code: str, reason: str
) -> None:
    from gbdraw.web_support.sequence_search import read_search_queries

    path = tmp_path / "queries.fasta"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValidationError) as error:
        read_search_queries(str(path))
    assert error.value.diagnostic is not None
    assert (error.value.diagnostic["code"], error.value.diagnostic.get("reason")) == (code, reason)


def test_a_search_set_keeps_its_highlight_overlay_through_the_request(tmp_path: Path) -> None:
    result = _search()
    assert result.annotation_set is not None and result.annotation_set.keep_highlight_overlay
    record = SeqRecord(Seq("ACGT" * 500), id="rec1", annotations={"molecule_type": "DNA"})
    searched = AnnotationSet("Search 1", tuple(
        row for row in result.annotation_set.annotations if row.target.record.record_id == "rec1"
    ), keep_highlight_overlay=True)
    plain = AnnotationSet("plain", (RegionAnnotation("p", CoordinateSpan(None, 10, 50), mark="band"),))
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(annotations=AnnotationOptions(sets=(searched, plain))),
    )
    encoded = encode_canonical_request(request)
    sets = encoded.payload["diagramOptions"]["annotations"]["sets"]
    # Written only on the search's set; every other set is written as before.
    assert sets[0]["keepHighlightOverlay"] is True and "keepHighlightOverlay" not in sets[1]
    resource_paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        target.write_bytes(resource.content or b"")
        resource_paths[resource.resource_id] = target
    decoded = decode_canonical_request(encoded.payload, resource_paths=resource_paths, output_directory=tmp_path / "out")
    assert [item.keep_highlight_overlay for item in decoded.options.annotations.sets] == [True, False]
    bad = json.loads(json.dumps(encoded.payload))
    bad["diagramOptions"]["annotations"]["sets"][0]["keepHighlightOverlay"] = "yes"
    with pytest.raises(Exception, match="keepHighlightOverlay must be a boolean"):
        decode_canonical_request(bad, resource_paths=resource_paths, output_directory=tmp_path / "out")


def _custom_stack_svg(mode: str, *, keep: bool) -> str:
    from gbdraw.api.options import CircularRequestTrackOptions, LinearRequestTrackOptions
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    record = SeqRecord(Seq("ACGT" * 2500), id="rec1", annotations={"molecule_type": "DNA", "topology": "circular"})
    hit = RegionAnnotation("q#1", CoordinateSpan(None, 1000, 4000), mark="highlight",
                           style=RegionAnnotationStyle(fill="#c8607e"))
    annotations = AnnotationOptions(sets=(AnnotationSet("Search 1", (hit,), keep_highlight_overlay=keep),))
    if mode == "circular":
        options = CircularDiagramOptions(annotations=annotations, tracks=CircularRequestTrackOptions(circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"), CircularTrackSlot(id="ticks", renderer="ticks"),
        ]))
        return build_circular_diagram(record, options=options).tostring()
    options = LinearDiagramOptions(annotations=annotations, tracks=LinearRequestTrackOptions(linear_track_slots=[
        LinearTrackSlot(id="features", renderer="features"),
    ]))
    return build_linear_diagram([record], options=options).tostring()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_custom_stack_draws_a_search_set_behind_its_features(mode: str) -> None:
    # PD-OI-116: with a custom track stack a search's set is drawn without a slot of its own.
    assert 'data-gbdraw-annotation-id="q#1"' in _custom_stack_svg(mode, keep=True)
    # Another set keeps following the stack: no slot names it, so it is not drawn.
    assert 'data-gbdraw-annotation-id="q#1"' not in _custom_stack_svg(mode, keep=False)


# Receipt web.sequence-search (D-49, D-50): each hit sits on its strand. With separated strands a plus-strand hit
# covers the plus-strand lane and a minus-strand hit the minus-strand lane; otherwise, and for a row without a
# strand, a hit covers the whole feature track.
def test_each_hit_row_carries_its_strand() -> None:
    result = _search()
    assert result.annotation_set is not None
    assert [row.strand for row in result.annotation_set.annotations] == ["+", "-", "+"]
    assert [hit.strand for hit in result.hits] == ["+", "-", "+"]


def test_a_row_strand_is_typed_and_round_trips_through_the_request(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="strand") as error:
        RegionAnnotation("h", CoordinateSpan(None, 1, 10), mark="highlight", strand="x")  # type: ignore[arg-type]
    assert error.value.diagnostic == {"code": "INPUT_INVALID", "reason": "FIELDS"}
    record = SeqRecord(Seq("ACGT" * 500), id="rec1", annotations={"molecule_type": "DNA"})
    rows = (
        RegionAnnotation("p", CoordinateSpan(None, 10, 50), mark="highlight", strand="+"),
        RegionAnnotation("m", CoordinateSpan(None, 60, 90), mark="highlight", strand="-"),
        RegionAnnotation("n", CoordinateSpan(None, 100, 150), mark="band"),
    )
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(annotations=AnnotationOptions(sets=(AnnotationSet("s", rows),))),
    )
    encoded = encode_canonical_request(request)
    items = encoded.payload["diagramOptions"]["annotations"]["sets"][0]["annotations"]
    # Written only when set, so a row without one is written as before.
    assert [item.get("strand") for item in items] == ["+", "-", None] and "strand" not in items[2]
    resource_paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        target.write_bytes(resource.content or b"")
        resource_paths[resource.resource_id] = target
    decoded = decode_canonical_request(encoded.payload, resource_paths=resource_paths, output_directory=tmp_path / "out")
    assert [row.strand for row in decoded.options.annotations.sets[0].annotations] == ["+", "-", None]
    bad = json.loads(json.dumps(encoded.payload))
    bad["diagramOptions"]["annotations"]["sets"][0]["annotations"][0]["strand"] = "plus"
    with pytest.raises(Exception, match="strand"):
        decode_canonical_request(bad, resource_paths=resource_paths, output_directory=tmp_path / "out")


@pytest.mark.parametrize(("strand", "source_step", "shown"), [("+", 1, "+"), ("-", 1, "-"), ("+", -1, "-"), ("-", -1, "+"), (None, -1, None)])
def test_a_row_strand_is_shown_on_the_record_as_drawn(strand: str | None, source_step: int, shown: str | None) -> None:
    # A reverse-complemented record draws its source's plus strand as its minus strand, as its features.
    annotation = ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=((100, 200),), midpoint_bp=150.0, span_bp=100, label="",
        mark="highlight", lane=None, style=RegionAnnotationStyle(), source_step=source_step,  # type: ignore[arg-type]
        strand=strand,  # type: ignore[arg-type]
    )
    assert annotation.display_strand == shown


_PLUS_CDS, _MINUS_CDS = (300, 1500), (5000, 6200)


def _strand_record(plus: tuple[int, int] = _PLUS_CDS, minus: tuple[int, int] = _MINUS_CDS) -> SeqRecord:
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    record = SeqRecord(Seq("ACGT" * 2500), id="rec1", annotations={"molecule_type": "DNA", "topology": "circular"})
    record.features = [
        SeqFeature(FeatureLocation(*plus, strand=1), type="CDS"),
        SeqFeature(FeatureLocation(*minus, strand=-1), type="CDS"),
    ]
    return record


def _strand_hits(plus: tuple[int, int] = _PLUS_CDS, minus: tuple[int, int] = _MINUS_CDS) -> AnnotationOptions:
    style = RegionAnnotationStyle(fill="#c8607e", stroke="#c8607e", stroke_width=0.001)
    rows = (
        RegionAnnotation("p", CoordinateSpan(None, plus[0] + 101, plus[1] - 100), mark="highlight",
                         style=style, strand="+"),
        RegionAnnotation("m", CoordinateSpan(None, minus[0] + 101, minus[1] - 100), mark="highlight",
                         style=style, strand="-"),
        RegionAnnotation("n", CoordinateSpan(None, 8000, 9000), mark="highlight", style=style),
    )
    return AnnotationOptions(sets=(AnnotationSet("Search 1", rows, keep_highlight_overlay=True),))


def _arc_radii(element) -> tuple[float, float]:
    radii = [float(value) for value in re.findall(r"A ?([0-9.]+),", element.attrib.get("d", ""))]
    return min(radii), max(radii)


def _circular_strand_drawing(track_type: str, stranded: bool, custom: bool):
    from xml.etree import ElementTree

    from gbdraw.api.options import CircularRequestTrackOptions

    tracks = (
        CircularRequestTrackOptions(circular_track_slots=("features:features", "ticks:ticks"))
        if custom else None
    )
    drawing = build_circular_diagram(_strand_record(), options=CircularDiagramOptions(
        annotations=_strand_hits(),
        config_overrides={"canvas.strandedness": stranded, "canvas.circular.track_type": track_type},
        **({"tracks": tracks} if tracks is not None else {}),
    ))
    root = ElementTree.fromstring(drawing.tostring())
    features = [element for element in root.iter() if "data-gbdraw-feature-id" in element.attrib]
    hits = {
        element.attrib["data-gbdraw-annotation-id"]: next(path for path in element.iter() if path.tag.endswith("path"))
        for element in root.iter()
        if element.attrib.get("data-gbdraw-annotation-id") in {"p", "m", "n"}
    }
    return [_arc_radii(element) for element in features], {key: _arc_radii(path) for key, path in hits.items()}


@pytest.mark.parametrize("custom", [False, True], ids=["automatic", "custom"])
@pytest.mark.parametrize("track_type", ["tuckin", "middle", "spreadout"])
@pytest.mark.parametrize("stranded", [True, False], ids=["separated", "shared"])
def test_a_circular_hit_covers_its_strands_lane(track_type: str, stranded: bool, custom: bool) -> None:
    lanes, hits = _circular_strand_drawing(track_type, stranded, custom)
    assert set(hits) == {"p", "m", "n"}
    whole = (min(lane[0] for lane in lanes), max(lane[1] for lane in lanes))
    assert hits["n"] == pytest.approx(whole, abs=0.01)
    if not stranded:
        assert hits["p"] == pytest.approx(whole, abs=0.01) and hits["m"] == pytest.approx(whole, abs=0.01)
        return
    plus, minus = lanes
    assert plus != pytest.approx(minus, abs=0.5)
    assert hits["p"] == pytest.approx(plus, abs=0.01)
    assert hits["m"] == pytest.approx(minus, abs=0.01)


def _absolute_y_ranges(root, elements) -> list[tuple[float, float]]:
    parents = {child: parent for parent in root.iter() for child in parent}

    def offset(element) -> float:
        total = 0.0
        while element is not None:
            total += sum(
                float(value)
                for value in re.findall(r"translate\([^,)]+,\s*([^)]+)\)", element.attrib.get("transform", ""))
            )
            element = parents.get(element)
        return total

    ranges = []
    for element in elements:
        if element.tag.endswith("rect"):
            top = float(element.attrib["y"])
            values = [top, top + float(element.attrib["height"])]
        else:
            values = [float(y) for _x, y in re.findall(r"(-?[0-9.]+)[ ,](-?[0-9.]+)", element.attrib.get("d", ""))]
        shift = offset(element)
        ranges.append((min(values) + shift, max(values) + shift))
    return ranges


def _linear_strand_drawing(track_layout: str, stranded: bool, custom: bool, **wrapped: object):
    from xml.etree import ElementTree

    from gbdraw.api.options import LinearRequestTrackOptions

    tracks = LinearRequestTrackOptions(linear_track_slots=("features:features",)) if custom else None
    spans = {key: wrapped.pop(key) for key in ("plus", "minus") if key in wrapped}
    drawing = build_linear_diagram([_strand_record(**spans)], options=LinearDiagramOptions(  # type: ignore[arg-type]
        annotations=_strand_hits(**spans),  # type: ignore[arg-type]
        config_overrides={"canvas.strandedness": stranded, "canvas.linear.track_layout": track_layout, **wrapped},
        **({"tracks": tracks} if tracks is not None else {}),
    ))
    root = ElementTree.fromstring(drawing.tostring())
    features = [
        element for element in root.iter()
        if "data-gbdraw-feature-id" in element.attrib and element.tag.endswith("path")
    ]
    hits = {
        element.attrib["data-gbdraw-annotation-id"]: next(rect for rect in element.iter() if rect.tag.endswith("rect"))
        for element in root.iter()
        if element.attrib.get("data-gbdraw-annotation-id") in {"p", "m", "n"}
    }
    names = list(hits)
    return _absolute_y_ranges(root, features), dict(zip(names, _absolute_y_ranges(root, [hits[key] for key in names])))


def _inside(inner: tuple[float, float], outer: tuple[float, float]) -> bool:
    return outer[0] - 1e-6 <= inner[0] and inner[1] <= outer[1] + 1e-6


def _apart(first: tuple[float, float], second: tuple[float, float]) -> bool:
    return first[1] <= second[0] + 1e-6 or second[1] <= first[0] + 1e-6


@pytest.mark.parametrize("custom", [False, True], ids=["automatic", "custom"])
@pytest.mark.parametrize("track_layout", ["middle", "above", "below"])
@pytest.mark.parametrize("stranded", [True, False], ids=["separated", "shared"])
def test_a_linear_hit_covers_its_strands_lane(track_layout: str, stranded: bool, custom: bool) -> None:
    lanes, hits = _linear_strand_drawing(track_layout, stranded, custom)
    assert set(hits) == {"p", "m", "n"}
    plus, minus = lanes
    for lane in lanes:
        assert _inside(lane, hits["n"])
    if not stranded:
        assert hits["p"] == pytest.approx(hits["n"]) and hits["m"] == pytest.approx(hits["n"])
        return
    assert _apart(plus, minus)
    assert _inside(plus, hits["p"]) and _apart(hits["p"], minus)
    assert _inside(minus, hits["m"]) and _apart(hits["m"], plus)


def test_a_wrapped_linear_record_draws_a_hit_on_its_strands_lane_of_its_row() -> None:
    # Both features and both hits on the second of three rows.
    lanes, hits = _linear_strand_drawing(
        "above", True, False, plus=(4300, 5000), minus=(5500, 6200), **{"canvas.linear.wrap_length": 4000}
    )
    plus, minus = lanes
    assert _apart(plus, minus)
    assert _inside(plus, hits["p"]) and _apart(hits["p"], minus)
    assert _inside(minus, hits["m"]) and _apart(hits["m"], plus)


# OV-2571, D-52: a circular record drawn whole is also searched over its origin window (its last w and first w
# bases, w the longest query), so a hit across the origin is one row with wraps_origin, its halves in the record's
# own subject are dropped, and "the thresholds removed N" counts each alignment once.
_MITO_LENGTH = 16569
_WINDOW_START = _MITO_LENGTH - 600 + 1
ORIGIN_SUBJECTS = (
    SearchSubject(f"NC_012920.1:1-{_MITO_LENGTH}", parse_record_selector("NC_012920.1"), "mito"),
    SearchSubject(
        f"NC_012920.1:{_WINDOW_START}-600", parse_record_selector("NC_012920.1"), "mito",
        source_start=_WINDOW_START, record_length=_MITO_LENGTH,
    ),
)
ORIGIN_QUERIES = (SearchQuery("origin600", 600), SearchQuery("origin50", 50), SearchQuery("inner", 100))
_MAIN, _WINDOW = (subject.subject_id for subject in ORIGIN_SUBJECTS)
# The rows LOSAT writes for the last and first 300 bases of NC_012920.1 (origin600), and for the last and first 25
# (origin50, whose halves score under the bitscore threshold of 50): two halves on the record, one row on the window.
ORIGIN_ROWS = (
    ("origin600", _MAIN, "100.000", "300", "0", "0", "301", "600", "1", "300", "1e-150", "555"),
    ("origin600", _MAIN, "100.000", "300", "0", "0", "1", "300", "16270", "16569", "1e-150", "555"),
    ("origin600", _WINDOW, "100.000", "600", "0", "0", "1", "600", "301", "900", "0.0", "1109"),
    ("origin50", _MAIN, "100.000", "25", "0", "0", "26", "50", "1", "25", "3e-06", "46.4"),
    ("origin50", _MAIN, "100.000", "25", "0", "0", "1", "25", "16545", "16569", "3e-06", "46.4"),
    ("origin50", _WINDOW, "100.000", "50", "0", "0", "1", "50", "576", "625", "1e-20", "92.7"),
    # A hit near the origin that does not cross it is found on the record and on the window: read once.
    ("inner", _MAIN, "100.000", "100", "0", "0", "1", "100", "16301", "16400", "1e-50", "185"),
    ("inner", _WINDOW, "100.000", "100", "0", "0", "1", "100", "332", "431", "1e-50", "185"),
)


def _origin_search(rows=ORIGIN_ROWS, thresholds=THRESHOLDS):
    return search_hit_annotation_set(
        read_comparison_table(io.StringIO("\n".join("\t".join(row) for row in rows))),
        set_id="Search 1", subjects=ORIGIN_SUBJECTS, thresholds=thresholds, program="blastn", queries=ORIGIN_QUERIES,
    )


def test_a_hit_across_the_origin_is_one_row_and_its_halves_are_dropped() -> None:
    result = _origin_search()
    assert result.annotation_set is not None
    targets = {row.id: row.target for row in result.annotation_set.annotations}
    assert set(targets) == {"origin600#1", "origin50#1", "inner#1"}
    assert (targets["origin600#1"].start, targets["origin600#1"].end, targets["origin600#1"].wraps_origin) == (
        16270, 300, True,
    )
    assert (targets["inner#1"].start, targets["inner#1"].end, targets["inner#1"].wraps_origin) == (16301, 16400, False)
    hit = next(hit for hit in result.hits if hit.query == "origin600")
    # The window row's values: the whole alignment.
    assert (hit.start, hit.end, hit.strand, hit.alignment_length, hit.bitscore, hit.evalue) == (
        16270, 300, "+", 600, 1109.0, 0.0,
    )
    assert [hit.query for hit in result.hits].count("origin600") == 1


def test_a_short_hit_across_the_origin_is_found_when_the_whole_alignment_passes() -> None:
    result = _origin_search()
    [hit] = [hit for hit in result.hits if hit.query == "origin50"]
    assert (hit.start, hit.end, hit.bitscore) == (16545, 25, 92.7)


def test_a_minus_strand_hit_across_the_origin_keeps_its_strand() -> None:
    minus = (
        ("origin600", _MAIN, "100.000", "300", "0", "0", "1", "300", "300", "1", "1e-150", "555"),
        ("origin600", _MAIN, "100.000", "300", "0", "0", "301", "600", "16569", "16270", "1e-150", "555"),
        ("origin600", _WINDOW, "100.000", "600", "0", "0", "1", "600", "900", "301", "0.0", "1109"),
    )
    [hit] = _origin_search(minus).hits
    assert (hit.start, hit.end, hit.strand) == (16270, 300, "-")


def test_the_thresholds_removed_count_counts_an_origin_alignment_once() -> None:
    # Every alignment once: the two crossing hits and the inner hit, not their halves or the inner hit's twin.
    assert (_origin_search().tally.rows_read, _origin_search().tally.rows_kept) == (3, 3)
    strict = ComparisonThresholds(evalue=1e-2, bitscore=5000.0, identity=70.0, alignment_length=0)
    result = _origin_search(thresholds=strict)
    assert result.annotation_set is None
    assert result.rows_removed_by_thresholds == 3


def test_a_half_that_passes_stands_for_a_crossing_alignment_that_does_not() -> None:
    # The whole alignment is under the identity threshold, one half is over it: the half is the hit, and the
    # alignment is counted once.
    rows = (
        ("origin600", _MAIN, "100.000", "300", "0", "0", "301", "600", "1", "300", "1e-150", "555"),
        ("origin600", _MAIN, "40.000", "300", "0", "0", "1", "300", "16270", "16569", "1e-5", "60"),
        ("origin600", _WINDOW, "69.000", "600", "0", "0", "1", "600", "301", "900", "1e-100", "600"),
    )
    result = _origin_search(rows)
    assert [(hit.start, hit.end) for hit in result.hits] == [(1, 300)]
    assert (result.tally.rows_read, result.tally.rows_kept) == (1, 1)


def test_an_origin_window_spans_the_last_and_first_bases_of_at_most_half_the_record() -> None:
    assert origin_window_span(_MITO_LENGTH, 600) == (_WINDOW_START, 600)
    assert origin_window_span(_MITO_LENGTH, 20000) == (_MITO_LENGTH - 8284 + 1, 8284)
    with pytest.raises(ValidationError):
        origin_window_span(1, 10)
    with pytest.raises(ValidationError):
        SearchSubject("w", parse_record_selector("r"), "r", source_start=1, record_length=100)


def test_the_worker_helper_reads_an_origin_window_subject() -> None:
    selector = {"kind": "recordId", "value": "NC_012920.1"}
    payload = {
        "table": "\n".join("\t".join(row) for row in ORIGIN_ROWS[:3]), "setId": "Search 1",
        "subjects": [
            {"subjectId": _MAIN, "record": selector, "recordKey": "mito", "sourceStart": 1,
             "wrapsOrigin": False, "recordLength": _MITO_LENGTH},
            {"subjectId": _WINDOW, "record": selector, "recordKey": "mito", "sourceStart": _WINDOW_START,
             "wrapsOrigin": True, "recordLength": _MITO_LENGTH},
        ],
        "thresholds": {"evalue": 1e-2, "bitscore": 50, "identity": 70, "alignmentLength": 0},
        "program": "blastn", "queries": [{"name": "origin600", "length": 600}],
        "colorMode": "query", "queryColors": {},
    }
    reply = json.loads(sequence_search_annotation_set_json(json.dumps(payload)))
    [row] = reply["annotationSet"]["annotations"]
    assert (row["target"]["start"], row["target"]["end"], row["target"]["wrapsOrigin"]) == (16270, 300, True)
    assert [(hit["start"], hit["end"]) for hit in reply["hits"]] == [(16270, 300)]
    assert (reply["rowsRead"], reply["rowsKept"]) == (1, 1)


# LOSAT 0.1.0 (blastn -task blastn -outfmt 6) on the subjects the Web sends for NC_012920.1 drawn whole with a 600 bp
# longest query: the record from base 1 and its origin window 15970-600 (prb/2d2/losat-origin-window.tsv).
LOSAT_ORIGIN_TABLE = (
    "origin600\tNC_012920.1:15970-600\t100.000\t600\t0\t0\t1\t600\t301\t900\t0.0\t1083\n"
    "origin600\tNC_012920.1:15970-600\t100.000\t12\t0\t0\t508\t519\t819\t808\t1.3\t22.9\n"
    "origin600\tNC_012920.1:15970-600\t100.000\t11\t0\t0\t39\t49\t185\t195\t4.5\t21.1\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t300\t0\t0\t301\t600\t1\t300\t5.85e-157\t542\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t300\t0\t0\t1\t300\t16270\t16569\t5.85e-157\t542\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t13\t0\t0\t310\t322\t2275\t2287\t0.37\t24.7\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t12\t0\t0\t508\t519\t219\t208\t1.3\t22.9\n"
    "origin600\tNC_012920.1:1-16569\t88.235\t17\t2\t0\t12\t28\t1991\t2007\t1.3\t22.9\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t12\t0\t0\t449\t460\t4310\t4299\t1.3\t22.9\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t12\t0\t0\t321\t332\t14786\t14797\t1.3\t22.9\n"
    "origin600\tNC_012920.1:1-16569\t93.750\t16\t0\t1\t502\t517\t14795\t14781\t4.5\t22.0\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t11\t0\t0\t343\t353\t6083\t6093\t4.5\t21.1\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t11\t0\t0\t261\t271\t13140\t13150\t4.5\t21.1\n"
    "origin600\tNC_012920.1:1-16569\t100.000\t11\t0\t0\t39\t49\t16154\t16164\t4.5\t21.1\n"
    "origin50\tNC_012920.1:15970-600\t100.000\t50\t0\t0\t1\t50\t576\t625\t1.93e-22\t91.5\n"
    "origin50\tNC_012920.1:1-16569\t100.000\t25\t0\t0\t26\t50\t1\t25\t7.22e-09\t46.4\n"
    "origin50\tNC_012920.1:1-16569\t100.000\t25\t0\t0\t1\t25\t16545\t16569\t7.22e-09\t46.4\n"
    "origin50\tNC_012920.1:1-16569\t100.000\t13\t0\t0\t35\t47\t2275\t2287\t0.024\t24.7\n"
)


def test_losat_rows_of_an_origin_window_give_each_alignment_once() -> None:
    result = _origin_search(tuple(tuple(line.split("\t")) for line in LOSAT_ORIGIN_TABLE.splitlines()))
    assert [(hit.query, hit.start, hit.end, hit.strand) for hit in result.hits] == [
        ("origin600", 16270, 300, "+"), ("origin50", 16545, 25, "+"),
    ]
    # The window rows that do not cross the origin (twins of record rows) and the four halves are not read; the ten
    # short record rows are, and the thresholds remove them.
    assert (result.tally.rows_read, result.tally.rows_kept) == (12, 2)
