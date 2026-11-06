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
            {"subjectId": "rec1:1-5000", "record": {"kind": "recordId", "value": "rec1"}, "recordKey": "record-1", "sourceStart": 1},
            {"subjectId": "rec2:1-3000", "record": {"kind": "recordId", "value": "rec2"}, "recordKey": "record-2", "sourceStart": 1001},
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
    steps = gradient_steps(annotation, gradient, 0, 100.0, 200.0, record_length=1000, length_px=40.0)
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
