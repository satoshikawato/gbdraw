"""One exact contact rule for a join drawn as one block over a circular record's origin (OV-3661, OV-3659, OV-3689).

The feature factory decides it once per feature and view (``origin_join_order``, ``FeatureObject.joins_over_origin``);
the Circular drawer, the label anchors and lane placement read that decision or the parts built from it.
"""
from __future__ import annotations

from collections import Counter
from xml.etree import ElementTree as ET

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import build_circular_multi_diagram
from gbdraw.api.options import CircularDiagramOptions, CircularMultiRecordOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest, GenBankInputSource, GffFastaInputSource, InMemoryRecordSource, RecordDisplayOptions,
    RecordInput,
)
from gbdraw.features.factory import FeatureLayerMode, create_feature_layers
from gbdraw.labels.circular_candidates import _candidate_segment
from gbdraw.labels.coordinates import display_label_segment
from gbdraw.layout.record_coordinates import CirclePlace

_LENGTH = 3000
_SEQUENCE = ("ACGTTGCA" * 400)[:_LENGTH]


def _record(*features: tuple[str, object], topology: str = "circular") -> SeqRecord:
    record = SeqRecord(Seq(_SEQUENCE), id="rec", name="rec", annotations={
        "topology": topology, "molecule_type": "DNA",
    })
    record.features = [
        SeqFeature(location, type="CDS", qualifiers={"locus_tag": [tag]}) for tag, location in features
    ]
    return record


def _join(*parts: tuple[int, int], strand: int = 1) -> CompoundLocation:
    return CompoundLocation([SimpleLocation(start, end, strand=strand) for start, end in parts])


def _plan(source, start: int | None):
    cfg = apply_config_overrides(None, {
        "labels.circular.scope": "none", "canvas.show_gc": False, "canvas.show_skew": False,
    })
    return plan_request(CircularDiagramRequest(
        records=(RecordInput(source, display=RecordDisplayOptions(start_coordinate=start)),),
        options=CircularDiagramOptions(config=cfg),
    ))


def _drawn_parts(source, start: int | None) -> Counter[str]:
    """How many paths of each kind (block, connector) the Circular drawing draws for its features."""
    result = _plan(source, start).build()
    svg = result.drawing if hasattr(result, "drawing") else result
    return Counter(
        node.get("data-gbdraw-feature-part") for node in ET.fromstring(svg.tostring()).iter()
        if node.get("data-gbdraw-feature-part")
    )


def _features(record: SeqRecord, mode: str, start: int | None = None):
    """The feature objects of ``record`` as a ``mode`` ("circular" or "linear") diagram builds them, lanes resolved."""
    plan = _plan(InMemoryRecordSource(record), start)
    layer_mode = (
        FeatureLayerMode(joins_over_origin=True) if mode == "circular" else FeatureLayerMode(past_end_as_join=True)
    )
    result = create_feature_layers(
        plan.records[0], {}, ["CDS"], {"CDS": "gray"}, False, True, {},
        layer_mode=layer_mode, record_transform=plan.transforms[0],
    )
    return {feature.qualifiers["locus_tag"][0]: feature for feature in result.foreground_features.values()}


_IN_ORDER = ("ends", _join((0, 300), (2700, 3000)))
_MID = ("mid", SimpleLocation(1200, 1800, strand=1))


@pytest.mark.parametrize("start", [None, 1001])
def test_a_join_written_in_order_over_both_ends_draws_one_block_over_the_origin_in_every_circular_view(start):
    # OV-3661: the whole view joined the two blocks; a display start drew them with a connector all round the ring.
    parts = _drawn_parts(InMemoryRecordSource(_record(_IN_ORDER)), start)
    assert parts["connector"] == 0
    # The whole view draws one block; a display start draws two that meet at the origin.
    assert parts["block"] == (1 if start is None else 2)


@pytest.mark.parametrize("start", [None, 1001])
def test_sorted_gff3_rows_of_a_gene_over_the_origin_draw_one_block_in_every_circular_view(tmp_path, start):
    gff = tmp_path / "rows.gff3"
    gff.write_text(
        "##gff-version 3\n##sequence-region rows 1 3000\n"
        "rows\t.\tregion\t1\t3000\t.\t+\t.\tID=rows;Is_circular=true\n"
        "rows\t.\tCDS\t1\t300\t.\t+\t0\tID=cds1;locus_tag=rowsorigin\n"
        "rows\t.\tCDS\t2701\t3000\t.\t+\t0\tID=cds1;locus_tag=rowsorigin\n"
    )
    fasta = tmp_path / "rows.fa"
    fasta.write_text(f">rows\n{_SEQUENCE}\n")
    parts = _drawn_parts(GffFastaInputSource(gff, fasta), start)
    assert parts["connector"] == 0


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("start", [None, 1001])
def test_a_join_over_the_origin_reads_one_anchor_and_ends_where_its_strand_does(strand, start):
    features = _features(_record(("ends", _join((0, 300), (2700, 3000), strand=strand))), "circular", start)
    feature = features["ends"]
    assert feature.joins_over_origin
    if start is None:
        start_bp, end_bp, _, span, label_strand = _candidate_segment(feature, _LENGTH)
        assert (start_bp, end_bp, span) == (2700, 300, 600)
        # OV-3868: the joined block's label reads the strand it is drawn on (dev read "undefined").
        assert label_strand == ("positive" if strand == 1 else "negative")
        # The arrowhead is at the 3' end over the origin: base 300 on the plus strand, base 2701 on the minus.
        last = next(part for part in feature.location if part.is_last)
        assert (last.start, last.end) == ((0, 300) if strand == 1 else (2700, 3000))
    else:
        segment = display_label_segment(feature, _LENGTH, circular=True)
        assert segment is not None and segment.span == 600
        assert not any(part.kind == "line" for part in feature.display_parts)


@pytest.mark.parametrize("start", [None, 1001])
def test_a_one_base_gap_at_the_first_base_is_drawn_with_its_line_in_every_circular_view(start):
    # OV-3659: join(2901..3000,2..100). The whole view joined the blocks over base 1 (a part from base 2 counted as
    # base 1); the display start drew the gap and its 1-base line.
    record = _record(("gapone", _join((2900, 3000), (1, 100))))
    parts = _drawn_parts(InMemoryRecordSource(record), start)
    assert parts["connector"] == 1
    assert parts["block"] == 2
    feature = _features(record, "circular", start)["gapone"]
    assert not feature.joins_over_origin
    # The label goes to the longest part, 2901..3000.
    if start is None:
        assert _candidate_segment(feature, _LENGTH)[:2] == (2900, 3000)
    else:
        segment = display_label_segment(feature, _LENGTH, circular=True)
        assert segment is not None and segment.source_span == (2900, 3000)


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_a_feature_under_the_connector_of_a_join_in_order_takes_its_own_lane_in_linear_mode(topology):
    # OV-3689: Linear mode draws join(1..300,2701..3000) as written, its line across 301..2700.
    features = _features(_record(_IN_ORDER, _MID, topology=topology), "linear")
    assert features["ends"].feature_track_id != features["mid"].feature_track_id
    assert not features["ends"].joins_over_origin


def test_a_circular_join_that_is_not_one_block_reserves_its_connector_line():
    # Three parts in order over both ends: the Circular view draws them as written, lines across the record.
    record = _record(("three", _join((0, 300), (1000, 1100), (2700, 3000))), ("mid", SimpleLocation(1300, 1800, 1)))
    features = _features(record, "circular")
    assert features["three"].feature_track_id != features["mid"].feature_track_id
    assert not features["three"].joins_over_origin


def test_a_join_over_the_origin_leaves_the_middle_of_the_record_to_other_features():
    features = _features(_record(_IN_ORDER, _MID), "circular")
    assert features["ends"].joins_over_origin
    assert features["ends"].feature_track_id == features["mid"].feature_track_id


def test_a_linear_record_drawn_circular_draws_the_join_as_written():
    # A linear record has no origin to join over: its join keeps the line across the record.
    features = _features(_record(_IN_ORDER, _MID, topology="linear"), "circular")
    assert not features["ends"].joins_over_origin
    assert any(part.kind == "line" for part in features["ends"].location)


def test_a_two_strand_join_at_both_ends_is_not_one_block():
    features = _features(_record(("mixed", CompoundLocation([
        SimpleLocation(0, 300, strand=-1), SimpleLocation(2700, 3000, strand=1),
    ]))), "circular")
    assert not features["mixed"].joins_over_origin


def test_a_genbank_join_over_the_origin_keeps_its_coordinates_as_written(tmp_path):
    # The drawing order is the factory's; the feature's coordinates (selectors, comparison hashes) stay as written.
    path = tmp_path / "ends.gb"
    SeqIO.write(_record(_IN_ORDER), path, "genbank")
    assert _drawn_parts(GenBankInputSource(path), None)["connector"] == 0
    feature = _features(_record(_IN_ORDER), "circular")["ends"]
    assert [(int(part.start), int(part.end)) for part in feature.coordinates] == [(0, 300), (2700, 3000)]


def _on_a_shared_circle(record: SeqRecord, start: int | None):
    """The feature objects of ``record`` drawn on a circle it shares with another record (one-circle layout)."""
    plan = _plan(InMemoryRecordSource(record), start)
    place = CirclePlace(offset=100, length=5200, scale_length=_LENGTH)
    result = create_feature_layers(
        plan.records[0], {}, ["CDS"], {"CDS": "gray"}, False, True, {},
        layer_mode=FeatureLayerMode(joins_over_origin=True),
        record_transform=plan.transforms[0].with_circle(place, _LENGTH),
    )
    return {feature.qualifiers["locus_tag"][0]: feature for feature in result.foreground_features.values()}


_OVER = ("over", _join((2700, 3000), (0, 300)))


def test_a_record_on_a_shared_circle_draws_a_join_in_order_as_written_inside_its_arc():
    # Review r1 HIGH 1: the record's ends lie next to the gaps to the other records, so its origin is no point of
    # the drawing. The join keeps its connector inside the record's arc, as on dev, and its label stands on a part.
    features = _on_a_shared_circle(_record(_IN_ORDER, _MID), None)
    ends = features["ends"]
    assert not ends.joins_over_origin
    assert [part.kind for part in ends.display_parts] == ["block", "line", "block"]
    segment = display_label_segment(ends, 5200, circular=True)
    assert segment is not None and segment.span == 300
    # Its connector line runs under "mid", which takes a lane of its own.
    assert ends.feature_track_id != features["mid"].feature_track_id


def test_the_label_of_a_join_written_over_the_origin_of_a_record_on_a_shared_circle_stands_on_a_part():
    # OV-3864: dev anchored the label at the middle of both parts, the record's arc end, under the next record.
    feature = _on_a_shared_circle(_record(_OVER), None)["over"]
    assert not feature.joins_over_origin
    segment = display_label_segment(feature, 5200, circular=True)
    assert segment is not None and segment.span == 300
    assert 100 < segment.middle < 3100


@pytest.mark.parametrize("location", [_IN_ORDER, _OVER])
def test_a_display_start_inside_a_record_on_a_shared_circle_joins_over_the_origin_it_draws(location):
    # From display start 1001 the record's origin lies inside its arc: the parts meet there, as on its own circle.
    feature = _on_a_shared_circle(_record(location), 1001)[location[0]]
    assert feature.joins_over_origin
    assert not any(part.kind == "line" for part in feature.display_parts)
    segment = display_label_segment(feature, 5200, circular=True)
    assert segment is not None and segment.span == 600


def test_a_one_circle_drawing_keeps_the_connector_of_a_join_in_order_inside_the_record_arc():
    other = SeqRecord(Seq(_SEQUENCE[:2000]), id="other", name="other", annotations={
        "topology": "circular", "molecule_type": "DNA",
    })
    drawing = build_circular_multi_diagram(
        [_record(_IN_ORDER), other], options=CircularDiagramOptions(),
        layout=CircularMultiRecordOptions(record_layout="one_circle"),
    )
    parts = Counter(
        node.get("data-gbdraw-feature-part") for node in ET.fromstring(drawing.tostring()).iter()
        if node.get("data-gbdraw-feature-part")
    )
    assert parts["connector"] == 1
