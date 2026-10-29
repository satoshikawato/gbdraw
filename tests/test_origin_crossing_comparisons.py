"""OV-2258: a CDS across a circular record's origin covers only its own arc in comparisons."""
from __future__ import annotations

from xml.etree import ElementTree as ET

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity import (
    LosslessCollinearityParameters,
    build_collinearity_blocks_from_hits,
    convert_collinearity_blocks_to_comparisons,
)
from gbdraw.analysis.collinearity_units import build_collinearity_unit_index
from gbdraw.analysis.origin_arcs import cover_arcs, location_origin_arc
from gbdraw.analysis.protein_colinearity import (
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS

LENGTH = 720


def _record(record_id: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(
        Seq("ATGAAATAG" * 80),
        id=record_id,
        annotations={"topology": "circular", "molecule_type": "DNA"},
    )
    record.features = features
    return record


def _cds(location, name: str) -> SeqFeature:
    return SeqFeature(
        location,
        type="CDS",
        qualifiers={"translation": ["MKK"], "locus_tag": [name], "protein_id": [name]},
    )


def _origin_cds(name: str, strand: int) -> SeqFeature:
    tail, head = FeatureLocation(660, LENGTH, strand), FeatureLocation(0, 30, strand)
    parts = [tail, head] if strand == 1 else [head, tail]
    return _cds(CompoundLocation(parts), name)


def _records(strand: int) -> list[SeqRecord]:
    query = _record(
        "query",
        [_cds(FeatureLocation(start, start + 90, 1), f"a{index}") for index, start in enumerate((100, 200, 300))],
    )
    # b0 crosses the origin; b1 and b2 follow it on the plus strand, or precede it on the minus strand.
    others = [
        _cds(FeatureLocation(40, 130, strand), "b1"),
        _cds(FeatureLocation(140, 230, strand), "b2"),
    ]
    subject = _record("subject", [_origin_cds("b0", strand), *others])
    return [query, subject]


def _hit(query: str, subject: str) -> dict[str, object]:
    return {
        "query": query, "subject": subject, "identity": 90.0, "alignment_length": 30,
        "mismatches": 0, "gap_opens": 0, "qstart": 1, "qend": 30, "sstart": 1,
        "send": 30, "evalue": 1e-30, "bitscore": 200.0,
    }


def _pairs(strand: int) -> list[tuple[str, str]]:
    if strand == 1:
        return [("a0", "b0"), ("a1", "b1"), ("a2", "b2")]
    return [("a0", "b2"), ("a1", "b1"), ("a2", "b0")]


def _blocks(strand: int):
    records = _records(strand)
    extraction = extract_cds_proteins(records)
    forward = pd.DataFrame.from_records([_hit(q, s) for q, s in _pairs(strand)], columns=COMPARISON_COLUMNS)
    reverse = pd.DataFrame.from_records([_hit(s, q) for q, s in _pairs(strand)], columns=COMPARISON_COLUMNS)
    result = build_collinearity_blocks_from_hits(
        [forward], extraction, records=records,
        params=LosslessCollinearityParameters(min_anchors=3, max_unit_gap=1),
        reverse_hits_by_pair=[reverse],
    )
    return records, extraction, result


def _spans(frame: pd.DataFrame) -> list[tuple[int, int, int, int]]:
    return [tuple(int(value) for value in row) for row in frame[["qstart", "qend", "sstart", "send"]].itertuples(index=False)]


@pytest.mark.parametrize("strand", [1, -1])
def test_protein_and_unit_cover_only_their_arc(strand: int) -> None:
    records = _records(strand)
    extraction = extract_cds_proteins(records)
    protein = extraction.protein_map["b0"]
    assert protein.origin_arc == (660, 30)
    # The linear extent stays the location's extent (identity and order keys read it).
    assert (protein.start, protein.end) == (0, LENGTH)
    unit = build_collinearity_unit_index(extraction, records=records).unit_by_protein_id["b0"]
    assert (unit.start, unit.end) == (660, 30)
    assert unit.order == 0


@pytest.mark.parametrize("strand", [1, -1])
def test_block_across_the_origin_is_cut_into_rows_inside_the_record(strand: int) -> None:
    records, _extraction, result = _blocks(strand)
    assert len(result.blocks) == 1
    assert result.blocks[0].orientation == ("plus" if strand == 1 else "minus")
    frame = convert_collinearity_blocks_to_comparisons(result, records=records)[0]
    if strand == 1:
        # Query 101..390 runs with subject 661..720 then 1..230.
        assert _spans(frame) == [(101, 160, 661, 720), (161, 390, 1, 230)]
    else:
        # Query 101..390 runs with subject 230..1 then 720..661.
        assert _spans(frame) == [(101, 330, 230, 1), (331, 390, 720, 661)]
    assert set(frame["collinearity_block_id"]) == {result.blocks[0].block_id}
    assert set(frame["collinearity_anchor_count"]) == {3}


@pytest.mark.parametrize("strand", [1, -1])
def test_pairwise_link_of_an_origin_cds_is_cut_at_the_origin(strand: int) -> None:
    extraction = extract_cds_proteins(_records(strand))
    hits = pd.DataFrame.from_records([_hit("a0", "b0")], columns=COMPARISON_COLUMNS)
    frame = convert_protein_hits_to_genomic_links(hits, extraction.protein_map)
    if strand == 1:
        assert _spans(frame) == [(101, 160, 661, 720), (161, 190, 1, 30)]
    else:
        assert _spans(frame) == [(101, 130, 30, 1), (131, 190, 720, 661)]
    assert list(frame["subject_protein_id"]) == ["b0", "b0"]


def test_locations_and_covers_without_an_origin_keep_their_linear_extent() -> None:
    assert location_origin_arc([(10, 20, 1), (30, 40, 1)], LENGTH) is None
    assert location_origin_arc([(30, 40, -1), (10, 20, -1)], LENGTH) is None
    assert location_origin_arc([(10, 20, 1), (30, 40, -1)], LENGTH) is None
    assert location_origin_arc([(660, LENGTH, 1), (0, 30, 1)], LENGTH) == (660, 30)
    assert location_origin_arc([(0, 30, -1), (660, LENGTH, -1)], LENGTH) == (660, 30)
    assert cover_arcs([(10, 20), (5, 40)]) == (5, 40)
    assert cover_arcs([(660, 30), (40, 130), (500, 600)]) == (500, 130)
    assert cover_arcs([(600, 300), (200, 650)]) is None


def _linear_svg(records, result, starts, reverses) -> ET.Element:
    config = apply_config_overrides(None, {
        "canvas.show_gc": False, "canvas.show_skew": False, "labels.linear.scope": "none",
    })
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(
                InMemoryRecordSource(record),
                display=RecordDisplayOptions(start_coordinate=start),
                presentation=RecordPresentation(reverse_complement=reverse),
            )
            for record, start, reverse in zip(records, starts, reverses, strict=True)
        ),
        options=LinearDiagramOptions(config=config, collinearity_blocks=result),
    )
    result_ = plan_request(request).build()
    drawing = result_.drawing if hasattr(result_, "drawing") else result_
    return ET.fromstring(drawing.tostring())


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("starts,reverses", [
    ((None, None), (False, False)),
    ((None, 361), (False, False)),
    ((None, None), (False, True)),
    ((51, 361), (True, True)),
])
def test_drawn_block_never_spans_the_whole_record(strand, starts, reverses) -> None:
    records, _extraction, result = _blocks(strand)
    root = _linear_svg(records, result, starts, reverses)
    paths = [node for node in root.iter() if node.get("data-collinearity-block-id")]
    assert paths
    for path in paths:
        for role in ("q", "s"):
            first, last = int(path.get(f"data-{role}start")), int(path.get(f"data-{role}end"))
            assert abs(last - first) < LENGTH / 2


# OV-2259: a Similarity alignment anchored on a CDS across the origin centres on
# the middle of the CDS's own arc, on every path that reads that centre.


def _alignment_collection(strand: int, start: int | None, reverse: bool):
    from gbdraw.api.record_planning import ResolvedRecordCollection

    record = _record("subject", [_origin_cds("b0", strand), _cds(FeatureLocation(40, 130, strand), "b1")])
    request = LinearDiagramRequest(
        records=(RecordInput(
            InMemoryRecordSource(record),
            display=RecordDisplayOptions(start_coordinate=start),
            presentation=RecordPresentation(reverse_complement=reverse),
        ),),
        options=LinearDiagramOptions(),
    )
    plan = plan_request(request)
    return ResolvedRecordCollection(plan.records, plan.provenance)


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("start,reverse,center", [
    # b0 covers 661..720 and 1..30, so its middle is source boundary 705.
    (None, False, 705.0),
    (361, False, 345.0),
    (None, True, 15.0),
    (361, True, 376.0),
])
def test_alignment_centre_of_an_origin_cds_is_the_middle_of_its_arc(strand, start, reverse, center) -> None:
    from gbdraw.analysis.protein_colinearity import OrthogroupMember
    from gbdraw.api.record_planning import (
        _similarity_candidate_from_orthogroup_member,
        project_similarity_alignment_anchor_fact,
    )
    from gbdraw.layout.similarity_alignment import AlignmentAnchorIdentity

    collection = _alignment_collection(strand, start, reverse)
    provenance = collection.provenance[0]
    feature = provenance.source_feature_catalog[0]
    anchor = AlignmentAnchorIdentity(
        record_key=provenance.record_key,
        biological_feature_id=feature.biological_feature_id,
        source_feature_index=feature.source_feature_index,
        stable_feature_svg_id=feature.stable_feature_id,
    )
    fact = project_similarity_alignment_anchor_fact(collection, anchor)
    assert fact.display_center == pytest.approx(center)
    # The member path resolves the same feature and reads the same centre.
    member = OrthogroupMember(
        orthogroup_id="og1", protein_id="b0", record_index=0,
        feature_index=feature.source_feature_index, record_id="subject", label="b0",
        start=0, end=LENGTH, strand=strand, feature_svg_id=feature.stable_feature_id,
        source_protein_id="b0",
    )
    candidate = _similarity_candidate_from_orthogroup_member(collection, member)
    assert candidate.display_center == pytest.approx(center)


def test_alignment_centre_without_an_origin_keeps_the_middle_of_the_extent() -> None:
    from gbdraw.api.record_planning import project_similarity_alignment_anchor_fact
    from gbdraw.layout.similarity_alignment import AlignmentAnchorIdentity

    collection = _alignment_collection(1, None, False)
    provenance = collection.provenance[0]
    feature = provenance.source_feature_catalog[1]
    anchor = AlignmentAnchorIdentity(
        record_key=provenance.record_key,
        biological_feature_id=feature.biological_feature_id,
        source_feature_index=feature.source_feature_index,
        stable_feature_svg_id=feature.stable_feature_id,
    )
    assert project_similarity_alignment_anchor_fact(collection, anchor).display_center == 85.0


# OV-2260: parts that step back once away from the origin (a minus-strand CDS
# listed in ascending order, as GFF3 input gives it, or a trans-spliced gene)
# keep their linear extent; only a step from the record's end to its start
# crosses the origin.


@pytest.mark.parametrize("topology", ["linear", "circular"])
@pytest.mark.parametrize("locations,extent", [
    ([FeatureLocation(100, 130, -1), FeatureLocation(210, 261, -1)], (100, 261)),
    ([FeatureLocation(500, 560, 1), FeatureLocation(100, 130, 1)], (100, 560)),
    ([FeatureLocation(660, 700, 1), FeatureLocation(0, 30, 1)], (0, 700)),
])
def test_parts_that_step_back_away_from_the_origin_keep_their_extent(topology, locations, extent) -> None:
    record = _record("subject", [_cds(CompoundLocation(locations), "c0")])
    record.annotations["topology"] = topology
    protein = extract_cds_proteins([record]).protein_map["c0"]
    assert protein.origin_arc is None
    assert (protein.start, protein.end) == extent
