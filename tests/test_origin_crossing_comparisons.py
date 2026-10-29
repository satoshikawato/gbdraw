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
    assert location_origin_arc([(10, 20, 1), (30, 40, 1)]) is None
    assert location_origin_arc([(30, 40, -1), (10, 20, -1)]) is None
    assert location_origin_arc([(10, 20, 1), (30, 40, -1)]) is None
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
