"""OV-2258: a CDS across a circular record's origin covers only its own arc in comparisons."""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity import (
    LosslessCollinearityParameters,
    build_collinearity_blocks_from_hits,
    build_orthogroup_collinearity_blocks_from_hits,
    convert_collinearity_blocks_to_comparisons,
)
from gbdraw.analysis.collinearity_units import build_collinearity_unit_index
from gbdraw.analysis.origin_arcs import cover_arcs, location_origin_arc
from gbdraw.analysis.protein_colinearity import (
    CdsProtein,
    convert_pair_protein_hits_to_genomic_links,
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
    select_rbh_orthogroup_edges_from_directional_hits,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import project_collinearity_frames
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
from gbdraw.io.regions import crop_then_reverse_records, parse_region_specs
from tests.utils.comparison_tables import converted_pair_rows

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
        infer_orthogroups=True,
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


# The Web path: python-helpers.js extracts each record alone, sends the protein
# map through the page, and converts the raw LOSATP tables from that map.


def _web_helpers() -> dict[str, object]:
    helpers_js = (Path(__file__).parents[1] / "gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(helpers_js.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    return namespace


def _web_records(strand: int, reverse: bool) -> list[SeqRecord]:
    """``_records``, whose subject a reversed view draws in the order that forms a block."""
    query, subject = _records(strand)
    return [query, crop_then_reverse_records([subject], [], [reverse])[0][0]]


def _web_conversion(
    tmp_path: Path, stage, mode: str, strand: int, *, reverse: bool = False, subject_region: str | None = None,
) -> list[tuple[int, int, int, int]]:
    """The rows the Web converter stores (search frame), from a real extraction of written files."""
    helpers = _web_helpers()
    locus_to_id: dict[str, str] = {}
    payload_records = []
    for index, record in enumerate(_web_records(strand, reverse)):
        path = tmp_path / f"{record.id}.gb"
        SeqIO.write([record], path, "genbank")
        extracted = json.loads(json.dumps(json.loads(helpers["extract_cds_protein_fasta"](
            str(path), "genbank",
            region_spec=subject_region if index == 1 else None,
            record_index=index, record_instance_key=f"r_{index + 1}",
        ))))
        assert "error" not in extracted
        for protein_id, entry in extracted["protein_map"].items():
            locus_to_id[str(entry["locus_tag"])] = protein_id
        payload_records.append({
            "recordIndex": index,
            "recordId": extracted["record_id"],
            "proteinMap": extracted["protein_map"],
            "runtimeBindingHash": extracted["runtime_binding_hash"],
            "displayBindingHash": extracted["display_binding_hash"],
            "viewTransform": {"length": extracted["record_length"], "reverse": reverse and index == 1},
        })

    def table(rows: list[tuple[str, str]]) -> str:
        return pd.DataFrame.from_records(
            [_hit(locus_to_id[q], locus_to_id[s]) for q, s in rows if q in locus_to_id and s in locus_to_id],
            columns=COMPARISON_COLUMNS,
        ).to_csv(sep="\t", header=False, index=False, lineterminator="\n")

    pairs = _pairs(strand)
    payload = {
        "records": payload_records,
        "pairs": [
            {"pairIndex": 0, "queryIndex": 0, "subjectIndex": 1, "displayPair": True,
             "cacheKey": "forward", "blastText": table(pairs)},
            # Pairwise searches each displayed pair in one direction.
            *([] if mode == "pairwise" else [
                {"pairIndex": 1, "queryIndex": 1, "subjectIndex": 0, "displayPair": False,
                 "cacheKey": "reverse", "blastText": table([(s, q) for q, s in pairs])},
            ]),
        ],
    }
    pairs_path, raw_path = stage(tmp_path, payload)
    settings = {
        "collinear_min_anchors": 3, "collinear_max_unit_gap": 1, "collinear_anchor_mode": "rbh",
        "collinear_merge_orientation": "strand", "collinear_infer_orthogroups": False,
    }
    result = json.loads(str(helpers["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path), str(raw_path), mode, 5, 50, "1e-2", 0, 0, **(settings if mode == "collinear" else {}),
    )))
    assert "error" not in result, result.get("error")
    return [
        (int(row["qstart"]), int(row["qend"]), int(row["sstart"]), int(row["send"]))
        for pair in result["pairs"]
        for row in converted_pair_rows(pair)
    ]


def _cli_conversion(mode: str, strand: int, *, reverse: bool = False, region: str | None = None) -> list[tuple[int, int, int, int]]:
    """The rows of the CLI's path for the same records, written in the search frame.

    The CLI extracts from the drawn record, so a reversed subject's rows are mapped
    back with x -> L + 1 - x, the rule of the planner's search-frame projection.
    """
    query, subject = _web_records(strand, reverse)
    subject = crop_then_reverse_records([subject], parse_region_specs([region] if region else []), [False])[0][0]
    drawn = crop_then_reverse_records([subject], [], [reverse])[0][0]
    records = [query, drawn]
    extraction = extract_cds_proteins(records)
    pairs = [(q, s) for q, s in _pairs(strand) if s in extraction.protein_map]
    tables = {
        (0, 1): pd.DataFrame.from_records([_hit(q, s) for q, s in pairs], columns=COMPARISON_COLUMNS),
        (1, 0): pd.DataFrame.from_records([_hit(s, q) for q, s in pairs], columns=COMPARISON_COLUMNS),
    }
    if mode == "collinear":
        result = build_orthogroup_collinearity_blocks_from_hits(
            tables, extraction,
            params=LosslessCollinearityParameters(min_anchors=3, max_unit_gap=1),
            edge_mode="rbh", search_scope="adjacent", infer_orthogroups=False, comparison_pairs=((0, 1),),
        )
        frames = project_collinearity_frames(
            result, records, color_mode="orientation", search_scope="adjacent", comparison_pairs=((0, 1),),
        ).pairs
        frame = frames[0].matches if frames else pd.DataFrame(columns=["qstart", "qend", "sstart", "send"])
    elif mode == "pairwise":
        frame = convert_pair_protein_hits_to_genomic_links(
            tables[(0, 1)], extraction.protein_map, extraction.protein_map,
        )
    else:
        selection = select_rbh_orthogroup_edges_from_directional_hits(
            tables, extraction.protein_map, comparison_pairs=((0, 1),),
        )
        frame = convert_pair_protein_hits_to_genomic_links(
            selection.adjacent_display_edges_by_pair[(0, 1)], extraction.protein_map, extraction.protein_map,
            orthogroups=selection.orthogroups,
        )
    rows = _spans(frame)
    if not reverse:
        return rows
    length = len(subject)
    return [(qstart, qend, length + 1 - sstart, length + 1 - send) for qstart, qend, sstart, send in rows]


MODES = ["collinear", "pairwise", "orthogroup"]


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_web_rows_of_an_origin_cds_equal_the_cli_rows(tmp_path, stage_web_losatp_transport, mode, strand, reverse) -> None:
    rows = _web_conversion(tmp_path, stage_web_losatp_transport, mode, strand, reverse=reverse)
    assert rows and sorted(rows) == sorted(_cli_conversion(mode, strand, reverse=reverse))
    # The origin CDS (61 + 30 bp arc) never reads as the whole 720 bp record.
    assert all(abs(sstart - send) < LENGTH / 2 for _, _, sstart, send in rows)


@pytest.mark.parametrize("mode", MODES)
def test_web_rows_of_a_crop_equal_the_cli_rows(tmp_path, stage_web_losatp_transport, mode) -> None:
    rows = _web_conversion(tmp_path, stage_web_losatp_transport, mode, 1, subject_region="subject:1-400")
    assert rows and sorted(rows) == sorted(_cli_conversion(mode, 1, region="subject:1-400"))


# ``sequence`` travels in the FASTA, read back by _web_ordered_proteins_from_fasta,
# so the map stays small. Every other field a comparison can read crosses in the map.
TRANSPORT_EXEMPT = {"sequence"}


@pytest.mark.parametrize("strand", [1, -1])
def test_every_cds_protein_field_crosses_the_web_transport(strand: int) -> None:
    helpers = _web_helpers()
    records = _records(strand)
    extraction = extract_cds_proteins(records)
    subject = records[1]
    raw_map = json.loads(json.dumps({
        protein.protein_id: helpers["_serialize_cds_protein"](protein, subject)
        for protein in extraction.proteins_by_record[1]
    }))
    rebuilt = helpers["_build_web_cds_protein_map"](raw_map)
    missing = {
        field.name
        for protein in extraction.proteins_by_record[1]
        for field in dataclasses.fields(CdsProtein)
        if field.name not in TRANSPORT_EXEMPT
        and getattr(rebuilt[protein.protein_id], field.name) != getattr(protein, field.name)
    }
    assert not missing
    assert any(protein.origin_arc is not None for protein in rebuilt.values())


@pytest.mark.parametrize("strand", [1, -1])
def test_web_display_map_flips_the_arc_of_a_reversed_record(strand: int) -> None:
    helpers = _web_helpers()
    records = _records(strand)
    flipped = crop_then_reverse_records([records[1]], [], [True])[0][0]
    source = extract_cds_proteins([records[1]]).proteins_by_record[0]
    expected = extract_cds_proteins([flipped]).proteins_by_record[0]
    raw_map = json.loads(json.dumps({
        protein.protein_id: helpers["_serialize_cds_protein"](protein, records[1]) for protein in source
    }))
    display = helpers["_build_display_web_cds_protein_map"](raw_map, {"length": LENGTH, "reverse": True})
    assert sorted((protein.arc, protein.strand) for protein in display.values()) == sorted(
        (protein.arc, protein.strand) for protein in expected
    )
