from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.analysis.collinearity as collinearity_module
import gbdraw.diagrams.linear.assemble as linear_assemble_module
from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityParameters,
    CollinearityResult,
    LosslessCollinearityParameters,
    build_collinearity_blocks_from_hits,
    calculate_collinearity_block_evalue,
    call_collinearity_blocks,
    cluster_lossless_collinearity_anchors,
    convert_collinearity_blocks_to_comparisons,
    deduplicate_unit_pair_anchors,
    iter_collinearity_search_pairs,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
    orthogroup_edges_to_lossless_collinearity_anchors,
)
from gbdraw.api import assemble_linear_diagram_from_records
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.protein_colinearity import (
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
    select_rbh_orthogroup_edges_from_directional_hits,
)
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators.blast import BlastMatchConfigurator
from gbdraw.core.color import interpolate_color
from gbdraw.core.text import calculate_bbox_dimensions
from gbdraw.io.collinearity import parse_native_collinearity_tsv, write_native_collinearity_tsv
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.legend.table import prepare_legend_table
from gbdraw.render.groups.linear.legend import LegendGroup
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup
from gbdraw.exceptions import ValidationError


def _record(record_id: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(Seq("ATGAAATAG" * 80), id=record_id)
    record.features = features
    return record


def _cds(start: int, end: int, qualifiers: dict[str, list[str]]) -> SeqFeature:
    merged = {"translation": ["MKK"]}
    merged.update(qualifiers)
    return SeqFeature(FeatureLocation(start, end, strand=1), type="CDS", qualifiers=merged)


def _hit_row(
    query: str,
    subject: str,
    bitscore: float = 200.0,
    *,
    alignment_length: int = 100,
    qend: int | None = None,
    send: int | None = None,
) -> dict[str, object]:
    qend = alignment_length if qend is None else qend
    send = alignment_length if send is None else send
    return {
        "query": query,
        "subject": subject,
        "identity": 90.0,
        "alignment_length": alignment_length,
        "mismatches": 0,
        "gap_opens": 0,
        "qstart": 1,
        "qend": qend,
        "sstart": 1,
        "send": send,
        "evalue": 1e-30,
        "bitscore": bitscore,
    }


def _first_fasta_id(fasta_text: str) -> str:
    for line in str(fasta_text).splitlines():
        if line.startswith(">"):
            return line[1:].split(None, 1)[0]
    return ""


def _hex_distance(color_a: str, color_b: str) -> float:
    rgb_a = tuple(int(color_a[index : index + 2], 16) for index in (1, 3, 5))
    rgb_b = tuple(int(color_b[index : index + 2], 16) for index in (1, 3, 5))
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(rgb_a, rgb_b)))


def _translate_xy(transform: str | None) -> tuple[float, float]:
    assert transform is not None
    assert transform.startswith("translate(")
    parts = transform.removeprefix("translate(").removesuffix(")").replace(",", " ").split()
    assert len(parts) >= 2
    return float(parts[0]), float(parts[1])


def _translate_xy_or_zero(transform: str | None) -> tuple[float, float]:
    if transform is None:
        return 0.0, 0.0
    return _translate_xy(transform)


def _build_collinearity_match_group() -> PairWiseMatchGroup:
    group = PairWiseMatchGroup.__new__(PairWiseMatchGroup)
    group.canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=1000,
    )
    group.records = [_record("record_a", []), _record("record_b", [])]
    group.comparison_count = 1
    group.comparison_height = 40
    group.query_offset_x = 0
    group.subject_offset_x = 0
    group.query_alignment_offset_x = 0
    group.subject_alignment_offset_x = 0
    group.min_identity = 0
    group.match_min_color = "#ffffff"
    group.match_max_color = "#000000"
    group.match_fill_opacity = 0.75
    group.match_stroke_color = "none"
    group.match_stroke_width = 0
    group.collinearity_orientation_colors = {"plus": "#112233", "minus": "#445566"}
    group.collinearity_orientation_min_colors = {"plus": "#eeeeee", "minus": "#ffeeee"}
    return group


def _web_protein_entry(
    protein_id: str,
    *,
    record_index: int,
    record_id: str,
    feature_index: int,
    start: int,
    end: int,
) -> dict[str, object]:
    return {
        "protein_id": protein_id,
        "record_index": record_index,
        "feature_index": feature_index,
        "record_id": record_id,
        "start": start,
        "end": end,
        "strand": 1,
        "label": protein_id,
        "protein_length": 30,
        "feature_svg_id": f"feature_{protein_id}",
    }


def _anchor(
    query_order: int,
    subject_order: int,
    *,
    bitscore: float = 200.0,
    query_strand: int | None = 1,
    subject_strand: int | None = 1,
    query_start: int | None = None,
    query_end: int | None = None,
    subject_start: int | None = None,
    subject_end: int | None = None,
) -> CollinearityAnchor:
    return CollinearityAnchor(
        query_protein_id=f"q{query_order}",
        subject_protein_id=f"s{subject_order}",
        query_record_index=0,
        subject_record_index=1,
        query_order=query_order,
        subject_order=subject_order,
        query_start=query_start if query_start is not None else query_order * 10 + 1,
        query_end=query_end if query_end is not None else query_order * 10 + 9,
        subject_start=subject_start if subject_start is not None else subject_order * 10 + 1,
        subject_end=subject_end if subject_end is not None else subject_order * 10 + 9,
        query_strand=query_strand,
        subject_strand=subject_strand,
        identity=90.0,
        evalue=1e-30,
        bitscore=bitscore,
        alignment_length=100,
        query_feature_svg_id=f"fq{query_order}",
        subject_feature_svg_id=f"fs{subject_order}",
        source="test",
        query_unit_id=f"qu{query_order}",
        subject_unit_id=f"su{subject_order}",
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=None,
        subject_locus_id=None,
        query_display_name=f"q{query_order}",
        subject_display_name=f"s{subject_order}",
    )


@pytest.mark.linear
def test_deduplicate_unit_pair_anchors_keeps_strongest_hit() -> None:
    weak = _anchor(0, 0, bitscore=100)
    strong = _anchor(0, 0, bitscore=300)

    deduplicated = deduplicate_unit_pair_anchors([weak, strong])

    assert deduplicated == (strong,)


@pytest.mark.linear
def test_collinearity_parameters_validates_block_evalue() -> None:
    CollinearityParameters(block_evalue=None).validate()
    CollinearityParameters(block_evalue=0).validate()
    CollinearityParameters(block_evalue=1e-3).validate()

    with pytest.raises(ValidationError, match="collinear_block_evalue"):
        CollinearityParameters(block_evalue=-1).validate()


@pytest.mark.linear
def test_collinearity_blocks_are_called_without_block_evalue_acceptance() -> None:
    anchors = [_anchor(index, index) for index in range(8)]
    result = call_collinearity_blocks(
        anchors,
        params=CollinearityParameters(min_anchors=5, max_gene_gap=1),
    )

    assert len(result.blocks) == 1
    assert result.blocks[0].kind == "syntenic"
    assert result.blocks[0].block_evalue is None


@pytest.mark.linear
def test_collinearity_block_evalue_uses_genomic_coordinates_not_unit_order() -> None:
    positions = [100, 110, 120, 130, 1100]
    anchors = [
        _anchor(
            index,
            index,
            query_start=position,
            query_end=position + 8,
            subject_start=position,
            subject_end=position + 8,
        )
        for index, position in enumerate(positions)
    ]

    observed = calculate_collinearity_block_evalue(
        SimpleNamespace(orientation="plus", score=250.0, anchors=tuple(anchors)),
        anchors,
    )
    m = len(positions)
    expected_log = math.log(2.0) + math.lgamma(m + 1)
    region_length = positions[-1] - positions[0]
    for previous, current in zip(positions, positions[1:]):
        expected_log += math.log(current - previous) * 2.0
    expected_log -= (m - 1) * (math.log(region_length) * 2.0)
    expected = math.exp(expected_log)
    old_unit_order_value = 2.0 * math.factorial(m) / ((m - 1) ** (2 * (m - 1)))

    assert observed == pytest.approx(expected)
    assert observed < 1e-5
    assert old_unit_order_value > 1e-5


@pytest.mark.linear
def test_pairwise_singleton_anchor_is_emitted_by_default() -> None:
    result = call_collinearity_blocks(
        [_anchor(0, 0)],
        params=CollinearityParameters(),
    )

    assert len(result.blocks) == 1
    assert result.blocks[0].kind == "singleton"
    assert result.blocks[0].anchors[0].query_protein_id == "q0"


@pytest.mark.linear
def test_min_anchors_drops_singleton_blocks() -> None:
    anchor = _anchor(0, 0)
    result = call_collinearity_blocks(
        [anchor],
        params=CollinearityParameters(min_anchors=2),
    )

    assert result.blocks == ()
    assert result.unblocked_anchors == (anchor,)


@pytest.mark.linear
def test_native_block_caller_detects_plus_and_minus_orientation() -> None:
    params = CollinearityParameters(min_anchors=3, max_gene_gap=1, block_evalue=None)

    plus = call_collinearity_blocks([_anchor(0, 0), _anchor(1, 1), _anchor(2, 2)], params=params)
    minus = call_collinearity_blocks([_anchor(0, 2), _anchor(1, 1), _anchor(2, 0)], params=params)

    assert len(plus.blocks) == 1
    assert plus.blocks[0].orientation == "plus"
    assert [anchor.query_order for anchor in plus.blocks[0].anchors] == [0, 1, 2]
    assert len(minus.blocks) == 1
    assert minus.blocks[0].orientation == "minus"


@pytest.mark.linear
def test_single_anchor_block_orientation_uses_feature_strands() -> None:
    params = CollinearityParameters(min_anchors=1, max_gene_gap=0, block_evalue=2)

    same_strand = call_collinearity_blocks(
        [_anchor(0, 0, query_strand=1, subject_strand=1)],
        params=params,
    )
    opposite_strand = call_collinearity_blocks(
        [_anchor(0, 0, query_strand=1, subject_strand=-1)],
        params=params,
    )

    assert same_strand.blocks[0].orientation == "plus"
    assert opposite_strand.blocks[0].orientation == "minus"


@pytest.mark.linear
def test_collinearity_color_mode_defaults_to_orientation_and_aliases_identity() -> None:
    assert normalize_collinearity_color_mode(None) == "orientation"
    assert normalize_collinearity_color_mode("identity") == "average_identity"
    assert normalize_collinearity_color_mode("orientation_identity") == "orientation_identity"
    assert normalize_collinearity_color_mode("orientation-identity") == "orientation_identity"


@pytest.mark.linear
def test_collinearity_anchor_mode_defaults_to_rbh_and_aliases_top1() -> None:
    assert normalize_collinearity_anchor_mode(None) == "rbh"
    assert normalize_collinearity_anchor_mode("top1") == "one_to_one"
    assert normalize_collinearity_anchor_mode("mutual-best") == "one_to_one"
    assert normalize_collinearity_anchor_mode("reciprocal-best") == "rbh"


@pytest.mark.linear
def test_collinearity_search_scope_defaults_to_adjacent() -> None:
    assert normalize_collinearity_search_scope(None) == "adjacent"
    assert normalize_collinearity_search_scope("all") == "all"

    with pytest.raises(ValidationError, match="collinear_search_scope"):
        normalize_collinearity_search_scope("pairwise")


@pytest.mark.linear
def test_iter_collinearity_search_pairs_adjacent_returns_neighboring_pairs() -> None:
    assert iter_collinearity_search_pairs(5, scope="adjacent") == (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
    )


@pytest.mark.linear
def test_iter_collinearity_search_pairs_all_returns_every_unordered_pair() -> None:
    assert iter_collinearity_search_pairs(5, scope="all") == (
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
    )


@pytest.mark.linear
def test_build_blocks_from_hits_one_to_one_removes_query_to_many_noise() -> None:
    records = [
        _record(
            "record_a",
            [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]}),
                _cds(12, 21, {"locus_tag": ["sb1"], "protein_id": ["sb1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame.from_records(
        [_hit_row("qa0", "sb0", bitscore=300), _hit_row("qa0", "sb1", bitscore=250)],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row("sb0", "qa0", bitscore=300), _hit_row("sb1", "qa0", bitscore=250)],
        columns=COMPARISON_COLUMNS,
    )
    params = CollinearityParameters(
        min_anchors=1,
        max_gene_gap=0,
        nearby_duplicate_window=0,
        block_evalue=2,
    )

    all_hits = build_collinearity_blocks_from_hits(
        [hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="all",
        reverse_hits_by_pair=[reverse_hits],
    )
    one_to_one = build_collinearity_blocks_from_hits(
        [hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="one_to_one",
        reverse_hits_by_pair=[reverse_hits],
    )

    assert len(all_hits.blocks) == 2
    assert all_hits.blocks[0].anchors[0].subject_protein_id == "sb0"
    assert all_hits.blocks[0].anchors[0].subject_orthogroup_member_count == 2
    assert len(one_to_one.blocks) == 1
    assert one_to_one.blocks[0].anchors[0].subject_protein_id == "sb0"


@pytest.mark.linear
def test_build_blocks_from_hits_rbh_uses_reverse_best_hits() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]}),
                _cds(12, 21, {"locus_tag": ["qa1"], "protein_id": ["qa1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]}),
                _cds(12, 21, {"locus_tag": ["sb1"], "protein_id": ["sb1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    forward_hits = pd.DataFrame.from_records(
        [_hit_row("qa0", "sb0", bitscore=300), _hit_row("qa1", "sb1", bitscore=250)],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row("sb0", "qa0", bitscore=300), _hit_row("sb1", "qa0", bitscore=400)],
        columns=COMPARISON_COLUMNS,
    )
    params = CollinearityParameters(
        min_anchors=1,
        max_gene_gap=0,
        nearby_duplicate_window=0,
        block_evalue=2,
    )

    one_to_one = build_collinearity_blocks_from_hits(
        [forward_hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="one_to_one",
    )
    rbh = build_collinearity_blocks_from_hits(
        [forward_hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="rbh",
        reverse_hits_by_pair=[reverse_hits],
    )

    assert sum(len(block.anchors) for block in one_to_one.blocks) == 2
    assert len(rbh.blocks) == 1
    assert rbh.blocks[0].anchors[0].query_protein_id == "qa0"
    assert rbh.blocks[0].anchors[0].subject_protein_id == "sb0"


@pytest.mark.linear
def test_lossless_collinearity_preserves_every_adjacent_orthogroup_edge() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(index * 12, index * 12 + 9, {"locus_tag": [protein_id], "protein_id": [protein_id]})
                for index, protein_id in enumerate(["a0", "a1", "a2", "a3", "a4"])
            ],
        ),
        _record(
            "record_b",
            [
                _cds(index * 12, index * 12 + 9, {"locus_tag": [protein_id], "protein_id": [protein_id]})
                for index, protein_id in enumerate(["b0", "b1", "b2", "b3", "b4", "b5"])
            ],
        ),
        _record(
            "record_c",
            [_cds(0, 9, {"locus_tag": ["c0"], "protein_id": ["c0"]})],
        ),
    ]
    extraction = extract_cds_proteins(records)

    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a0", "b0"),
                _hit_row("a1", "b1"),
                _hit_row("a2", "b4"),
                _hit_row("a3", "b5"),
                _hit_row("a4", "b2"),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("b0", "a0"),
                _hit_row("b1", "a1"),
                _hit_row("b4", "a2"),
                _hit_row("b5", "a3"),
                _hit_row("b2", "a4"),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (0, 2): pd.DataFrame.from_records([_hit_row("a3", "c0")], columns=COMPARISON_COLUMNS),
        (2, 0): pd.DataFrame.from_records([_hit_row("c0", "a3")], columns=COMPARISON_COLUMNS),
        (1, 2): pd.DataFrame.from_records([_hit_row("b4", "c0")], columns=COMPARISON_COLUMNS),
        (2, 1): pd.DataFrame.from_records([_hit_row("c0", "b4")], columns=COMPARISON_COLUMNS),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
    )
    expected_edge_ids = {
        (str(row.query), str(row.subject), query_index, subject_index)
        for (query_index, subject_index), edges in edge_selection.adjacent_display_edges_by_pair.items()
        for row in edges.itertuples(index=False)
    }
    anchors = orthogroup_edges_to_lossless_collinearity_anchors(
        edge_selection.adjacent_display_edges_by_pair,
        extraction.protein_map,
        edge_selection.orthogroups,
    )

    assert {
        (
            anchor.query_protein_id,
            anchor.subject_protein_id,
            anchor.query_record_index,
            anchor.subject_record_index,
        )
        for anchor in anchors
    } == expected_edge_ids

    result = cluster_lossless_collinearity_anchors(
        anchors,
        params=LosslessCollinearityParameters(max_unit_gap=0, max_diagonal_drift=0),
    )
    rendered_frames = convert_collinearity_blocks_to_comparisons(
        result,
        records=records,
        color_mode="orientation",
    )
    rendered_edge_ids = set()
    for frame in rendered_frames:
        for row in frame.itertuples(index=False):
            query_ids = str(row.query_protein_id).split(";")
            subject_ids = str(row.subject_protein_id).split(";")
            rendered_edge_ids.update(zip(query_ids, subject_ids))

    assert {(query_id, subject_id) for query_id, subject_id, *_ in expected_edge_ids} <= rendered_edge_ids
    assert any(block.kind == "singleton" for block in result.blocks)
    assert any(block.kind == "cluster" and len(block.anchors) > 1 for block in result.blocks)


@pytest.mark.linear
def test_lossless_min_anchors_filters_blocks_after_clustering() -> None:
    anchors = [_anchor(0, 0), _anchor(1, 1), _anchor(4, 4)]

    min_one = cluster_lossless_collinearity_anchors(
        anchors,
        params=LosslessCollinearityParameters(min_anchors=1, max_unit_gap=0),
    )
    min_two = cluster_lossless_collinearity_anchors(
        anchors,
        params=LosslessCollinearityParameters(min_anchors=2, max_unit_gap=0),
    )
    min_three = cluster_lossless_collinearity_anchors(
        anchors,
        params=LosslessCollinearityParameters(min_anchors=3, max_unit_gap=0),
    )

    assert [len(block.anchors) for block in min_one.blocks] == [2, 1]
    assert [len(block.anchors) for block in min_two.blocks] == [2]
    assert min_two.unblocked_anchors == (anchors[2],)
    assert min_three.blocks == ()
    assert min_three.unblocked_anchors == tuple(anchors)


@pytest.mark.linear
def test_build_blocks_from_hits_maps_proteins_to_units() -> None:
    records = [
        _record(
            "record_a",
            [_cds(index * 12, index * 12 + 9, {"locus_tag": [f"qa{index}"], "protein_id": [f"qa{index}"]}) for index in range(3)],
        ),
        _record(
            "record_b",
            [_cds(index * 12, index * 12 + 9, {"locus_tag": [f"sb{index}"], "protein_id": [f"sb{index}"]}) for index in range(3)],
        ),
    ]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame.from_records(
        [_hit_row(f"qa{index}", f"sb{index}") for index in range(3)],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row(f"sb{index}", f"qa{index}") for index in range(3)],
        columns=COMPARISON_COLUMNS,
    )

    result = build_collinearity_blocks_from_hits(
        [hits],
        extraction,
        records=records,
        params=CollinearityParameters(min_anchors=3, max_gene_gap=1, block_evalue=None),
        anchor_mode="one_to_one",
        reverse_hits_by_pair=[reverse_hits],
    )
    comparisons = convert_collinearity_blocks_to_comparisons(result, records=records, color_mode="orientation")

    assert result.blocks[0].block_id == "block_0001"
    assert result.blocks[0].kind == "cluster"
    assert comparisons[0].shape[0] == 1
    assert comparisons[0].iloc[0]["collinearity_block_id"] == "block_0001"
    assert comparisons[0].iloc[0]["collinearity_block_kind"] == "cluster"
    assert comparisons[0].iloc[0]["collinearity_anchor_count"] == 3
    assert comparisons[0].iloc[0]["collinearity_color_mode"] == "orientation"
    assert comparisons[0].iloc[0]["orthogroup_id"] == "og_1;og_2;og_3"
    assert comparisons[0].iloc[0]["query_locus_id"] == "qa0;qa1;qa2"


@pytest.mark.linear
def test_true_global_singleton_orthogroup_does_not_create_link() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]}),
                _cds(12, 21, {"locus_tag": ["qa_single"], "protein_id": ["qa_single"]}),
            ],
        ),
        _record(
            "record_b",
            [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})],
        ),
    ]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame.from_records([_hit_row("qa0", "sb0")], columns=COMPARISON_COLUMNS)
    reverse_hits = pd.DataFrame.from_records([_hit_row("sb0", "qa0")], columns=COMPARISON_COLUMNS)

    result = build_collinearity_blocks_from_hits(
        [hits],
        extraction,
        records=records,
        params=LosslessCollinearityParameters(min_anchors=1),
        anchor_mode="one_to_one",
        reverse_hits_by_pair=[reverse_hits],
    )

    assert len(result.blocks) == 1
    assert result.blocks[0].kind == "singleton"
    assert result.blocks[0].anchors[0].orthogroup_id == "og_1"
    assert result.orthogroups is not None
    assert len(result.orthogroups.orthogroups) == 1
    assert all(
        "qa_single" not in [member.protein_id for member in members]
        for members in result.orthogroups.orthogroups.values()
    )


@pytest.mark.linear
def test_anchor_core_membership_adds_strong_non_rbh_member_without_new_anchor() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]})]),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", bitscore=300), _hit_row("a0", "b1", bitscore=250)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", bitscore=300), _hit_row("b1", "a0", bitscore=250)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        directional_hits,
        extraction,
        records=records,
        params=LosslessCollinearityParameters(min_anchors=1),
        edge_mode="rbh",
        orthogroup_membership_mode="anchor_core_v1",
        orthogroup_member_max_hits=2,
    )

    assert result.orthogroups is not None
    member = result.orthogroups.member_by_protein_id["b1"]
    assert member.orthogroup_id == "og_1"
    assert member.representative is False
    assert member.role == "coortholog"
    assert member.confidence == "high"
    assert result.blocks[0].anchors[0].subject_protein_id == "b0"
    assert {
        edge.edge_kind
        for edge in result.orthogroups.ortholog_edges_by_orthogroup_id["og_1"]
    } == {"rbh", "coortholog"}

@pytest.mark.linear
def test_anchor_core_unions_connected_near_reciprocal_groups() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a0", "b0", bitscore=300),
                _hit_row("a0", "b1", bitscore=280),
                _hit_row("a1", "b1", bitscore=310),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("b0", "a0", bitscore=300),
                _hit_row("b1", "a1", bitscore=310),
                _hit_row("b1", "a0", bitscore=280),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        directional_hits,
        extraction,
        records=records,
        params=LosslessCollinearityParameters(min_anchors=1),
        edge_mode="rbh",
        orthogroup_membership_mode="anchor_core_v1",
        orthogroup_member_max_hits=2,
    )

    assert result.orthogroups is not None
    assert set(result.orthogroups.orthogroups) == {"og_1"}
    assert {member.protein_id for member in result.orthogroups.orthogroups["og_1"]} == {
        "a0",
        "a1",
        "b0",
        "b1",
    }
    assert result.orthogroups.rbh_orthogroups == {
        "og_1": ("a0", "a1", "b0", "b1"),
    }
    path_sets = {
        tuple(path.protein_ids)
        for path in result.orthogroups.ortholog_paths_by_orthogroup_id["og_1"]
    }
    assert ("a0", "b0") in path_sets
    assert ("a1", "b1") in path_sets


@pytest.mark.linear
def test_orthogroup_display_edges_skip_cross_group_related_hits() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a0", "b0", bitscore=300),
                _hit_row("a1", "b1", bitscore=300),
                _hit_row("a0", "b1", bitscore=120),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("b0", "a0", bitscore=300),
                _hit_row("b1", "a1", bitscore=300),
                _hit_row("b1", "a0", bitscore=120),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
    )
    display_edges = edge_selection.adjacent_display_edges_by_pair[(0, 1)]
    links = convert_protein_hits_to_genomic_links(
        display_edges,
        extraction.protein_map,
        edge_selection.orthogroups,
    )

    assert {
        (str(row.query), str(row.subject))
        for row in display_edges.itertuples(index=False)
    } == {("a0", "b0"), ("a1", "b1")}
    assert set(links["orthogroup_id"]) == {"og_1", "og_2"}


@pytest.mark.linear
def test_anchor_core_same_record_hits_do_not_merge_distinct_cross_record_cores() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a0", "a1", bitscore=1000), _hit_row("a1", "a0", bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 1): pd.DataFrame.from_records(
            [_hit_row("b0", "b1", bitscore=1000), _hit_row("b1", "b0", bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", bitscore=300), _hit_row("a1", "b1", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", bitscore=300), _hit_row("b1", "a1", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    member_sets = {
        frozenset(member.protein_id for member in members)
        for members in result.orthogroups.values()
    }
    assert member_sets == {frozenset({"a0", "b0"}), frozenset({"a1", "b1"})}


@pytest.mark.linear
def test_anchor_core_same_record_hits_assign_inparalogs_to_existing_core() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a1", "a0", bitscore=250), _hit_row("a0", "a1", bitscore=250)],
            columns=COMPARISON_COLUMNS,
        ),
        (0, 1): pd.DataFrame.from_records([_hit_row("a0", "b0", bitscore=300)], columns=COMPARISON_COLUMNS),
        (1, 0): pd.DataFrame.from_records([_hit_row("b0", "a0", bitscore=300)], columns=COMPARISON_COLUMNS),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert set(result.orthogroups) == {"og_1"}
    assert {member.protein_id for member in result.orthogroups["og_1"]} == {"a0", "a1", "b0"}
    member = result.member_by_protein_id["a1"]
    assert member.role == "inparalog"
    assert member.confidence == "high"
    assert {
        edge.edge_kind
        for edge in result.ortholog_edges_by_orthogroup_id["og_1"]
    } == {"rbh", "same_record_inparalog"}


@pytest.mark.linear
def test_anchor_core_record_local_paralog_cluster_forms_species_specific_orthogroup() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a0", "a1", bitscore=300), _hit_row("a1", "a0", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert set(result.orthogroups) == {"og_1"}
    assert {member.protein_id for member in result.orthogroups["og_1"]} == {"a0", "a1"}
    assert result.scope_by_orthogroup_id["og_1"] == "record_local"
    assert result.source_record_index_by_orthogroup_id["og_1"] == 0
    assert {member.role for member in result.orthogroups["og_1"]} == {"local_paralog"}
    assert {member.confidence for member in result.orthogroups["og_1"]} == {"high"}
    assert {
        edge.edge_kind
        for edge in result.ortholog_edges_by_orthogroup_id["og_1"]
    } == {"record_local_paralog"}
    assert len({member.record_index for member in result.orthogroups["og_1"]}) == 1


@pytest.mark.linear
def test_anchor_core_record_local_singleton_remains_unassigned_with_singletons_requested() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]})]),
    ]
    extraction = extract_cds_proteins(records)

    result = select_rbh_orthogroup_edges_from_directional_hits(
        {},
        extraction.protein_map,
        record_count=len(records),
        include_singletons=True,
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert result.orthogroups == {}
    assert result.member_by_protein_id == {}


@pytest.mark.linear
def test_anchor_core_record_local_ignores_already_assigned_members() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a0", "a1", bitscore=280), _hit_row("a1", "a0", bitscore=280)],
            columns=COMPARISON_COLUMNS,
        ),
        (0, 1): pd.DataFrame.from_records([_hit_row("a0", "b0", bitscore=300)], columns=COMPARISON_COLUMNS),
        (1, 0): pd.DataFrame.from_records([_hit_row("b0", "a0", bitscore=300)], columns=COMPARISON_COLUMNS),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert set(result.orthogroups) == {"og_1"}
    assert {member.protein_id for member in result.orthogroups["og_1"]} == {"a0", "a1", "b0"}
    assert result.member_by_protein_id["a1"].role == "inparalog"
    assert set(result.scope_by_orthogroup_id.values()) == {"cross_record"}


@pytest.mark.linear
def test_anchor_core_record_local_rejects_domain_only_component() -> None:
    long_translation = "M" * 100
    records = [
        _record(
            "record_a",
            [
                _cds(0, 300, {"locus_tag": ["a0"], "protein_id": ["a0"], "translation": [long_translation]}),
                _cds(303, 603, {"locus_tag": ["a1"], "protein_id": ["a1"], "translation": [long_translation]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [
                _hit_row("a0", "a1", bitscore=300, alignment_length=20, qend=20, send=20),
                _hit_row("a1", "a0", bitscore=300, alignment_length=20, qend=20, send=20),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert result.orthogroups == {}
    assert "a0" not in result.member_by_protein_id
    assert "a1" not in result.member_by_protein_id


@pytest.mark.linear
def test_anchor_core_record_local_ids_are_deterministic_after_cross_record_groups() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
                _cds(24, 33, {"locus_tag": ["a2"], "protein_id": ["a2"]}),
            ],
        ),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a1", "a2", bitscore=300), _hit_row("a2", "a1", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
        (0, 1): pd.DataFrame.from_records([_hit_row("a0", "b0", bitscore=300)], columns=COMPARISON_COLUMNS),
        (1, 0): pd.DataFrame.from_records([_hit_row("b0", "a0", bitscore=300)], columns=COMPARISON_COLUMNS),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert {member.protein_id for member in result.orthogroups["og_1"]} == {"a0", "b0"}
    assert {member.protein_id for member in result.orthogroups["og_2"]} == {"a1", "a2"}
    assert result.scope_by_orthogroup_id == {"og_1": "cross_record", "og_2": "record_local"}


@pytest.mark.linear
def test_anchor_core_record_local_member_hit_limit_does_not_change_membership() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a0", "a1", bitscore=300), _hit_row("a1", "a0", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    max_one = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_member_max_hits=1,
    ).orthogroups
    max_five = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_member_max_hits=5,
    ).orthogroups

    assert {
        frozenset(member.protein_id for member in members)
        for members in max_one.orthogroups.values()
    } == {frozenset({"a0", "a1"})}
    assert {
        frozenset(member.protein_id for member in members)
        for members in max_five.orthogroups.values()
    } == {frozenset({"a0", "a1"})}


@pytest.mark.linear
def test_anchor_core_record_local_membership_edge_not_duplicated_as_related() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a0", "a1", bitscore=300), _hit_row("a1", "a0", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
    ).orthogroups

    membership_pairs = {
        frozenset((edge.query_protein_id, edge.subject_protein_id))
        for edge in result.ortholog_edges_by_orthogroup_id["og_1"]
    }
    related_pairs = {
        frozenset((edge.query_protein_id, edge.subject_protein_id))
        for edges in result.related_edges_by_orthogroup_id.values()
        for edge in edges
    }
    assert membership_pairs == {frozenset({"a0", "a1"})}
    assert not membership_pairs.intersection(related_pairs)


@pytest.mark.linear
def test_anchor_core_record_local_edges_do_not_create_collinearity_anchors() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 0): pd.DataFrame.from_records(
            [_hit_row("a0", "a1", bitscore=300), _hit_row("a1", "a0", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
    )
    anchors = orthogroup_edges_to_lossless_collinearity_anchors(
        edge_selection.adjacent_display_edges_by_pair,
        extraction.protein_map,
        edge_selection.orthogroups,
    )

    assert edge_selection.orthogroups.scope_by_orthogroup_id["og_1"] == "record_local"
    assert edge_selection.adjacent_display_edges_by_pair[(0, 1)].empty
    assert anchors == ()


@pytest.mark.linear
def test_anchor_core_member_hit_limit_does_not_change_inferred_membership() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]})]),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", bitscore=300), _hit_row("a0", "b1", bitscore=250)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", bitscore=300), _hit_row("b1", "a0", bitscore=250)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    max_one = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
        orthogroup_member_max_hits=1,
    ).orthogroups
    max_five = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
        orthogroup_member_max_hits=5,
    ).orthogroups

    one_sets = {
        frozenset(member.protein_id for member in members)
        for members in max_one.orthogroups.values()
    }
    five_sets = {
        frozenset(member.protein_id for member in members)
        for members in max_five.orthogroups.values()
    }
    assert one_sets == five_sets == {frozenset({"a0", "b0", "b1"})}


@pytest.mark.linear
def test_anchor_core_records_ambiguous_paralog_without_assigning_member() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["c0"], "protein_id": ["c0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", bitscore=300), _hit_row("a1", "b1", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", bitscore=300), _hit_row("b1", "a1", bitscore=300)],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 0): pd.DataFrame.from_records(
            [_hit_row("c0", "a0", bitscore=200), _hit_row("c0", "a1", bitscore=190)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert "c0" not in result.member_by_protein_id
    assert any(
        edge.edge_kind == "ambiguous_paralog"
        for edges in result.related_edges_by_orthogroup_id.values()
        for edge in edges
    )


@pytest.mark.linear
def test_shared_node_paths_are_serialized_without_collapsing_to_one_to_one() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]}),
                _cds(12, 21, {"locus_tag": ["a1"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]}),
                _cds(12, 21, {"locus_tag": ["b1"], "protein_id": ["b1"]}),
            ],
        ),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["c0"], "protein_id": ["c0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", bitscore=300), _hit_row("a1", "b1", bitscore=240)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", bitscore=300), _hit_row("b1", "a1", bitscore=240)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 2): pd.DataFrame.from_records(
            [_hit_row("b0", "c0", bitscore=300), _hit_row("b1", "c0", bitscore=260)],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 1): pd.DataFrame.from_records(
            [_hit_row("c0", "b0", bitscore=300), _hit_row("c0", "b1", bitscore=260)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        directional_hits,
        extraction,
        records=records,
        params=LosslessCollinearityParameters(min_anchors=1),
        edge_mode="rbh",
        orthogroup_membership_mode="anchor_core_v1",
        orthogroup_member_max_hits=2,
    )

    assert result.orthogroups is not None
    paths = result.orthogroups.ortholog_paths_by_orthogroup_id["og_1"]
    path_sets = {tuple(path.protein_ids): tuple(path.shared_protein_ids) for path in paths}
    assert path_sets[("a0", "b0", "c0")] == ("c0",)
    assert path_sets[("a1", "b1", "c0")] == ("c0",)


@pytest.mark.linear
def test_build_blocks_from_hits_adjacent_scope_ignores_non_adjacent_tables() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["c0"], "protein_id": ["c0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 2): pd.DataFrame.from_records([_hit_row("a0", "c0")], columns=COMPARISON_COLUMNS),
    }

    result = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        directional_hits,
        extraction,
        records=records,
        params=LosslessCollinearityParameters(),
        edge_mode="one_to_one",
        search_scope="adjacent",
    )

    assert result.blocks == ()
    assert result.orthogroups is not None
    assert result.orthogroups.orthogroups == {}


@pytest.mark.linear
def test_build_blocks_from_hits_all_scope_uses_non_adjacent_connectivity_only_as_evidence() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["c0"], "protein_id": ["c0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records([_hit_row("a0", "b0")], columns=COMPARISON_COLUMNS),
        (1, 0): pd.DataFrame.from_records([_hit_row("b0", "a0")], columns=COMPARISON_COLUMNS),
        (0, 2): pd.DataFrame.from_records([_hit_row("a0", "c0")], columns=COMPARISON_COLUMNS),
        (2, 0): pd.DataFrame.from_records([_hit_row("c0", "a0")], columns=COMPARISON_COLUMNS),
    }

    result = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        directional_hits,
        extraction,
        records=records,
        params=LosslessCollinearityParameters(),
        edge_mode="one_to_one",
        search_scope="all",
    )
    comparisons = convert_collinearity_blocks_to_comparisons(
        result,
        records=records,
        color_mode="orientation",
    )

    assert result.orthogroups is not None
    assert any(
        {member.protein_id for member in members} == {"a0", "b0", "c0"}
        for members in result.orthogroups.orthogroups.values()
    )
    assert all(
        block.subject_record_index == block.query_record_index + 1
        for block in result.blocks
    )
    assert comparisons[0].shape[0] == 1
    assert comparisons[1].shape[0] == 0


@pytest.mark.linear
def test_orthogroup_collinearity_rbh_search_defaults_to_top_candidate(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_limits: list[int | None] = []
    observed_pairs: list[tuple[str, str]] = []

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert losatp_threads is None
        assert max_hsps_per_subject is None
        observed_limits.append(candidate_limit)
        observed_pairs.append((_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta)))
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(collinearity_module, "_run_losatp_search", fake_search)

    result = collinearity_module.build_orthogroup_collinearity_blocks(
        records,
        params=LosslessCollinearityParameters(),
        edge_mode="rbh",
    )

    assert observed_limits == [None, None, None, None, None, None, None]
    assert observed_pairs == [
        ("qa0", "qa0"),
        ("sb0", "sb0"),
        ("tc0", "tc0"),
        ("qa0", "sb0"),
        ("sb0", "qa0"),
        ("sb0", "tc0"),
        ("tc0", "sb0"),
    ]
    assert len(result.blocks) == 2
    assert {block.query_record_index for block in result.blocks} == {0, 1}


@pytest.mark.linear
def test_orthogroup_collinearity_anchor_core_ignores_member_hit_depth(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_limits: list[int | None] = []

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert max_hsps_per_subject is None
        observed_limits.append(candidate_limit)
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(collinearity_module, "_run_losatp_search", fake_search)

    collinearity_module.build_orthogroup_collinearity_blocks(
        records,
        params=LosslessCollinearityParameters(),
        edge_mode="rbh",
        orthogroup_membership_mode="anchor_core_v1",
        orthogroup_member_max_hits=3,
    )

    assert observed_limits == [None, None, None, None, None, None, None]


@pytest.mark.linear
def test_orthogroup_collinearity_all_hits_keeps_uncapped_search(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_limits: list[int | None] = []
    observed_pairs: list[tuple[str, str]] = []

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert losatp_threads is None
        assert max_hsps_per_subject is None
        observed_limits.append(candidate_limit)
        observed_pairs.append((_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta)))
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(collinearity_module, "_run_losatp_search", fake_search)

    result = collinearity_module.build_orthogroup_collinearity_blocks(
        records,
        params=LosslessCollinearityParameters(),
        edge_mode="all",
        search_scope="all",
    )

    assert observed_limits == [None, None, None, None, None, None, None, None, None]
    assert observed_pairs == [
        ("qa0", "qa0"),
        ("sb0", "sb0"),
        ("tc0", "tc0"),
        ("qa0", "sb0"),
        ("sb0", "qa0"),
        ("qa0", "tc0"),
        ("tc0", "qa0"),
        ("sb0", "tc0"),
        ("tc0", "sb0"),
    ]
    assert all(
        block.subject_record_index == block.query_record_index + 1
        for block in result.blocks
    )


@pytest.mark.linear
def test_orthogroup_collinearity_all_scope_rbh_searches_every_direction(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_pairs: list[tuple[str, str]] = []

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert losatp_threads is None
        assert max_hsps_per_subject is None
        observed_pairs.append((_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta)))
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(collinearity_module, "_run_losatp_search", fake_search)

    collinearity_module.build_orthogroup_collinearity_blocks(
        records,
        params=LosslessCollinearityParameters(),
        edge_mode="rbh",
        search_scope="all",
    )

    assert observed_pairs == [
        ("qa0", "qa0"),
        ("sb0", "sb0"),
        ("tc0", "tc0"),
        ("qa0", "sb0"),
        ("sb0", "qa0"),
        ("qa0", "tc0"),
        ("tc0", "qa0"),
        ("sb0", "tc0"),
        ("tc0", "sb0"),
    ]


@pytest.mark.linear
def test_native_tsv_parser_resolves_records_and_units_and_writer_round_trips() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]}),
                _cds(12, 21, {"locus_tag": ["qa1"], "protein_id": ["qa1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]}),
                _cds(12, 21, {"locus_tag": ["sb1"], "protein_id": ["sb1"]}),
            ],
        ),
    ]
    text = "\n".join(
        [
            "block_id\tanchor_index\tquery_record\tquery_unit\tsubject_record\tsubject_unit\torientation\tidentity\tevalue\tbitscore\talignment_length\tscore\tblock_evalue",
            "block_a\t1\t#1\tqa0\t#2\tsb0\tplus\t91\t1e-20\t200\t100\t50\t1e-8",
            "block_a\t2\trecord_a\tqa1\trecord_b\tsb1\tplus\t92\t1e-25\t210\t100\t50\t0.00000001",
            "",
        ]
    )

    result = parse_native_collinearity_tsv(
        text,
        records,
        params=CollinearityParameters(min_anchors=2),
    )
    written = write_native_collinearity_tsv(result)
    reparsed = parse_native_collinearity_tsv(
        StringIO(written).getvalue(),
        records,
        params=CollinearityParameters(min_anchors=2),
    )

    assert result.blocks[0].block_id == "block_a"
    assert result.blocks[0].kind == "syntenic"
    assert result.blocks[0].block_evalue == pytest.approx(1e-8)
    assert [anchor.query_display_name for anchor in result.blocks[0].anchors] == ["qa0", "qa1"]
    assert "block_kind" in written.splitlines()[0].split("\t")
    assert "orthogroup_id" in written.splitlines()[0].split("\t")
    assert "block_evalue" in written.splitlines()[0].split("\t")
    assert reparsed.blocks[0].anchors[1].subject_locus_id == "sb1"
    assert reparsed.blocks[0].block_evalue == pytest.approx(1e-8)


@pytest.mark.linear
def test_native_tsv_parser_round_trips_singleton_block_kind() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
    ]
    text = "\n".join(
        [
            "block_id\tblock_kind\tanchor_index\torthogroup_id\tquery_record\tquery_unit\tsubject_record\tsubject_unit\torientation",
            "singleton_a\tsingleton\t1\tog_x\t#1\tqa0\t#2\tsb0\tplus",
            "",
        ]
    )

    result = parse_native_collinearity_tsv(
        text,
        records,
        params=LosslessCollinearityParameters(min_anchors=1),
    )
    written = write_native_collinearity_tsv(result)

    assert result.blocks[0].kind == "singleton"
    assert result.blocks[0].anchors[0].orthogroup_id == "og_x"
    assert "\tsingleton\t" in written


@pytest.mark.linear
def test_native_tsv_parser_min_anchors_drops_small_blocks() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
    ]
    text = "\n".join(
        [
            "block_id\tblock_kind\tanchor_index\torthogroup_id\tquery_record\tquery_unit\tsubject_record\tsubject_unit\torientation",
            "singleton_a\tsingleton\t1\tog_x\t#1\tqa0\t#2\tsb0\tplus",
            "",
        ]
    )

    result = parse_native_collinearity_tsv(
        text,
        records,
        params=LosslessCollinearityParameters(min_anchors=2),
    )

    assert result.blocks == ()
    assert len(result.unblocked_anchors) == 1
    assert result.unblocked_anchors[0].orthogroup_id == "og_x"


@pytest.mark.linear
def test_native_tsv_parser_rejects_conflicting_block_evalues() -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]}),
                _cds(12, 21, {"locus_tag": ["qa1"], "protein_id": ["qa1"]}),
            ],
        ),
        _record(
            "record_b",
            [
                _cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]}),
                _cds(12, 21, {"locus_tag": ["sb1"], "protein_id": ["sb1"]}),
            ],
        ),
    ]
    text = "\n".join(
        [
            "block_id\tanchor_index\tquery_record\tquery_unit\tsubject_record\tsubject_unit\torientation\tscore\tblock_evalue",
            "block_a\t1\t#1\tqa0\t#2\tsb0\tplus\t50\t1e-8",
            "block_a\t2\t#1\tqa1\t#2\tsb1\tplus\t50\t1e-7",
            "",
        ]
    )

    with pytest.raises(ValidationError, match="conflicting block_evalue"):
        parse_native_collinearity_tsv(
            text,
            records,
            params=CollinearityParameters(min_anchors=2),
        )


@pytest.mark.linear
def test_linear_cli_builds_and_saves_native_collinearity(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"]})]), _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"]})])]
    captured: dict[str, object] = {}
    save_path = tmp_path / "blocks.collinear.tsv"

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_build(*_args, **kwargs):
        captured["collinearity_params"] = kwargs["params"]
        captured["unit_mode"] = kwargs["unit_mode"]
        captured["edge_mode"] = kwargs["edge_mode"]
        captured["search_scope"] = kwargs["search_scope"]
        captured["orthogroup_membership_mode"] = kwargs["orthogroup_membership_mode"]
        captured["orthogroup_member_max_hits"] = kwargs["orthogroup_member_max_hits"]
        captured["max_paralog_links_per_orthogroup"] = kwargs["max_paralog_links_per_orthogroup"]
        anchor = _anchor(0, 0)
        return CollinearityResult(
            blocks=(
                CollinearityBlock(
                    block_id="block_0001",
                    query_record_index=0,
                    subject_record_index=1,
                    orientation="plus",
                    score=50.0,
                    anchors=(anchor,),
                ),
            )
        )

    def fake_assemble(*_args, **kwargs):
        captured["protein_comparisons"] = kwargs.get("protein_comparisons")
        captured["protein_blastp_mode"] = kwargs.get("protein_blastp_mode")
        return kwargs.get("canvas") or object()

    monkeypatch.setattr(linear_cli_module, "build_orthogroup_collinearity_blocks", fake_build)
    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "collinear",
            "--collinear_min_anchors",
            "1",
            "--collinear_unit_mode",
            "cds",
            "--collinear_anchor_mode",
            "all",
            "--collinear_search_scope",
            "all",
            "--orthogroup_membership_mode",
            "family_merge",
            "--orthogroup_member_max_hits",
            "3",
            "--collinear_block_merge_gap",
            "12",
            "--collinear_singleton_merge_gap",
            "7",
            "--collinear_color_mode",
            "orientation",
            "--save_collinear_blocks",
            str(save_path),
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    params = captured["collinearity_params"]
    assert isinstance(params, LosslessCollinearityParameters)
    assert params.min_anchors == 1
    assert params.max_unit_gap == 0
    assert params.max_diagonal_drift == 0
    assert captured["unit_mode"] == "cds"
    assert captured["edge_mode"] == "rbh"
    assert captured["search_scope"] == "all"
    assert captured["orthogroup_membership_mode"] == "anchor_core_v1"
    assert captured["orthogroup_member_max_hits"] == 3
    assert captured["max_paralog_links_per_orthogroup"] == 2
    assert captured["protein_blastp_mode"] == "none"
    comparisons = captured["protein_comparisons"]
    assert isinstance(comparisons, list)
    assert comparisons[0].iloc[0]["collinearity_color_mode"] == "orientation"
    assert "block_0001" in save_path.read_text(encoding="utf-8")


@pytest.mark.linear
def test_linear_cli_parses_collinear_block_evalue_none() -> None:
    args = linear_cli_module._get_args(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "collinear",
            "--collinear_block_evalue",
            "none",
        ]
    )

    assert args.collinear_block_evalue is None


@pytest.mark.linear
def test_linear_cli_accepts_orientation_identity_collinear_color_mode() -> None:
    args = linear_cli_module._get_args(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "collinear",
            "--collinear_color_mode",
            "orientation_identity",
        ]
    )

    assert args.collinear_color_mode == "orientation_identity"


@pytest.mark.linear
def test_linear_cli_parses_collinear_min_anchors() -> None:
    default_args = linear_cli_module._get_args(["--gbk", "a.gb", "b.gb"])
    explicit_args = linear_cli_module._get_args(
        ["--gbk", "a.gb", "b.gb", "--collinear_min_anchors", "2"]
    )

    assert default_args.collinear_min_anchors == 1
    assert explicit_args.collinear_min_anchors == 2


@pytest.mark.linear
def test_linear_cli_rejects_nonpositive_collinear_min_anchors() -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(
            ["--gbk", "a.gb", "b.gb", "--collinear_min_anchors", "0"]
        )


@pytest.mark.linear
def test_linear_cli_rejects_negative_collinear_block_evalue() -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(
            [
                "--gbk",
                "a.gb",
                "b.gb",
                "--protein_blastp_mode",
                "collinear",
                "--collinear_block_evalue",
                "-1",
            ]
        )


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_returns_collinear_rows() -> None:
    class JsNull:
        def __str__(self) -> str:
            return "null"

    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    hits = pd.DataFrame.from_records(
        [_hit_row(f"qa{index}", f"sb{index}") for index in range(8)],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row(f"sb{index}", f"qa{index}") for index in range(8)],
        columns=COMPARISON_COLUMNS,
    )
    query_map = {
        f"qa{index}": _web_protein_entry(
            f"qa{index}",
            record_index=0,
            record_id="record_a",
            feature_index=index,
            start=index * 12,
            end=index * 12 + 9,
        )
        for index in range(8)
    }
    subject_map = {
        f"sb{index}": _web_protein_entry(
            f"sb{index}",
            record_index=1,
            record_id="record_b",
            feature_index=index,
            start=index * 12,
            end=index * 12 + 9,
        )
        for index in range(8)
    }
    payload = {
        "records": [
            {
                "recordIndex": 0,
                "recordId": "record_a",
                "proteinMap": query_map,
                "proteinCacheKey": "record-a-cache",
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": subject_map,
                "proteinCacheKey": "record-b-cache",
            },
        ],
        "pairs": [
            {
                "pairIndex": 0,
                "queryIndex": 0,
                "subjectIndex": 1,
                "cacheKey": "pair-a-b",
                "blastText": hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
            {
                "pairIndex": 0,
                "queryIndex": 1,
                "subjectIndex": 0,
                "cacheKey": "pair-b-a",
                "blastText": reverse_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            }
        ],
    }

    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        json.dumps(payload),
        "collinear",
        5,
        50,
        "1e-5",
        0,
        0,
        3,
        1,
        "cds",
        "orientation",
        "one_to_one",
        50,
        25,
        25,
        1,
        2,
        "adjacent",
        "rbh",
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    assert result["collinearityBlocks"][0]["id"] == "block_0001"
    rows = result["pairs"][0]["rows"]
    assert len(rows) == 1
    assert rows[0]["collinearity_block_id"] == "block_0001"
    assert rows[0]["collinearity_block_kind"] == "cluster"
    assert rows[0]["collinearity_anchor_count"] == 8
    assert rows[0]["collinearity_color_mode"] == "orientation"
    assert rows[0]["collinearity_block_evalue"] == ""
    assert rows[0]["group_kind"] == "collinear_gene_group"
    assert rows[0]["group_scope"] == "adjacent_local"
    assert rows[0]["collinear_group_scope"] == "adjacent_local"
    assert result["orthogroups"] == []
    assert result["collinearGroupScope"] == "adjacent_local"
    assert len(result["collinearGroups"]) == 8
    assert {group["scope"] for group in result["collinearGroups"]} == {"adjacent_local"}


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_uses_rbh_collinear_anchor_mode() -> None:
    class JsNull:
        def __str__(self) -> str:
            return "null"

    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    query_map = {
        f"qa{index}": _web_protein_entry(
            f"qa{index}",
            record_index=0,
            record_id="record_a",
            feature_index=index,
            start=index * 12,
            end=index * 12 + 9,
        )
        for index in range(2)
    }
    subject_map = {
        f"sb{index}": _web_protein_entry(
            f"sb{index}",
            record_index=1,
            record_id="record_b",
            feature_index=index,
            start=index * 12,
            end=index * 12 + 9,
        )
        for index in range(2)
    }
    forward_hits = pd.DataFrame.from_records(
        [_hit_row("qa0", "sb0", bitscore=300), _hit_row("qa1", "sb1", bitscore=250)],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row("sb0", "qa0", bitscore=300), _hit_row("sb1", "qa0", bitscore=400)],
        columns=COMPARISON_COLUMNS,
    )
    payload = {
        "records": [
            {
                "recordIndex": 0,
                "recordId": "record_a",
                "proteinMap": query_map,
                "proteinCacheKey": "record-a-cache",
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": subject_map,
                "proteinCacheKey": "record-b-cache",
            },
        ],
        "pairs": [
            {
                "pairIndex": 0,
                "queryIndex": 0,
                "subjectIndex": 1,
                "cacheKey": "pair-a-b",
                "blastText": forward_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
            {
                "pairIndex": 0,
                "queryIndex": 1,
                "subjectIndex": 0,
                "cacheKey": "pair-b-a",
                "blastText": reverse_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
        ],
    }

    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        json.dumps(payload),
        "collinear",
        5,
        50,
        "1e-5",
        0,
        0,
        1,
        0,
        "cds",
        "orientation",
        "rbh",
        50,
        25,
        25,
        1,
        2,
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    rows = result["pairs"][0]["rows"]
    assert len(rows) == 1
    assert rows[0]["collinearity_anchor_count"] == 1
    assert rows[0]["collinearity_block_kind"] == "singleton"
    assert rows[0]["query_protein_id"] == "qa0"
    assert rows[0]["subject_protein_id"] == "sb0"

    raw_min_two = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        json.dumps(payload),
        "collinear",
        5,
        50,
        "1e-5",
        0,
        0,
        2,
        0,
        "cds",
        "orientation",
        "rbh",
        50,
        25,
        25,
        1,
        2,
    )
    min_two = json.loads(str(raw_min_two))

    assert "error" not in min_two
    assert min_two["pairs"][0]["rows"] == []


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_applies_collinear_search_scope() -> None:
    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    a_map = {
        "a0": _web_protein_entry(
            "a0",
            record_index=0,
            record_id="record_a",
            feature_index=0,
            start=0,
            end=9,
        )
    }
    b_map = {
        "b0": _web_protein_entry(
            "b0",
            record_index=1,
            record_id="record_b",
            feature_index=0,
            start=0,
            end=9,
        )
    }
    c_map = {
        "c0": _web_protein_entry(
            "c0",
            record_index=2,
            record_id="record_c",
            feature_index=0,
            start=0,
            end=9,
        )
    }
    adjacent_hits = pd.DataFrame.from_records([_hit_row("a0", "b0")], columns=COMPARISON_COLUMNS)
    adjacent_reverse_hits = pd.DataFrame.from_records([_hit_row("b0", "a0")], columns=COMPARISON_COLUMNS)
    non_adjacent_hits = pd.DataFrame.from_records([_hit_row("a0", "c0")], columns=COMPARISON_COLUMNS)
    non_adjacent_reverse_hits = pd.DataFrame.from_records([_hit_row("c0", "a0")], columns=COMPARISON_COLUMNS)
    payload = {
        "records": [
            {
                "recordIndex": 0,
                "recordId": "record_a",
                "proteinMap": a_map,
                "proteinCacheKey": "record-a-cache",
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": b_map,
                "proteinCacheKey": "record-b-cache",
            },
            {
                "recordIndex": 2,
                "recordId": "record_c",
                "proteinMap": c_map,
                "proteinCacheKey": "record-c-cache",
            },
        ],
        "pairs": [
            {
                "pairIndex": 0,
                "queryIndex": 0,
                "subjectIndex": 1,
                "cacheKey": "pair-a-b",
                "blastText": adjacent_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
            {
                "pairIndex": 0,
                "queryIndex": 1,
                "subjectIndex": 0,
                "cacheKey": "pair-b-a",
                "blastText": adjacent_reverse_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
            {
                "pairIndex": 0,
                "queryIndex": 0,
                "subjectIndex": 2,
                "cacheKey": "pair-a-c",
                "blastText": non_adjacent_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
            {
                "pairIndex": 0,
                "queryIndex": 2,
                "subjectIndex": 0,
                "cacheKey": "pair-c-a",
                "blastText": non_adjacent_reverse_hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            },
        ],
    }

    def convert(scope: str) -> dict[str, object]:
        raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
            json.dumps(payload),
            "collinear",
            5,
            50,
            "1e-5",
            0,
            0,
            1,
            0,
            "cds",
            "orientation",
            "one_to_one",
            50,
            25,
            25,
            1,
            2,
            scope,
            "rbh",
        )
        return json.loads(str(raw_result))

    adjacent = convert("adjacent")
    all_records = convert("all")

    assert "error" not in adjacent
    assert "error" not in all_records
    adjacent_member_sets = [
        {member["proteinId"] for member in group["members"]}
        for group in adjacent["collinearGroups"]
    ]
    all_member_sets = [
        {member["proteinId"] for member in group["members"]}
        for group in all_records["orthogroups"]
    ]
    assert adjacent["orthogroups"] == []
    assert adjacent["collinearGroupScope"] == "adjacent_local"
    assert all_records["collinearGroups"] == []
    assert all_records["collinearGroupScope"] == "global_collinear"
    assert {group["scope"] for group in all_records["orthogroups"]} == {"global_collinear"}
    assert {"a0", "b0"} in adjacent_member_sets
    assert {"a0", "b0", "c0"} in all_member_sets
    assert all(
        block["subjectRecordIndex"] == block["queryRecordIndex"] + 1
        for block in all_records["collinearityBlocks"]
    )


@pytest.mark.linear
def test_collinearity_comparison_rows_use_block_spans() -> None:
    result = CollinearityResult(
        blocks=(
            CollinearityBlock(
                block_id="block_plus",
                query_record_index=0,
                subject_record_index=1,
                orientation="plus",
                score=150.0,
                block_evalue=1e-9,
                anchors=(_anchor(0, 0), _anchor(1, 1), _anchor(2, 2)),
            ),
            CollinearityBlock(
                block_id="block_minus",
                query_record_index=0,
                subject_record_index=1,
                orientation="minus",
                score=150.0,
                block_evalue=2e-9,
                anchors=(_anchor(0, 2), _anchor(1, 1), _anchor(2, 0)),
            ),
        )
    )

    comparisons = convert_collinearity_blocks_to_comparisons(
        result,
        record_ids=["record_a", "record_b"],
        color_mode="orientation",
    )

    rows = comparisons[0].set_index("collinearity_block_id")
    assert rows.loc["block_plus", "qstart"] == 1
    assert rows.loc["block_plus", "qend"] == 29
    assert rows.loc["block_plus", "sstart"] == 1
    assert rows.loc["block_plus", "send"] == 29
    assert rows.loc["block_minus", "qstart"] == 1
    assert rows.loc["block_minus", "qend"] == 29
    assert rows.loc["block_minus", "sstart"] == 29
    assert rows.loc["block_minus", "send"] == 1
    assert rows.loc["block_plus", "collinearity_block_evalue"] == pytest.approx(1e-9)


@pytest.mark.linear
def test_collinearity_match_path_with_metadata_serializes() -> None:
    group = PairWiseMatchGroup.__new__(PairWiseMatchGroup)
    group.canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=1000,
    )
    group.records = [_record("record_a", []), _record("record_b", [])]
    group.comparison_count = 1
    group.comparison_height = 40
    group.query_offset_x = 0
    group.subject_offset_x = 0
    group.query_alignment_offset_x = 0
    group.subject_alignment_offset_x = 0
    group.min_identity = 0
    group.match_min_color = "#ffffff"
    group.match_max_color = "#000000"
    group.match_fill_opacity = 0.75
    group.match_stroke_color = "none"
    group.match_stroke_width = 0
    group.collinearity_orientation_colors = {"plus": "#112233", "minus": "#445566"}

    row = SimpleNamespace(
        identity=95,
        qstart=1,
        qend=100,
        sstart=10,
        send=110,
        collinearity_block_id="block_0001",
        collinearity_block_kind="singleton",
        collinearity_orientation="plus",
        collinearity_block_score=42.5,
        collinearity_block_evalue=1e-9,
        collinearity_anchor_index=1,
        collinearity_anchor_count=1,
        collinearity_color_mode="orientation",
        orthogroup_id="og_1",
        query_protein_id="qa0",
        subject_protein_id="sb0",
        query_feature_svg_id="feature_qa0",
        subject_feature_svg_id="feature_sb0",
        query_unit_id="qa0",
        subject_unit_id="sb0",
        query_locus_id="locus_a",
        subject_locus_id="locus_b",
        query_display_name="geneA",
        subject_display_name="geneB",
    )

    drawing = Drawing(debug=False)
    drawing.add(group.generate_linear_match_path(row))

    svg_text = drawing.tostring()
    assert 'data-gbdraw-pairwise-match-id="comparison1_match1"' in svg_text
    assert 'data-match-kind="collinear"' in svg_text
    assert 'data-query-record-index="0"' in svg_text
    assert 'data-subject-record-index="1"' in svg_text
    assert 'data-query-record-id="record_a"' in svg_text
    assert 'data-subject-record-id="record_b"' in svg_text
    assert 'data-qstart="1"' in svg_text
    assert 'data-qend="100"' in svg_text
    assert 'data-sstart="10"' in svg_text
    assert 'data-send="110"' in svg_text
    assert 'data-identity="95"' in svg_text
    assert 'data-alignment-length=" "' in svg_text
    assert 'data-evalue=" "' in svg_text
    assert 'data-bitscore=" "' in svg_text
    assert 'data-mismatches=" "' in svg_text
    assert 'data-gap-opens=" "' in svg_text
    assert "data-collinearity-block-id" in drawing.tostring()
    assert 'data-collinearity-block-kind="singleton"' in svg_text
    assert 'data-orthogroup-id="og_1"' in svg_text
    assert 'data-collinearity-block-score="42.5"' in svg_text
    assert 'data-collinearity-block-evalue="1e-09"' in svg_text
    assert 'data-collinearity-anchor-index="1"' in svg_text
    assert 'data-collinearity-anchor-count="1"' in svg_text
    assert 'data-collinearity-color-mode="orientation"' in svg_text
    assert 'data-query-locus-id="locus_a"' in svg_text
    assert 'data-subject-locus-id="locus_b"' in svg_text
    assert 'data-query-display-name="geneA"' in svg_text
    assert 'data-subject-display-name="geneB"' in svg_text
    assert 'fill="#112233"' in svg_text


@pytest.mark.linear
def test_pairwise_match_path_emits_plain_required_metadata() -> None:
    group = _build_collinearity_match_group()
    row = SimpleNamespace(
        query="query_record",
        subject="subject_record",
        identity=87.5,
        alignment_length=120,
        mismatches=3,
        gap_opens=1,
        qstart=5,
        qend=125,
        sstart=40,
        send=160,
        evalue=2e-30,
        bitscore=240.5,
    )

    drawing = Drawing(debug=False)
    drawing.add(group.generate_linear_match_path(row, match_index=7))
    svg_text = drawing.tostring()

    assert 'data-gbdraw-pairwise-match-id="comparison1_match7"' in svg_text
    assert 'data-match-kind="pairwise"' in svg_text
    assert 'data-query-record-id="query_record"' in svg_text
    assert 'data-subject-record-id="subject_record"' in svg_text
    assert 'data-identity="87.5"' in svg_text
    assert 'data-alignment-length="120"' in svg_text
    assert 'data-evalue="2e-30"' in svg_text
    assert 'data-bitscore="240.5"' in svg_text
    assert 'data-mismatches="3"' in svg_text
    assert 'data-gap-opens="1"' in svg_text
    assert "data-orthogroup-id" not in svg_text
    assert "data-collinearity-block-id" not in svg_text


@pytest.mark.linear
def test_pairwise_match_path_marks_orthogroup_kind() -> None:
    group = _build_collinearity_match_group()
    row = SimpleNamespace(
        identity=91,
        alignment_length=95,
        mismatches=0,
        gap_opens=0,
        qstart=10,
        qend=105,
        sstart=20,
        send=115,
        evalue=1e-20,
        bitscore=180,
        orthogroup_id="og_2",
        query_protein_id="qa2",
        subject_protein_id="sb2",
        query_feature_svg_id="feature_qa2",
        subject_feature_svg_id="feature_sb2",
    )

    drawing = Drawing(debug=False)
    drawing.add(group.generate_linear_match_path(row, match_index=2))
    svg_text = drawing.tostring()

    assert 'data-match-kind="orthogroup"' in svg_text
    assert 'data-orthogroup-id="og_2"' in svg_text
    assert 'data-query-protein-id="qa2"' in svg_text
    assert 'data-subject-protein-id="sb2"' in svg_text
    assert 'data-query-feature-svg-id="feature_qa2"' in svg_text
    assert 'data-subject-feature-svg-id="feature_sb2"' in svg_text


@pytest.mark.linear
def test_blast_configurator_reads_collinearity_colors_from_default_colors() -> None:
    default_colors = pd.DataFrame(
        [
            ("pairwise_match_min", "#ffffff"),
            ("pairwise_match_max", "#000000"),
            ("pairwise_match", "#cccccc"),
            ("collinear_block_plus_min", "#ddeeff"),
            ("collinear_block_minus", "#445566"),
            ("collinear_block_minus_min", "#ffdddd"),
        ],
        columns=["feature_type", "color"],
    )

    config = BlastMatchConfigurator(
        evalue=1e-5,
        bitscore=50,
        identity=0,
        alignment_length=0,
        sequence_length_dict={},
        config_dict=load_config_toml("gbdraw.data", "config.toml"),
        default_colors_df=default_colors,
    )

    assert config.collinearity_orientation_colors == {"plus": "#8b9cc1", "minus": "#445566"}
    assert config.collinearity_orientation_min_colors == {"plus": "#ddeeff", "minus": "#ffdddd"}


@pytest.mark.linear
def test_blast_configurator_accepts_plus_max_collinearity_color_alias() -> None:
    default_colors = pd.DataFrame(
        [
            ("pairwise_match_min", "#ffffff"),
            ("pairwise_match_max", "#000000"),
            ("pairwise_match", "#cccccc"),
            ("collinear_block_plus_min", "#ddeeff"),
            ("collinear_block_plus_max", "#123456"),
        ],
        columns=["feature_type", "color"],
    )

    config = BlastMatchConfigurator(
        evalue=1e-5,
        bitscore=50,
        identity=0,
        alignment_length=0,
        sequence_length_dict={},
        config_dict=load_config_toml("gbdraw.data", "config.toml"),
        default_colors_df=default_colors,
    )

    assert config.collinearity_orientation_colors["plus"] == "#123456"


@pytest.mark.linear
def test_orientation_collinearity_suppresses_pairwise_identity_legend() -> None:
    feature_config = SimpleNamespace(
        color_table=None,
        default_colors=pd.DataFrame([("CDS", "#54bcf8"), ("default", "#d3d3d3")], columns=["feature_type", "color"]),
        block_stroke_color="none",
        block_stroke_width=0,
    )
    gc_config = SimpleNamespace(
        show_gc=False,
        dinucleotide="GC",
        stroke_color="none",
        stroke_width=0,
        high_fill_color="#ffffff",
        low_fill_color="#000000",
    )
    skew_config = SimpleNamespace(
        show_skew=False,
        stroke_color="none",
        stroke_width=0,
        high_fill_color="#ffffff",
        low_fill_color="#000000",
    )
    blast_config = SimpleNamespace(
        min_color="#ffffff",
        max_color="#000000",
        identity=0,
        hide_pairwise_identity_legend=True,
    )

    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        [],
        blast_config=blast_config,
        has_blast=True,
    )

    assert "Pairwise match identity" not in legend_table


@pytest.mark.linear
def test_average_identity_collinearity_legend_uses_average_identity_label() -> None:
    feature_config = SimpleNamespace(
        color_table=None,
        default_colors=pd.DataFrame([("CDS", "#54bcf8"), ("default", "#d3d3d3")], columns=["feature_type", "color"]),
        block_stroke_color="none",
        block_stroke_width=0,
    )
    gc_config = SimpleNamespace(
        show_gc=False,
        dinucleotide="GC",
        stroke_color="none",
        stroke_width=0,
        high_fill_color="#ffffff",
        low_fill_color="#000000",
    )
    skew_config = SimpleNamespace(
        show_skew=False,
        stroke_color="none",
        stroke_width=0,
        high_fill_color="#ffffff",
        low_fill_color="#000000",
    )
    blast_config = SimpleNamespace(
        min_color="#ffffff",
        max_color="#000000",
        identity=0,
        hide_pairwise_identity_legend=False,
        pairwise_identity_legend_label="Average identity",
    )

    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        [],
        blast_config=blast_config,
        has_blast=True,
    )

    assert "Average identity" in legend_table
    assert "Pairwise match identity" not in legend_table


@pytest.mark.linear
def test_orientation_identity_collinearity_legend_uses_two_orientation_gradients() -> None:
    feature_config = SimpleNamespace(
        color_table=None,
        default_colors=pd.DataFrame([("CDS", "#54bcf8"), ("default", "#d3d3d3")], columns=["feature_type", "color"]),
        block_stroke_color="none",
        block_stroke_width=0,
    )
    gc_config = SimpleNamespace(
        show_gc=False,
        dinucleotide="GC",
        stroke_color="none",
        stroke_width=0,
        high_fill_color="#ffffff",
        low_fill_color="#000000",
    )
    skew_config = SimpleNamespace(
        show_skew=False,
        stroke_color="none",
        stroke_width=0,
        high_fill_color="#ffffff",
        low_fill_color="#000000",
    )
    blast_config = SimpleNamespace(
        min_color="#ffffff",
        max_color="#000000",
        identity=0,
        hide_pairwise_identity_legend=False,
        pairwise_identity_legend_entries=[
            {
                "label": "Collinear",
                "min_color": "#eeeeee",
                "max_color": "#112233",
            },
            {
                "label": "Inverted",
                "min_color": "#ffeeee",
                "max_color": "#445566",
            },
        ],
    )

    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        [],
        blast_config=blast_config,
        has_blast=True,
    )

    assert "Collinear" in legend_table
    assert "Inverted" in legend_table
    assert "Pairwise match identity" not in legend_table
    assert legend_table["Collinear"]["min_color"] == "#eeeeee"
    assert legend_table["Inverted"]["max_color"] == "#445566"


@pytest.mark.linear
def test_orientation_identity_pairwise_legend_renders_collinear_above_inverted() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = SimpleNamespace(legend_position="bottom", total_width=800)
    legend_config = SimpleNamespace(
        font_family="'Liberation Sans', 'Arial', sans-serif",
        font_weight="normal",
        font_size=16.0,
        color_rect_size=20.0,
        num_of_columns=1,
        has_gradient=True,
        total_feature_legend_width=0.0,
        pairwise_legend_width=160.0,
    )
    legend_table = {
        "Collinear": {
            "type": "gradient",
            "min_color": "#e7ffff",
            "max_color": "#57e1df",
            "stroke": "none",
            "width": 0,
            "min_value": 0,
        },
        "Inverted": {
            "type": "gradient",
            "min_color": "#ffeeee",
            "max_color": "#e15759",
            "stroke": "none",
            "width": 0,
            "min_value": 0,
        },
    }

    drawing = Drawing(debug=False)
    legend_group = LegendGroup(config_dict, canvas_config, legend_config, legend_table)
    drawing.add(legend_group.get_group())
    svg_text = drawing.tostring()

    assert svg_text.index('data-legend-key="Collinear"') < svg_text.index(
        'data-legend-key="Inverted"'
    )
    assert svg_text.count("<linearGradient") == 4

    root = ET.fromstring(svg_text)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    collinear_entry = root.find(".//svg:g[@data-legend-key='Collinear']", namespace)
    assert collinear_entry is not None
    collinear_bar = next(
        path
        for path in collinear_entry.findall("svg:path", namespace)
        if str(path.get("fill", "")).startswith("url(")
    )
    bar_x, _ = _translate_xy(collinear_bar.get("transform"))
    label_width, _ = calculate_bbox_dimensions(
        "Collinear", legend_config.font_family, legend_config.font_size, legend_group.dpi
    )
    assert bar_x == pytest.approx(label_width + 0.2 * legend_config.color_rect_size)


@pytest.mark.linear
def test_vertical_linear_pairwise_legend_aligns_to_feature_color_left_edge() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = SimpleNamespace(legend_position="right", total_width=800)
    legend_config = SimpleNamespace(
        font_family="'Liberation Sans', 'Arial', sans-serif",
        font_weight="normal",
        font_size=10.0,
        color_rect_size=12.0,
        num_of_columns=1,
        has_gradient=True,
        total_feature_legend_width=0.0,
        pairwise_legend_width=120.0,
    )
    feature_key = "tyrosine recombinase"
    gradient_key = "Pairwise match identity"
    legend_table = {
        feature_key: {
            "type": "solid",
            "fill": "#54bcf8",
            "stroke": "none",
            "width": 0,
        },
        gradient_key: {
            "type": "gradient",
            "min_color": "#ffffff",
            "max_color": "#000000",
            "stroke": "none",
            "width": 0,
            "min_value": 20,
        },
    }

    drawing = Drawing(debug=False)
    legend_group = LegendGroup(config_dict, canvas_config, legend_config, legend_table)
    drawing.add(legend_group.get_group())
    root = ET.fromstring(drawing.tostring())
    namespace = {"svg": "http://www.w3.org/2000/svg"}

    vertical_legend = root.find(".//svg:g[@id='legend_vertical']", namespace)
    assert vertical_legend is not None
    feature_legend = vertical_legend.find("svg:g[@id='feature_legend_v']", namespace)
    pairwise_legend = vertical_legend.find("svg:g[@id='pairwise_legend']", namespace)
    assert feature_legend is not None
    assert pairwise_legend is not None

    feature_group_x, _ = _translate_xy_or_zero(feature_legend.get("transform"))
    feature_entry = feature_legend.find(f"svg:g[@data-legend-key='{feature_key}']", namespace)
    assert feature_entry is not None
    feature_rect = feature_entry.find("svg:path", namespace)
    feature_text = feature_entry.find("svg:text", namespace)
    assert feature_rect is not None
    assert feature_text is not None
    feature_rect_x, _ = _translate_xy(feature_rect.get("transform"))

    pairwise_group_x, _ = _translate_xy_or_zero(pairwise_legend.get("transform"))
    pairwise_entry = pairwise_legend.find(
        f"svg:g[@data-legend-key='{gradient_key}']", namespace
    )
    assert pairwise_entry is not None
    pairwise_bar = next(
        path
        for path in pairwise_entry.findall("svg:path", namespace)
        if str(path.get("fill", "")).startswith("url(")
    )
    pairwise_bar_x, _ = _translate_xy(pairwise_bar.get("transform"))
    title = pairwise_entry.find("svg:text", namespace)
    assert title is not None
    title_x, _ = _translate_xy(title.get("transform"))

    text_width, _ = calculate_bbox_dimensions(
        feature_key, legend_config.font_family, legend_config.font_size, legend_group.dpi
    )
    text_x_offset = (22 / 14) * legend_config.color_rect_size
    feature_width = text_x_offset + text_width
    pairwise_alignment_width = legend_config.pairwise_legend_width

    assert feature_width > pairwise_alignment_width
    feature_left = feature_group_x + feature_rect_x
    assert pairwise_group_x == pytest.approx(feature_left)
    assert pairwise_group_x + pairwise_bar_x == pytest.approx(feature_left)
    assert pairwise_bar_x == pytest.approx(0)
    assert title_x == pytest.approx(pairwise_bar_x + pairwise_alignment_width / 2)
    assert feature_group_x == pytest.approx(0)
    assert title.get("text-anchor") == "middle"


@pytest.mark.linear
def test_vertical_linear_pairwise_legend_centers_when_wider_than_features() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = SimpleNamespace(legend_position="right", total_width=800)
    legend_config = SimpleNamespace(
        font_family="'Liberation Sans', 'Arial', sans-serif",
        font_weight="normal",
        font_size=10.0,
        color_rect_size=12.0,
        num_of_columns=1,
        has_gradient=True,
        total_feature_legend_width=0.0,
        pairwise_legend_width=220.0,
    )
    feature_key = "CDS"
    gradient_key = "Pairwise match identity"
    legend_table = {
        feature_key: {
            "type": "solid",
            "fill": "#54bcf8",
            "stroke": "none",
            "width": 0,
        },
        gradient_key: {
            "type": "gradient",
            "min_color": "#ffffff",
            "max_color": "#000000",
            "stroke": "none",
            "width": 0,
            "min_value": 20,
        },
    }

    drawing = Drawing(debug=False)
    legend_group = LegendGroup(config_dict, canvas_config, legend_config, legend_table)
    drawing.add(legend_group.get_group())
    root = ET.fromstring(drawing.tostring())
    namespace = {"svg": "http://www.w3.org/2000/svg"}

    vertical_legend = root.find(".//svg:g[@id='legend_vertical']", namespace)
    assert vertical_legend is not None
    feature_legend = vertical_legend.find("svg:g[@id='feature_legend_v']", namespace)
    pairwise_legend = vertical_legend.find("svg:g[@id='pairwise_legend']", namespace)
    assert feature_legend is not None
    assert pairwise_legend is not None
    pairwise_entry = pairwise_legend.find(
        f"svg:g[@data-legend-key='{gradient_key}']", namespace
    )
    assert pairwise_entry is not None
    title = pairwise_entry.find("svg:text", namespace)
    assert title is not None

    feature_group_x, _ = _translate_xy_or_zero(feature_legend.get("transform"))
    pairwise_group_x, _ = _translate_xy_or_zero(pairwise_legend.get("transform"))

    text_width, _ = calculate_bbox_dimensions(
        feature_key, legend_config.font_family, legend_config.font_size, legend_group.dpi
    )
    text_x_offset = (22 / 14) * legend_config.color_rect_size
    feature_width = text_x_offset + text_width
    pairwise_alignment_width = legend_config.pairwise_legend_width
    expected_feature_offset = (pairwise_alignment_width - feature_width) / 2
    title_x, _ = _translate_xy(title.get("transform"))

    assert pairwise_group_x == pytest.approx(0)
    assert feature_group_x == pytest.approx(expected_feature_offset)
    assert title_x == pytest.approx(pairwise_alignment_width / 2)
    assert title.get("text-anchor") == "middle"


@pytest.mark.linear
def test_blast_file_collinearity_metadata_drives_orientation_identity_legend(monkeypatch) -> None:
    frame = pd.DataFrame(
        [
            {
                **_hit_row("record_a", "record_b"),
                "collinearity_block_id": "block_1",
                "collinearity_orientation": "plus",
                "collinearity_color_mode": "orientation_identity",
            }
        ]
    )

    def fake_load_comparisons(_blast_files, _blast_config):
        return [frame]

    monkeypatch.setattr(linear_assemble_module, "load_comparisons", fake_load_comparisons)

    svg_text = assemble_linear_diagram_from_records(
        [_record("record_a", []), _record("record_b", [])],
        blast_files=["/virtual/blast_0.txt"],
        legend="bottom",
        config_overrides={"show_gc": False, "show_skew": False},
        bitscore=0,
        identity=0,
        evalue=1,
        alignment_length=0,
    ).tostring()

    assert 'data-legend-key="Collinear"' in svg_text
    assert 'data-legend-key="Inverted"' in svg_text
    assert 'data-legend-key="Pairwise match identity"' not in svg_text


@pytest.mark.linear
def test_collinearity_orientation_colors_use_default_color_overrides() -> None:
    group = PairWiseMatchGroup.__new__(PairWiseMatchGroup)
    group.canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=1000,
    )
    group.records = [_record("record_a", []), _record("record_b", [])]
    group.comparison_count = 1
    group.comparison_height = 40
    group.query_offset_x = 0
    group.subject_offset_x = 0
    group.query_alignment_offset_x = 0
    group.subject_alignment_offset_x = 0
    group.min_identity = 0
    group.match_min_color = "#ffffff"
    group.match_max_color = "#000000"
    group.match_fill_opacity = 0.75
    group.match_stroke_color = "none"
    group.match_stroke_width = 0
    group.collinearity_orientation_colors = {"plus": "#112233", "minus": "#445566"}

    row = SimpleNamespace(
        identity=95,
        qstart=1,
        qend=100,
        sstart=10,
        send=110,
        collinearity_block_id="block_0001",
        collinearity_orientation="minus",
        collinearity_color_mode="orientation",
    )

    drawing = Drawing(debug=False)
    drawing.add(group.generate_linear_match_path(row))

    assert 'fill="#445566"' in drawing.tostring()


@pytest.mark.linear
def test_orientation_identity_collinearity_uses_orientation_specific_identity_ramp() -> None:
    group = _build_collinearity_match_group()
    row = SimpleNamespace(
        identity=95,
        qstart=1,
        qend=100,
        sstart=10,
        send=110,
        collinearity_block_id="block_0001",
        collinearity_orientation="minus",
        collinearity_color_mode="orientation_identity",
    )

    path = group.generate_linear_match_path(row)

    assert path.attribs["fill"] == interpolate_color("#ffeeee", "#445566", 0.95)
    assert path.attribs["fill"] != "#445566"
    assert path.attribs["data-collinearity-color-mode"] == "orientation_identity"
    assert path.attribs["data-collinearity-orientation"] == "minus"
    assert float(path.attribs["data-identity-factor"]) == pytest.approx(0.95)


@pytest.mark.linear
def test_orientation_identity_collinearity_low_identity_is_paler_than_high_identity() -> None:
    group = _build_collinearity_match_group()
    high = group.generate_linear_match_path(
        SimpleNamespace(
            identity=95,
            qstart=1,
            qend=100,
            sstart=10,
            send=110,
            collinearity_block_id="block_high",
            collinearity_orientation="minus",
            collinearity_color_mode="orientation_identity",
        )
    )
    low = group.generate_linear_match_path(
        SimpleNamespace(
            identity=25,
            qstart=1,
            qend=100,
            sstart=10,
            send=110,
            collinearity_block_id="block_low",
            collinearity_orientation="minus",
            collinearity_color_mode="orientation_identity",
        )
    )

    assert _hex_distance(str(high.attribs["fill"]), "#445566") < _hex_distance(
        str(low.attribs["fill"]),
        "#445566",
    )
    assert _hex_distance(str(low.attribs["fill"]), "#ffeeee") < _hex_distance(
        str(high.attribs["fill"]),
        "#ffeeee",
    )


@pytest.mark.linear
def test_orientation_identity_collinearity_uses_two_orientation_gradients() -> None:
    group = _build_collinearity_match_group()
    plus = group.generate_linear_match_path(
        SimpleNamespace(
            identity=50,
            qstart=1,
            qend=100,
            sstart=10,
            send=110,
            collinearity_block_id="block_plus",
            collinearity_orientation="plus",
            collinearity_color_mode="orientation_identity",
        )
    )
    minus = group.generate_linear_match_path(
        SimpleNamespace(
            identity=50,
            qstart=1,
            qend=100,
            sstart=10,
            send=110,
            collinearity_block_id="block_minus",
            collinearity_orientation="minus",
            collinearity_color_mode="orientation_identity",
        )
    )

    assert plus.attribs["fill"] == interpolate_color("#eeeeee", "#112233", 0.5)
    assert minus.attribs["fill"] == interpolate_color("#ffeeee", "#445566", 0.5)
    assert plus.attribs["fill"] != minus.attribs["fill"]


@pytest.mark.linear
def test_average_identity_collinearity_ignores_orientation_colors() -> None:
    group = _build_collinearity_match_group()
    row = SimpleNamespace(
        identity=95,
        qstart=1,
        qend=100,
        sstart=10,
        send=110,
        collinearity_block_id="block_0001",
        collinearity_orientation="minus",
        collinearity_color_mode="average_identity",
    )

    path = group.generate_linear_match_path(row)

    assert path.attribs["fill"] == interpolate_color("#ffffff", "#000000", 0.95)
    assert path.attribs["fill"] != "#445566"


@pytest.mark.linear
def test_non_collinearity_match_path_does_not_add_metadata_attributes() -> None:
    group = PairWiseMatchGroup.__new__(PairWiseMatchGroup)
    group.canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=1000,
    )
    group.records = [_record("record_a", []), _record("record_b", [])]
    group.comparison_count = 1
    group.comparison_height = 40
    group.query_offset_x = 0
    group.subject_offset_x = 0
    group.query_alignment_offset_x = 0
    group.subject_alignment_offset_x = 0
    group.min_identity = 0
    group.match_min_color = "#ffffff"
    group.match_max_color = "#000000"
    group.match_fill_opacity = 0.75
    group.match_stroke_color = "none"
    group.match_stroke_width = 0

    row = SimpleNamespace(
        identity=95,
        qstart=1,
        qend=100,
        sstart=10,
        send=110,
        query_protein_id="qa0",
        subject_protein_id="sb0",
    )

    drawing = Drawing(debug=False)
    drawing.add(group.generate_linear_match_path(row))
    svg = drawing.tostring()

    assert "data-query-protein-id" not in svg
    assert "data-subject-protein-id" not in svg


@pytest.mark.linear
def test_collinearity_match_group_draws_inversions_above_plus_blocks() -> None:
    rows = [
        {
            **_hit_row("record_a", "record_b"),
            "qstart": 80,
            "qend": 120,
            "sstart": 120,
            "send": 80,
            "collinearity_block_id": "block_minus",
            "collinearity_orientation": "minus",
            "collinearity_color_mode": "orientation",
        },
        {
            **_hit_row("record_a", "record_b"),
            "qstart": 1,
            "qend": 200,
            "sstart": 1,
            "send": 200,
            "collinearity_block_id": "block_plus",
            "collinearity_orientation": "plus",
            "collinearity_color_mode": "orientation",
        },
    ]
    comparison_df = pd.DataFrame.from_records(rows)
    canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=1000,
        align_center=False,
    )
    blast_config = SimpleNamespace(
        fill_color="#d3d3d3",
        min_color="#ffffff",
        max_color="#000000",
        fill_opacity=1.0,
        stroke_color="none",
        stroke_width=0,
        identity=0,
        sequence_length_dict={},
    )
    match_group = PairWiseMatchGroup(
        canvas_config,
        blast_config.sequence_length_dict,
        comparison_df,
        40,
        1,
        blast_config,
        [_record("record_a", []), _record("record_b", [])],
    ).get_group()

    drawing = Drawing(debug=False)
    drawing.add(match_group)
    svg = drawing.tostring()

    assert svg.index('data-collinearity-block-id="block_plus"') < svg.index(
        'data-collinearity-block-id="block_minus"'
    )
