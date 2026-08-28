from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.analysis.collinearity as collinearity_module
import gbdraw.api.request_render as request_render_module
import gbdraw.diagrams.linear.assemble as linear_assemble_module
from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
    LosslessCollinearityParameters,
    build_collinearity_blocks_from_hits,
    calculate_collinearity_block_evalue,
    call_collinearity_blocks,
    cluster_lossless_collinearity_anchors,
    convert_collinearity_blocks_to_comparisons,
    convert_collinearity_blocks_to_pair_comparisons,
    iter_collinearity_search_pairs,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
    orthogroup_edges_to_lossless_collinearity_anchors,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.requests import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.protein_colinearity import (
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
    select_rbh_orthogroup_edges_from_directional_hits,
)
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import (
    LegendDrawingConfigurator,
    LegendMeasurement,
)
from gbdraw.configurators.blast import BlastMatchConfigurator
from gbdraw.core.color import interpolate_color
from gbdraw.core.text import calculate_bbox_dimensions
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.legend.table import prepare_legend_table
from gbdraw.render.groups.linear.legend import LegendGroup
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup
from gbdraw.exceptions import ValidationError


def _linear_legend_measurement(
    legend_position: str,
    legend_table: dict[str, dict[str, object]],
) -> tuple[LinearCanvasConfigurator, LegendMeasurement, GbdrawConfig]:
    cfg = GbdrawConfig.from_dict(
        load_config_toml("gbdraw.data", "config.toml")
    )
    profile = LinearRenderProfile(cfg)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1_000,
        profile=profile,
        legend=legend_position,
    )
    configurator = LegendDrawingConfigurator(
        color_table=None,
        default_colors=None,
        selected_features_set=[],
        profile=profile,
        gc_config=None,
        skew_config=None,
        feature_config=None,
        canvas_config=canvas_config,
    )
    return (
        canvas_config,
        configurator.measure_legend(
            legend_table,
            placement=canvas_config.legend_position,
            wrap_width=canvas_config.total_width,
        ),
        cfg,
    )


def _record(record_id: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(Seq("ATGAAATAG" * 80), id=record_id)
    record.features = features
    return record


def _cds(
    start: int,
    end: int,
    qualifiers: dict[str, list[str]],
    *,
    strand: int = 1,
) -> SeqFeature:
    merged = {"translation": ["MKK"]}
    merged.update(qualifiers)
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers=merged,
    )


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


def _fasta_ids(fasta_text: str) -> list[str]:
    return [
        line[1:].split(None, 1)[0]
        for line in str(fasta_text).splitlines()
        if line.startswith(">")
    ]


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
def test_lossless_collinearity_parameters_validate_domains() -> None:
    LosslessCollinearityParameters().validate()

    with pytest.raises(ValidationError, match="collinear_max_unit_gap"):
        LosslessCollinearityParameters(max_unit_gap=-1).validate()
    with pytest.raises(ValidationError, match="collinear_merge_orientation"):
        LosslessCollinearityParameters(merge_orientation="invalid").validate()  # type: ignore[arg-type]


@pytest.mark.linear
def test_collinearity_blocks_are_called_with_lossless_parameters() -> None:
    anchors = [_anchor(index, index) for index in range(8)]
    result = call_collinearity_blocks(
        anchors,
        params=LosslessCollinearityParameters(min_anchors=5, max_unit_gap=1),
    )

    assert len(result.blocks) == 1
    assert result.blocks[0].kind == "cluster"
    assert result.blocks[0].block_evalue is None


@pytest.mark.linear
def test_collinearity_block_evalue_uses_genomic_coordinates_not_unit_order() -> None:
    assert "calculate_collinearity_block_evalue" in collinearity_module.__all__
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
        params=LosslessCollinearityParameters(),
    )

    assert len(result.blocks) == 1
    assert result.blocks[0].kind == "singleton"
    assert result.blocks[0].anchors[0].query_protein_id == "q0"


@pytest.mark.linear
def test_min_anchors_drops_singleton_blocks() -> None:
    anchor = _anchor(0, 0)
    result = call_collinearity_blocks(
        [anchor],
        params=LosslessCollinearityParameters(min_anchors=2),
    )

    assert result.blocks == ()
    assert result.unblocked_anchors == (anchor,)


@pytest.mark.linear
def test_native_block_caller_detects_plus_and_minus_orientation() -> None:
    params = LosslessCollinearityParameters(
        min_anchors=3,
        max_unit_gap=1,
        merge_orientation="order",
    )

    plus = call_collinearity_blocks([_anchor(0, 0), _anchor(1, 1), _anchor(2, 2)], params=params)
    minus = call_collinearity_blocks([_anchor(0, 2), _anchor(1, 1), _anchor(2, 0)], params=params)

    assert len(plus.blocks) == 1
    assert plus.blocks[0].orientation == "plus"
    assert [anchor.query_order for anchor in plus.blocks[0].anchors] == [0, 1, 2]
    assert len(minus.blocks) == 1
    assert minus.blocks[0].orientation == "minus"


@pytest.mark.linear
def test_single_anchor_block_orientation_uses_feature_strands() -> None:
    params = LosslessCollinearityParameters(min_anchors=1, max_unit_gap=0)

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
    params = LosslessCollinearityParameters(
        min_anchors=1,
        max_unit_gap=0,
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
    params = LosslessCollinearityParameters(
        min_anchors=1,
        max_unit_gap=0,
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
@pytest.mark.parametrize(
    ("anchors", "max_unit_gap", "expected_zero", "expected_one"),
    [
        (
            [
                _anchor(0, 0),
                _anchor(1, 1),
                _anchor(2, 5, subject_strand=-1),
                _anchor(3, 3),
                _anchor(4, 4),
            ],
            1,
            [[0, 1, 3, 4], [2]],
            [[0, 1, 3, 4], [2]],
        ),
        (
            [
                _anchor(0, 0),
                _anchor(1, 1),
                _anchor(2, 2, subject_strand=-1),
                _anchor(3, 3),
                _anchor(4, 4),
            ],
            1,
            [[0, 1], [2], [3, 4]],
            [[0, 1, 3, 4], [2]],
        ),
        (
            [
                _anchor(0, 0),
                _anchor(1, 1),
                _anchor(2, 3, subject_strand=-1),
                _anchor(3, 2),
                _anchor(4, 4),
                _anchor(5, 5),
            ],
            2,
            [[0, 1], [2], [3], [4, 5]],
            [[0, 1], [2], [3], [4, 5]],
        ),
    ],
    ids=(
        "zero-interior-conflicts",
        "one-interior-conflict",
        "two-interior-conflicts",
    ),
)
def test_lossless_max_conflicts_threshold_controls(
    anchors: list[CollinearityAnchor],
    max_unit_gap: int,
    expected_zero: list[list[int]],
    expected_one: list[list[int]],
) -> None:
    input_members = {
        (anchor.query_order, anchor.subject_order) for anchor in anchors
    }
    for max_conflicts, expected in enumerate((expected_zero, expected_one)):
        result = cluster_lossless_collinearity_anchors(
            anchors,
            params=LosslessCollinearityParameters(
                max_unit_gap=max_unit_gap,
                max_conflicts=max_conflicts,
                merge_orientation="strand",
            ),
        )

        assert [
            [anchor.query_order for anchor in block.anchors]
            for block in result.blocks
        ] == expected
        assert sum(len(block.anchors) for block in result.blocks) == len(anchors)
        assert {
            (anchor.query_order, anchor.subject_order)
            for block in result.blocks
            for anchor in block.anchors
        } == input_members


@pytest.mark.linear
def test_typed_request_max_conflicts_reaches_real_collinearity_consumer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _record(
            "record_a",
            [
                _cds(
                    index * 12,
                    index * 12 + 9,
                    {"locus_tag": [f"q{index}"], "protein_id": [f"q{index}"]},
                )
                for index in range(5)
            ],
        ),
        _record(
            "record_b",
            [
                _cds(
                    index * 12,
                    index * 12 + 9,
                    {"locus_tag": [f"s{index}"], "protein_id": [f"s{index}"]},
                    strand=-1 if index == 2 else 1,
                )
                for index in range(5)
            ],
        ),
    ]

    def fake_search(query_fasta, subject_fasta, **_kwargs):
        query_ids = _fasta_ids(query_fasta)
        subject_ids = _fasta_ids(subject_fasta)
        rows = [
            _hit_row(query_id, subject_id)
            for query_id, subject_id in zip(query_ids, subject_ids, strict=True)
        ]
        return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)

    monkeypatch.setattr(collinearity_module, "_run_losatp_search", fake_search)

    results: dict[int, CollinearityResult] = {}
    for max_conflicts in (0, 1):
        params = LosslessCollinearityParameters(
            max_unit_gap=1,
            max_conflicts=max_conflicts,
            merge_orientation="strand",
        )
        request = LinearDiagramRequest(
            records=tuple(
                RecordInput(source=InMemoryRecordSource(record)) for record in records
            ),
            options=LinearDiagramOptions(
                protein_blastp_mode="collinear",
                collinearity_unit_mode="cds",
                collinearity_params=params,
            ),
        )
        prepared = request_render_module.build_request_diagram(request)
        assert prepared.request.options.collinearity_params is params
        assert prepared.linear_metadata is not None
        result = prepared.linear_metadata.collinearity_result
        assert result is not None
        results[max_conflicts] = result

    assert [
        [anchor.query_order for anchor in block.anchors]
        for block in results[0].blocks
    ] == [[0, 1], [2], [3, 4]]
    assert [
        [anchor.query_order for anchor in block.anchors]
        for block in results[1].blocks
    ] == [[0, 1, 3, 4], [2]]


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
        params=LosslessCollinearityParameters(min_anchors=3, max_unit_gap=1),
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
def test_build_blocks_from_hits_adjacent_scope_uses_explicit_adjacent_row_pairs() -> None:
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
        params=LosslessCollinearityParameters(min_anchors=1),
        edge_mode="one_to_one",
        search_scope="adjacent",
        comparison_pairs=((0, 2),),
    )

    assert len(result.blocks) == 1
    assert result.blocks[0].query_record_index == 0
    assert result.blocks[0].subject_record_index == 2


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
def test_build_blocks_from_hits_can_render_selected_non_adjacent_pairs() -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["a0"], "protein_id": ["a0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["b0"], "protein_id": ["b0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["c0"], "protein_id": ["c0"]})]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
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
        comparison_pairs=((0, 2),),
    )
    comparisons = convert_collinearity_blocks_to_pair_comparisons(
        result,
        records=records,
        color_mode="orientation",
    )

    assert {(block.query_record_index, block.subject_record_index) for block in result.blocks} == {(0, 2)}
    assert comparisons[(0, 2)].shape[0] == 1


@pytest.mark.linear
def test_orthogroup_collinearity_rbh_search_defaults_to_top_candidate(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_limits: list[int | None] = []
    observed_pairs: list[tuple[str, str]] = []

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, ncbi_blastp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert ncbi_blastp_bin is None
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

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, ncbi_blastp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert ncbi_blastp_bin is None
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

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, ncbi_blastp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert ncbi_blastp_bin is None
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

    def fake_search(query_fasta, subject_fasta, *, losatp_bin, ncbi_blastp_bin, losatp_threads, candidate_limit, max_hsps_per_subject, runner):
        assert ncbi_blastp_bin is None
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
def test_linear_cli_forwards_collinearity_options(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"]})]), _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"]})])]
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _path: None)

    def fake_render_request(request, **_kwargs):
        resolved = request_render_module.resolve_request(request)
        captured["canonical_request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render_request)

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
            "--collinear_search_scope",
            "all",
            "--collinear_max_conflicts_in_merge_gap",
            "3",
            "--collinear_color_mode",
            "orientation",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    options = captured["canonical_request"].options
    params = options.collinearity_params
    assert isinstance(params, LosslessCollinearityParameters)
    assert params.min_anchors == 1
    assert params.max_unit_gap == 0
    assert params.max_diagonal_drift == 0
    assert params.max_conflicts == 3
    assert options.collinearity_unit_mode == "cds"
    assert options.collinearity_anchor_mode == "rbh"
    assert options.collinearity_search_scope == "all"
    assert options.orthogroup_membership_mode == "anchor_core_v1"
    assert options.orthogroup_member_max_hits == 5
    assert options.collinear_max_paralog_links_per_orthogroup == 2
    assert options.protein_blastp_mode == "collinear"
    assert options.protein_comparisons is None
    assert options.collinearity_color_mode == "orientation"


@pytest.mark.linear
def test_linear_cli_help_omits_removed_collinearity_options(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        linear_cli_module._get_args(["--help"])

    assert exc_info.value.code == 0
    help_text = capsys.readouterr().out
    removed_option_names = [
        "orthogroup_membership_mode",
        "orthogroup-membership-mode",
        "orthogroup_member_max_hits",
        "orthogroup-member-max-hits",
        "collinear_anchor_mode",
        "collinear-anchor-mode",
        "collinear_orthogroup_edge_mode",
        "collinear-orthogroup-edge-mode",
        "collinear_blocks",
        "collinear-blocks",
        "save_collinear_blocks",
        "save-collinear-blocks",
        "collinear_block_merge_gap",
        "collinear-block-merge-gap",
        "collinear_singleton_merge_gap",
        "collinear-singleton-merge-gap",
        "collinear_gap_penalty",
        "collinear-gap-penalty",
        "collinear_nearby_duplicate_window",
        "collinear-nearby-duplicate-window",
        "collinear_score_mode",
        "collinear-score-mode",
        "collinear_constant_anchor_score",
        "collinear-constant-anchor-score",
        "collinear_min_block_score",
        "collinear-min-block-score",
        "collinear_block_evalue",
        "collinear-block-evalue",
    ]
    for option_name in removed_option_names:
        assert f"--{option_name}" not in help_text


@pytest.mark.linear
def test_linear_cli_help_uses_underscore_option_aliases(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        linear_cli_module._get_args(["--help"])

    assert exc_info.value.code == 0
    help_text = capsys.readouterr().out
    visible_option_names = [
        "losatp_bin",
        "losatp_threads",
        "protein_blastp_mode",
        "collinear_min_anchors",
        "collinear_max_unit_gap",
        "collinear_color_mode",
        "keep_definition_left_aligned",
        "pairwise_match_style",
        "save_session",
        "session_output",
    ]
    hidden_option_names = [
        "losatp-bin",
        "losatp-threads",
        "protein-blastp-mode",
        "collinear-min-anchors",
        "collinear-max-unit-gap",
        "collinear-color-mode",
        "keep-definition-left-aligned",
        "pairwise-match-style",
        "save-session",
        "session-output",
    ]
    for option_name in visible_option_names:
        assert f"--{option_name}" in help_text
    for option_name in hidden_option_names:
        assert f"--{option_name}" not in help_text


@pytest.mark.linear
@pytest.mark.parametrize(
    "option_args",
    [
        ["--losatp-bin", "losatp"],
        ["--losatp-threads", "2"],
        ["--protein-blastp-mode", "pairwise"],
        ["--collinear-min-anchors", "1"],
        ["--collinear-max-unit-gap", "0"],
        ["--collinear-color-mode", "orientation"],
        ["--keep-definition-left-aligned"],
        ["--pairwise-match-style", "ribbon"],
        ["--save-session"],
        ["--session-output", "session.json"],
        ["--collinear_max_gene_gap", "3"],
    ],
)
def test_linear_cli_rejects_removed_hidden_aliases(
    option_args: list[str],
) -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(["--gbk", "dummy.gb", *option_args])


@pytest.mark.linear
@pytest.mark.parametrize(
    "option_name",
    [
        "--collinear_blocks",
        "--save_collinear_blocks",
        "--collinear-blocks",
        "--save-collinear-blocks",
        "--collinear_block_merge_gap",
        "--collinear_singleton_merge_gap",
        "--collinear_gap_penalty",
        "--collinear_nearby_duplicate_window",
        "--collinear_score_mode",
        "--collinear_constant_anchor_score",
        "--collinear_min_block_score",
        "--collinear_block_evalue",
        "--collinear-block-merge-gap",
        "--collinear-singleton-merge-gap",
        "--collinear-gap-penalty",
        "--collinear-nearby-duplicate-window",
        "--collinear-score-mode",
        "--collinear-constant-anchor-score",
        "--collinear-min-block-score",
        "--collinear-block-evalue",
    ],
)
def test_linear_cli_rejects_removed_collinear_block_options(option_name: str) -> None:
    with pytest.raises(SystemExit) as exc_info:
        linear_cli_module._get_args(["--gbk", "a.gb", "b.gb", option_name, "blocks.tsv"])

    assert exc_info.value.code == 2


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
def test_web_losatp_blastp_payload_helper_returns_collinear_rows(tmp_path: Path) -> None:
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

    pairs_path = tmp_path / "losatp-pairs.json"
    pairs_path.write_text(json.dumps(payload), encoding="utf-8")
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
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
        25,
        1,
        2,
        "adjacent",
        "rbh",
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    assert result["collinearityResult"]["schema"] == 2
    assert result["collinearityResult"]["kind"] == "result"
    assert result["collinearityResult"]["value"]["type"] == "CollinearityResult"
    typed_fields = result["collinearityResult"]["value"]["fields"]
    assert typed_fields["blocks"][0]["fields"]["blockId"] == "block_0001"
    assert len(typed_fields["orthogroups"]["fields"]["orthogroups"]) == 8
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
    assert rows[0]["query_feature_index"] == ";".join(map(str, range(8)))
    assert rows[0]["subject_feature_index"] == ";".join(map(str, range(8)))
    assert rows[0]["query_feature_svg_id"] == ";".join(
        f"feature_qa{index}" for index in range(8)
    )
    assert rows[0]["query_view_feature_svg_id"] == rows[0][
        "query_feature_svg_id"
    ]
    assert set(result) == {"pairs", "collinearityResult", "cache"}


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_uses_rbh_collinear_anchor_mode(
    tmp_path: Path,
) -> None:
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

    pairs_path = tmp_path / "losatp-pairs.json"
    pairs_path.write_text(json.dumps(payload), encoding="utf-8")
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
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
        str(pairs_path),
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
        25,
        1,
        2,
    )
    min_two = json.loads(str(raw_min_two))

    assert "error" not in min_two
    assert min_two["pairs"][0]["rows"] == []


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_applies_collinear_search_scope(
    tmp_path: Path,
) -> None:
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

    pairs_path = tmp_path / "losatp-pairs.json"
    pairs_path.write_text(json.dumps(payload), encoding="utf-8")

    def convert(scope: str, input_payload: Path = pairs_path) -> dict[str, object]:
        raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
            str(input_payload),
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
            25,
            1,
            2,
            scope,
            "rbh",
        )
        return json.loads(str(raw_result))

    adjacent = convert("adjacent")
    all_records = convert("all")
    multi_row_payload = json.loads(json.dumps(payload))
    for pair in multi_row_payload["pairs"]:
        pair["displayPair"] = (
            pair["queryIndex"] == 0 and pair["subjectIndex"] == 2
        )
    pairs_path = tmp_path / "losatp-pairs.json"
    pairs_path.write_text(json.dumps(multi_row_payload), encoding="utf-8")
    multi_row_adjacent = convert("adjacent", pairs_path)

    assert "error" not in adjacent
    assert "error" not in all_records
    assert "error" not in multi_row_adjacent
    def member_sets(result: dict[str, object]) -> list[set[str]]:
        groups = result["collinearityResult"]["value"]["fields"][
            "orthogroups"
        ]["fields"]["orthogroups"]
        return [
            {member["fields"]["proteinId"] for member in members}
            for members in groups.values()
        ]

    adjacent_member_sets = member_sets(adjacent)
    all_member_sets = member_sets(all_records)
    assert set(adjacent) == {"pairs", "collinearityResult", "cache"}
    assert set(all_records) == {"pairs", "collinearityResult", "cache"}
    assert {"a0", "b0"} in adjacent_member_sets
    assert {"a0", "b0", "c0"} in all_member_sets
    assert {
        (
            block["fields"]["queryRecordIndex"],
            block["fields"]["subjectRecordIndex"],
        )
        for block in multi_row_adjacent["collinearityResult"]["value"]["fields"]["blocks"]
    } == {(0, 2)}
    assert all(
        block["fields"]["subjectRecordIndex"]
        == block["fields"]["queryRecordIndex"] + 1
        for block in all_records["collinearityResult"]["value"]["fields"]["blocks"]
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
        search_scope="adjacent",
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
    assert set(rows["group_kind"]) == {"collinear_gene_group"}
    assert set(rows["group_scope"]) == {"adjacent_local"}
    assert set(rows["collinear_group_scope"]) == {"adjacent_local"}


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
    assert 'data-query-feature-svg-id="feature_qa2_record_1"' in svg_text
    assert 'data-subject-feature-svg-id="feature_sb2_record_2"' in svg_text
    assert 'data-query-stable-feature-svg-id="feature_qa2"' in svg_text
    assert 'data-subject-stable-feature-svg-id="feature_sb2"' in svg_text


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
        profile=LinearRenderProfile(
            GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
        ),
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
        profile=LinearRenderProfile(
            GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
        ),
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
        show_gc=False,
        show_skew=False,
        show_depth=False,
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
        show_gc=False,
        show_skew=False,
        show_depth=False,
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
        show_gc=False,
        show_skew=False,
        show_depth=False,
    )

    assert "Collinear" in legend_table
    assert "Inverted" in legend_table
    assert "Pairwise match identity" not in legend_table
    assert legend_table["Collinear"]["min_color"] == "#eeeeee"
    assert legend_table["Inverted"]["max_color"] == "#445566"


@pytest.mark.linear
def test_orientation_identity_pairwise_legend_renders_collinear_above_inverted() -> None:
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

    canvas_config, legend_measurement, cfg = _linear_legend_measurement(
        "bottom",
        legend_table,
    )
    drawing = Drawing(debug=False)
    legend_group = LegendGroup(
        canvas_config,
        legend_measurement,
        legend_table,
        cfg=cfg,
    )
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
        "Collinear",
        legend_measurement.font_family,
        legend_measurement.font_size,
        legend_group.dpi,
    )
    assert bar_x == pytest.approx(
        label_width + 0.2 * legend_measurement.color_rect_size
    )


@pytest.mark.linear
def test_vertical_linear_pairwise_legend_aligns_to_feature_color_left_edge() -> None:
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

    canvas_config, legend_measurement, cfg = _linear_legend_measurement(
        "right",
        legend_table,
    )
    drawing = Drawing(debug=False)
    legend_group = LegendGroup(
        canvas_config,
        legend_measurement,
        legend_table,
        cfg=cfg,
    )
    drawing.add(legend_group.get_group())
    root = ET.fromstring(drawing.tostring())
    namespace = {"svg": "http://www.w3.org/2000/svg"}

    vertical_legend = root.find(".//svg:g[@id='legend_vertical']", namespace)
    assert vertical_legend is not None
    feature_legend = vertical_legend.find("svg:g[@id='feature_legend_v']", namespace)
    pairwise_legend = vertical_legend.find("svg:g[@id='pairwise_legend_v']", namespace)
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
        feature_key,
        legend_measurement.font_family,
        legend_measurement.font_size,
        legend_group.dpi,
    )
    text_x_offset = (22 / 14) * legend_measurement.color_rect_size
    feature_width = text_x_offset + text_width
    assert legend_measurement.linear_layout is not None
    gradient_layout = legend_measurement.linear_layout.vertical.gradient
    assert gradient_layout is not None
    pairwise_alignment_width = gradient_layout.width

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

    canvas_config, legend_measurement, cfg = _linear_legend_measurement(
        "right",
        legend_table,
    )
    drawing = Drawing(debug=False)
    legend_group = LegendGroup(
        canvas_config,
        legend_measurement,
        legend_table,
        cfg=cfg,
    )
    drawing.add(legend_group.get_group())
    root = ET.fromstring(drawing.tostring())
    namespace = {"svg": "http://www.w3.org/2000/svg"}

    vertical_legend = root.find(".//svg:g[@id='legend_vertical']", namespace)
    assert vertical_legend is not None
    feature_legend = vertical_legend.find("svg:g[@id='feature_legend_v']", namespace)
    pairwise_legend = vertical_legend.find("svg:g[@id='pairwise_legend_v']", namespace)
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
        feature_key,
        legend_measurement.font_family,
        legend_measurement.font_size,
        legend_group.dpi,
    )
    text_x_offset = (22 / 14) * legend_measurement.color_rect_size
    feature_width = text_x_offset + text_width
    assert legend_measurement.linear_layout is not None
    gradient_layout = legend_measurement.linear_layout.vertical.gradient
    assert gradient_layout is not None
    pairwise_alignment_width = gradient_layout.width
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
        cfg=apply_config_overrides(
            None,
            {"canvas.show_gc": False, "canvas.show_skew": False},
        ),
        blast_files=["/virtual/blast_0.txt"],
        legend="bottom",
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
