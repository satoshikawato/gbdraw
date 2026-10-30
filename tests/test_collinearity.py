from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.analysis.collinearity as collinearity_module
import gbdraw.analysis.protein_colinearity as protein_colinearity_module
import gbdraw.api.request_render as request_render_module
import gbdraw.diagrams.linear.assemble as linear_assemble_module
import gbdraw.io.comparisons as comparisons_io_module
from gbdraw.legend.sources import ComparisonLegendSource
from tests.utils.comparison_tables import converted_pair_rows, ribbon_group
from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
    LosslessCollinearityParameters,
    build_collinearity_blocks_from_hits,
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
from gbdraw.api.options import LinearDiagramOptions, LosatSearchOptions
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
from gbdraw.configurators.legend import measure_legend
from gbdraw.configurators.blast import BlastMatchConfigurator
from gbdraw.core.color import interpolate_color
from gbdraw.core.text import drawn_text_box
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.legend.table import (
    PairwiseIdentityGradient,
    PairwiseIdentityLegend,
    prepare_legend_table,
)
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
        profile=profile,
        canvas_config=canvas_config,
    )
    return (
        canvas_config,
        measure_legend(
            configurator,
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


def _record_local_mixed_fixture():
    """Issue #460 evidence shared by analysis, catalog, session and browser checks."""
    records = [
        _record(record_id, [
            _cds(start, start + 299, {
                "locus_tag": [protein_id], "protein_id": [protein_id],
                "translation": ["M" * 100],
            })
            for protein_id, start in features
        ])
        for record_id, features in (
            ("record_a", [("a0", 100), ("a1", 700), ("a3", 1200), ("a4", 1600), ("a2", 2100)]),
            ("record_b", [("b0", 120), ("b1", 900)]),
        )
    ]
    for record in records:
        record.seq = Seq("ATG" * 900)
        record.annotations["molecule_type"] = "DNA"
    hits = {
        pair: pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)
        for pair, rows in {
            (0, 0): [
                _hit_row("a0", "a1", 280), _hit_row("a1", "a0", 280),
                _hit_row("a3", "a4", 300), _hit_row("a4", "a3", 300),
            ],
            (0, 1): [_hit_row("a0", "b0", 300)],
            (1, 0): [_hit_row("b0", "a0", 300)],
        }.items()
    }
    return records, extract_cds_proteins(records), hits


@pytest.mark.parametrize("search_scope", ["adjacent", "all"])
@pytest.mark.parametrize("edge_mode", ["rbh", "all_hits"])
def test_collinear_record_local_mixed_fixture(search_scope, edge_mode):
    records, extraction, hits = _record_local_mixed_fixture()
    result = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        hits, extraction, records=records, search_scope=search_scope, edge_mode=edge_mode,
    )
    groups = result.orthogroups
    assert groups is not None
    assert {
        group_id: {member.protein_id for member in members}
        for group_id, members in groups.orthogroups.items()
    } == {"og_1": {"a0", "a1", "b0"}, "og_2": {"a3", "a4"}}
    assert groups.scope_by_orthogroup_id == {"og_1": "cross_record", "og_2": "record_local"}
    assert groups.source_record_index_by_orthogroup_id == {"og_2": 0}
    assert {member.record_id for member in groups.orthogroups["og_2"]} == {"record_a"}
    assert groups.member_by_protein_id["a1"].role == "inparalog"
    assert {member.role for member in groups.orthogroups["og_2"]} == {"local_paralog"}
    assert {"a2", "b1"}.isdisjoint(groups.member_by_protein_id)

    # Removing only local evidence gives the pre-fix adjacent geometry baseline.
    geometry_only = collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        {pair: table for pair, table in hits.items() if pair[0] != pair[1]},
        extraction, records=records, search_scope=search_scope, edge_mode=edge_mode,
    )
    assert len(result.blocks) == len(geometry_only.blocks) == 1
    anchor = result.blocks[0].anchors[0]
    assert anchor.query_orthogroup_member_count == 2
    assert (replace(result.blocks[0], anchors=(
        replace(anchor, query_orthogroup_member_count=1),
    )),) == geometry_only.blocks
    assert result.unblocked_anchors == geometry_only.unblocked_anchors == ()
    assert {anchor.orthogroup_id for block in result.blocks for anchor in block.anchors} == {"og_1"}
    comparisons = convert_collinearity_blocks_to_comparisons(result, records=records)
    baseline = convert_collinearity_blocks_to_comparisons(geometry_only, records=records)
    assert len(comparisons) == len(baseline) == 1
    # Complete membership is semantic metadata, not a geometry coordinate.
    assert comparisons[0]["query_orthogroup_member_count"].tolist() == ["2"]
    assert comparisons[0]["subject_orthogroup_member_count"].tolist() == ["1"]
    geometry_columns = comparisons[0].columns.difference([
        "query_orthogroup_member_count", "subject_orthogroup_member_count",
    ])
    pd.testing.assert_frame_equal(
        comparisons[0][geometry_columns], baseline[0][geometry_columns],
    )


@pytest.mark.parametrize("search_scope", ["adjacent", "all"])
@pytest.mark.parametrize("edge_mode", ["all_hits"])
def test_collinear_scope_partitions_table_identity(monkeypatch, search_scope, edge_mode):
    records, extraction, hits = _record_local_mixed_fixture()
    records.append(_record("record_c", [_cds(0, 9, {"protein_id": ["c0"]})]))
    extraction = extract_cds_proteins(records)
    for pair, query, subject in (
        ((1, 1), "b0", "b1"), ((2, 2), "c0", "c0"),
        ((1, 2), "b0", "c0"), ((2, 1), "c0", "b0"),
        ((0, 2), "a2", "c0"), ((2, 0), "c0", "a2"),
    ):
        hits[pair] = pd.DataFrame.from_records([_hit_row(query, subject)], columns=COMPARISON_COLUMNS)
    captured = {}
    infer = collinearity_module.select_rbh_orthogroup_edges_from_directional_hits
    select_geometry = collinearity_module._select_orthogroup_edges

    def capture_semantic(tables, *args, **kwargs):
        captured["semantic"] = tables
        return infer(tables, *args, **kwargs)

    def capture_geometry(tables, **kwargs):
        captured["geometry"] = tables
        return select_geometry(tables, **kwargs)

    monkeypatch.setattr(collinearity_module, "select_rbh_orthogroup_edges_from_directional_hits", capture_semantic)
    monkeypatch.setattr(collinearity_module, "_select_orthogroup_edges", capture_geometry)
    collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        hits, extraction, records=records, search_scope=search_scope, edge_mode=edge_mode,
    )
    cross_pairs = {(0, 1), (1, 0), (1, 2), (2, 1)}
    if search_scope == "all":
        cross_pairs |= {(0, 2), (2, 0)}
    assert set(captured["semantic"]) == cross_pairs | {(0, 0), (1, 1), (2, 2)}
    assert set(captured["geometry"]) == cross_pairs
    for tables in captured.values():
        assert all(table is hits[pair] for pair, table in tables.items())


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
    return ribbon_group(
        [_record("record_a", []), _record("record_b", [])],
        collinearity_orientation_colors={"plus": "#112233", "minus": "#445566"},
        collinearity_orientation_min_colors={"plus": "#eeeeee", "minus": "#ffeeee"},
    )


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
def test_collinearity_anchor_mode_defaults_to_rbh_and_refuses_one_to_one() -> None:
    assert normalize_collinearity_anchor_mode(None) == "rbh"
    assert normalize_collinearity_anchor_mode("reciprocal-best") == "rbh"
    # The removed anchor mode one_to_one and its spellings (OV-2303); a stored
    # value is read as rbh by the request and Session readers instead.
    for retired in ("one_to_one", "top1", "mutual-best"):
        with pytest.raises(ValidationError):
            normalize_collinearity_anchor_mode(retired)


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
def test_build_blocks_from_hits_rbh_removes_query_to_many_noise() -> None:
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
    rbh = build_collinearity_blocks_from_hits(
        [hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="rbh",
        reverse_hits_by_pair=[reverse_hits],
    )

    assert len(all_hits.blocks) == 2
    assert all_hits.blocks[0].anchors[0].subject_protein_id == "sb0"
    assert all_hits.blocks[0].anchors[0].subject_orthogroup_member_count == 2
    assert len(rbh.blocks) == 1
    assert rbh.blocks[0].anchors[0].subject_protein_id == "sb0"


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

    all_hits = build_collinearity_blocks_from_hits(
        [forward_hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="all",
    )
    rbh = build_collinearity_blocks_from_hits(
        [forward_hits],
        extraction,
        records=records,
        params=params,
        anchor_mode="rbh",
        reverse_hits_by_pair=[reverse_hits],
    )

    assert sum(len(block.anchors) for block in all_hits.blocks) == 2
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
@pytest.mark.parametrize(
    ("unit_mode", "expected_unit_kind"),
    (("auto", "locus"), ("cds", "cds"), ("locus", "locus")),
)
def test_typed_request_max_conflicts_reaches_real_collinearity_consumer(
    monkeypatch: pytest.MonkeyPatch,
    unit_mode: str,
    expected_unit_kind: str,
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

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

    params_by_case = {
        "explicit-zero": LosslessCollinearityParameters(
            max_unit_gap=1,
            max_conflicts=0,
            merge_orientation="strand",
        ),
        "default": LosslessCollinearityParameters(
            max_unit_gap=1,
            merge_orientation="strand",
        ),
    }
    results: dict[str, CollinearityResult] = {}
    for case, params in params_by_case.items():
        request = LinearDiagramRequest(
            records=tuple(
                RecordInput(source=InMemoryRecordSource(record)) for record in records
            ),
            options=LinearDiagramOptions(
                losat_search=LosatSearchOptions(
                    program="losatp",
                    losatp_mode="collinear",
                ),
                collinearity_unit_mode=unit_mode,
                collinearity_params=params,
            ),
        )
        prepared = request_render_module.build_request_diagram(request)
        assert prepared.request.options.collinearity_params is params
        assert prepared.request.options.collinearity_unit_mode == unit_mode
        assert prepared.linear_metadata is not None
        result = prepared.linear_metadata.collinearity_result
        assert result is not None
        assert {
            anchor.query_unit_kind
            for block in result.blocks
            for anchor in block.anchors
        } == {expected_unit_kind}
        results[case] = result

    assert [
        [anchor.query_order for anchor in block.anchors]
        for block in results["explicit-zero"].blocks
    ] == [[0, 1], [2], [3, 4]]
    assert [
        [anchor.query_order for anchor in block.anchors]
        for block in results["default"].blocks
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


def _block_orders(result: CollinearityResult) -> list[list[tuple[int, int]]]:
    return [
        [(anchor.query_order, anchor.subject_order) for anchor in block.anchors]
        for block in result.blocks
    ]


@pytest.mark.linear
def test_lossless_chaining_keeps_a_run_when_each_query_unit_has_a_second_anchor() -> None:
    # OV-2294: a second anchor per query unit used to break the run into 12 singletons.
    anchors = [_anchor(index, index) for index in range(6)] + [
        _anchor(index, index + 50) for index in range(6)
    ]

    result = cluster_lossless_collinearity_anchors(anchors)

    assert _block_orders(result) == [
        [(index, index) for index in range(6)],
        [(index, index + 50) for index in range(6)],
    ]
    assert [block.kind for block in result.blocks] == ["cluster", "cluster"]


@pytest.mark.linear
def test_lossless_chaining_keeps_a_run_past_a_lone_second_anchor() -> None:
    # OV-2294: one unrelated second anchor in the run is its own block, not a split.
    result = cluster_lossless_collinearity_anchors(
        [_anchor(0, 0), _anchor(1, 1), _anchor(1, 5), _anchor(2, 2)]
    )

    assert _block_orders(result) == [[(0, 0), (1, 1), (2, 2)], [(1, 5)]]
    assert [block.kind for block in result.blocks] == ["cluster", "singleton"]


@pytest.mark.linear
def test_lossless_chaining_bridges_an_unrelated_anchor_within_the_unit_gap() -> None:
    # OV-2295: an anchor off the run between two anchors within max_unit_gap.
    params = LosslessCollinearityParameters(max_unit_gap=1)

    assert _block_orders(
        cluster_lossless_collinearity_anchors([_anchor(0, 0), _anchor(2, 2)], params=params)
    ) == [[(0, 0), (2, 2)]]
    assert _block_orders(
        cluster_lossless_collinearity_anchors(
            [_anchor(0, 0), _anchor(1, 9), _anchor(2, 2)], params=params
        )
    ) == [[(0, 0), (2, 2)], [(1, 9)]]


@pytest.mark.linear
def test_lossless_chaining_joins_runs_across_an_off_diagonal_run() -> None:
    # OV-2296: an off-diagonal run with no anchor in the gap box no longer blocks the join.
    params = LosslessCollinearityParameters(max_unit_gap=2)
    anchors = [
        _anchor(0, 0),
        _anchor(1, 1),
        _anchor(2, 40),
        _anchor(3, 41),
        _anchor(4, 4),
        _anchor(5, 5),
    ]

    result = cluster_lossless_collinearity_anchors(anchors, params=params)

    assert _block_orders(result) == [
        [(0, 0), (1, 1), (4, 4), (5, 5)],
        [(2, 40), (3, 41)],
    ]


@pytest.mark.linear
def test_lossless_chaining_keeps_runs_when_inference_adds_coortholog_anchors() -> None:
    # OV-2294: the CLI default (inference on) gave 8 singletons for a duplicated subject region.
    records = [
        _record(
            f"record_{record}",
            [
                _cds(index * 60, index * 60 + 30, {"locus_tag": [f"r{record}g{index}"], "protein_id": [f"r{record}g{index}"]})
                for index in range(count)
            ],
        )
        for record, count in ((0, 4), (1, 8))
    ]
    extraction = extract_cds_proteins(records)
    forward = [_hit_row(f"r0g{i}", f"r1g{i}", 400.0) for i in range(4)] + [
        _hit_row(f"r0g{i}", f"r1g{i + 4}", 390.0) for i in range(4)
    ]
    reverse = [_hit_row(f"r1g{i}", f"r0g{i}", 400.0) for i in range(4)] + [
        _hit_row(f"r1g{i + 4}", f"r0g{i}", 390.0) for i in range(4)
    ]
    tables = {
        (0, 1): pd.DataFrame.from_records(forward, columns=COMPARISON_COLUMNS),
        (1, 0): pd.DataFrame.from_records(reverse, columns=COMPARISON_COLUMNS),
        (0, 0): pd.DataFrame(columns=COMPARISON_COLUMNS),
        (1, 1): pd.DataFrame(columns=COMPARISON_COLUMNS),
    }

    blocks_by_inference = {
        infer: _block_orders(
            collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
                tables,
                extraction,
                records=records,
                infer_orthogroups=infer,
                edge_mode="rbh",
            )
        )
        for infer in (False, True)
    }

    assert blocks_by_inference[False] == [[(index, index) for index in range(4)]]
    assert blocks_by_inference[True][0] == [(index, index) for index in range(4)]
    assert all(len(block) == 4 for block in blocks_by_inference[True])


def _random_anchors(rng, count: int) -> list[CollinearityAnchor]:
    anchors = []
    seen = set()
    for _ in range(count):
        q, s = rng.randrange(20), rng.randrange(20)
        if (q, s) in seen:
            continue
        seen.add((q, s))
        anchors.append(
            _anchor(
                q,
                s,
                bitscore=float(rng.choice((50, 100, 200))),
                query_strand=rng.choice((1, -1, None)),
                subject_strand=rng.choice((1, -1, None)),
            )
        )
    return anchors


@pytest.mark.linear
@pytest.mark.parametrize("seed", range(20))
def test_lossless_chaining_is_deterministic_under_input_order(seed: int) -> None:
    import random

    rng = random.Random(seed)
    anchors = [
        replace(anchor, query_record_index=pair[0], subject_record_index=pair[1])
        for pair in ((0, 1), (1, 2))
        for anchor in _random_anchors(rng, 30)
    ]
    params = LosslessCollinearityParameters(
        max_unit_gap=seed % 4,
        max_diagonal_drift=seed % 3,
        max_conflicts=seed % 2,
        merge_orientation=("strand", "order", "either")[seed % 3],
    )
    first = cluster_lossless_collinearity_anchors(anchors, params=params)
    shuffled = list(anchors)
    rng.shuffle(shuffled)

    assert cluster_lossless_collinearity_anchors(anchors, params=params) == first
    assert cluster_lossless_collinearity_anchors(shuffled, params=params) == first
    assert sorted(
        (anchor for block in first.blocks for anchor in block.anchors),
        key=collinearity_module._anchor_sort_key,
    ) == sorted(anchors, key=collinearity_module._anchor_sort_key)


def _rechained_from_scratch(
    anchors: list[CollinearityAnchor],
    params: LosslessCollinearityParameters,
) -> list[tuple[str, list[tuple[int, int]]]]:
    """The chaining rule restated without incremental updates: every round
    recomputes each orientation's best chains over all remaining anchors."""

    subjects = collinearity_module._pair_order_index(anchors)
    remaining = list(range(len(anchors)))
    chains = []
    while True:
        best = None
        for orientation_rank, orientation in enumerate(("plus", "minus")):
            nodes = sorted(
                (
                    index
                    for index in remaining
                    if collinearity_module._lossless_anchor_supports_orientation(
                        anchors[index], orientation, params
                    )
                ),
                key=lambda index: collinearity_module._path_sort_key(anchors[index], orientation),
            )
            value: dict[int, tuple[int, float, int | None]] = {}
            for rank, node in enumerate(nodes):
                count, score, previous = 1, float(anchors[node].bitscore), None
                for previous_node in nodes[:rank]:
                    if not collinearity_module._lossless_anchors_are_compatible(
                        anchors[previous_node], anchors[node], orientation=orientation, params=params
                    ) or collinearity_module._lossless_step_conflicts_exceed(
                        anchors[previous_node], anchors[node], subjects, int(params.max_conflicts)
                    ):
                        continue
                    candidate_count = value[previous_node][0] + 1
                    candidate_score = value[previous_node][1] + float(anchors[node].bitscore)
                    if (candidate_count, candidate_score) > (count, score):
                        count, score, previous = candidate_count, candidate_score, previous_node
                value[node] = (count, score, previous)
                key = (-count, -score, orientation_rank, rank)
                if count > 1 and (best is None or key < best[0]):
                    best = (key, orientation, node, value)
        if best is None:
            break
        _key, orientation, node, value = best
        taken = []
        while node is not None:
            taken.append(node)
            node = value[node][2]
        remaining = [index for index in remaining if index not in taken]
        chains.append(
            (orientation, sorted((anchors[i].query_order, anchors[i].subject_order) for i in taken))
        )
    return chains + [
        (
            collinearity_module._lossless_singleton_orientation(anchors[index]),
            [(anchors[index].query_order, anchors[index].subject_order)],
        )
        for index in remaining
    ]


@pytest.mark.linear
@pytest.mark.parametrize("seed", range(40))
def test_lossless_chaining_matches_the_rule_recomputed_from_scratch(seed: int) -> None:
    import random

    rng = random.Random(seed)
    anchors = _random_anchors(rng, 40)
    params = LosslessCollinearityParameters(
        max_unit_gap=seed % 4,
        max_diagonal_drift=seed % 3,
        max_conflicts=seed % 3,
        merge_orientation=("strand", "order", "either")[seed % 3],
    )

    result = cluster_lossless_collinearity_anchors(anchors, params=params)

    assert sorted(
        (block.orientation, sorted((a.query_order, a.subject_order) for a in block.anchors))
        for block in result.blocks
    ) == sorted(_rechained_from_scratch(anchors, params))


@pytest.mark.linear
def test_lossless_chaining_with_a_huge_unit_gap_visits_only_existing_units() -> None:
    # Review of OV-2294: the step window walked every unit index up to max_unit_gap;
    # it now visits only the query orders that hold an anchor (a huge gap stays fast).
    anchors = [_anchor(3 * index, 3 * index) for index in range(200)] + [
        _anchor(3 * index + 1, 3 * index + 5000) for index in range(200)
    ]

    def blocks(max_unit_gap: int):
        return _block_orders(
            cluster_lossless_collinearity_anchors(
                anchors,
                params=LosslessCollinearityParameters(
                    max_unit_gap=max_unit_gap, max_diagonal_drift=5
                ),
            )
        )

    assert blocks(10**9) == blocks(10)
    assert [len(block) for block in blocks(10)] == [200, 200]

@pytest.mark.linear
@pytest.mark.parametrize("seed", range(30))
def test_lossless_settings_without_effect_keep_their_contract(seed: int) -> None:
    import random

    rng = random.Random(1000 + seed)
    anchors = _random_anchors(rng, 25)

    def blocks(**kwargs):
        return _block_orders(
            cluster_lossless_collinearity_anchors(
                anchors, params=LosslessCollinearityParameters(**kwargs)
            )
        )

    # Diagonal drift and merge conflicts have no effect at max_unit_gap 0 (OV-2302).
    assert blocks(max_diagonal_drift=0) == blocks(max_diagonal_drift=7)
    assert blocks(max_conflicts=0) == blocks(max_conflicts=5)
    # "either" is the old spelling of "strand" (OV-2301).
    gap = seed % 4
    assert blocks(max_unit_gap=gap, merge_orientation="strand") == blocks(
        max_unit_gap=gap, merge_orientation="either"
    )


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
        anchor_mode="rbh",
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
        anchor_mode="rbh",
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
        for path in tuple(result.orthogroups.path_indexes_by_orthogroup_id["og_1"].iter_paths())
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
def test_anchor_core_record_local_membership_survives_when_hits_fit_limit() -> None:
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
def test_anchor_core_membership_survives_reciprocal_candidates_within_limit() -> None:
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
    paths = tuple(result.orthogroups.path_indexes_by_orthogroup_id["og_1"].iter_paths())
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
        edge_mode="all",
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
        edge_mode="all",
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
        edge_mode="all",
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
        edge_mode="all",
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
def test_orthogroup_collinearity_rbh_search_defaults_to_uncapped_candidates(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_limits: list[int | None] = []
    observed_pairs: list[tuple[str, str]] = []

    def fake_search(
        query_fasta,
        subject_fasta,
        *,
        ncbi_blastp_bin,
        max_hits,
        max_hsps_per_subject,
        **_kwargs,
    ):
        assert ncbi_blastp_bin is None
        assert max_hsps_per_subject is None
        observed_limits.append(max_hits)
        observed_pairs.append((_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta)))
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

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
def test_orthogroup_collinearity_member_hit_depth_is_post_search(monkeypatch) -> None:
    records = [
        _record("record_a", [_cds(0, 9, {"locus_tag": ["qa0"], "protein_id": ["qa0"]})]),
        _record("record_b", [_cds(0, 9, {"locus_tag": ["sb0"], "protein_id": ["sb0"]})]),
        _record("record_c", [_cds(0, 9, {"locus_tag": ["tc0"], "protein_id": ["tc0"]})]),
    ]
    observed_limits: list[int | None] = []

    def fake_search(
        query_fasta,
        subject_fasta,
        *,
        ncbi_blastp_bin,
        max_hits,
        max_hsps_per_subject,
        **_kwargs,
    ):
        assert ncbi_blastp_bin is None
        assert max_hsps_per_subject is None
        observed_limits.append(max_hits)
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

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

    def fake_search(
        query_fasta,
        subject_fasta,
        *,
        ncbi_blastp_bin,
        max_hits,
        max_hsps_per_subject,
        **_kwargs,
    ):
        assert ncbi_blastp_bin is None
        assert max_hsps_per_subject is None
        observed_limits.append(max_hits)
        observed_pairs.append((_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta)))
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

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

    def fake_search(
        query_fasta,
        subject_fasta,
        *,
        ncbi_blastp_bin,
        max_hsps_per_subject,
        **_kwargs,
    ):
        assert ncbi_blastp_bin is None
        assert max_hsps_per_subject is None
        observed_pairs.append((_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta)))
        return pd.DataFrame.from_records(
            [_hit_row(_first_fasta_id(query_fasta), _first_fasta_id(subject_fasta))],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

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
            annotation_warnings=(),
            feature_identity_notices=(),
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render_request)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--losat",
            "losatp",
            "--losatp_mode",
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
    assert options.losat_search is not None
    assert options.losat_search.losatp_member_max_hits is None
    assert options.collinear_max_paralog_links_per_orthogroup == 2
    assert options.losat_search.program == "losatp"
    assert options.losat_search.losatp_mode == "collinear"
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
        "losat_bin",
        "losat_threads",
        "losatp_mode",
        "collinear_min_anchors",
        "collinear_max_unit_gap",
        "collinear_color_mode",
        "keep_definition_left_aligned",
        "pairwise_match_style",
        "save_session",
        "session_output",
    ]
    hidden_option_names = [
        "losat-bin",
        "losat-threads",
        "losatp-mode",
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
        ["--losat-bin", "losatp"],
        ["--losat-threads", "2"],
        ["--losatp-mode", "pairwise"],
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
            "--losat",
            "losatp",
            "--losatp_mode",
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
def test_web_losatp_blastp_payload_helper_returns_collinear_rows(
    tmp_path: Path,
    stage_web_losatp_transport,
) -> None:
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
                "runtimeBindingHash": "record-a-cache-runtime",
                "displayBindingHash": "record-a-cache-display",
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": subject_map,
                "runtimeBindingHash": "record-b-cache-runtime",
                "displayBindingHash": "record-b-cache-display",
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

    pairs_path, raw_tsv_path = stage_web_losatp_transport(tmp_path, payload)
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
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
        "rbh",
        25,
        1,
        2,
        "adjacent",
        "rbh",
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    assert result["collinearityResult"]["schema"] == 3
    assert result["collinearityResult"]["kind"] == "result"
    assert result["collinearityResult"]["value"]["type"] == "CollinearityResult"
    typed_fields = result["collinearityResult"]["value"]["fields"]
    assert typed_fields["blocks"][0]["fields"]["blockId"] == "block_0001"
    assert len(typed_fields["orthogroups"]["fields"]["orthogroups"]) == 8
    rows = converted_pair_rows(result["pairs"][0])
    assert len(rows) == 1
    assert rows[0]["collinearity_block_id"] == "block_0001"
    assert rows[0]["collinearity_block_kind"] == "cluster"
    assert rows[0]["collinearity_anchor_count"] == 8
    assert rows[0]["collinearity_color_mode"] == "orientation"
    assert "collinearity_block_evalue" not in rows[0]  # OV-2300: the block E-value is removed
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
    assert set(result) == {"pairs", "collinearityResult", "identity", "identityKey", "provenance", "cache"}
    assert result["provenance"]["upstreamRawKeys"] == ["pair-a-b", "pair-b-a"]
    assert result["provenance"]["orthogroup"] == {
        "membershipMode": "anchor_core_v1",
        "memberMaxHits": None,
        "maxParalogLinksPerOrthogroup": 2,
    }
    assert result["provenance"]["collinear"]["unitMode"] == {
        "requested": "cds",
        "effectiveKinds": ["cds"],
    }
    assert result["provenance"]["collinear"]["anchorMode"] == "rbh"


@pytest.mark.linear
@pytest.mark.parametrize(
    ("anchor_mode", "merge_orientation", "expected_anchor_count"),
    (
        ("all", "strand", 2),
        ("all", "order", 2),
        # "either" is the old spelling of "strand" (OV-2301).
        ("rbh", "either", 1),
    ),
)
def test_web_losatp_blastp_payload_helper_uses_collinear_option_domains(
    tmp_path: Path,
    stage_web_losatp_transport,
    anchor_mode: str,
    merge_orientation: str,
    expected_anchor_count: int,
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
                "runtimeBindingHash": "record-a-cache-runtime",
                "displayBindingHash": "record-a-cache-display",
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": subject_map,
                "runtimeBindingHash": "record-b-cache-runtime",
                "displayBindingHash": "record-b-cache-display",
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

    pairs_path, raw_tsv_path = stage_web_losatp_transport(tmp_path, payload)
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
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
        anchor_mode,
        25,
        1,
        2,
        "adjacent",
        "anchor_core_v1",
        5,
        merge_orientation,
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    rows = converted_pair_rows(result["pairs"][0])
    assert len(rows) == 1
    assert rows[0]["collinearity_anchor_count"] == expected_anchor_count
    assert result["provenance"]["collinear"]["anchorMode"] == anchor_mode
    assert (
        result["provenance"]["collinear"]["parameters"]["mergeOrientation"]
        == {"either": "strand"}.get(merge_orientation, merge_orientation)
    )

    if anchor_mode == "rbh":
        raw_min_two = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
            str(pairs_path),
            str(raw_tsv_path),
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
        assert converted_pair_rows(min_two["pairs"][0]) == []


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_applies_collinear_search_scope(
    tmp_path: Path,
    stage_web_losatp_transport,
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
                "runtimeBindingHash": "record-a-cache-runtime",
                "displayBindingHash": "record-a-cache-display",
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": b_map,
                "runtimeBindingHash": "record-b-cache-runtime",
                "displayBindingHash": "record-b-cache-display",
            },
            {
                "recordIndex": 2,
                "recordId": "record_c",
                "proteinMap": c_map,
                "runtimeBindingHash": "record-c-cache-runtime",
                "displayBindingHash": "record-c-cache-display",
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

    pairs_path, raw_tsv_path = stage_web_losatp_transport(tmp_path, payload)

    def convert(
        scope: str,
        input_payload: Path = pairs_path,
        input_raw_tsv: Path = raw_tsv_path,
        explicit_display_pairs: bool = False,
    ) -> dict[str, object]:
        raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
            str(input_payload),
            str(input_raw_tsv),
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
            "all",
            25,
            1,
            2,
            scope,
            "rbh",
            None,
            "strand",
            True,
            None,
            explicit_display_pairs,
        )
        return json.loads(str(raw_result))

    adjacent = convert("adjacent")
    all_records = convert("all")
    invalid_scope = convert("global")
    invalid_mode = json.loads(
        str(
            namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
                str(pairs_path),
                str(raw_tsv_path),
                "similarity",
            )
        )
    )
    multi_row_payload = json.loads(json.dumps(payload))
    for pair in multi_row_payload["pairs"]:
        pair["displayPair"] = (
            pair["queryIndex"] == 0 and pair["subjectIndex"] == 2
        )
    multi_row_pairs_path, multi_row_raw_tsv_path = stage_web_losatp_transport(
        tmp_path,
        multi_row_payload,
        "losatp-multi-row-pairs",
    )
    multi_row_adjacent = convert(
        "adjacent",
        multi_row_pairs_path,
        multi_row_raw_tsv_path,
    )
    multi_row_all = convert(
        "all",
        multi_row_pairs_path,
        multi_row_raw_tsv_path,
        explicit_display_pairs=True,
    )

    assert "error" not in adjacent
    assert "error" not in all_records
    assert "error" not in multi_row_adjacent
    assert "error" not in multi_row_all
    assert invalid_scope["error"] == {
        "code": "INPUT_INVALID", "operation": "convertLosatpPairsToGenomicPayload",
        "stage": "helper", "context": {"field": "collinear_search_scope", "reason": "ADJACENT_ALL"},
    }
    assert invalid_mode["error"]["code"] == "INPUT_INVALID"
    assert invalid_mode["error"]["context"] == {"field": "protein_blastp_mode", "reason": "BLASTP_MODE"}

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
    assert set(adjacent) == {"pairs", "collinearityResult", "identity", "identityKey", "provenance", "cache"}
    assert set(all_records) == {"pairs", "collinearityResult", "identity", "identityKey", "provenance", "cache"}
    assert {"a0", "b0"} in adjacent_member_sets
    assert {"a0", "b0", "c0"} in all_member_sets
    assert {
        (
            block["fields"]["queryRecordIndex"],
            block["fields"]["subjectRecordIndex"],
        )
        for block in multi_row_adjacent["collinearityResult"]["value"]["fields"]["blocks"]
    } == {(0, 2)}
    assert {
        (block["fields"]["queryRecordIndex"], block["fields"]["subjectRecordIndex"])
        for block in multi_row_all["collinearityResult"]["value"]["fields"]["blocks"]
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
                # The rows name a group scope only when an anchor has a group (OV-2302).
                anchors=(replace(_anchor(0, 0), orthogroup_id="og_1"), _anchor(1, 1), _anchor(2, 2)),
            ),
            CollinearityBlock(
                block_id="block_minus",
                query_record_index=0,
                subject_record_index=1,
                orientation="minus",
                score=150.0,
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
    # OV-2300: a caller-supplied block E-value no longer reaches the rows.
    assert "collinearity_block_evalue" not in rows.columns
    assert set(rows["group_kind"]) == {"collinear_gene_group"}
    assert set(rows["group_scope"]) == {"adjacent_local"}
    assert set(rows["collinear_group_scope"]) == {"adjacent_local"}


@pytest.mark.linear
def test_collinearity_match_path_with_metadata_serializes() -> None:
    group = ribbon_group(
        [_record("record_a", []), _record("record_b", [])],
        collinearity_orientation_colors={"plus": "#112233", "minus": "#445566"},
    )

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
    assert 'data-gbdraw-match-id="comparison1_match1"' in svg_text
    assert "data-gbdraw-pairwise-match-id" not in svg_text
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
    assert "data-collinearity-block-evalue" not in svg_text  # OV-2300
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

    assert 'data-gbdraw-match-id="comparison1_match7"' in svg_text
    assert 'data-match-kind="pairwise"' in svg_text
    # Record identity comes from the endpoint records, never from table IDs (CO-06).
    assert 'data-query-record-id="record_a"' in svg_text
    assert 'data-subject-record-id="record_b"' in svg_text
    assert "query_record" not in svg_text
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
    )

    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        [],
        comparison=ComparisonLegendSource.of(blast_config),
        show_gc=False,
        show_skew=False,
        show_depth=False,
        pairwise_identity_legend=PairwiseIdentityLegend(hidden=True),
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
    )

    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        [],
        comparison=ComparisonLegendSource.of(blast_config),
        show_gc=False,
        show_skew=False,
        show_depth=False,
        pairwise_identity_legend=PairwiseIdentityLegend(label="Average identity"),
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
    )

    legend_table = prepare_legend_table(
        gc_config,
        skew_config,
        feature_config,
        [],
        comparison=ComparisonLegendSource.of(blast_config),
        show_gc=False,
        show_skew=False,
        show_depth=False,
        pairwise_identity_legend=PairwiseIdentityLegend(
            gradients=(
                PairwiseIdentityGradient("Collinear", "#eeeeee", "#112233"),
                PairwiseIdentityGradient("Inverted", "#ffeeee", "#445566"),
            ),
        ),
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
    legend_group = LegendGroup(legend_measurement, legend_table)
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
    label_width, _ = drawn_text_box(
        "Collinear",
        legend_measurement.font_family,
        legend_measurement.font_size,
    )
    assert bar_x == pytest.approx(
        label_width + 0.2 * legend_measurement.color_rect_size
    )


@pytest.mark.linear
def test_vertical_linear_pairwise_legend_aligns_to_feature_color_left_edge() -> None:
    feature_key = "tyrosine recombinase family protein"
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
    legend_group = LegendGroup(legend_measurement, legend_table)
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

    text_width, _ = drawn_text_box(
        feature_key,
        legend_measurement.font_family,
        legend_measurement.font_size,
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
    legend_group = LegendGroup(legend_measurement, legend_table)
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

    text_width, _ = drawn_text_box(
        feature_key,
        legend_measurement.font_family,
        legend_measurement.font_size,
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

    def fake_read_comparison_table(_blast_file):
        return frame

    # Linear assembly reads the table through read_comparison_resource, which parses with this reader.
    monkeypatch.setattr(comparisons_io_module, "read_comparison_table", fake_read_comparison_table)

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
    group = ribbon_group(
        [_record("record_a", []), _record("record_b", [])],
        collinearity_orientation_colors={"plus": "#112233", "minus": "#445566"},
    )

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
    group = ribbon_group([_record("record_a", []), _record("record_b", [])])

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
    match_group = ribbon_group(
        [_record("record_a", []), _record("record_b", [])], comparison_df, fill_opacity=1.0
    ).get_group()

    drawing = Drawing(debug=False)
    drawing.add(match_group)
    svg = drawing.tostring()

    assert svg.index('data-collinearity-block-id="block_plus"') < svg.index(
        'data-collinearity-block-id="block_minus"'
    )


@pytest.mark.parametrize("scope,pair_count", [("adjacent", 4), ("all", 6)])
@pytest.mark.parametrize("anchor_mode", ["rbh", "all"])
def test_collinear_without_orthogroups_skips_self_search_and_inference(
    monkeypatch, scope, pair_count, anchor_mode,
):
    records = [
        _record(name, [_cds(0, 9, {"protein_id": [name]})])
        for name in ["a", "b", "c"]
    ]
    observed = []

    def search(query, subject, **kwargs):
        pair = (_first_fasta_id(query), _first_fasta_id(subject))
        assert pair[0] != pair[1]
        assert kwargs["max_hits"] == 5
        observed.append(pair)
        return pd.DataFrame([_hit_row(*pair)], columns=COMPARISON_COLUMNS)

    def forbidden_inference(*args, **kwargs):
        pytest.fail("Collinear OFF must not infer orthogroups")

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", search)
    monkeypatch.setattr(collinearity_module, "select_rbh_orthogroup_edges_from_directional_hits", forbidden_inference)
    result = collinearity_module.build_orthogroup_collinearity_blocks(
        records, infer_orthogroups=False, candidate_limit=5,
        orthogroup_member_max_hits=5, search_scope=scope, edge_mode=anchor_mode,
    )
    assert len(observed) == pair_count
    assert result.orthogroups is None
    assert len(result.blocks) == 2
    assert all(anchor.source == "protein_hit" and not anchor.orthogroup_id
               for block in result.blocks for anchor in block.anchors)


def test_collinear_without_inference_applies_member_limit_after_raw_filters():
    records = [_record(name, [
        _cds(index * 15, index * 15 + 9, {"protein_id": [f"{name}{index}"]})
        for index in range(2)
    ]) for name in ["a", "b"]]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame([
        _hit_row("a0", "b0", bitscore=200), _hit_row("a0", "b1", bitscore=100),
    ], columns=COMPARISON_COLUMNS)
    results = [collinearity_module.build_orthogroup_collinearity_blocks_from_hits(
        {(0, 1): hits}, extraction, records=records, infer_orthogroups=False,
        edge_mode="all", orthogroup_member_max_hits=limit,
    ) for limit in [1, None]]
    assert [sum(len(block.anchors) for block in result.blocks) for result in results] == [1, 2]


def _left_to_right_sum(values, /, start=0):
    """``sum()`` as Python 3.10 and 3.11 add floats (3.12 compensates them)."""

    total = start
    for value in values:
        total = total + value
    return total


def _protein_float_reductions():
    anchors = tuple(
        replace(_anchor(index, index, bitscore=0.1), identity=0.1, alignment_length=1)
        for index in range(10)
    )
    block = CollinearityBlock(
        block_id="block",
        query_record_index=0,
        subject_record_index=1,
        orientation="plus",
        score=collinearity_module._anchor_bitscore_total(anchors),
        anchors=anchors,
    )
    (row,) = convert_collinearity_blocks_to_pair_comparisons(
        CollinearityResult(blocks=(block,)), record_ids=["a", "b"]
    )[(0, 1)][["bitscore", "identity"]].itertuples(index=False)
    size = 12
    hits = pd.DataFrame(
        {
            "query": [f"q{index % 3}" for index in range(size)],
            "subject": [f"s{index % 2}" for index in range(size)],
            "length_product": [float(100 + (index % 7) * 17) for index in range(size)],
            "bitscore": [float(10 + (index % 4) * 3) for index in range(size)],
        }
    )
    return (
        block.score,
        collinearity_module._weighted_identity(anchors),
        tuple(row),
        protein_colinearity_module._fit_expected_bitscore_model(hits, top_fraction=1.0),
    )


def test_protein_float_reductions_do_not_depend_on_the_python_version(monkeypatch) -> None:
    """Stored Gallery protein results equal a rebuild on every Python (OV-1237).

    Python 3.12 changed ``sum()`` of floats to compensated summation, so a
    result built with ``sum()`` on 3.10 or 3.11 differed in its last bits from
    the 3.13 Sessions and the Web's Pyodide.
    """

    current = _protein_float_reductions()
    with monkeypatch.context() as patched:
        patched.setattr("builtins.sum", _left_to_right_sum)
        left_to_right = _protein_float_reductions()

    assert current == left_to_right
    score, identity, (bitscore, row_identity), model = current
    assert score == bitscore == 1.0
    assert identity == row_identity == 0.1
    # CPython 3.12+ sum() values; left-to-right addition gives (0.2345015637431571, 0.645579356667975).
    assert model == (0.23450156374315714, 0.6455793566679746)


@pytest.mark.linear
def test_a_comparison_interpolates_each_gradient_color_once(monkeypatch: pytest.MonkeyPatch) -> None:
    """The match colors of one comparison are interpolated once per distinct
    identity; every match still gets its own color."""
    import gbdraw.render.groups.linear.pairwise_match as pairwise_match

    calls: list[float] = []
    real = pairwise_match.interpolate_default_colors

    def counting(min_color: tuple[str, str], max_color: tuple[str, str], factor: float) -> str:
        calls.append(factor)
        return real(min_color, max_color, factor)

    monkeypatch.setattr(pairwise_match, "interpolate_default_colors", counting)
    identities = [90.0, 95.0, 90.0, 100.0, 95.0, 90.0]
    comparison = pd.DataFrame(
        {
            "identity": identities,
            "qstart": [10 * i + 1 for i in range(6)],
            "qend": [10 * i + 9 for i in range(6)],
            "sstart": [10 * i + 1 for i in range(6)],
            "send": [10 * i + 9 for i in range(6)],
        }
    )

    group = ribbon_group([_record("record_a", []), _record("record_b", [])], comparison)

    paths = group.get_group().elements
    assert len(paths) == 6 and len(calls) == 3
    fills = {path.attribs["data-identity"]: path.attribs["fill"] for path in paths}
    assert len(set(fills.values())) == 3
