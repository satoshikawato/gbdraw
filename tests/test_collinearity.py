from __future__ import annotations

import json
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityParameters,
    CollinearityResult,
    build_collinearity_blocks_from_hits,
    call_collinearity_blocks,
    convert_collinearity_blocks_to_comparisons,
    deduplicate_unit_pair_anchors,
    normalize_collinearity_color_mode,
)
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.protein_colinearity import extract_cds_proteins
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators.blast import BlastMatchConfigurator
from gbdraw.io.collinearity import parse_native_collinearity_tsv, write_native_collinearity_tsv
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.legend.table import prepare_legend_table
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup


def _record(record_id: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(Seq("ATGAAATAG" * 80), id=record_id)
    record.features = features
    return record


def _cds(start: int, end: int, qualifiers: dict[str, list[str]]) -> SeqFeature:
    merged = {"translation": ["MKK"]}
    merged.update(qualifiers)
    return SeqFeature(FeatureLocation(start, end, strand=1), type="CDS", qualifiers=merged)


def _hit_row(query: str, subject: str, bitscore: float = 200.0) -> dict[str, object]:
    return {
        "query": query,
        "subject": subject,
        "identity": 90.0,
        "alignment_length": 100,
        "mismatches": 0,
        "gap_opens": 0,
        "qstart": 1,
        "qend": 100,
        "sstart": 1,
        "send": 100,
        "evalue": 1e-30,
        "bitscore": bitscore,
    }


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


def _anchor(query_order: int, subject_order: int, *, bitscore: float = 200.0) -> CollinearityAnchor:
    return CollinearityAnchor(
        query_protein_id=f"q{query_order}",
        subject_protein_id=f"s{subject_order}",
        query_record_index=0,
        subject_record_index=1,
        query_order=query_order,
        subject_order=subject_order,
        query_start=query_order * 10 + 1,
        query_end=query_order * 10 + 9,
        subject_start=subject_order * 10 + 1,
        subject_end=subject_order * 10 + 9,
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
def test_native_block_caller_detects_plus_and_minus_orientation() -> None:
    params = CollinearityParameters(min_anchors=3, max_gene_gap=1)

    plus = call_collinearity_blocks([_anchor(0, 0), _anchor(1, 1), _anchor(2, 2)], params=params)
    minus = call_collinearity_blocks([_anchor(0, 2), _anchor(1, 1), _anchor(2, 0)], params=params)

    assert len(plus.blocks) == 1
    assert plus.blocks[0].orientation == "plus"
    assert [anchor.query_order for anchor in plus.blocks[0].anchors] == [0, 1, 2]
    assert len(minus.blocks) == 1
    assert minus.blocks[0].orientation == "minus"


@pytest.mark.linear
def test_collinearity_color_mode_defaults_to_orientation_and_aliases_identity() -> None:
    assert normalize_collinearity_color_mode(None) == "orientation"
    assert normalize_collinearity_color_mode("identity") == "average_identity"


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

    result = build_collinearity_blocks_from_hits(
        [hits],
        extraction,
        records=records,
        params=CollinearityParameters(min_anchors=3, max_gene_gap=1),
    )
    comparisons = convert_collinearity_blocks_to_comparisons(result, records=records, color_mode="orientation")

    assert result.blocks[0].block_id == "block_0001"
    assert comparisons[0].shape[0] == 1
    assert comparisons[0].iloc[0]["collinearity_block_id"] == "block_0001"
    assert comparisons[0].iloc[0]["collinearity_anchor_count"] == 3
    assert comparisons[0].iloc[0]["collinearity_color_mode"] == "orientation"
    assert comparisons[0].iloc[0]["query_locus_id"] == "qa0;qa1;qa2"


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
            "block_id\tanchor_index\tquery_record\tquery_unit\tsubject_record\tsubject_unit\torientation\tidentity\tevalue\tbitscore\talignment_length\tscore",
            "block_a\t1\t#1\tqa0\t#2\tsb0\tplus\t91\t1e-20\t200\t100\t50",
            "block_a\t2\trecord_a\tqa1\trecord_b\tsb1\tplus\t92\t1e-25\t210\t100\t50",
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
    assert [anchor.query_display_name for anchor in result.blocks[0].anchors] == ["qa0", "qa1"]
    assert reparsed.blocks[0].anchors[1].subject_locus_id == "sb1"


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

    monkeypatch.setattr(linear_cli_module, "build_native_collinearity_blocks", fake_build)
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
    assert isinstance(params, CollinearityParameters)
    assert params.min_anchors == 1
    assert captured["unit_mode"] == "cds"
    assert captured["protein_blastp_mode"] == "none"
    comparisons = captured["protein_comparisons"]
    assert isinstance(comparisons, list)
    assert comparisons[0].iloc[0]["collinearity_color_mode"] == "orientation"
    assert "block_0001" in save_path.read_text(encoding="utf-8")


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
        [_hit_row(f"qa{index}", f"sb{index}") for index in range(3)],
        columns=COMPARISON_COLUMNS,
    )
    payload = [
        {
            "pairIndex": 0,
            "queryIndex": 0,
            "subjectIndex": 1,
            "blastText": hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            "queryProteinMap": {
                f"qa{index}": _web_protein_entry(
                    f"qa{index}",
                    record_index=0,
                    record_id="record_a",
                    feature_index=index,
                    start=index * 12,
                    end=index * 12 + 9,
                )
                for index in range(3)
            },
            "subjectProteinMap": {
                f"sb{index}": _web_protein_entry(
                    f"sb{index}",
                    record_index=1,
                    record_id="record_b",
                    feature_index=index,
                    start=index * 12,
                    end=index * 12 + 9,
                )
                for index in range(3)
            },
        }
    ]

    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        json.dumps(payload),
        "collinear",
        5,
        50,
        "1e-5",
        0,
        0,
        3,
        JsNull(),
        JsNull(),
        JsNull(),
        JsNull(),
        JsNull(),
        JsNull(),
        "cds",
        "orientation",
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    assert result["collinearityBlocks"][0]["id"] == "block_0001"
    rows = result["pairs"][0]["rows"]
    assert len(rows) == 1
    assert rows[0]["collinearity_block_id"] == "block_0001"
    assert rows[0]["collinearity_anchor_count"] == 3
    assert rows[0]["collinearity_color_mode"] == "orientation"


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
                anchors=(_anchor(0, 0), _anchor(1, 1), _anchor(2, 2)),
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
        collinearity_orientation="plus",
        collinearity_color_mode="orientation",
        query_protein_id="qa0",
        subject_protein_id="sb0",
        query_feature_svg_id="feature_qa0",
        subject_feature_svg_id="feature_sb0",
        query_unit_id="qa0",
        subject_unit_id="sb0",
    )

    drawing = Drawing(debug=False)
    drawing.add(group.generate_linear_match_path(row))

    assert "data-collinearity-block-id" in drawing.tostring()
    assert 'data-collinearity-color-mode="orientation"' in drawing.tostring()
    assert 'fill="#112233"' in drawing.tostring()


@pytest.mark.linear
def test_blast_configurator_reads_collinearity_colors_from_default_colors() -> None:
    default_colors = pd.DataFrame(
        [
            ("pairwise_match_min", "#ffffff"),
            ("pairwise_match_max", "#000000"),
            ("pairwise_match", "#cccccc"),
            ("collinear_block_minus", "#445566"),
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

    assert config.collinearity_orientation_colors == {"plus": "#d3d3d3", "minus": "#445566"}


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
