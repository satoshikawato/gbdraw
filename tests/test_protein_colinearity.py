from __future__ import annotations

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.diagram as api_diagram_module
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.protein_colinearity import (
    build_orthogroups_from_protein_hits,
    build_protein_colinearity_comparisons,
    cap_hits_per_query,
    convert_pair_protein_hits_to_genomic_links,
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
    parse_losatp_outfmt6,
    proteins_to_fasta,
)
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import DiagramOptions
from gbdraw.diagrams.linear.orthogroup_alignment import calculate_orthogroup_alignment_offsets
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup


def _record(
    record_id: str,
    sequence: str = "ATGAAATAG" * 20,
    features: list[SeqFeature] | None = None,
) -> SeqRecord:
    record = SeqRecord(Seq(sequence), id=record_id)
    record.features = list(features or [])
    return record


def _cds(
    start: int,
    end: int,
    *,
    strand: int = 1,
    qualifiers: dict[str, list[str]] | None = None,
) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers=qualifiers or {"translation": ["MK*"]},
    )


def _hit_row(
    query: str,
    subject: str,
    *,
    identity: float = 90.0,
    alignment_length: int = 100,
    evalue: float = 1e-30,
    bitscore: float = 200.0,
) -> dict[str, object]:
    return {
        "query": query,
        "subject": subject,
        "identity": identity,
        "alignment_length": alignment_length,
        "mismatches": 0,
        "gap_opens": 0,
        "qstart": 1,
        "qend": 100,
        "sstart": 1,
        "send": 100,
        "evalue": evalue,
        "bitscore": bitscore,
    }


@pytest.mark.linear
def test_extract_cds_proteins_uses_translation_and_stable_synthetic_ids() -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MKT*"],
                        "protein_id": ["duplicate"],
                        "locus_tag": ["gene_a"],
                    },
                )
            ],
        ),
        _record(
            "record_b",
            features=[
                _cds(
                    9,
                    18,
                    qualifiers={
                        "translation": ["MGG*"],
                        "protein_id": ["duplicate"],
                    },
                )
            ],
        ),
    ]

    result = extract_cds_proteins(records)

    assert [protein.protein_id for protein in result.protein_map.values()] == [
        "gbd_r0001_cds000001",
        "gbd_r0002_cds000001",
    ]
    first = result.proteins_by_record[0][0]
    assert first.sequence == "MKT"
    assert first.label == "gene_a"
    assert first.start == 0
    assert first.end == 9


@pytest.mark.linear
def test_extract_cds_proteins_prefers_unique_source_protein_id() -> None:
    record = _record(
        "record_a",
        features=[
            _cds(
                0,
                9,
                qualifiers={
                    "translation": ["MKT*"],
                    "protein_id": ["WP_123456789.1"],
                },
            )
        ],
    )

    result = extract_cds_proteins([record])

    protein = result.proteins_by_record[0][0]
    assert protein.protein_id == "WP_123456789.1"
    assert protein.source_protein_id == "WP_123456789.1"
    assert ">WP_123456789.1" in proteins_to_fasta([protein])


@pytest.mark.linear
def test_extract_cds_proteins_accepts_record_index_offset() -> None:
    record = _record("record_b", features=[_cds(9, 18)])

    result = extract_cds_proteins([record], record_index_offset=1)

    protein = result.proteins_by_record[0][0]
    assert protein.protein_id == "gbd_r0002_cds000001"
    assert protein.record_index == 1


@pytest.mark.linear
def test_extract_cds_proteins_translates_when_translation_missing() -> None:
    record = _record(
        "record_a",
        sequence="ATGAAATAG",
        features=[_cds(0, 9, qualifiers={"locus_tag": ["fallback"]})],
    )

    result = extract_cds_proteins([record])

    assert result.proteins_by_record[0][0].sequence == "MK"


@pytest.mark.linear
def test_extract_cds_proteins_handles_compound_location_span() -> None:
    feature = SeqFeature(
        CompoundLocation(
            [
                FeatureLocation(0, 6, strand=1),
                FeatureLocation(12, 18, strand=1),
            ]
        ),
        type="CDS",
        qualifiers={"translation": ["MKM"]},
    )
    record = _record("record_a", features=[feature])

    result = extract_cds_proteins([record])

    protein = result.proteins_by_record[0][0]
    assert protein.start == 0
    assert protein.end == 18
    assert protein.strand == 1


@pytest.mark.linear
def test_cap_hits_per_query_keeps_top_five_distinct_subjects() -> None:
    rows = [
        _hit_row("q1", "s1", bitscore=100),
        _hit_row("q1", "s1", bitscore=90),
        _hit_row("q1", "s2", bitscore=80),
        _hit_row("q1", "s3", bitscore=70),
        _hit_row("q1", "s4", bitscore=60),
        _hit_row("q1", "s5", bitscore=50),
        _hit_row("q1", "s6", bitscore=40),
        _hit_row("q2", "s7", bitscore=10),
    ]
    hits = pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)

    capped = cap_hits_per_query(hits, max_hits=5)

    assert capped[capped["query"] == "q1"]["subject"].tolist() == [
        "s1",
        "s2",
        "s3",
        "s4",
        "s5",
    ]
    assert capped[capped["query"] == "q2"]["subject"].tolist() == ["s7"]


@pytest.mark.linear
def test_convert_protein_hits_to_genomic_links_uses_cds_spans_and_strand() -> None:
    records = [
        _record("query_record", features=[_cds(0, 9, strand=1)]),
        _record("subject_record", features=[_cds(20, 50, strand=-1)]),
    ]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame.from_records(
        [
            _hit_row(
                "gbd_r0001_cds000001",
                "gbd_r0002_cds000001",
                alignment_length=20,
            )
        ],
        columns=COMPARISON_COLUMNS,
    )

    converted = convert_protein_hits_to_genomic_links(hits, extraction.protein_map)

    row = converted.iloc[0]
    assert row["query"] == "query_record"
    assert row["subject"] == "subject_record"
    assert row["qstart"] == 1
    assert row["qend"] == 9
    assert row["sstart"] == 50
    assert row["send"] == 21
    assert row["alignment_length"] == 20


@pytest.mark.linear
def test_separately_extracted_offset_protein_maps_convert_without_collision() -> None:
    query_record = _record("query_record", features=[_cds(0, 9, strand=1)])
    subject_record = _record("subject_record", features=[_cds(50, 80, strand=1)])
    query_extraction = extract_cds_proteins([query_record], record_index_offset=0)
    subject_extraction = extract_cds_proteins([subject_record], record_index_offset=1)
    protein_map = {
        **query_extraction.protein_map,
        **subject_extraction.protein_map,
    }
    hits = pd.DataFrame.from_records(
        [
            _hit_row(
                "gbd_r0001_cds000001",
                "gbd_r0002_cds000001",
                alignment_length=20,
            )
        ],
        columns=COMPARISON_COLUMNS,
    )

    converted = convert_protein_hits_to_genomic_links(hits, protein_map)

    row = converted.iloc[0]
    assert row["query"] == "query_record"
    assert row["subject"] == "subject_record"
    assert row["qstart"] == 1
    assert row["qend"] == 9
    assert row["sstart"] == 51
    assert row["send"] == 80


@pytest.mark.linear
def test_pair_conversion_disambiguates_same_source_protein_id_between_records() -> None:
    query_record = _record(
        "query_record",
        features=[_cds(0, 9, strand=1, qualifiers={"translation": ["MK*"], "protein_id": ["WP_same.1"]})],
    )
    subject_record = _record(
        "subject_record",
        features=[_cds(50, 80, strand=1, qualifiers={"translation": ["MK*"], "protein_id": ["WP_same.1"]})],
    )
    query_extraction = extract_cds_proteins([query_record], record_index_offset=0)
    subject_extraction = extract_cds_proteins([subject_record], record_index_offset=1)
    hits = pd.DataFrame.from_records(
        [_hit_row("WP_same.1", "WP_same.1", alignment_length=20)],
        columns=COMPARISON_COLUMNS,
    )

    converted = convert_pair_protein_hits_to_genomic_links(
        hits,
        query_extraction.protein_map,
        subject_extraction.protein_map,
    )

    row = converted.iloc[0]
    assert row["query"] == "query_record"
    assert row["subject"] == "subject_record"
    assert row["qstart"] == 1
    assert row["qend"] == 9
    assert row["sstart"] == 51
    assert row["send"] == 80


@pytest.mark.linear
def test_parse_losatp_outfmt6_returns_standard_columns() -> None:
    parsed = parse_losatp_outfmt6(
        "# comment\nq1\ts1\t91.5\t42\t1\t0\t1\t42\t2\t43\t1e-20\t150\n"
    )

    assert parsed.columns.tolist() == list(COMPARISON_COLUMNS)
    assert parsed.iloc[0]["identity"] == pytest.approx(91.5)
    assert parsed.iloc[0]["bitscore"] == pytest.approx(150)


@pytest.mark.linear
def test_build_protein_colinearity_comparisons_accepts_test_runner() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]

    def fake_runner(query_fasta: str, subject_fasta: str) -> pd.DataFrame:
        assert ">gbd_r0001_cds000001" in query_fasta
        assert ">gbd_r0002_cds000001" in subject_fasta
        return pd.DataFrame.from_records(
            [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
            columns=COMPARISON_COLUMNS,
        )

    comparisons = build_protein_colinearity_comparisons(
        records,
        runner=fake_runner,
    )

    assert len(comparisons) == 1
    assert comparisons[0].iloc[0]["query"] == "record_a"
    assert comparisons[0].iloc[0]["subject"] == "record_b"
    assert comparisons[0].iloc[0]["orthogroup_id"] == "og_1"
    assert comparisons[0].iloc[0]["query_protein_id"] == "gbd_r0001_cds000001"
    assert comparisons[0].iloc[0]["subject_protein_id"] == "gbd_r0002_cds000001"


@pytest.mark.linear
def test_build_orthogroups_selects_record_representatives_with_paralogs() -> None:
    records = [
        _record("record_a", features=[_cds(0, 30)]),
        _record("record_b", features=[_cds(100, 130), _cds(200, 230)]),
        _record("record_c", features=[_cds(400, 430)]),
    ]
    extraction = extract_cds_proteins(records)
    hits_ab = pd.DataFrame.from_records(
        [
            _hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001", bitscore=120, evalue=1e-20),
            _hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000002", bitscore=250, evalue=1e-30),
        ],
        columns=COMPARISON_COLUMNS,
    )
    hits_bc = pd.DataFrame.from_records(
        [_hit_row("gbd_r0002_cds000001", "gbd_r0003_cds000001", bitscore=110, evalue=1e-10)],
        columns=COMPARISON_COLUMNS,
    )

    orthogroups = build_orthogroups_from_protein_hits(
        [hits_ab, hits_bc],
        extraction.protein_map,
    )

    members = orthogroups.orthogroups["og_1"]
    assert {member.protein_id for member in members} == {
        "gbd_r0001_cds000001",
        "gbd_r0002_cds000001",
        "gbd_r0002_cds000002",
        "gbd_r0003_cds000001",
    }
    record_b_reps = [
        member.protein_id
        for member in members
        if member.record_index == 1 and member.representative
    ]
    assert record_b_reps == ["gbd_r0002_cds000002"]


@pytest.mark.linear
def test_orthogroup_alignment_offsets_align_selected_member_to_representatives() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    comparison = pd.DataFrame.from_records(
        [
            {
                **_hit_row("record_a", "record_b", bitscore=200),
                "qstart": 100,
                "qend": 200,
                "sstart": 400,
                "send": 500,
                "query_protein_id": "prot_a",
                "subject_protein_id": "prot_b",
                "query_source_protein_id": "",
                "subject_source_protein_id": "",
                "query_record_index": 0,
                "subject_record_index": 1,
                "query_feature_index": 0,
                "subject_feature_index": 0,
                "query_feature_svg_id": "fanchor",
                "subject_feature_svg_id": "fsubject",
                "orthogroup_id": "og_1",
                "query_orthogroup_representative": True,
                "subject_orthogroup_representative": True,
            }
        ]
    )
    canvas_config = type(
        "CanvasConfig",
        (),
        {
            "normalize_length": False,
            "align_center": False,
            "alignment_width": 1000.0,
            "longest_genome": 1000,
        },
    )()

    offsets = calculate_orthogroup_alignment_offsets(
        records,
        [comparison],
        canvas_config,
        "fanchor",
    )

    assert offsets[0] == pytest.approx(0.0)
    assert offsets[1] == pytest.approx(-300.0)


@pytest.mark.linear
def test_pairwise_match_group_applies_record_specific_alignment_offsets() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    comparison = pd.DataFrame.from_records(
        [
            {
                **_hit_row("record_a", "record_b", identity=90),
                "qstart": 100,
                "qend": 200,
                "sstart": 400,
                "send": 500,
            }
        ]
    )
    canvas_config = type(
        "CanvasConfig",
        (),
        {
            "normalize_length": False,
            "align_center": False,
            "longest_genome": 1000,
            "alignment_width": 1000.0,
        },
    )()
    blast_config = type(
        "BlastConfig",
        (),
        {
            "fill_color": "#cccccc",
            "identity": 0,
            "min_color": "#eeeeee",
            "max_color": "#111111",
            "fill_opacity": 0.5,
            "stroke_color": "#000000",
            "stroke_width": 0.1,
        },
    )()

    group = PairWiseMatchGroup(
        canvas_config,
        {"record_a": 1000, "record_b": 1000},
        comparison,
        100.0,
        1,
        blast_config,
        records,
        record_offsets_x={1: -300.0},
    ).get_group()

    path = group.elements[0]
    assert 'd="M 100.0,0L200.0,0 L200.0,100.0L100.0,100.0 z"' in path.tostring()


@pytest.mark.linear
def test_assemble_linear_diagram_accepts_precomputed_protein_comparisons(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured["comparison_dataframes"] = kwargs.get("comparison_dataframes")
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram", fake_assemble)

    records = [_record("record_a"), _record("record_b")]
    comparison = pd.DataFrame.from_records(
        [_hit_row("record_a", "record_b")],
        columns=COMPARISON_COLUMNS,
    )
    canvas = assemble_linear_diagram_from_records(
        records,
        protein_comparisons=[comparison],
    )

    assert isinstance(canvas, Drawing)
    frames = captured["comparison_dataframes"]
    assert isinstance(frames, list)
    assert len(frames) == 1
    assert frames[0].equals(comparison)


@pytest.mark.linear
def test_build_linear_diagram_forwards_protein_colinearity_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_record("record_a"), _record("record_b")],
        options=DiagramOptions(
            protein_colinearity=True,
            losatp_bin="custom-losat",
            losatp_max_hits=7,
        ),
    )

    assert isinstance(canvas, Drawing)
    assert captured["protein_colinearity"] is True
    assert captured["losatp_bin"] == "custom-losat"
    assert captured["losatp_max_hits"] == 7
    assert captured["align_orthogroup_feature"] is None


@pytest.mark.linear
def test_build_linear_diagram_forwards_orthogroup_alignment_option(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_record("record_a"), _record("record_b")],
        options=DiagramOptions(align_orthogroup_feature="fanchor"),
    )

    assert isinstance(canvas, Drawing)
    assert captured["align_orthogroup_feature"] == "fanchor"


@pytest.mark.linear
def test_linear_cli_rejects_blast_with_protein_colinearity() -> None:
    with pytest.raises(SystemExit):
        linear_cli_module.linear_main(
            [
                "--gbk",
                "a.gb",
                "b.gb",
                "-b",
                "a_b.tsv",
                "--protein_colinearity",
            ]
        )


@pytest.mark.linear
def test_linear_cli_requires_two_records_for_protein_colinearity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [_record("only")])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)

    with pytest.raises(ValidationError, match="requires at least two"):
        linear_cli_module.linear_main(
            [
                "--gbk",
                "dummy.gb",
                "--protein_colinearity",
                "--format",
                "svg",
            ]
        )


@pytest.mark.linear
def test_linear_cli_forwards_protein_colinearity_options(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a"), _record("record_b")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_colinearity",
            "--losatp_bin",
            "custom-losat",
            "--losatp_max_hits",
            "9",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["protein_colinearity"] is True
    assert captured["losatp_bin"] == "custom-losat"
    assert captured["losatp_max_hits"] == 9
    assert captured["align_orthogroup_feature"] is None


@pytest.mark.linear
def test_linear_cli_forwards_orthogroup_alignment_option(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a"), _record("record_b")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--align_orthogroup_feature",
            "fanchor",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["align_orthogroup_feature"] == "fanchor"
