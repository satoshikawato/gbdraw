from __future__ import annotations

import json
from pathlib import Path

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
    build_pairwise_protein_blastp_comparisons,
    build_rbh_orthogroup_protein_blastp_comparisons,
    cap_hits_per_query,
    convert_pair_protein_hits_to_genomic_links,
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
    filter_protein_hits_by_thresholds,
    parse_losatp_outfmt6,
    proteins_to_fasta,
    select_best_hits_per_query,
    select_reciprocal_best_hit_edges,
    select_reciprocal_best_hits,
    select_top_hits_per_query,
)
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import DiagramOptions
from gbdraw.diagrams.linear.orthogroup_alignment import (
    calculate_orthogroup_alignment_canvas_adjustment,
    calculate_orthogroup_alignment_offsets,
)
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


def _web_protein_entry(
    protein_id: str,
    *,
    record_index: int,
    record_id: str,
    feature_index: int = 0,
    start: int = 0,
    end: int = 90,
    gene: str | None = None,
    product: str | None = None,
    note: str | None = None,
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
        "gene": gene,
        "product": product,
        "note": note,
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
def test_extract_cds_proteins_carries_annotation_fields() -> None:
    record = _record(
        "record_a",
        features=[
            _cds(
                0,
                9,
                qualifiers={
                    "translation": ["MKT*"],
                    "gene": ["rpoB"],
                    "product": ["DNA-directed RNA polymerase beta subunit"],
                    "note": ["core polymerase subunit"],
                },
            )
        ],
    )

    result = extract_cds_proteins([record])

    protein = result.proteins_by_record[0][0]
    assert protein.gene == "rpoB"
    assert protein.product == "DNA-directed RNA polymerase beta subunit"
    assert protein.note == "core polymerase subunit"


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
def test_filter_protein_hits_by_thresholds_removes_low_confidence_bridge() -> None:
    raw_hits = pd.DataFrame.from_records(
        [
            _hit_row("a", "b"),
            _hit_row("b", "c", bitscore=20, evalue=1e-1),
        ],
        columns=COMPARISON_COLUMNS,
    )
    filtered_hits = filter_protein_hits_by_thresholds(
        raw_hits,
        bitscore=50,
        evalue=1e-5,
        identity=0,
        alignment_length=0,
    )
    assert filtered_hits[["query", "subject"]].to_records(index=False).tolist() == [
        ("a", "b")
    ]

    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
        _record("record_c", features=[_cds(18, 27)]),
    ]
    calls = 0

    def fake_runner(_query_fasta: str, _subject_fasta: str) -> pd.DataFrame:
        nonlocal calls
        calls += 1
        if calls == 1:
            rows = [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")]
        else:
            rows = [
                _hit_row(
                    "gbd_r0002_cds000001",
                    "gbd_r0003_cds000001",
                    bitscore=20,
                    evalue=1e-1,
                )
            ]
        return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)

    result = build_pairwise_protein_blastp_comparisons(
        records,
        runner=fake_runner,
        bitscore=50,
        evalue=1e-5,
        identity=0,
        alignment_length=0,
    )

    assert result.orthogroups is None
    assert result.comparisons[0].iloc[0]["orthogroup_id"] == ""
    assert result.comparisons[1].empty


@pytest.mark.linear
def test_select_best_hits_per_query_avoids_secondary_paralog_merge() -> None:
    hits = pd.DataFrame.from_records(
        [
            _hit_row("a", "b1", bitscore=300),
            _hit_row("a", "b2", bitscore=120),
        ],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_best_hits_per_query(hits)

    assert selected[["query", "subject"]].to_records(index=False).tolist() == [
        ("a", "b1")
    ]


@pytest.mark.linear
def test_select_reciprocal_best_hits_rejects_many_to_one_nonreciprocal_hit() -> None:
    hits = pd.DataFrame.from_records(
        [
            _hit_row("a1", "b", bitscore=300),
            _hit_row("a2", "b", bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_reciprocal_best_hits(hits)

    assert selected[["query", "subject"]].to_records(index=False).tolist() == [
        ("a1", "b")
    ]


@pytest.mark.linear
def test_select_reciprocal_best_hit_edges_requires_directional_reciprocity() -> None:
    forward = pd.DataFrame.from_records(
        [
            _hit_row("a1", "b", bitscore=300),
            _hit_row("a2", "b", bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )
    reverse = pd.DataFrame.from_records(
        [_hit_row("b", "a3", bitscore=400)],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_reciprocal_best_hit_edges(forward, reverse)

    assert selected.empty


@pytest.mark.linear
def test_select_top_hits_per_query_is_deterministic_on_subject_tie() -> None:
    hits = pd.DataFrame.from_records(
        [
            _hit_row("a", "b2", bitscore=300),
            _hit_row("a", "b1", bitscore=300),
        ],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_top_hits_per_query(hits, max_hits=1)

    assert selected.iloc[0]["subject"] == "b1"


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
def test_build_pairwise_protein_blastp_comparisons_accepts_test_runner_without_orthogroups() -> None:
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

    result = build_pairwise_protein_blastp_comparisons(
        records,
        runner=fake_runner,
    )

    comparisons = result.comparisons
    assert result.orthogroups is None
    assert len(comparisons) == 1
    assert comparisons[0].iloc[0]["query"] == "record_a"
    assert comparisons[0].iloc[0]["subject"] == "record_b"
    assert comparisons[0].iloc[0]["orthogroup_id"] == ""
    assert comparisons[0].iloc[0]["query_protein_id"] == "gbd_r0001_cds000001"
    assert comparisons[0].iloc[0]["subject_protein_id"] == "gbd_r0002_cds000001"


@pytest.mark.linear
def test_build_rbh_orthogroup_protein_blastp_comparisons_keeps_transitive_all_vs_all_grouping() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
        _record("record_c", features=[_cds(18, 27)]),
    ]
    calls: list[tuple[str, str]] = []
    rows_by_call = [
        [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
        [_hit_row("gbd_r0002_cds000001", "gbd_r0001_cds000001")],
        [],
        [],
        [_hit_row("gbd_r0002_cds000001", "gbd_r0003_cds000001")],
        [_hit_row("gbd_r0003_cds000001", "gbd_r0002_cds000001")],
    ]

    def fake_runner(query_fasta: str, subject_fasta: str) -> pd.DataFrame:
        calls.append((query_fasta, subject_fasta))
        rows = rows_by_call[len(calls) - 1]
        return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)

    result = build_rbh_orthogroup_protein_blastp_comparisons(
        records,
        runner=fake_runner,
        identity=0,
    )

    comparisons = result.comparisons
    assert len(calls) == 6
    assert result.orthogroups is not None
    assert set(result.orthogroups.member_by_protein_id) == {
        "gbd_r0001_cds000001",
        "gbd_r0002_cds000001",
        "gbd_r0003_cds000001",
    }
    assert comparisons[0].iloc[0]["orthogroup_id"] == "og_1"
    assert comparisons[1].iloc[0]["orthogroup_id"] == "og_1"


@pytest.mark.linear
def test_build_orthogroups_suggests_names_from_cds_annotations() -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MKT*"],
                        "gene": ["rpoB"],
                        "product": ["DNA-directed RNA polymerase beta subunit"],
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
                        "translation": ["MKT*"],
                        "gene": ["rpoB"],
                        "product": ["DNA-directed RNA polymerase beta subunit"],
                    },
                )
            ],
        ),
        _record(
            "record_c",
            features=[
                _cds(
                    18,
                    27,
                    qualifiers={
                        "translation": ["MKT*"],
                        "product": ["hypothetical protein"],
                        "note": ["product: DNA-directed RNA polymerase beta subunit"],
                    },
                )
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame.from_records(
        [
            _hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001"),
            _hit_row("gbd_r0002_cds000001", "gbd_r0003_cds000001"),
        ],
        columns=COMPARISON_COLUMNS,
    )

    orthogroups = build_orthogroups_from_protein_hits(
        [hits],
        extraction.protein_map,
    )

    assert orthogroups.names_by_orthogroup_id["og_1"] == "DNA-directed RNA polymerase beta subunit"
    assert orthogroups.confidence_by_orthogroup_id["og_1"] == "high"
    assert orthogroups.descriptions_by_orthogroup_id["og_1"] == (
        "Suggested from product annotations in 3 of 3 records."
    )
    candidates = orthogroups.name_candidates_by_orthogroup_id["og_1"]
    assert candidates[0].source == "product"
    assert candidates[0].record_coverage_count == 3
    assert orthogroups.orthogroups["og_1"][0].product == "DNA-directed RNA polymerase beta subunit"


@pytest.mark.linear
def test_convert_protein_hits_to_genomic_links_only_sets_matching_orthogroup_id() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
        _record("record_c", features=[_cds(18, 27)]),
    ]
    extraction = extract_cds_proteins(records)
    orthogroup_edges = pd.DataFrame.from_records(
        [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
        columns=COMPARISON_COLUMNS,
    )
    display_hits = pd.DataFrame.from_records(
        [
            _hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001"),
            _hit_row("gbd_r0001_cds000001", "gbd_r0003_cds000001"),
        ],
        columns=COMPARISON_COLUMNS,
    )
    orthogroups = build_orthogroups_from_protein_hits(
        [orthogroup_edges],
        extraction.protein_map,
    )

    converted = convert_protein_hits_to_genomic_links(
        display_hits,
        extraction.protein_map,
        orthogroups=orthogroups,
    )

    assert converted.iloc[0]["orthogroup_id"] == "og_1"
    assert converted.iloc[1]["orthogroup_id"] == ""


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_uses_rbh_edges_for_orthogroups() -> None:
    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    forward_hits = pd.DataFrame.from_records(
        [
            _hit_row("a1", "b", bitscore=300),
            _hit_row("a2", "b", bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row("b", "a1", bitscore=300)],
        columns=COMPARISON_COLUMNS,
    )
    payload = [
        {
            "pairIndex": 0,
            "queryIndex": 0,
            "subjectIndex": 1,
            "blastText": forward_hits.to_csv(
                sep="\t",
                header=False,
                index=False,
                lineterminator="\n",
            ),
            "queryProteinMap": {
                "a1": _web_protein_entry(
                    "a1",
                    record_index=0,
                    record_id="record_a",
                    feature_index=0,
                    gene="rpoB",
                    product="DNA-directed RNA polymerase beta subunit",
                ),
                "a2": _web_protein_entry(
                    "a2",
                    record_index=0,
                    record_id="record_a",
                    feature_index=1,
                    start=100,
                    end=190,
                ),
            },
            "subjectProteinMap": {
                "b": _web_protein_entry(
                    "b",
                    record_index=1,
                    record_id="record_b",
                    gene="rpoB",
                    product="DNA-directed RNA polymerase beta subunit",
                ),
            },
        },
        {
            "pairIndex": 0,
            "queryIndex": 1,
            "subjectIndex": 0,
            "blastText": reverse_hits.to_csv(
                sep="\t",
                header=False,
                index=False,
                lineterminator="\n",
            ),
            "queryProteinMap": {
                "b": _web_protein_entry(
                    "b",
                    record_index=1,
                    record_id="record_b",
                    gene="rpoB",
                    product="DNA-directed RNA polymerase beta subunit",
                ),
            },
            "subjectProteinMap": {
                "a1": _web_protein_entry(
                    "a1",
                    record_index=0,
                    record_id="record_a",
                    feature_index=0,
                    gene="rpoB",
                    product="DNA-directed RNA polymerase beta subunit",
                ),
                "a2": _web_protein_entry(
                    "a2",
                    record_index=0,
                    record_id="record_a",
                    feature_index=1,
                    start=100,
                    end=190,
                ),
            },
        },
    ]

    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        json.dumps(payload),
        "orthogroup",
        2,
        50,
        "1e-5",
        0,
        0,
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    assert result["orthogroups"][0]["member_count"] == 2
    assert result["orthogroups"][0]["name"] == "DNA-directed RNA polymerase beta subunit"
    assert result["orthogroups"][0]["nameConfidence"] == "high"
    assert result["orthogroups"][0]["nameCandidates"][0]["recordCoverageCount"] == 2
    assert result["orthogroups"][0]["members"][0]["product"] == "DNA-directed RNA polymerase beta subunit"
    rows = result["pairs"][0]["rows"]
    assert rows[0]["orthogroup_id"] == "og_1"
    assert len(rows) == 1


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
def test_orthogroup_alignment_canvas_adjustment_fits_negative_record_offsets() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
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

    shift_x, width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        {1: -300.0},
    )

    assert shift_x == pytest.approx(300.0)
    assert width_extension == pytest.approx(300.0)


@pytest.mark.linear
def test_orthogroup_alignment_canvas_adjustment_extends_positive_record_offsets() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
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

    shift_x, width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        {1: 250.0},
    )

    assert shift_x == pytest.approx(0.0)
    assert width_extension == pytest.approx(250.0)


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
def test_build_linear_diagram_forwards_protein_blastp_options(
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
            protein_blastp_mode="orthogroup",
            losatp_bin="custom-losat",
            protein_blastp_max_hits=7,
            protein_blastp_candidate_limit=99,
        ),
    )

    assert isinstance(canvas, Drawing)
    assert captured["protein_blastp_mode"] == "orthogroup"
    assert captured["losatp_bin"] == "custom-losat"
    assert captured["protein_blastp_max_hits"] == 7
    assert captured["protein_blastp_candidate_limit"] == 99
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
        options=DiagramOptions(
            protein_blastp_mode="orthogroup",
            align_orthogroup_feature="fanchor",
        ),
    )

    assert isinstance(canvas, Drawing)
    assert captured["align_orthogroup_feature"] == "fanchor"


@pytest.mark.linear
def test_linear_cli_rejects_blast_with_protein_blastp_mode() -> None:
    with pytest.raises(SystemExit):
        linear_cli_module.linear_main(
            [
                "--gbk",
                "a.gb",
                "b.gb",
                "-b",
                "a_b.tsv",
                "--protein_blastp_mode",
                "pairwise",
            ]
        )


@pytest.mark.linear
def test_linear_cli_requires_two_records_for_protein_blastp_mode(
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
                "--protein_blastp_mode",
                "pairwise",
                "--format",
                "svg",
            ]
        )


@pytest.mark.linear
def test_linear_cli_forwards_protein_blastp_options(
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
            "--protein_blastp_mode",
            "orthogroup",
            "--losatp_bin",
            "custom-losat",
            "--protein_blastp_max_hits",
            "9",
            "--protein_blastp_candidate_limit",
            "123",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["protein_blastp_mode"] == "orthogroup"
    assert captured["losatp_bin"] == "custom-losat"
    assert captured["protein_blastp_max_hits"] == 9
    assert captured["protein_blastp_candidate_limit"] == 123
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
            "--protein_blastp_mode",
            "orthogroup",
            "--align_orthogroup_feature",
            "fanchor",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["align_orthogroup_feature"] == "fanchor"
