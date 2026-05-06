from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import assemble_linear_diagram_from_rows
from gbdraw.exceptions import ValidationError
from gbdraw.io.genome import load_gbk_rows, load_gbks
from gbdraw.layout.linear_rows import build_linear_layout_index, clone_rows_for_rendering, make_linear_input_row


def _record(record_id: str, length: int) -> SeqRecord:
    return SeqRecord(
        Seq("A" * length),
        id=record_id,
        name=record_id,
        description=record_id,
        annotations={"molecule_type": "DNA"},
    )


def _write_gbk(path: Path, records: list[SeqRecord]) -> None:
    SeqIO.write(records, path, "genbank")


def _comparison(query: str = "chr2", subject: str = "chr1") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "query": query,
                "subject": subject,
                "identity": 90.0,
                "alignment_length": 100,
                "mismatches": 0,
                "gap_opens": 0,
                "qstart": 10,
                "qend": 100,
                "sstart": 20,
                "send": 120,
                "evalue": 1e-20,
                "bitscore": 200,
            }
        ]
    )


def test_gbk_row_loader_preserves_all_records_with_comparison_inputs(tmp_path: Path) -> None:
    first = tmp_path / "first.gbk"
    second = tmp_path / "second.gbk"
    _write_gbk(first, [_record("a1", 100), _record("a2", 80)])
    _write_gbk(second, [_record("b1", 90), _record("b2", 70)])

    legacy = load_gbks([str(first), str(second)], "linear", load_comparison=True)
    rows = load_gbk_rows([str(first), str(second)])

    assert [record.id for record in legacy] == ["a1", "b1"]
    assert [[record.id for record in row.records] for row in rows] == [["a1", "a2"], ["b1", "b2"]]


def test_linear_layout_index_resolves_duplicate_ids_locally() -> None:
    rows = [
        make_linear_input_row(
            row_index=0,
            source_id="a",
            source_path="a.gbk",
            label="A",
            records=[_record("shared", 100), _record("other", 50)],
        ),
        make_linear_input_row(
            row_index=1,
            source_id="b",
            source_path="b.gbk",
            label="B",
            records=[_record("shared", 80)],
        ),
    ]
    layout = build_linear_layout_index(
        rows,
        alignment_width=1000,
        record_gap=20,
        normalize_length=False,
        align_center=False,
    )

    assert layout.resolve_record(0, "shared").ref.source_id == "a"
    assert layout.resolve_record(1, "shared").ref.source_id == "b"


def test_source_row_layout_uses_original_record_id_as_segment_label() -> None:
    rows = [
        make_linear_input_row(
            row_index=0,
            source_id="seq_0",
            source_path="seq_0.gb",
            label="seq_0.gb",
            records=[_record("NC_002204.1", 100), _record("NC_002211.1", 80)],
        )
    ]
    render_rows, _id_map = clone_rows_for_rendering(rows)
    layout = build_linear_layout_index(
        render_rows,
        alignment_width=1000,
        record_gap=20,
        normalize_length=False,
        align_center=False,
    )

    assert [segment.ref.label for segment in layout.rows[0].records] == ["NC_002204.1", "NC_002211.1"]
    assert [segment.ref.unique_key for segment in layout.rows[0].records] == [
        "seq_0:1:NC_002204.1",
        "seq_0:2:NC_002211.1",
    ]


def test_source_row_assembly_links_any_segment_in_adjacent_rows() -> None:
    rows = [
        make_linear_input_row(
            row_index=0,
            source_id="srcA",
            source_path="a.gbk",
            label="A",
            records=[_record("chr1", 1000), _record("chr2", 500)],
        ),
        make_linear_input_row(
            row_index=1,
            source_id="srcB",
            source_path="b.gbk",
            label="B",
            records=[_record("chr1", 700)],
        ),
    ]

    canvas = assemble_linear_diagram_from_rows(
        rows,
        comparison_dataframes=[_comparison()],
        legend="none",
        config_overrides={"show_gc": False, "show_skew": False, "show_labels": "none"},
    )
    svg = canvas.tostring()

    assert "comparison1" in svg
    assert "data-query-record-key=\"srcA:2:chr2\"" in svg
    assert "data-subject-record-key=\"srcB:1:chr1\"" in svg


def test_source_row_assembly_rejects_ambiguous_ids_within_row() -> None:
    rows = [
        make_linear_input_row(
            row_index=0,
            source_id="srcA",
            source_path="a.gbk",
            label="A",
            records=[_record("dup", 200), _record("dup", 200)],
        ),
        make_linear_input_row(
            row_index=1,
            source_id="srcB",
            source_path="b.gbk",
            label="B",
            records=[_record("target", 300)],
        ),
    ]

    with pytest.raises(ValidationError, match="ambiguous"):
        assemble_linear_diagram_from_rows(
            rows,
            comparison_dataframes=[_comparison(query="dup", subject="target")],
            legend="none",
            config_overrides={"show_gc": False, "show_skew": False, "show_labels": "none"},
        )
