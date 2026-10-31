"""Automatic orientation of Linear records from nucleotide hits (``--auto_orient``).

The planner reverse-complements a record when its reverse-strand hits against
the row above, as that row is drawn, weigh more than half of its hits
(``alignment_length``). An explicit reverse wins, protein comparisons are not
read, and the figure equals the one an explicit reverse of the same record
draws.
"""
from __future__ import annotations

import csv
import logging
import os
import subprocess
import sys
from dataclasses import fields
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

from gbdraw import LinearComparisonOptions
from gbdraw import linear as linear_cli
from gbdraw.api import (
    GenBankInputSource,
    InMemoryRecordSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    plan_linear_request,
)
from gbdraw.api.record_planning import record_is_reverse_complemented
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear_comparison import LinearComparison, auto_orient_reversed_records

_EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
_QUERY = _EXAMPLES / "MjeNMV.gb"
_SUBJECT = _EXAMPLES / "MelaMJNV.gb"
_TABLE = _EXAMPLES / "MjeNMV.MelaMJNV.tblastx.out"


def _rows(*rows: tuple[int, int, int, int, int], bound: bool = False) -> pd.DataFrame:
    """Rows ``(qstart, qend, sstart, send, alignment_length)``."""

    frame = pd.DataFrame(
        [("Q", "S", 99.0, length, 0, 0, qs, qe, ss, se, 1e-20, 50.0) for qs, qe, ss, se, length in rows],
        columns=list(COMPARISON_COLUMNS),
    )
    if bound:
        for role in ("query", "subject"):
            frame[f"{role}_feature_index"] = "0"
            frame[f"{role}_feature_svg_id"] = f"{role}-feature"
    return frame


@pytest.mark.linear
def test_reverse_hits_over_half_the_aligned_length_reverse_the_record() -> None:
    reverse = _rows((1, 60, 60, 1, 60), (100, 140, 100, 140, 40))
    forward = _rows((1, 40, 40, 1, 40), (100, 160, 100, 160, 60))
    assert auto_orient_reversed_records([LinearComparison(0, 1, reverse)], (0, 1), (False, False)) == {1}
    assert auto_orient_reversed_records([LinearComparison(0, 1, forward)], (0, 1), (False, False)) == set()
    # Exactly half, or no hits, keeps the record as it is.
    tie = _rows((1, 50, 50, 1, 50), (100, 149, 100, 149, 50))
    assert auto_orient_reversed_records([LinearComparison(0, 1, tie)], (0, 1), (False, False)) == set()
    assert auto_orient_reversed_records([], (0, 1), (False, False)) == set()


@pytest.mark.linear
def test_each_record_is_oriented_against_the_row_above_as_drawn() -> None:
    # 0-1 is reverse, so record 1 is reversed; 1-2 is forward in the search
    # frame, which is reverse against record 1 as drawn: record 2 follows it.
    comparisons = [
        LinearComparison(0, 1, _rows((1, 60, 60, 1, 60))),
        LinearComparison(2, 1, _rows((1, 60, 1, 60, 60))),
    ]
    assert auto_orient_reversed_records(comparisons, (0, 1, 2), (False, False, False)) == {1, 2}
    # An explicitly reversed record keeps its orientation and orients the next one.
    assert auto_orient_reversed_records(comparisons, (0, 1, 2), (False, True, False)) == {2}
    # Records of the first row are never reversed; rows bound to features
    # (protein results) are not read.
    assert auto_orient_reversed_records(
        [LinearComparison(1, 0, _rows((1, 60, 60, 1, 60)))], (1, 0), (False, False)
    ) == {0}
    assert auto_orient_reversed_records(
        [LinearComparison(0, 1, _rows((1, 60, 60, 1, 60), bound=True))], (0, 1), (False, False)
    ) == set()


def _reversed_subject_inputs(directory: Path) -> tuple[Path, Path]:
    """MelaMJNV reverse-complemented with its features, and the table in its frame."""

    record = SeqIO.read(_SUBJECT, "genbank")
    reversed_record = record.reverse_complement(
        id=True, name=True, description=True, annotations=True, dbxrefs=True
    )
    genbank = directory / "MelaMJNV_rc.gb"
    SeqIO.write(reversed_record, genbank, "genbank")
    length = len(record)
    table = directory / "MjeNMV.MelaMJNV_rc.tsv"
    with _TABLE.open(encoding="utf-8") as source, table.open("w", encoding="utf-8", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        for row in csv.reader(source, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            row[8], row[9] = str(length + 1 - int(row[8])), str(length + 1 - int(row[9]))
            writer.writerow(row)
    return genbank, table


def _render(directory: Path, prefix: str, subject: Path, table: Path, *extra: str) -> str:
    linear_cli.linear_main([
        "--gbk", str(_QUERY), str(subject), "-b", str(table), "-f", "svg",
        "-o", str(directory / prefix), *extra,
    ])
    return (directory / f"{prefix}.svg").read_text(encoding="utf-8")


@pytest.mark.linear
def test_auto_orient_draws_the_figure_of_an_explicit_reverse(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    subject, table = _reversed_subject_inputs(tmp_path)
    with caplog.at_level(logging.INFO):
        auto = _render(tmp_path, "auto", subject, table, "--auto_orient")
    assert "INFO: auto_orient reverse-complemented record #2" in caplog.text
    explicit = _render(
        tmp_path, "explicit", subject, table, "--reverse_complement", "0", "--reverse_complement", "1"
    )
    unoriented = _render(tmp_path, "unoriented", subject, table)
    assert auto == explicit
    assert auto != unoriented


@pytest.mark.linear
def test_an_explicit_reverse_wins_over_auto_orient(tmp_path: Path) -> None:
    flags = ("--reverse_complement", "0", "--reverse_complement", "1")
    explicit = _render(tmp_path, "explicit", _SUBJECT, _TABLE, *flags)
    both = _render(tmp_path, "both", _SUBJECT, _TABLE, *flags, "--auto_orient")
    assert both == explicit
    # Mostly forward hits keep the record as it is.
    assert _render(tmp_path, "auto", _SUBJECT, _TABLE, "--auto_orient") == _render(
        tmp_path, "plain", _SUBJECT, _TABLE
    )


def _plan_subject_reversed(table: pd.DataFrame, *, protein: bool) -> bool:
    comparisons = (
        {"protein_comparisons": (table,)}
        if protein
        else {"linear_comparisons": (LinearComparison(0, 1, table),)}
    )
    request = LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(_QUERY)), RecordInput(GenBankInputSource(_SUBJECT))),
        options=LinearDiagramOptions(auto_orient=True, **comparisons),
    )
    return record_is_reverse_complemented(plan_linear_request(request).records[1])


@pytest.mark.linear
def test_protein_comparisons_leave_the_orientation_unchanged() -> None:
    reverse = _rows((1, 600, 6000, 5401, 600))
    assert _plan_subject_reversed(reverse, protein=False) is True
    assert _plan_subject_reversed(reverse, protein=True) is False


@pytest.mark.linear
def test_a_session_does_not_store_auto_orient_yet(tmp_path: Path) -> None:
    from gbdraw.exceptions import ValidationError
    from gbdraw.session_request_codec import encode_canonical_request

    request = LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(_QUERY)), RecordInput(GenBankInputSource(_SUBJECT))),
        options=LinearDiagramOptions(auto_orient=True),
    )
    with pytest.raises(ValidationError, match="auto_orient"):
        encode_canonical_request(request)
    out = tmp_path / "out"
    out.mkdir()
    completed = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear",
         "--gbk", str(_QUERY), str(_SUBJECT), "-b", str(_TABLE), "-f", "svg",
         "-o", str(out / "saved"), "--auto_orient", "--save_session"],
        capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
    )
    assert completed.returncode == 1
    assert "A Session does not store auto_orient yet; save the Session without it." in completed.stderr
    # Refused before anything is drawn.
    assert list(out.iterdir()) == []


@pytest.mark.linear
def test_the_plan_request_of_in_memory_records_carries_the_reversal() -> None:
    # Exact-one in-memory records normally keep the request as given; after a
    # reversal the plan's request must name the reversed record it draws.
    records = (SeqIO.read(_QUERY, "genbank"), SeqIO.read(_SUBJECT, "genbank"))
    request = LinearDiagramRequest(
        records=tuple(RecordInput(InMemoryRecordSource(record)) for record in records),
        options=LinearDiagramOptions(
            auto_orient=True,
            linear_comparisons=(LinearComparison(0, 1, _rows((1, 600, 6000, 5401, 600))),),
        ),
    )
    plan = plan_linear_request(request)
    assert record_is_reverse_complemented(plan.records[1])
    planned = plan.request.records[1].source
    assert isinstance(planned, InMemoryRecordSource)
    assert record_is_reverse_complemented(planned.record)


@pytest.mark.linear
def test_auto_orient_reads_a_comparisons_table_on_a_multi_record_row(tmp_path: Path) -> None:
    # Row 1: MjeNMV; row 2: the reversed MelaMJNV and PeseMJNV, which has no hits.
    subject, table = _reversed_subject_inputs(tmp_path)
    manifest = tmp_path / "comparisons.tsv"
    manifest.write_text(f"blast\tquery\tsubject\n{table}\t#1\t#2\n", encoding="utf-8")

    def render(prefix: str, *extra: str) -> str:
        linear_cli.linear_main([
            "--gbk", str(_QUERY), str(subject), str(_EXAMPLES / "PeseMJNV.gb"),
            "--comparisons_table", str(manifest),
            "--multi_record_position", "#1@1", "--multi_record_position", "#2@2",
            "--multi_record_position", "#3@2", "-f", "svg", "-o", str(tmp_path / prefix), *extra,
        ])
        return (tmp_path / f"{prefix}.svg").read_text(encoding="utf-8")

    explicit = render(
        "explicit", "--reverse_complement", "0", "--reverse_complement", "1", "--reverse_complement", "0"
    )
    auto = render("auto", "--auto_orient")
    assert auto == explicit
    assert auto != render("plain")


def test_auto_orient_is_the_last_positional_comparison_field() -> None:
    # LinearComparisonOptions is public; released positional fields keep their places.
    assert [item.name for item in fields(LinearComparisonOptions)][-2:] == [
        "similarity_alignment", "auto_orient",
    ]

