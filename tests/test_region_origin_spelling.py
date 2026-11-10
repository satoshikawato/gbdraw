"""The command-line spelling of a window across the origin: ``START-END:origin`` (PD-OI-122, P0-Q15=A).

A span whose start is greater than its end and has no ``:origin`` stays the reversed span it was in v0.14.5.
``:origin`` comes after the span and before ``:rc``; it is refused, naming the record, on a linear record and on a
span whose start is not greater than its end. Every writer of a region string writes it for a window across the
origin, so a round trip keeps the window.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from gbdraw.api.record_planning import record_input_manifest_from_paths, record_input_manifest_from_table
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import (
    GenBankInputSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import RecordWindow, parse_region_spec
from tests.utils.session_files import read_session_payload

REPO = Path(__file__).resolve().parents[1]
MITO = REPO / "tests" / "test_inputs" / "HmmtDNA.gbk"  # circular, 16,569 bp


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("16000-1000:origin", (None, 16000, 1000, False)),
        ("16000-1000:origin:rc", (None, 16000, 1000, True)),
        ("16000..1000:ORIGIN:minus", (None, 16000, 1000, True)),
        ("NC_012920:16000-1000:origin", ("NC_012920", 16000, 1000, False)),
        # Without ``:origin`` the larger number first is the reversed span, as in v0.14.5.
        ("16000-1000", (None, 1000, 16000, True)),
        ("16000-1000:rc", (None, 1000, 16000, True)),
    ],
)
def test_origin_marks_a_span_across_the_origin_and_its_absence_keeps_the_reversed_span(text, expected) -> None:
    spec = parse_region_spec(text)
    assert (spec.record_id, spec.start, spec.end, spec.reverse_complement) == expected


@pytest.mark.parametrize("text", ["NC_012920:100-200:origin", "NC_012920:200-200:origin:rc"])
def test_origin_on_a_span_whose_start_is_not_greater_than_its_end_is_refused_naming_the_record(text) -> None:
    with pytest.raises(ValueError, match=r"NC_012920.*not greater than its end"):
        parse_region_spec(text)


def test_origin_after_rc_is_not_a_region() -> None:
    with pytest.raises(ValueError, match=r"start-end\[:origin\]\[:rc\]"):
        parse_region_spec("16000-1000:rc:origin")


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_record_input_region_reads_origin_as_region_does_and_writes_it_back(reverse: bool) -> None:
    text = "16000-1000:origin" + (":rc" if reverse else "")
    record = RecordInput(GenBankInputSource(MITO), region=parse_region_spec(text))
    assert record.windows == (RecordWindow(16000, 1000, reverse),)
    assert record.region == parse_region_spec(text)
    assert record.region is not None and record.region.raw == text
    assert RecordWindow(16000, 1000, reverse).region().raw == text


def _cli_manifest(regions: list[str]):
    return record_input_manifest_from_paths(
        gbk_paths=[str(MITO)],
        gff_paths=None,
        fasta_paths=None,
        cardinalities=RecordCardinality.EXACTLY_ONE,
        selectors=None,
        reverse_flags=None,
        labels=None,
        subtitles=None,
        regions=regions,
    )


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_a_command_line_region_with_origin_draws_the_window_across_the_origin(reverse: bool) -> None:
    manifest = _cli_manifest(["16000-1000:origin" + (":rc" if reverse else "")])
    collection = resolve_request_records(
        LinearDiagramRequest(records=manifest.records, record_options=manifest.record_options)
    )
    windowed = resolve_request_records(
        LinearDiagramRequest(records=(RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(16000, 1000, reverse),)),))
    )
    [unit] = collection.records
    assert len(unit) == 570 + 1000
    assert str(unit.seq) == str(windowed.records[0].seq)


def _table(tmp_path: Path, header: str, rows: list[str]) -> Path:
    table = tmp_path / "records.tsv"
    table.write_text(header + "\n" + "".join(f"{row}\n" for row in rows), encoding="utf-8")
    return table


def test_a_records_table_region_with_origin_is_a_window_across_the_origin(tmp_path: Path) -> None:
    table = _table(
        tmp_path,
        "gbk\trow\tcolumn\tregion\trecord_group",
        [f"{MITO}\t1\t1\t16000-1000:origin:rc\tm", f"{MITO}\t1\t2\t5000-5600\tm", f"{MITO}\t2\t1\t16000-1000\t"],
    )
    manifest = record_input_manifest_from_table(str(table))
    assert manifest.records[0].windows == (RecordWindow(16000, 1000, True, 1), RecordWindow(5000, 5600, False, 2))
    # Without ``:origin`` the row stays the reversed span.
    assert manifest.records[1].windows == (RecordWindow(1000, 16000, True),)


def test_origin_on_a_linear_record_is_refused_naming_the_record(tmp_path: Path) -> None:
    table = _table(tmp_path, "gbk\ttopology\tregion", [f"{MITO}\tlinear\t16000-1000:origin"])
    manifest = record_input_manifest_from_table(str(table))
    with pytest.raises(ValidationError, match=r"NC_012920.*linear record"):
        resolve_request_records(LinearDiagramRequest(records=manifest.records))


def test_a_records_table_origin_whose_start_is_not_greater_than_its_end_names_the_row_and_record(
    tmp_path: Path,
) -> None:
    table = _table(tmp_path, "gbk\trecord_id\tregion", [f"{MITO}\tNC_012920\t100-200:origin"])
    with pytest.raises(ValidationError, match=r"row 2.*NC_012920.*not greater than its end"):
        record_input_manifest_from_table(str(table))


def _cli(monkeypatch: pytest.MonkeyPatch, arguments: list[str]) -> None:
    from gbdraw.cli import main

    monkeypatch.setattr(sys, "argv", ["gbdraw", *arguments])
    main()


def test_a_cli_session_of_a_window_across_the_origin_replays_the_figure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    _cli(monkeypatch, ["linear", "--gbk", str(MITO), "-o", "w", "-f", "svg", "--save_session",
                       "--region", "16000-1000:origin:rc", "--region", "5000-5600"])
    document = read_session_payload(tmp_path / "w.gbdraw")
    [record] = document["drawings"][0]["renderRequest"]["records"]
    assert [(w["start"], w["end"], w["reverseComplement"]) for w in record["windows"]] == [
        (16000, 1000, True), (5000, 5600, False)
    ]
    _cli(monkeypatch, ["render", "--session", "w.gbdraw", "-o", "replay", "-f", "svg"])
    assert (tmp_path / "replay.svg").read_bytes() == (tmp_path / "w.svg").read_bytes()
