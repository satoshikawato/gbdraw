"""A FASTA given without GFF3 is a set of records with sequence and no features (PD-OI-068 rev 2).

One loader (``gbdraw.io.genome.load_fasta``) owns the records and the
nucleotide refusal for the CLI, the Python API, the typed request, and the
Web's Worker helper. A FASTA in a GFF3 + FASTA pair keeps its own reader.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
from pathlib import Path

import pytest

import gbdraw
from gbdraw.api import FastaInputSource
from gbdraw.api.record_planning import record_input_manifest_from_table, resolve_record_inputs
from gbdraw.circular import circular_main
from gbdraw.cli_utils.common import add_input_args, validate_input_args
from gbdraw.exceptions import ValidationError
from gbdraw.io.genome import load_fasta, load_gff_fasta
from gbdraw.layout.record_transform import detected_topology
from gbdraw.linear import linear_main
from gbdraw.web_support.error_adapter import serialize_web_error

def _two_records(path: Path) -> Path:
    path.write_text(">first some description\nACGTACGTNN\nacgtu\n>second\n" + "GC" * 30 + "\n")
    return path


def test_load_fasta_reads_every_entry_as_a_record_without_features(tmp_path: Path) -> None:
    records = load_fasta([str(_two_records(tmp_path / "two.fasta"))])

    assert [(record.id, len(record), len(record.features)) for record in records] == [
        ("first", 15, 0), ("second", 60, 0),
    ]
    assert records[0].annotations["gbdraw_source_file"] == str(tmp_path / "two.fasta")
    # The file makes no topology claim.
    assert [detected_topology(record, "fasta") for record in records] == ["unknown", "unknown"]


def test_load_fasta_reads_gzip_and_applies_selector_and_reverse(tmp_path: Path) -> None:
    plain = _two_records(tmp_path / "two.fasta")
    packed = tmp_path / "two.fasta.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()))

    (record,) = load_fasta([str(packed)], record_selectors=["first"], reverse_flags=[True])

    assert record.id == "first"
    assert str(record.seq) == str(load_fasta([str(plain)])[0].seq.reverse_complement())
    assert record.annotations["gbdraw_coord_step"] == -1


@pytest.mark.parametrize(
    ("text", "reason", "letter"),
    [
        (">dna\nACGT\n>LvMJNV_ORF1 protein\nMKLVE\n", "NOT_NUCLEOTIDE", "L"),
        (">dna\nACGT\n>LvMJNV_ORF1\n\n", "EMPTY_SEQUENCE", None),
        (">LvMJNV_ORF1\nACGT-ACGT\n", "NOT_NUCLEOTIDE", "-"),
    ],
)
def test_load_fasta_refuses_a_record_that_is_not_nucleotide(
    tmp_path: Path, text: str, reason: str, letter: str | None
) -> None:
    path = tmp_path / "input.fasta"
    path.write_text(text)

    with pytest.raises(ValidationError) as caught:
        load_fasta([str(path)])

    position = 1 if text.startswith(">LvMJNV_ORF1") else 2
    expected = {"code": "INPUT_INVALID", "reason": reason, "recordName": "LvMJNV_ORF1", "recordPosition": position}
    if letter is not None:
        expected["codepoint"] = ord(letter)
    assert caught.value.diagnostic == expected
    assert "LvMJNV_ORF1" in str(caught.value)
    # The Web receives the record's name (not its position, as the name can be shown) and the letter, not the message.
    payload = serialize_web_error(caught.value, operation="listSequenceRecords", stage="helper")
    assert payload["code"] == "INPUT_INVALID"
    assert payload["context"] == {
        key: value for key, value in expected.items() if key not in {"code", "recordPosition"}
    }


@pytest.mark.parametrize(
    ("header", "shown"),
    [
        ("\u30bf\u30f3\u30d1\u30af", "\u30bf\u30f3\u30d1\u30af"),  # a non-ASCII ID: the CLI shows it, the Web gets the position
        ("x" * 81, "x" * 81),  # an ID over the bound
        ("", "#2"),  # ">" alone: no ID
    ],
)
def test_a_refused_record_whose_id_cannot_be_shown_is_named_by_its_position(
    tmp_path: Path, header: str, shown: str
) -> None:
    path = tmp_path / "input.fasta"
    path.write_text(f">dna\nACGT\n>{header}\nMKLVE\n", encoding="utf-8")

    with pytest.raises(ValidationError) as caught:
        load_fasta([str(path)])

    assert f"FASTA record {shown} in input.fasta holds 'L'" in str(caught.value)
    payload = serialize_web_error(caught.value, operation="listSequenceRecords", stage="helper")
    assert payload["context"] == {"reason": "NOT_NUCLEOTIDE", "recordPosition": 2, "codepoint": ord("L")}


def test_the_web_wording_reads_the_record_context_keys_a_fasta_refusal_sends() -> None:
    # The Web wording owner reads the context keys Python sends (the vocabulary guard in
    # test_web_error_producer_coverage.py does not list these two).
    text = (Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "js" / "utils" / "error-normalization.js").read_text(
        encoding="utf-8"
    )
    block = text[text.index("const contextFor") : text.index("export const normalizeUserFacingError")]
    web_context_keys = set(re.findall(r"'([A-Za-z]+)'", block)) | set(re.findall(r"value\.([A-Za-z]+)", block))
    assert {"recordName", "recordPosition"} <= web_context_keys


def test_load_fasta_warns_on_repeated_ids(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    path = tmp_path / "repeat.fasta"
    path.write_text(">same\nACGT\n>same\nGGCC\n")

    with caplog.at_level(logging.WARNING):
        records = load_fasta([str(path)])

    assert [str(record.seq) for record in records] == ["ACGT", "GGCC"]
    assert "Record same seems to have been already loaded" in caplog.text


def test_gff_fasta_pair_still_ignores_entries_without_gff_rows(tmp_path: Path) -> None:
    fasta = _two_records(tmp_path / "two.fasta")
    gff = tmp_path / "one.gff3"
    gff.write_text("##gff-version 3\nsecond\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n")
    # The pair does not apply the FASTA-only nucleotide check either.
    fasta.write_text(fasta.read_text() + ">protein\nMKLVE\n")

    assert [record.id for record in load_gff_fasta([str(gff)], [str(fasta)])] == ["second"]


def test_read_fasta_is_the_beginner_reader(tmp_path: Path) -> None:
    assert "read_fasta" in gbdraw.__all__
    records = gbdraw.read_fasta(_two_records(tmp_path / "two.fasta"))

    assert [record.id for record in records] == ["first", "second"]


def _parse_input_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    add_input_args(parser)
    parser.add_argument("--records_table", default=None)
    args = parser.parse_args(argv)
    validate_input_args(parser, args)
    return args


def test_cli_accepts_fasta_alone_and_still_refuses_gff_alone(capsys: pytest.CaptureFixture[str]) -> None:
    assert _parse_input_args(["--fasta", "a.fasta"]).fasta == ["a.fasta"]
    with pytest.raises(SystemExit):
        _parse_input_args(["--gff", "a.gff3"])
    assert "--gff requires --fasta" in capsys.readouterr().err


def test_cli_draws_a_fasta_alone_in_both_modes(tmp_path: Path) -> None:
    fasta = _two_records(tmp_path / "two.fasta")

    circular_main(["--fasta", str(fasta), "-o", str(tmp_path / "circular"), "-f", "svg"])
    linear_main(["--fasta", str(fasta), "-o", str(tmp_path / "linear"), "-f", "svg"])

    # Circular draws each record of a multi-record input to its own file; Linear stacks them.
    assert "first" in (tmp_path / "circular_1.svg").read_text()
    assert "second" in (tmp_path / "circular_2.svg").read_text()
    linear = (tmp_path / "linear.svg").read_text()
    assert "first" in linear and "second" in linear


def test_cli_refuses_a_protein_fasta_alone(tmp_path: Path) -> None:
    fasta = tmp_path / "protein.fasta"
    fasta.write_text(">LvMJNV_ORF1\nMKLVE\n")

    with pytest.raises(ValidationError, match="LvMJNV_ORF1"):
        circular_main(["--fasta", str(fasta), "-o", str(tmp_path / "circular"), "-f", "svg"])


def test_records_table_row_with_only_a_fasta_path_is_a_fasta_source(tmp_path: Path) -> None:
    _two_records(tmp_path / "two.fasta")
    table = tmp_path / "records.tsv"
    table.write_text("fasta\trecord_id\ntwo.fasta\tsecond\n")

    manifest = record_input_manifest_from_table(str(table))

    assert manifest.records[0].source == FastaInputSource(tmp_path / "two.fasta")
    resolved = resolve_record_inputs(manifest.records, gff_candidate_features=None, gff_keep_all_features=True)
    assert [record.id for record in resolved.records] == ["second"]
    assert resolved.provenance[0].source_kind == "fasta"


def test_records_table_refuses_mixing_fasta_rows_with_gff_rows(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text("gff\tfasta\na.gff3\ta.fasta\n\tb.fasta\n")

    with pytest.raises(ValidationError, match="cannot mix GFF3/FASTA rows with FASTA rows"):
        record_input_manifest_from_table(str(table))


def test_worker_helper_lists_a_fasta_alone_through_the_fasta_only_loader(tmp_path: Path) -> None:
    source = (Path(__file__).resolve().parents[1] / "gbdraw/web/js/app/python-helpers.js").read_text()
    helpers: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0], helpers)
    list_records = helpers["list_sequence_records"]
    assert callable(list_records)

    listed = json.loads(list_records(str(_two_records(tmp_path / "two.fasta")), "fasta"))
    assert [(record["record_id"], record["record_length"], record["topology"]) for record in listed["records"]] == [
        ("first", 15, "unknown"), ("second", 60, "unknown"),
    ]
    protein = tmp_path / "protein.fasta"
    protein.write_text(">LvMJNV_ORF1\nMKLVE\n")
    refused = json.loads(list_records(str(protein), "fasta"))
    assert refused["error"] == {
        "code": "INPUT_INVALID", "operation": "listSequenceRecords", "stage": "helper",
        "context": {"reason": "NOT_NUCLEOTIDE", "recordName": "LvMJNV_ORF1", "codepoint": ord("L")},
    }


FASTA_SESSION = Path(__file__).resolve().parent / "fixtures/sessions/cli-circular-fasta.v51.gbdraw"


def test_a_session_51_with_a_fasta_source_loads_and_replays(tmp_path: Path) -> None:
    from gbdraw.api import load_session_document
    from gbdraw.session import session_drawing_artifacts

    session = load_session_document(FASTA_SESSION)
    drawing = session_drawing_artifacts(session, "circular").fields
    assert drawing["renderRequest"]["records"][0]["source"] == {"kind": "fasta", "fastaResourceId": "record-1-fasta"}

    circular_main(["--session", str(FASTA_SESSION), "-o", str(tmp_path / "replay"), "-f", "svg"])

    assert "NC_012920.1" in (tmp_path / "replay.svg").read_text()
