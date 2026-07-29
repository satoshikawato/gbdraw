from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.circular import circular_main
from gbdraw.core.sequence import determine_output_file_prefixes


def _record(record_id: str) -> SeqRecord:
    return SeqRecord(
        Seq("ATGC" * 25),
        id=record_id,
        annotations={"molecule_type": "DNA"},
    )


def test_implicit_circular_batch_prefixes_avoid_existing_suffix_collisions() -> None:
    records = [_record("record"), _record("record_2"), _record("record")]

    assert determine_output_file_prefixes(records, None) == [
        "record",
        "record_2",
        "record_3",
    ]
    assert determine_output_file_prefixes(records, "diagram") == [
        "diagram_1",
        "diagram_2",
        "diagram_3",
    ]


def test_circular_cli_duplicate_record_ids_do_not_overwrite(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_paths = [tmp_path / "first.gb", tmp_path / "second.gb"]
    for path in input_paths:
        SeqIO.write(_record("duplicate"), path, "genbank")
    monkeypatch.chdir(tmp_path)

    circular_main(
        [
            "--gbk",
            *(str(path) for path in input_paths),
            "--format",
            "svg",
        ]
    )

    assert (tmp_path / "duplicate.svg").is_file()
    assert (tmp_path / "duplicate_2.svg").is_file()


def test_circular_cli_batch_splits_dotted_path_prefix(
    tmp_path: Path,
) -> None:
    input_paths = [tmp_path / "first.gb", tmp_path / "second.gb"]
    for index, path in enumerate(input_paths, start=1):
        SeqIO.write(_record(f"record-{index}"), path, "genbank")
    output_directory = tmp_path / "nested"
    output_directory.mkdir()
    output_prefix = output_directory / "diagram.v1"

    circular_main(
        [
            "--gbk",
            *(str(path) for path in input_paths),
            "--format",
            "svg",
            "-o",
            str(output_prefix),
        ]
    )

    assert (output_directory / "diagram.v1_1.svg").is_file()
    assert (output_directory / "diagram.v1_2.svg").is_file()


def test_circular_batch_session_writes_all_outputs_and_sidecar(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_paths = [tmp_path / "first.gb", tmp_path / "second.gb"]
    for index, path in enumerate(input_paths, start=1):
        SeqIO.write(_record(f"record-{index}"), path, "genbank")
    monkeypatch.chdir(tmp_path)

    circular_main(
        [
            "--gbk",
            *(str(path) for path in input_paths),
            "--format",
            "svg",
            "--save_session",
        ]
    )

    assert (tmp_path / "record-1.svg").is_file()
    assert (tmp_path / "record-2.svg").is_file()
    assert (tmp_path / "gbdraw.gbdraw-session.json").is_file()
