from __future__ import annotations

import inspect
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.io import load_gbks as api_load_gbks
from gbdraw.api.io import load_gff_fasta as api_load_gff_fasta
from gbdraw.cli_utils.common import load_records_table_records
from gbdraw.io.cli_tables import RecordsTable
from gbdraw.io.cli_tables import RecordsTableRow
from gbdraw.io.genome import load_gbks
from gbdraw.io.genome import load_gff_fasta


def _write_genbank(path: Path, records: list[SeqRecord]) -> None:
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    SeqIO.write(records, path, "genbank")


def test_record_loader_signatures_are_mode_neutral() -> None:
    for loader in (load_gbks, load_gff_fasta, api_load_gbks, api_load_gff_fasta):
        parameters = inspect.signature(loader).parameters
        assert "mode" not in parameters
        assert "load_comparison" not in parameters


def test_load_gbks_preserves_all_records_without_topology_policy(
    tmp_path: Path,
    caplog,
) -> None:
    path = tmp_path / "records.gb"
    records = [
        SeqRecord(Seq("AAGC"), id="linear_record"),
        SeqRecord(Seq("TTAA"), id="circular_record"),
    ]
    records[0].annotations["topology"] = "linear"
    records[1].annotations["topology"] = "circular"
    _write_genbank(path, records)

    caplog.set_level(logging.WARNING, logger="gbdraw.io.genome")
    loaded = load_gbks([str(path)])

    assert [record.id for record in loaded] == ["linear_record", "circular_record"]
    assert all(
        record.annotations["gbdraw_source_file"] == str(path) for record in loaded
    )
    assert "visualize it as circular" not in caplog.text


def test_load_gbks_keeps_neutral_selector_and_reverse_controls(
    tmp_path: Path,
) -> None:
    path = tmp_path / "records.gb"
    _write_genbank(
        path,
        [
            SeqRecord(Seq("AAAA"), id="first"),
            SeqRecord(Seq("AAGC"), id="second"),
        ],
    )

    loaded = api_load_gbks(
        [str(path)],
        record_selectors=["#2"],
        reverse_flags=[True],
    )

    assert [record.id for record in loaded] == ["second"]
    assert str(loaded[0].seq) == "GCTT"
    assert loaded[0].annotations["gbdraw_coord_base"] == 4
    assert loaded[0].annotations["gbdraw_coord_step"] == -1


def test_records_table_loader_forwards_only_neutral_record_controls() -> None:
    row = RecordsTableRow(
        row_index=0,
        row_number=2,
        gbk="/inputs/record.gb",
        gff="",
        fasta="",
        record_label="",
        record_subtitle="",
        record_id="#2",
        region="",
        reverse_complement=True,
        order=None,
        row=None,
        column=None,
    )
    table = RecordsTable(
        table_path="/inputs/records.tsv",
        input_kind="gbk",
        rows=(row,),
        path_dependencies=(),
    )
    record = SeqRecord(Seq("AAGC"), id="second")
    captured: dict[str, object] = {}

    def fake_loader(paths, **kwargs):
        captured["paths"] = paths
        captured["kwargs"] = kwargs
        return [record]

    loaded = load_records_table_records(
        table,
        selected_features_set={"CDS"},
        color_table=None,
        feature_table=None,
        gbk_loader=fake_loader,
    )

    assert loaded == [record]
    assert captured == {
        "paths": ["/inputs/record.gb"],
        "kwargs": {
            "record_selectors": ["#2"],
            "reverse_flags": [True],
        },
    }
