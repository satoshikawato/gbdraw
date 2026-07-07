from __future__ import annotations

from pathlib import Path

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    read_circular_track_table,
    read_conservation_table,
    read_records_table,
)


def test_conservation_table_resolves_relative_paths(tmp_path: Path) -> None:
    blast = tmp_path / "blast.tsv"
    blast.write_text("", encoding="utf-8")
    table_dir = tmp_path / "tables"
    table_dir.mkdir()
    table = table_dir / "conservation.tsv"
    table.write_text("blast\tlabel\tcolor\n../blast.tsv\tReference A\t#E15759\n", encoding="utf-8")

    parsed = read_conservation_table(str(table))

    assert parsed.conservation_blast_files == [str(blast.resolve())]
    assert parsed.labels == ["Reference A"]
    assert parsed.colors == ["#E15759"]
    assert parsed.path_dependencies[0].column == "blast"
    assert parsed.path_dependencies[0].row_index == 0


def test_conservation_table_rejects_unknown_columns(tmp_path: Path) -> None:
    table = tmp_path / "bad.tsv"
    table.write_text("blast\tcolour\nx.tsv\tred\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="unknown table column 'colour'"):
        read_conservation_table(str(table))


def test_circular_track_table_defaults_first_features_row_to_axis(tmp_path: Path) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text(
        "id\trenderer\tside\tw\tparams\n"
        "features\tfeatures\t\t\t\n"
        "gc_content\tdinucleotide_content\tinside\t0.1\tnt=GC\n",
        encoding="utf-8",
    )

    parsed = read_circular_track_table(str(table))

    assert parsed.axis_index == 0
    assert parsed.slot_specs[0] == "features:features@side=overlay,lane_direction=split"
    assert parsed.slot_specs[1] == "gc_content:dinucleotide_content@w=0.1,side=inside,nt=GC"


def test_circular_track_table_rejects_duplicate_axis_rows(tmp_path: Path) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text(
        "id\trenderer\tside\n"
        "features\tfeatures\taxis\n"
        "features2\tfeatures\taxis\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="only one side=axis"):
        read_circular_track_table(str(table))


def test_circular_track_table_rejects_non_feature_axis_row(tmp_path: Path) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text("id\trenderer\tside\ngc\tdinucleotide_content\taxis\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="side=axis is only supported"):
        read_circular_track_table(str(table))


def test_records_table_sorts_by_order_and_generates_positions(tmp_path: Path) -> None:
    first = tmp_path / "a.gbk"
    second = tmp_path / "b.gbk"
    first.write_text("", encoding="utf-8")
    second.write_text("", encoding="utf-8")
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_label\torder\trow\tcolumn\tregion\treverse_complement\n"
        "b.gbk\tB\t2\t1\t2\t200-300\t0\n"
        "a.gbk\tA\t1\t1\t1\t100-150:rc\ttrue\n",
        encoding="utf-8",
    )

    parsed = read_records_table(str(table))

    assert parsed.gbk_files == [str(first.resolve()), str(second.resolve())]
    assert parsed.record_labels == ["A", "B"]
    assert parsed.reverse_flags == [True, False]
    assert parsed.row_scoped_region_specs() == ["#1:100-150:rc", "#2:200-300"]
    assert parsed.multi_record_positions() == ["#1@1", "#2@1"]


def test_records_table_rejects_invalid_boolean(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text("gbk\treverse_complement\na.gbk\tmaybe\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="expected a boolean value"):
        read_records_table(str(table))


def test_records_table_rejects_mixed_input_kinds(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\tgff\tfasta\n"
        "a.gbk\t\t\n"
        "\tb.gff3\tb.fna\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="cannot mix GenBank rows"):
        read_records_table(str(table))


def test_records_table_rejects_duplicate_grid_cells(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trow\tcolumn\n"
        "a.gbk\t1\t1\n"
        "b.gbk\t1\t1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="duplicate records table placement"):
        read_records_table(str(table))
