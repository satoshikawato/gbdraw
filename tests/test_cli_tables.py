from __future__ import annotations

from pathlib import Path

import pytest

from gbdraw.api import (
    CircularTrackTable,
    ConservationTable,
    RecordsTable,
    read_circular_track_table as api_read_circular_track_table,
    read_conservation_table as api_read_conservation_table,
    read_label_override_table,
    read_label_whitelist_table,
    read_qualifier_priority_table,
    read_records_table as api_read_records_table,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    read_circular_track_table,
    read_conservation_table,
    read_records_table,
)


def test_public_table_readers_share_cli_models_and_validation(tmp_path: Path) -> None:
    records_path = tmp_path / "records.tsv"
    records_path.write_text("gbk\nrecord.gbk\n", encoding="utf-8")
    conservation_path = tmp_path / "conservation.tsv"
    conservation_path.write_text("blast\nhits.tsv\n", encoding="utf-8")
    tracks_path = tmp_path / "tracks.tsv"
    tracks_path.write_text("id\trenderer\nfeatures\tfeatures\n", encoding="utf-8")

    records = api_read_records_table(str(records_path))
    conservation = api_read_conservation_table(str(conservation_path))
    tracks = api_read_circular_track_table(str(tracks_path))

    assert isinstance(records, RecordsTable)
    assert isinstance(conservation, ConservationTable)
    assert isinstance(tracks, CircularTrackTable)
    assert records == read_records_table(str(records_path))
    assert conservation == read_conservation_table(str(conservation_path))
    assert tracks == read_circular_track_table(str(tracks_path))


def test_public_label_table_readers(tmp_path: Path) -> None:
    whitelist = tmp_path / "whitelist.tsv"
    whitelist.write_text("CDS\tproduct\tpolymerase\n", encoding="utf-8")
    priority = tmp_path / "priority.tsv"
    priority.write_text("CDS\tgene,product\n", encoding="utf-8")
    override = tmp_path / "override.tsv"
    override.write_text("rec1\tCDS\tgene\tpol\tpolymerase\n", encoding="utf-8")

    assert list(read_label_whitelist_table(str(whitelist)).columns) == [
        "feature_type",
        "qualifier",
        "keyword",
    ]
    assert list(read_qualifier_priority_table(str(priority)).columns) == [
        "feature_type",
        "priorities",
    ]
    assert list(read_label_override_table(str(override)).columns) == [
        "record_id",
        "feature_type",
        "qualifier",
        "value",
        "label_text",
    ]


def test_conservation_table_resolves_relative_paths(tmp_path: Path) -> None:
    blast = tmp_path / "blast.tsv"
    blast.write_text("", encoding="utf-8")
    comparison_fasta = tmp_path / "comparison.fna"
    comparison_fasta.write_text(">comparison\nACGT\n", encoding="utf-8")
    table_dir = tmp_path / "tables"
    table_dir.mkdir()
    table = table_dir / "conservation.tsv"
    table.write_text(
        "blast\tlabel\tcolor\tcomparison_fasta\n"
        "../blast.tsv\tReference A\t#E15759\t../comparison.fna\n",
        encoding="utf-8",
    )

    parsed = read_conservation_table(str(table))

    assert parsed.conservation_blast_files == [str(blast.resolve())]
    assert parsed.labels == ["Reference A"]
    assert parsed.colors == ["#E15759"]
    assert parsed.comparison_fasta_files == [str(comparison_fasta.resolve())]
    assert parsed.path_dependencies[0].column == "blast"
    assert parsed.path_dependencies[0].row_index == 0
    assert parsed.path_dependencies[1].column == "comparison_fasta"


def test_conservation_table_rejects_unknown_columns(tmp_path: Path) -> None:
    table = tmp_path / "bad.tsv"
    table.write_text("blast\tcolour\nx.tsv\tred\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="unknown table column 'colour'"):
        read_conservation_table(str(table))


@pytest.mark.parametrize(
    ("filename", "content", "reader"),
    [
        ("records.tsv", "gbk\na.gbk\n", read_records_table),
        ("conservation.tsv", "blast\na.tsv\n", read_conservation_table),
        ("tracks.tsv", "id\trenderer\nfeatures\tfeatures\n", read_circular_track_table),
    ],
)
def test_cli_tables_accept_utf8_bom(
    tmp_path: Path,
    filename: str,
    content: str,
    reader,
) -> None:
    table = tmp_path / filename
    table.write_text(content, encoding="utf-8-sig")

    parsed = reader(str(table))

    assert parsed.table_path == str(table)


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


@pytest.mark.parametrize(
    "key",
    [
        "id",
        "renderer",
        "type",
        "side",
        "r",
        "radius",
        "w",
        "width",
        "spacing",
        "inner_gap_px",
        "outer_gap_px",
        "z",
        "z_index",
        "zindex",
        "enabled",
        "show",
        "visible",
        "strict",
        "compress",
        "reserve",
    ],
)
def test_circular_track_table_rejects_structural_params(
    tmp_path: Path,
    key: str,
) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text(
        f"id\trenderer\tparams\noriginal\tticks\t{key}=replacement\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValidationError,
        match=rf"{table}.*row 2, column 'params'.*{key}",
    ):
        read_circular_track_table(str(table))


@pytest.mark.parametrize("key", ["lane_direction", "lanes"])
def test_circular_track_table_rejects_feature_lane_params(
    tmp_path: Path,
    key: str,
) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text(
        f"id\trenderer\tparams\nfeatures\tfeatures\t{key}=outside\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match=rf"column 'params'.*{key}"):
        read_circular_track_table(str(table))


def test_circular_track_table_accepts_renderer_specific_params(tmp_path: Path) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text(
        "id\trenderer\tside\tparams\n"
        "ticks\tticks\toutside\ttick_label_layout=tick_only\n"
        "features\tfeatures\taxis\tlegend_label=Genes\n"
        "skew\tdinucleotide_skew\tinside\tnt=AT,positive_color=red,negative_color=blue\n",
        encoding="utf-8",
    )

    parsed = read_circular_track_table(str(table))

    assert parsed.axis_index == 1
    assert "tick_label_layout=tick_only" in parsed.slot_specs[0]
    assert "legend_label=Genes" in parsed.slot_specs[1]
    assert "nt=AT" in parsed.slot_specs[2]


def test_circular_track_table_runs_axis_aware_normalization(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    table = tmp_path / "tracks.tsv"
    table.write_text(
        "id\trenderer\tside\noutside\tticks\toutside\nfeatures\tfeatures\taxis\n",
        encoding="utf-8",
    )
    called: dict[str, int] = {}

    def reject_with_axis(_slots, axis_index):
        called["axis_index"] = axis_index
        raise ValueError("axis-aware failure")

    monkeypatch.setattr(
        "gbdraw.io.cli_tables.normalize_circular_track_slots_with_axis",
        reject_with_axis,
    )

    with pytest.raises(ValidationError, match="invalid circular track table.*axis-aware failure"):
        read_circular_track_table(str(table))
    assert called == {"axis_index": 1}


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


@pytest.mark.parametrize(
    "region",
    [
        "#2:1-10",
        "record_b:1-10",
        "b.gbk:#1:1-10",
    ],
)
def test_records_table_rejects_qualified_regions(tmp_path: Path, region: str) -> None:
    table = tmp_path / "records.tsv"
    table.write_text(f"gbk\tregion\na.gbk\t{region}\n", encoding="utf-8")

    with pytest.raises(
        ValidationError,
        match=rf"{table}.*row 2, column 'region'.*must not include a record or file selector",
    ):
        read_records_table(str(table))


def test_records_table_qualifies_regions_after_order_sorting(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\tregion\torder\n"
        "b.gbk\t200..300\t2\n"
        "a.gbk\t100-150:rc\t1\n",
        encoding="utf-8",
    )

    parsed = read_records_table(str(table))

    assert [Path(path).name for path in parsed.gbk_files] == ["a.gbk", "b.gbk"]
    assert parsed.row_scoped_region_specs() == ["#1:100-150:rc", "#2:200..300"]


def test_records_table_explicit_large_order_sorts_before_missing(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\torder\nmissing.gbk\t\nexplicit.gbk\t1000000001\n",
        encoding="utf-8",
    )

    parsed = read_records_table(str(table))

    assert [Path(path).name for path in parsed.gbk_files] == ["explicit.gbk", "missing.gbk"]


def test_records_table_order_sort_is_stable_for_equal_and_missing_values(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\torder\n"
        "equal_a.gbk\t7\n"
        "missing_a.gbk\t\n"
        "equal_b.gbk\t7\n"
        "missing_b.gbk\t\n",
        encoding="utf-8",
    )

    parsed = read_records_table(str(table))

    assert [Path(path).name for path in parsed.gbk_files] == [
        "equal_a.gbk",
        "equal_b.gbk",
        "missing_a.gbk",
        "missing_b.gbk",
    ]
