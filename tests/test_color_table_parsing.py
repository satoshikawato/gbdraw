from __future__ import annotations

from pathlib import Path

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import read_color_table


def test_read_color_table_accepts_four_column_rows(tmp_path: Path) -> None:
    color_table_path = tmp_path / "specific_colors.tsv"
    color_table_path.write_text("CDS\tproduct\tATPase\t#ff0000\n", encoding="utf-8")

    df = read_color_table(str(color_table_path))

    assert df is not None
    assert df.to_dict(orient="records") == [
        {
            "feature_type": "CDS",
            "qualifier_key": "product",
            "value": "ATPase",
            "color": "#ff0000",
            "caption": "",
        }
    ]


def test_read_color_table_normalizes_empty_caption_column(tmp_path: Path) -> None:
    color_table_path = tmp_path / "specific_colors.tsv"
    color_table_path.write_text(
        "CDS\tproduct\tATPase\t#ff0000\t\n",
        encoding="utf-8",
    )

    df = read_color_table(str(color_table_path))

    assert df is not None
    assert df.iloc[0].to_dict() == {
        "feature_type": "CDS",
        "qualifier_key": "product",
        "value": "ATPase",
        "color": "#ff0000",
        "caption": "",
    }


def test_read_color_table_rejects_missing_required_columns(tmp_path: Path) -> None:
    color_table_path = tmp_path / "specific_colors.tsv"
    color_table_path.write_text("CDS\tproduct\tATPase\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="Missing values"):
        read_color_table(str(color_table_path))


@pytest.mark.regression
@pytest.mark.circular
def test_circular_cli_accepts_existing_four_column_color_table(
    gbdraw_runner,
    examples_dir: Path,
    temp_output_dir: Path,
) -> None:
    gbk_file = examples_dir / "MjeNMV.gb"
    color_table = examples_dir / "feature_specific_color_table.txt"

    returncode, output, svg_path = gbdraw_runner.run_circular(
        [gbk_file],
        "circular_four_column_color_table",
        temp_output_dir,
        extra_args=["-t", str(color_table), "--legend", "right"],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    assert svg_path.exists()
    assert svg_path.stat().st_size > 0
