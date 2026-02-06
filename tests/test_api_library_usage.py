from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO

from gbdraw.api import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.config.toml import load_config_toml
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.io.genome import load_gbks
from gbdraw.render.export import save_figure


SELECTED_FEATURES = [
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "repeat_region",
]


@pytest.mark.circular
def test_api_circular_minimal(examples_dir: Path, temp_output_dir: Path) -> None:
    record_path = examples_dir / "MjeNMV.gb"
    record = next(SeqIO.parse(str(record_path), "genbank"))

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    default_colors = load_default_colors("", palette="default")
    color_table = read_color_table("")

    output_prefix = temp_output_dir / "api_circular_minimal"
    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        output_prefix=str(output_prefix),
        legend="right",
    )

    save_figure(canvas, ["svg"])
    output_svg = output_prefix.with_suffix(".svg")

    assert output_svg.exists()
    assert output_svg.stat().st_size > 0


@pytest.mark.linear
def test_api_linear_minimal(examples_dir: Path, temp_output_dir: Path) -> None:
    gbk_files = [
        str(examples_dir / "MjeNMV.gb"),
        str(examples_dir / "MelaMJNV.gb"),
    ]
    records = load_gbks(gbk_files, mode="linear", load_comparison=False)

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    default_colors = load_default_colors("", palette="default", load_comparison=False)
    color_table = read_color_table("")

    output_prefix = temp_output_dir / "api_linear_minimal"
    canvas = assemble_linear_diagram_from_records(
        records,
        blast_files=None,
        config_dict=config_dict,
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        output_prefix=str(output_prefix),
        legend="right",
    )

    save_figure(canvas, ["svg"])
    output_svg = output_prefix.with_suffix(".svg")

    assert output_svg.exists()
    assert output_svg.stat().st_size > 0
