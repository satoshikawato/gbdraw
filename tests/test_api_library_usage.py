from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

import gbdraw.api.diagram as api_diagram_module
from gbdraw.api import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
    load_gff_fasta,
)
from gbdraw.api.options import DiagramOptions
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
def test_documented_python_api_example_runs(
    examples_dir: Path,
    temp_output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    documentation = (Path(__file__).parents[1] / "docs" / "PYTHON_API.md").read_text(
        encoding="utf-8"
    )
    match = re.search(r"```python\n(.*?)\n```", documentation, re.DOTALL)
    assert match is not None
    monkeypatch.setenv("GBDRAW_EXAMPLE_GBK", str(examples_dir / "MjeNMV.gb"))
    monkeypatch.setenv("GBDRAW_API_OUTPUT_DIR", str(temp_output_dir))

    exec(compile(match.group(1), "docs/PYTHON_API.md", "exec"), {})

    assert (temp_output_dir / "api_circular.svg").exists()


@pytest.mark.linear
def test_documented_gff3_fasta_fixture_preserves_records_and_cds() -> None:
    fixture_dir = Path(__file__).parents[1] / "examples" / "gff3_lambda"
    records = load_gff_fasta(
        [str(fixture_dir / "lambda_two_contigs.gff3")],
        [str(fixture_dir / "lambda_two_contigs.fna")],
        mode="linear",
        selected_features_set=["CDS", "gene"],
    )

    assert [record.id for record in records] == ["lambda_left", "lambda_right"]
    cds_features = [
        feature for record in records for feature in record.features if feature.type == "CDS"
    ]
    assert len(cds_features) == 45
    assert {feature.location.strand for feature in cds_features} == {1, -1}
    assert all(feature.qualifiers.get("translation") for feature in cds_features)


@pytest.mark.linear
def test_api_diagram_options_forward_collinearity_search_scope(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(records, **kwargs):
        captured["records"] = records
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    api_diagram_module.build_linear_diagram(
        [],
        options=DiagramOptions(
            protein_blastp_mode="collinear",
            collinearity_anchor_mode="all",
            collinearity_search_scope="all",
        ),
    )

    assert captured["collinearity_anchor_mode"] == "rbh"
    assert captured["collinearity_search_scope"] == "all"


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
def test_api_linear_gene_specific_rule_uses_default_fallback_for_legend(
    test_inputs_dir: Path,
    temp_output_dir: Path,
) -> None:
    record_path = test_inputs_dir / "HmmtDNA.gbk"
    record = next(SeqIO.parse(str(record_path), "genbank"))

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    default_colors = load_default_colors("", palette="default", load_comparison=False)
    color_table = pd.DataFrame(
        [["gene", "gene", "^TRNF$", "#b56576", "TRNF gene"]],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )

    output_prefix = temp_output_dir / "api_linear_gene_legend_fallback"
    canvas = assemble_linear_diagram_from_records(
        [record],
        blast_files=None,
        config_dict=config_dict,
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=["gene"],
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
