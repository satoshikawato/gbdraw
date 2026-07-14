from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.api.diagram as api_diagram_module
from gbdraw.api import (
    CircularMultiRecordOptions,
    apply_config_overrides,
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
    load_gff_fasta,
)
from gbdraw.api.options import DiagramOptions
from gbdraw.analysis.collinearity import CollinearityResult
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
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
    blocks = re.findall(r"```python\n(.*?)\n```", documentation, re.DOTALL)
    assert blocks
    monkeypatch.setenv("GBDRAW_EXAMPLE_GBK", str(examples_dir / "MjeNMV.gb"))
    monkeypatch.setenv("GBDRAW_EXAMPLES_DIR", str(examples_dir))
    monkeypatch.setenv(
        "GBDRAW_TEST_INPUTS_DIR",
        str(Path(__file__).parent / "test_inputs"),
    )
    monkeypatch.setenv("GBDRAW_API_OUTPUT_DIR", str(temp_output_dir))

    namespace: dict[str, object] = {}
    for index, block in enumerate(blocks, start=1):
        exec(compile(block, f"docs/PYTHON_API.md:block-{index}", "exec"), namespace)

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

    assert captured["collinearity_anchor_mode"] == "all"
    assert captured["collinearity_search_scope"] == "all"


@pytest.mark.parametrize(
    ("requested", "expected"),
    [
        ("all", "all"),
        ("one_to_one", "one_to_one"),
        ("top1", "one_to_one"),
        ("rbh", "rbh"),
    ],
)
def test_linear_assembler_forwards_normalized_collinearity_anchor_mode(
    requested: str,
    expected: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_build(*_args, **kwargs):
        captured.update(kwargs)
        return CollinearityResult(blocks=())

    monkeypatch.setattr(
        api_diagram_module,
        "build_orthogroup_collinearity_blocks",
        fake_build,
    )
    records = [
        SeqRecord(Seq("ATGC" * 25), id="rec1", description="record 1"),
        SeqRecord(Seq("ATGC" * 25), id="rec2", description="record 2"),
    ]

    api_diagram_module.assemble_linear_diagram_from_records(
        records,
        protein_blastp_mode="collinear",
        collinearity_anchor_mode=requested,
        legend="none",
    )

    assert captured["edge_mode"] == expected


def test_typed_config_override_preserves_label_filtering_dataframes() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    whitelist = pd.DataFrame(
        [["CDS", "product", "polymerase"]],
        columns=["feature_type", "qualifier", "keyword"],
    )
    qualifier_priority = pd.DataFrame(
        [["CDS", "gene,product"]],
        columns=["feature_type", "priorities"],
    )
    label_override = pd.DataFrame(
        [["rec1", "CDS", "gene", "pol", "polymerase"]],
        columns=["record_id", "feature_type", "qualifier", "value", "label_text"],
    )
    filtering = config_dict["labels"]["filtering"]
    filtering["whitelist_df"] = whitelist
    filtering["qualifier_priority_df"] = qualifier_priority
    filtering["label_override_df"] = label_override
    filtering["extension_key"] = {"keep": True}

    typed = GbdrawConfig.from_dict(config_dict)
    updated = apply_config_overrides(typed, {"show_labels": True})
    updated_filtering = updated.labels.filtering.as_dict()

    pd.testing.assert_frame_equal(updated_filtering["whitelist_df"], whitelist)
    pd.testing.assert_frame_equal(
        updated_filtering["qualifier_priority_df"], qualifier_priority
    )
    pd.testing.assert_frame_equal(updated_filtering["label_override_df"], label_override)
    assert updated_filtering["extension_key"] == {"keep": True}


@pytest.mark.parametrize("builder_name", ["circular", "linear"])
def test_diagram_options_attach_explicit_label_tables(
    builder_name: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    target = (
        "assemble_circular_diagram_from_record"
        if builder_name == "circular"
        else "assemble_linear_diagram_from_records"
    )
    monkeypatch.setattr(api_diagram_module, target, fake_assemble)
    whitelist = pd.DataFrame(
        [["CDS", "product", "polymerase"]],
        columns=["feature_type", "qualifier", "keyword"],
    )
    priority = pd.DataFrame(
        [["CDS", "gene,product"]],
        columns=["feature_type", "priorities"],
    )
    override = pd.DataFrame(
        [["rec1", "CDS", "gene", "pol", "polymerase"]],
        columns=["record_id", "feature_type", "qualifier", "value", "label_text"],
    )
    options = DiagramOptions(
        label_whitelist_table=whitelist,
        qualifier_priority_table=priority,
        label_override_table=override,
    )

    if builder_name == "circular":
        api_diagram_module.build_circular_diagram(
            SeqRecord(Seq("ATGC"), id="rec1"), options=options
        )
    else:
        api_diagram_module.build_linear_diagram([], options=options)

    cfg = captured["cfg"]
    filtering = cfg.labels.filtering.as_dict()
    pd.testing.assert_frame_equal(filtering["whitelist_df"], whitelist)
    pd.testing.assert_frame_equal(filtering["qualifier_priority_df"], priority)
    pd.testing.assert_frame_equal(filtering["label_override_df"], override)


def test_diagram_options_reject_label_table_and_file_together() -> None:
    table = pd.DataFrame(
        [["CDS", "product", "polymerase"]],
        columns=["feature_type", "qualifier", "keyword"],
    )

    with pytest.raises(ValidationError, match="label_whitelist"):
        api_diagram_module.build_linear_diagram(
            [],
            options=DiagramOptions(
                label_whitelist_table=table,
                label_whitelist_file="whitelist.tsv",
            ),
        )


def test_circular_multi_builder_forwards_layout_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}
    records = [
        SeqRecord(Seq("ATGC"), id="rec1"),
        SeqRecord(Seq("ATGCATGC"), id="rec2"),
    ]

    def fake_assemble(received_records, **kwargs):
        captured["records"] = received_records
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(
        api_diagram_module,
        "assemble_circular_diagram_from_records",
        fake_assemble,
    )

    api_diagram_module.build_circular_multi_diagram(
        records,
        options=DiagramOptions(plot_title="Two records"),
        layout=CircularMultiRecordOptions(
            multi_record_size_mode="equal",
            multi_record_min_radius_ratio=0.7,
            multi_record_column_gap_ratio=0.2,
            multi_record_row_gap_ratio=0.1,
            multi_record_positions=["#2@1", "#1@2"],
        ),
    )

    assert captured["records"] == records
    assert captured["multi_record_size_mode"] == "equal"
    assert captured["multi_record_min_radius_ratio"] == 0.7
    assert captured["multi_record_column_gap_ratio"] == 0.2
    assert captured["multi_record_row_gap_ratio"] == 0.1
    assert captured["multi_record_positions"] == ["#2@1", "#1@2"]
    assert captured["plot_title"] == "Two records"


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
