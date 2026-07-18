from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw import CircularOptions
import gbdraw.api.diagram as api_diagram_module
from gbdraw.api import (
    CircularMultiRecordOptions,
    apply_config_overrides,
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
    load_gff_fasta,
)
from gbdraw.api.options import ColorOptions, DiagramOptions, OutputOptions, TrackOptions
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

    documented_record = namespace["record"]
    documented_options = namespace["options"]
    assert isinstance(documented_record, SeqRecord)
    assert isinstance(documented_options, CircularOptions)
    assert set(documented_options.features.types) <= {
        feature.type for feature in documented_record.features
    }
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
            collinearity_unit_mode="nt",
            collinear_max_paralog_links_per_orthogroup=7,
        ),
    )

    assert captured["collinearity_anchor_mode"] == "all"
    assert captured["collinearity_search_scope"] == "all"
    assert captured["collinearity_unit_mode"] == "nt"
    assert captured["collinear_max_paralog_links_per_orthogroup"] == 7


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


@pytest.mark.parametrize("builder_name", ["circular", "circular_multi", "linear"])
def test_diagram_options_attach_explicit_label_tables(
    builder_name: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    target = {
        "circular": "assemble_circular_diagram_from_record",
        "circular_multi": "assemble_circular_diagram_from_records",
        "linear": "assemble_linear_diagram_from_records",
    }[builder_name]
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
    elif builder_name == "circular_multi":
        api_diagram_module.build_circular_multi_diagram(
            [SeqRecord(Seq("ATGC"), id="rec1")], options=options
        )
    else:
        api_diagram_module.build_linear_diagram([], options=options)

    cfg = captured["cfg"]
    filtering = cfg.labels.filtering.as_dict()
    pd.testing.assert_frame_equal(filtering["whitelist_df"], whitelist)
    pd.testing.assert_frame_equal(filtering["qualifier_priority_df"], priority)
    pd.testing.assert_frame_equal(filtering["label_override_df"], override)


def test_diagram_options_attach_explicit_label_files(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}
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

    monkeypatch.setattr(api_diagram_module, "read_filter_list_file", lambda _path: whitelist)
    monkeypatch.setattr(
        api_diagram_module,
        "read_qualifier_priority_file",
        lambda _path: priority,
    )
    monkeypatch.setattr(
        api_diagram_module,
        "read_label_override_file",
        lambda _path: override,
    )

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(
        api_diagram_module,
        "assemble_linear_diagram_from_records",
        fake_assemble,
    )
    api_diagram_module.build_linear_diagram(
        [],
        options=DiagramOptions(
            label_whitelist_file="whitelist.tsv",
            qualifier_priority_file="priority.tsv",
            label_override_file="override.tsv",
        ),
    )

    filtering = captured["cfg"].labels.filtering.as_dict()
    pd.testing.assert_frame_equal(filtering["whitelist_df"], whitelist)
    pd.testing.assert_frame_equal(filtering["qualifier_priority_df"], priority)
    pd.testing.assert_frame_equal(filtering["label_override_df"], override)


@pytest.mark.parametrize(
    ("table_field", "file_field"),
    [
        ("label_whitelist_table", "label_whitelist_file"),
        ("qualifier_priority_table", "qualifier_priority_file"),
        ("label_override_table", "label_override_file"),
    ],
)
def test_diagram_options_reject_label_table_and_file_together(
    table_field: str,
    file_field: str,
) -> None:
    table = pd.DataFrame(
        [["CDS", "product", "polymerase"]],
        columns=["feature_type", "qualifier", "keyword"],
    )

    with pytest.raises(ValidationError, match=table_field.removesuffix("_table")):
        api_diagram_module.build_linear_diagram(
            [],
            options=DiagramOptions(**{table_field: table, file_field: "labels.tsv"}),
        )


_FORWARDING_TABLE = pd.DataFrame({"value": [1]})
_SHARED_FORWARDING_CASES = [
    ("config", {"extension": {"enabled": True}}, "config_dict", {"extension": {"enabled": True}}),
    ("config_overrides", {"show_labels": True}, "config_overrides", {"show_labels": True}),
    ("selected_features_set", ["gene"], "selected_features_set", ["gene"]),
    ("feature_table", _FORWARDING_TABLE, "feature_table", _FORWARDING_TABLE),
    ("feature_table_file", "feature.tsv", "feature_table_file", "feature.tsv"),
    ("feature_visibility_table", _FORWARDING_TABLE, "feature_visibility_table", _FORWARDING_TABLE),
    (
        "feature_visibility_table_file",
        "visibility.tsv",
        "feature_visibility_table_file",
        "visibility.tsv",
    ),
    ("feature_shapes", {"gene": "rectangle"}, "feature_shapes", {"gene": "rectangle"}),
    ("dinucleotide", "AT", "dinucleotide", "AT"),
    ("window", 111, "window", 111),
    ("step", 17, "step", 17),
    ("depth_window", 91, "depth_window", 91),
    ("depth_step", 13, "depth_step", 13),
    ("depth_table", _FORWARDING_TABLE, "depth_table", _FORWARDING_TABLE),
    ("depth_file", "depth.tsv", "depth_file", "depth.tsv"),
    (
        "depth_track_tables",
        [[_FORWARDING_TABLE]],
        "depth_track_tables",
        [[_FORWARDING_TABLE]],
    ),
    ("depth_track_files", [["depth-1.tsv"]], "depth_track_files", [["depth-1.tsv"]]),
    ("depth_track_labels", ["coverage"], "depth_track_labels", ["coverage"]),
    ("depth_track_colors", ["#123456"], "depth_track_colors", ["#123456"]),
    (
        "depth_track_large_tick_intervals",
        [25],
        "depth_track_large_tick_intervals",
        [25],
    ),
    (
        "depth_track_small_tick_intervals",
        [5],
        "depth_track_small_tick_intervals",
        [5],
    ),
    ("depth_track_tick_font_sizes", [9], "depth_track_tick_font_sizes", [9]),
    ("plot_title", "Forwarded title", "plot_title", "Forwarded title"),
    ("plot_title_font_size", 27.0, "plot_title_font_size", 27.0),
    ("evalue", 1e-12, "evalue", 1e-12),
    ("bitscore", 123.0, "bitscore", 123.0),
    ("identity", 88.0, "identity", 88.0),
    ("alignment_length", 42, "alignment_length", 42),
]


def _call_high_level_builder(builder_name: str, options: DiagramOptions) -> object:
    record = SeqRecord(Seq("ATGC"), id="rec1")
    if builder_name == "circular":
        return api_diagram_module.build_circular_diagram(record, options=options)
    if builder_name == "circular_multi":
        return api_diagram_module.build_circular_multi_diagram([record], options=options)
    return api_diagram_module.build_linear_diagram([record], options=options)


def _assert_forwarded_value(actual: object, expected: object) -> None:
    if isinstance(expected, pd.DataFrame):
        assert actual is expected
    elif isinstance(expected, list):
        assert isinstance(actual, list)
        assert len(actual) == len(expected)
        for actual_item, expected_item in zip(actual, expected):
            _assert_forwarded_value(actual_item, expected_item)
    else:
        assert actual == expected


@pytest.mark.parametrize("builder_name", ["circular", "circular_multi", "linear"])
@pytest.mark.parametrize(
    ("field_name", "value", "assembler_name", "expected"),
    _SHARED_FORWARDING_CASES,
    ids=[case[0] for case in _SHARED_FORWARDING_CASES],
)
def test_diagram_options_forward_shared_non_default_values(
    builder_name: str,
    field_name: str,
    value: object,
    assembler_name: str,
    expected: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    target = {
        "circular": "assemble_circular_diagram_from_record",
        "circular_multi": "assemble_circular_diagram_from_records",
        "linear": "assemble_linear_diagram_from_records",
    }[builder_name]
    monkeypatch.setattr(api_diagram_module, target, fake_assemble)

    _call_high_level_builder(
        builder_name,
        DiagramOptions(**{field_name: value}),
    )

    _assert_forwarded_value(captured[assembler_name], expected)


@pytest.mark.parametrize("builder_name", ["circular", "circular_multi", "linear"])
def test_diagram_option_bundles_forward_non_default_values(
    builder_name: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}
    color_table = pd.DataFrame({"color": ["#123456"]})
    default_colors = pd.DataFrame({"color": ["#abcdef"]})

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    target = {
        "circular": "assemble_circular_diagram_from_record",
        "circular_multi": "assemble_circular_diagram_from_records",
        "linear": "assemble_linear_diagram_from_records",
    }[builder_name]
    monkeypatch.setattr(api_diagram_module, target, fake_assemble)

    options = DiagramOptions(
        colors=ColorOptions(
            color_table=color_table,
            default_colors=default_colors,
            default_colors_palette="ajisai",
        ),
        tracks=TrackOptions(
            circular_track_slots=["features:features"],
            circular_track_axis_index=1,
            linear_track_slots=["features:features"],
            linear_track_axis_index=1,
            center_reserved_radius=0.2,
        ),
        output=OutputOptions(
            output_prefix="forwarded",
            legend="left",
            plot_title_position="top",
        ),
    )
    _call_high_level_builder(builder_name, options)

    assert captured["color_table"] is color_table
    assert captured["default_colors"] is default_colors
    assert captured["default_colors_palette"] == "ajisai"
    assert captured["output_prefix"] == "forwarded"
    assert captured["legend"] == "left"
    assert captured["plot_title_position"] == "top"
    if builder_name == "linear":
        assert captured["linear_track_slots"] == ["features:features"]
        assert captured["linear_track_axis_index"] == 1
    else:
        assert captured["circular_track_slots"] == ["features:features"]
        assert captured["circular_track_axis_index"] == 1
        assert captured["center_reserved_radius"] == 0.2


@pytest.mark.parametrize("builder_name", ["circular", "circular_multi", "linear"])
@pytest.mark.parametrize(
    ("field_name", "value", "assembler_name", "expected"),
    [
        ("depth_tables", [_FORWARDING_TABLE], "depth_tables", [_FORWARDING_TABLE]),
        ("depth_files", ["depth.tsv"], "depth_files", ["depth.tsv"]),
    ],
)
def test_plural_depth_options_forward_with_mode_specific_shape(
    builder_name: str,
    field_name: str,
    value: object,
    assembler_name: str,
    expected: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    target = {
        "circular": "assemble_circular_diagram_from_record",
        "circular_multi": "assemble_circular_diagram_from_records",
        "linear": "assemble_linear_diagram_from_records",
    }[builder_name]
    monkeypatch.setattr(api_diagram_module, target, fake_assemble)
    _call_high_level_builder(builder_name, DiagramOptions(**{field_name: value}))

    if builder_name == "circular":
        singular_name = "depth_table" if field_name == "depth_tables" else "depth_file"
        singular_value = expected[0]
        _assert_forwarded_value(captured[singular_name], singular_value)
    else:
        _assert_forwarded_value(captured[assembler_name], expected)


_CIRCULAR_ONLY_FORWARDING_CASES = [
    ("conservation_blast_files", ["conservation.tsv"]),
    ("conservation_dataframes", [_FORWARDING_TABLE]),
    ("conservation_reference", "subject"),
    ("conservation_labels", ["reference"]),
    ("conservation_colors", ["#123456"]),
    ("conservation_ring_width", 12.0),
    ("conservation_ring_gap", 3.0),
    ("keep_full_definition_with_plot_title", True),
    ("species", "Example species"),
    ("strain", "Example strain"),
]


@pytest.mark.parametrize("builder_name", ["circular", "circular_multi"])
@pytest.mark.parametrize(
    ("field_name", "value"),
    _CIRCULAR_ONLY_FORWARDING_CASES,
    ids=[case[0] for case in _CIRCULAR_ONLY_FORWARDING_CASES],
)
def test_diagram_options_forward_circular_only_non_default_values(
    builder_name: str,
    field_name: str,
    value: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    target = (
        "assemble_circular_diagram_from_record"
        if builder_name == "circular"
        else "assemble_circular_diagram_from_records"
    )
    monkeypatch.setattr(api_diagram_module, target, fake_assemble)
    _call_high_level_builder(builder_name, DiagramOptions(**{field_name: value}))

    _assert_forwarded_value(captured[field_name], value)


_LINEAR_ONLY_FORWARDING_CASES = [
    ("depth_track_heights", [12, 28], "depth_track_heights", [12, 28]),
    ("blast_files", ["comparison.tsv"], "blast_files", ["comparison.tsv"]),
    ("protein_comparisons", [_FORWARDING_TABLE], "protein_comparisons", [_FORWARDING_TABLE]),
    ("orthogroups", object(), "orthogroups", None),
    ("protein_blastp_mode", "pairwise", "protein_blastp_mode", "pairwise"),
    ("pairwise_match_style", "curve", "pairwise_match_style", "curve"),
    ("collinearity_blocks", object(), "collinearity_blocks", None),
    ("collinearity_params", object(), "collinearity_params", None),
    ("collinearity_unit_mode", "nt", "collinearity_unit_mode", "nt"),
    ("collinearity_anchor_mode", "top1", "collinearity_anchor_mode", "one_to_one"),
    ("collinearity_search_scope", "all", "collinearity_search_scope", "all"),
    ("collinearity_color_mode", "identity", "collinearity_color_mode", "identity"),
    ("losatp_bin", "custom-losat", "losatp_bin", "custom-losat"),
    ("ncbi_blastp_bin", "custom-blastp", "ncbi_blastp_bin", "custom-blastp"),
    ("losatp_threads", 3, "losatp_threads", 3),
    ("protein_blastp_max_hits", 8, "protein_blastp_max_hits", 8),
    ("protein_blastp_candidate_limit", 21, "protein_blastp_candidate_limit", 21),
    (
        "orthogroup_membership_mode",
        "anchor_core_v1",
        "orthogroup_membership_mode",
        "anchor_core_v1",
    ),
    ("orthogroup_member_max_hits", 9, "orthogroup_member_max_hits", 9),
    (
        "collinear_max_paralog_links_per_orthogroup",
        4,
        "collinear_max_paralog_links_per_orthogroup",
        4,
    ),
    ("align_orthogroup_feature", "anchor", "align_orthogroup_feature", "anchor"),
]


@pytest.mark.parametrize(
    ("field_name", "value", "assembler_name", "expected"),
    _LINEAR_ONLY_FORWARDING_CASES,
    ids=[case[0] for case in _LINEAR_ONLY_FORWARDING_CASES],
)
def test_diagram_options_forward_linear_only_non_default_values(
    field_name: str,
    value: object,
    assembler_name: str,
    expected: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(
        api_diagram_module,
        "assemble_linear_diagram_from_records",
        fake_assemble,
    )
    _call_high_level_builder("linear", DiagramOptions(**{field_name: value}))

    if expected is None:
        assert captured[assembler_name] is value
    else:
        _assert_forwarded_value(captured[assembler_name], expected)


_LINEAR_ONLY_WRONG_MODE_CASES = [
    (
        field_name,
        "unsupported" if field_name == "orthogroup_membership_mode" else value,
    )
    for field_name, value, _assembler_name, _expected in _LINEAR_ONLY_FORWARDING_CASES
]


@pytest.mark.parametrize("builder_name", ["circular", "circular_multi"])
@pytest.mark.parametrize(
    ("field_name", "value"),
    _LINEAR_ONLY_WRONG_MODE_CASES,
    ids=[case[0] for case in _LINEAR_ONLY_WRONG_MODE_CASES],
)
def test_circular_builders_reject_linear_only_non_default_options(
    builder_name: str,
    field_name: str,
    value: object,
) -> None:
    with pytest.raises(ValidationError, match=rf"{builder_name}.*{field_name}"):
        _call_high_level_builder(builder_name, DiagramOptions(**{field_name: value}))


@pytest.mark.parametrize(
    ("field_name", "value"),
    _CIRCULAR_ONLY_FORWARDING_CASES,
    ids=[case[0] for case in _CIRCULAR_ONLY_FORWARDING_CASES],
)
def test_linear_builder_rejects_circular_only_non_default_options(
    field_name: str,
    value: object,
) -> None:
    with pytest.raises(ValidationError, match=rf"linear.*{field_name}"):
        _call_high_level_builder("linear", DiagramOptions(**{field_name: value}))


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (DiagramOptions(depth_tables=[]), "one depth table"),
        (DiagramOptions(depth_tables=[_FORWARDING_TABLE, _FORWARDING_TABLE]), "one depth table"),
        (DiagramOptions(depth_files=[]), "one depth file"),
        (DiagramOptions(depth_files=["one.tsv", "two.tsv"]), "one depth file"),
        (
            DiagramOptions(depth_table=_FORWARDING_TABLE, depth_tables=[_FORWARDING_TABLE]),
            "singular.*plural",
        ),
        (DiagramOptions(depth_file="one.tsv", depth_files=["two.tsv"]), "singular.*plural"),
        (
            DiagramOptions(depth_table=_FORWARDING_TABLE, depth_files=["two.tsv"]),
            "singular.*plural",
        ),
        (
            DiagramOptions(depth_tables=[_FORWARDING_TABLE], depth_files=["one.tsv"]),
            "table.*file",
        ),
        (DiagramOptions(depth_table=_FORWARDING_TABLE, depth_file="one.tsv"), "table.*file"),
    ],
)
def test_single_circular_builder_rejects_ambiguous_or_lossy_depth_options(
    options: DiagramOptions,
    message: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assembly_called = False

    def fake_assemble(*_args, **_kwargs):
        nonlocal assembly_called
        assembly_called = True
        return object()

    monkeypatch.setattr(
        api_diagram_module,
        "assemble_circular_diagram_from_record",
        fake_assemble,
    )

    with pytest.raises(ValidationError, match=message):
        _call_high_level_builder("circular", options)

    assert not assembly_called


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
