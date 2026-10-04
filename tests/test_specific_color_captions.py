"""Caption allocation and real native rules/legend/request parity."""

import json
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import (
    CircularDiagramOptions,
    LinearDiagramOptions,
    CircularMultiRecordOptions,
    ColorOptions,
)
from gbdraw.api.prepared import resolve_feature_inputs
from gbdraw.api.request_render import render_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    CircularBatchRequest,
    LinearDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.features.colors import (
    normalize_specific_color_captions,
    get_color_with_info,
)
from gbdraw.web_support.rule_matching import evaluate_rules_json


def table(rows):
    return pd.DataFrame(
        rows, columns=["feature_type", "qualifier_key", "value", "color", "caption"]
    )


@pytest.mark.parametrize(
    "colors,captions",
    [
        (["#112233", "#112233"], ["Shared", "Shared"]),
        (["red", "#FF0000"], ["Shared", "Shared"]),
        (["#123", "#445566"], ["Shared [#112233]", "Shared [#445566]"]),
    ],
)
def test_shared_caption_color_domain(colors, captions):
    source = table(
        [
            ["CDS", "gene", gene, color, "Shared"]
            for gene, color in zip(["a", "b"], colors)
        ]
    )
    result = normalize_specific_color_captions(source)
    assert list(result.caption) == captions
    assert source.caption.tolist() == ["Shared", "Shared"]
    pd.testing.assert_frame_equal(normalize_specific_color_captions(result), result)


def test_literal_reservations_reorder_blank_and_source_identity():
    source = table(
        [
            ["CDS", "gene", "a", "#123", "Shared"],
            ["CDS", "gene", "b", "#445566", "Shared"],
            ["CDS", "gene", "c", "#abcdef", "Shared [#112233]"],
            ["CDS", "gene", "d", "#abcdef", "Shared [#112233] (2)"],
            ["CDS", "gene", "empty", "#445566", ""],
            ["CDS", "gene", "blank", "#123", ""],
        ]
    )
    result = normalize_specific_color_captions(source)
    assert result.caption.tolist() == [
        "Shared [#112233] (3)",
        "Shared [#445566]",
        "Shared [#112233]",
        "Shared [#112233] (2)",
        "",
        "",
    ]
    assert result.drop(columns="caption").equals(source.drop(columns="caption"))
    shuffled = normalize_specific_color_captions(source.iloc[::-1])
    assert dict(zip(shuffled.value, shuffled.caption)) == dict(
        zip(result.value, result.caption)
    )
    pd.testing.assert_frame_equal(normalize_specific_color_captions(result), result)
    browser = [
        dict(
            feat=row.feature_type,
            qual=row.qualifier_key,
            val=row.value,
            color=row.color,
            cap=row.caption,
            fromFile=True,
        )
        for row in source.itertuples()
    ]
    normalized = json.loads(
        evaluate_rules_json("[]", json.dumps(browser), "color-captions")
    )["rules"]
    assert [r["cap"] for r in normalized] == result.caption.tolist()
    assert all(r["fromFile"] for r in normalized)


@pytest.mark.parametrize("mode", ["circular", "grid", "batch", "linear"])
def test_real_render_request_legend_and_colors_use_canonical_rows(mode, tmp_path):
    record = SeqRecord(
        Seq("ATGC" * 250),
        id="record",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    record.features = [
        SeqFeature(
            SimpleLocation(i * 200 + 10, i * 200 + 160, strand=1),
            type="CDS",
            qualifiers={"gene": [gene]},
        )
        for i, gene in enumerate(["a", "b", "c", "d"])
    ]
    source = table(
        [
            ["CDS", "gene", "^a$", "#112233", "Shared"],
            ["CDS", "gene", "^b$", "#445566", "Shared"],
            ["CDS", "gene", "^c$", "#abcdef", ""],
            ["CDS", "gene", "^absent$", "#778899", "Shared"],
        ]
    )
    path = tmp_path / "original.tsv"
    source.to_csv(path, sep="\t", index=False, header=False)
    original = path.read_bytes()
    options_type = LinearDiagramOptions if mode == "linear" else CircularDiagramOptions
    options = options_type(colors=ColorOptions(color_table_file=str(path)))
    records = (RecordInput(source=InMemoryRecordSource(record)),)
    output = RenderOutputRequest(
        formats=("interactive_svg",), output_directory=tmp_path, output_prefix=mode
    )
    if mode == "batch":
        request = CircularBatchRequest(
            records=records, options=options, outputs=(output,)
        )
    elif mode == "grid":
        request = CircularDiagramRequest(
            records=records * 2,
            options=options,
            output=output,
            layout=CircularMultiRecordOptions(),
        )
    else:
        request_type = (
            LinearDiagramRequest if mode == "linear" else CircularDiagramRequest
        )
        request = request_type(records=records, options=options, output=output)
    result = render_request(request, include_feature_catalog=True)
    canonical = result.request.options.colors.color_table
    assert canonical.caption.tolist() == [
        "Shared [#112233]",
        "Shared [#445566]",
        "",
        "Shared [#778899]",
    ]
    assert result.request.options.colors.color_table_file is None
    assert path.read_bytes() == original
    svg = ET.fromstring(result.output_paths[0].read_text())
    entries = {
        node.get("data-legend-key"): node
        for node in svg.iter()
        if node.get("data-legend-key")
    }
    assert set(entries) >= {"Shared [#112233]", "Shared [#445566]"}
    assert "Shared [#778899]" not in entries
    assert "" not in entries
    for caption, color in [
        ("Shared [#112233]", "#112233"),
        ("Shared [#445566]", "#445566"),
    ]:
        assert any(node.get("fill") == color for node in entries[caption].iter())
    prepared = resolve_feature_inputs(
        color_table=source,
        default_colors=pd.DataFrame(
            [["CDS", "#cccccc"]], columns=["feature_type", "color"]
        ),
        feature_visibility_table=None,
    )
    assert get_color_with_info(
        record.features[0],
        prepared.specific_color_rules,
        prepared.default_color_map,
        "record",
    ) == ("#112233", "Shared [#112233]")
    assert get_color_with_info(
        record.features[1],
        prepared.specific_color_rules,
        prepared.default_color_map,
        "record",
    ) == ("#445566", "Shared [#445566]")


@pytest.mark.parametrize("caption", ["CDS", "GC content", "GC skew (+)"])
def test_normalized_specific_rows_preserve_default_and_numeric_legend_names(
    caption, tmp_path
):
    record = SeqRecord(
        Seq("ATGC" * 250),
        id="legend",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    record.features = [
        SeqFeature(
            SimpleLocation(i * 200 + 10, i * 200 + 160, strand=1),
            type=feature_type,
            qualifiers={"gene": [gene]},
        )
        for i, (gene, feature_type) in enumerate(
            [("a", "tRNA"), ("b", "tRNA"), ("c", "CDS")]
        )
    ]
    rows = table(
        [
            ["tRNA", "gene", gene, color, caption]
            for gene, color in [("a", "#112233"), ("b", "#445566")]
        ]
    )
    result = render_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(colors=ColorOptions(color_table=rows)),
            output=RenderOutputRequest(
                formats=("interactive_svg",),
                output_directory=tmp_path,
                output_prefix="legend",
            ),
        )
    )
    svg = ET.fromstring(result.output_paths[0].read_text())
    entries = {
        node.get("data-legend-key"): node
        for node in svg.iter()
        if node.get("data-legend-key")
    }
    for name in ["CDS", "GC content", "GC skew (+)"]:
        assert name in entries
    for color in ["#112233", "#445566"]:
        assert any(
            node.get("fill") == color for node in entries[f"{caption} [{color}]"].iter()
        )
    assert result.request.options.colors.color_table.caption.tolist() == [
        f"{caption} [#112233]",
        f"{caption} [#445566]",
    ]


@pytest.mark.parametrize("rule_color", ["#ff0000", "red"])
def test_rule_caption_equal_to_a_default_type_keeps_the_rule_color(rule_color, tmp_path):
    """N-06: a rule caption equal to a default legend name keeps its own row."""
    from gbdraw.io.colors import load_default_colors

    record = SeqRecord(
        Seq("ATGC" * 250),
        id="legend",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    record.features = [
        SeqFeature(
            SimpleLocation(i * 200 + 10, i * 200 + 90, strand=1),
            type=feature_type,
            qualifiers={"product": [product]},
        )
        for i, (feature_type, product) in enumerate(
            [("tRNA", "tRNA-Phe"), ("rRNA", "12S ribosomal RNA"), ("tRNA", "tRNA-Val")]
        )
    ]
    rows = table([["tRNA", "product", ".*", rule_color, "rRNA"]])
    result = render_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                colors=ColorOptions(color_table=rows),
                selected_features_set=("tRNA", "rRNA"),
            ),
            output=RenderOutputRequest(
                formats=("svg",), output_directory=tmp_path, output_prefix="legend"
            ),
        )
    )
    svg = ET.fromstring(result.output_paths[0].read_text())
    fills = {
        node.get("data-legend-key"): next(
            (
                child.get("fill").lower()
                for child in node.iter()
                if child.get("fill") not in (None, "none")
            ),
            None,
        )
        for node in svg.iter()
        if node.get("data-legend-key")
    }
    defaults = load_default_colors("", palette="default")
    default_rrna = str(
        defaults.loc[defaults["feature_type"] == "rRNA", "color"].iloc[0]
    ).lower()

    assert fills["rRNA"] == default_rrna
    assert fills["rRNA [#ff0000]"] == ("#ff0000" if rule_color.startswith("#") else "red")
    assert list(fills).count("rRNA") == 1


def test_rule_caption_equal_to_a_default_type_with_the_same_color_shares_one_row():
    from gbdraw.legend.table import prepare_legend_table

    feature_config = SimpleNamespace(
        color_table=table([["tRNA", "product", ".*", "#71ee7d", "rRNA"]]),
        default_colors=pd.DataFrame(
            [["tRNA", "#ff00ff"], ["rRNA", "#71EE7D"]], columns=["feature_type", "color"]
        ),
        block_stroke_color="black",
        block_stroke_width=1.0,
    )
    numeric = SimpleNamespace(
        stroke_color="none", stroke_width=0, high_fill_color="#111111",
        low_fill_color="#222222", dinucleotide="GC",
    )
    legend = prepare_legend_table(
        numeric,
        numeric,
        feature_config,
        ["tRNA", "rRNA"],
        used_color_rules={("rRNA", "#71ee7d")},
        default_used_features={"rRNA"},
        show_gc=False,
        show_skew=False,
        show_depth=False,
    )

    assert list(legend) == ["rRNA"]
