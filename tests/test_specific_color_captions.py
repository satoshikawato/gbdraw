"""Color-rule captions: kept as written, with row keys and drawn names derived from them."""

import json
import xml.etree.ElementTree as ET
from pathlib import Path
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
from gbdraw.features.colors import get_color_with_info
from gbdraw.legend.table import _legend_fill_identity, rule_legend_keys


def table(rows):
    return pd.DataFrame(
        rows, columns=["feature_type", "qualifier_key", "value", "color", "caption"]
    )


RULE_KEYS = json.loads(
    (Path(__file__).parent / "fixtures" / "legend_rule_keys.json").read_text(encoding="utf-8")
)["cases"]


def legend_rows(svg_path):
    """Each legend row of an SVG by key: its name, the text it draws, and its swatch fill."""
    rows = {}
    for node in ET.fromstring(svg_path.read_text()).iter():
        key = node.get("data-legend-key")
        if key is None:
            continue
        text = "".join(
            child.text or "" for child in node.iter() if child.tag.endswith("text")
        )
        fill = next(
            (child.get("fill") for child in node.iter() if child.tag.endswith("path")),
            None,
        )
        rows[key] = (node.get("data-legend-name"), text, fill)
    return rows


@pytest.mark.parametrize("case", RULE_KEYS, ids=[case["name"] for case in RULE_KEYS])
def test_rule_row_keys_answer_the_shared_fixture(case):
    """The Web app allocates the same keys (tests/web/specific-color-actions.test.mjs)."""
    rules = [(rule["cap"], rule["color"]) for rule in case["rules"]]
    keys = rule_legend_keys(rules, dict(case["renderer_rows"]))
    assert [
        keys.get((caption, _legend_fill_identity(color)), "") for caption, color in rules
    ] == case["keys"]
    assert rule_legend_keys(rules[::-1], dict(case["renderer_rows"])) == keys


@pytest.mark.parametrize("mode", ["circular", "grid", "batch", "linear"])
def test_a_caption_in_two_colors_is_kept_and_its_rows_are_numbered(mode, tmp_path):
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
    assert canonical.caption.tolist() == ["Shared", "Shared", "", "Shared"]
    assert path.read_bytes() == original
    rows = legend_rows(result.output_paths[0])
    # The key keeps the color ID as an identity; the name and the text never carry it.
    assert rows["Shared [#112233]"] == ("Shared", "Shared (1)", "#112233")
    assert rows["Shared [#445566]"] == ("Shared", "Shared (2)", "#445566")
    assert "Shared [#778899]" not in rows
    assert "" not in rows
    assert not any("[#" in text or "[#" in name for name, text, _ in rows.values())
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
    ) == ("#112233", "Shared")
    assert get_color_with_info(
        record.features[1],
        prepared.specific_color_rules,
        prepared.default_color_map,
        "record",
    ) == ("#445566", "Shared")


@pytest.mark.parametrize("caption", ["CDS", "GC content", "GC skew (+)"])
def test_a_caption_in_two_colors_naming_a_renderer_row_numbers_all_three_rows(
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
    rows = legend_rows(result.output_paths[0])
    for name in ["CDS", "GC content", "GC skew (+)"]:
        assert rows[name][0] == name
    for color in ["#112233", "#445566"]:
        assert rows[f"{caption} [{color}]"][::2] == (caption, color)
    assert sorted(text for name, text, _ in rows.values() if name == caption) == [
        f"{caption} (1)",
        f"{caption} (2)",
        f"{caption} (3)",
    ]
    assert all(text == name for name, text, _ in rows.values() if name != caption)
    assert result.request.options.colors.color_table.caption.tolist() == [caption, caption]


@pytest.mark.parametrize("rule_color", ["#ff0000", "red"])
def test_rule_caption_equal_to_a_default_type_keeps_the_rule_color(rule_color, tmp_path):
    """A rule named like a renderer row of another color keeps its own row; both are numbered."""
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
    rows = legend_rows(result.output_paths[0])
    defaults = load_default_colors("", palette="default")
    default_rrna = str(
        defaults.loc[defaults["feature_type"] == "rRNA", "color"].iloc[0]
    )

    # The rule's row comes first in the legend; both rows of the name are numbered.
    assert rows["rRNA [#ff0000]"] == ("rRNA", "rRNA (1)", rule_color)
    assert rows["rRNA"] == ("rRNA", "rRNA (2)", default_rrna)


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
    assert legend.names()["rRNA"] == ("rRNA", "rRNA")
