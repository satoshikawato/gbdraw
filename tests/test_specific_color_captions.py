"""Color-rule captions: kept as written, with row keys and drawn names derived from them."""

import json
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularTrackOptions,
    LinearDiagramOptions,
    CircularMultiRecordOptions,
    ColorOptions,
)
from gbdraw.api.diagram import resolve_feature_inputs
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
from gbdraw.io.colors import read_color_table
from gbdraw.legend.table import (
    _legend_fill_identity,
    legend_swatch,
    rule_legend_keys,
    saved_rule_caption,
)
from gbdraw.session import (
    build_session_document,
    load_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.web_support.request_render import render_canonical_web_request


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
        # A row without a name attribute is named by its key.
        rows[key] = (node.get("data-legend-name", key), text, fill)
    return rows


@pytest.mark.parametrize("case", RULE_KEYS, ids=[case["name"] for case in RULE_KEYS])
def test_rule_row_keys_answer_the_shared_fixture(case):
    """The Web app allocates the same keys (tests/web/specific-color-actions.test.mjs)."""
    def swatch(fill, stroke="none", width=0):
        return legend_swatch({"type": "solid", "fill": fill, "stroke": stroke, "width": width})

    rules = [(rule["cap"], rule["color"]) for rule in case["rules"]]
    renderer = {row[0]: swatch(*row[1:]) for row in case["renderer_rows"]}
    rule_stroke = case.get("rule_stroke", [])

    def allocate(table):
        return rule_legend_keys(table, renderer, lambda fill: swatch(fill, *rule_stroke))

    keys = allocate(rules)
    assert [
        keys.get((caption, _legend_fill_identity(color)), "") for caption, color in rules
    ] == case["keys"]
    assert allocate(rules[::-1]) == keys


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
        instance=None,
    ) == ("#112233", "Shared")
    assert get_color_with_info(
        record.features[1],
        prepared.specific_color_rules,
        prepared.default_color_map,
        "record",
        instance=None,
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


@pytest.mark.parametrize(
    ("block_stroke_width", "expected"),
    [
        # D-13: a swatch is its fill and its line, so the track's fill with the block line is another row.
        (1.0, {"GC content [#a1a1a1]": ("GC content", "GC content (1)"), "GC content": ("GC content", "GC content (2)")}),
        # With no visible line on either swatch the rule and the track share the track's row.
        (0.0, {"GC content": ("GC content", "GC content")}),
    ],
)
def test_a_rule_named_like_a_track_row_with_its_fill_is_another_row_only_with_another_line(
    block_stroke_width, expected, tmp_path
):
    record = SeqRecord(
        Seq("ATGC" * 250), id="record", annotations={"molecule_type": "DNA", "topology": "circular"}
    )
    # #a1a1a1 is the GC content track's default fill.
    record.features = [
        SeqFeature(SimpleLocation(10, 160, strand=1), type="CDS", qualifiers={"gene": ["a"]})
    ]
    path = tmp_path / "rules.tsv"
    table([["CDS", "gene", "^a$", "#a1a1a1", "GC content"]]).to_csv(
        path, sep="\t", index=False, header=False
    )
    result = render_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                colors=ColorOptions(color_table_file=str(path)),
                config_overrides={
                    "objects.features.block_stroke_width.short": block_stroke_width,
                    "objects.features.block_stroke_color": "black",
                },
            ),
            output=RenderOutputRequest(
                formats=("interactive_svg",), output_directory=tmp_path, output_prefix="track"
            ),
        )
    )
    rows = legend_rows(result.output_paths[0])
    assert {
        key: names[:2] for key, names in rows.items() if names[0] == "GC content"
    } == expected


@pytest.mark.parametrize(
    "slots",
    [
        # No slot draws the GC content track.
        ["features:features", "gc_skew:dinucleotide_skew"],
        # The GC content slot has its own label, so no slot re-adds the default row.
        ["features:features", "gc_content:dinucleotide_content@legend_label=Composition"],
    ],
)
def test_a_rule_row_that_shares_a_default_track_row_survives_explicit_track_slots(slots, tmp_path):
    """OV-588: the slot sync removes the default track's row, not the rule row that shared it."""
    record = SeqRecord(
        Seq("ATGC" * 250), id="record", annotations={"molecule_type": "DNA", "topology": "circular"}
    )
    record.features = [
        SeqFeature(SimpleLocation(10, 160, strand=1), type="CDS", qualifiers={"gene": ["a"]})
    ]
    path = tmp_path / "rules.tsv"
    # #a1a1a1 is the GC content track's default fill; with no block line the two swatches are equal.
    table([["CDS", "gene", "^a$", "#a1a1a1", "GC content"]]).to_csv(
        path, sep="\t", index=False, header=False
    )
    result = render_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(
                colors=ColorOptions(color_table_file=str(path)),
                tracks=CircularTrackOptions(circular_track_slots=slots),
                config_overrides={"objects.features.block_stroke_width.short": 0.0},
            ),
            output=RenderOutputRequest(
                formats=("interactive_svg",), output_directory=tmp_path, output_prefix="slots"
            ),
        )
    )
    rows = legend_rows(result.output_paths[0])
    assert rows["GC content"][:2] == ("GC content", "GC content")


# The compatibility reading of PD-OI-042 scenario revision 2 (accepted residual
# risk): a caption that ends in its own rule's color ID is read without it from
# a color table file and from a saved Session. The Web app's answer to the same
# cases is services/specific-color-rules.js::savedRuleCaption.
SAVED_CAPTIONS = json.loads(
    (Path(__file__).parent / "fixtures" / "saved_rule_captions.json").read_text(encoding="utf-8")
)["cases"]
# What the Web app of first-parent `main` 9839de8e downloaded and saved for one
# caption in two colors (tests/fixtures/sessions/rule-captions-with-color-id.provenance.json).
SAVED_COLOR_TABLE = (
    Path(__file__).parent / "fixtures" / "sessions" / "rule-captions-with-color-id.specific-colors.tsv"
)
SAVED_SESSION = SAVED_COLOR_TABLE.with_name("rule-captions-with-color-id.v46.gbdraw-session.json.gz")


@pytest.mark.parametrize("case", SAVED_CAPTIONS, ids=[case["name"] for case in SAVED_CAPTIONS])
def test_a_saved_caption_answers_the_shared_fixture(case):
    assert saved_rule_caption(case["caption"], case["color"]) == case["read"]


def test_a_color_table_file_is_read_without_the_color_ids_of_its_own_rules(tmp_path):
    path = tmp_path / "saved.tsv"
    path.write_text(
        "".join(
            f"CDS\tgene\t^g{index}$\t{case['color']}\t{case['caption']}\n"
            for index, case in enumerate(SAVED_CAPTIONS)
        ),
        encoding="utf-8",
    )
    original = path.read_bytes()
    assert read_color_table(str(path)).caption.tolist() == [case["read"] for case in SAVED_CAPTIONS]
    assert path.read_bytes() == original


def _saved_caption_request(tmp_path, colors):
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
        for i, gene in enumerate(["COX1", "COX2", "ND1"])
    ]
    return CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record), record_key="record"),),
        options=CircularDiagramOptions(colors=colors),
        output=RenderOutputRequest(
            formats=("interactive_svg",), output_directory=tmp_path, output_prefix="saved", overwrite=True
        ),
    )


def test_the_downloaded_color_table_of_an_older_build_draws_numbered_names(tmp_path):
    saved = SAVED_COLOR_TABLE.read_text(encoding="utf-8")
    assert "Complex IV [#1f77b4]" in saved and "Complex IV [#ffe700]" in saved
    request = _saved_caption_request(tmp_path, ColorOptions(color_table_file=str(SAVED_COLOR_TABLE)))
    result = render_request(request)
    assert result.request.options.colors.color_table.caption.tolist() == [
        "Complex IV",
        "Complex IV",
        "Complex I",
    ]
    rows = legend_rows(result.output_paths[0])
    # The keys are the texts the older build drew, so Legend edits saved under them stay on their rows.
    assert rows["Complex IV [#ffe700]"] == ("Complex IV", "Complex IV (1)", "#ffe700")
    assert rows["Complex IV [#1f77b4]"] == ("Complex IV", "Complex IV (2)", "#1f77b4")
    assert rows["Complex I"] == ("Complex I", "Complex I", "#2ca02c")


def test_the_session_of_an_older_build_replays_with_numbered_names(tmp_path):
    with materialize_session(load_session_document(SAVED_SESSION), output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        # The Web app saves its rule list as a color table file.
        assert request.options.colors.color_table_file is not None
        output = RenderOutputRequest(
            formats=("interactive_svg",), output_directory=tmp_path, output_prefix="replay", overwrite=True
        )
        result = render_request(replace(request, output=output))
        rows = legend_rows(result.output_paths[0])
    assert result.request.options.colors.color_table.caption.tolist() == [
        "Complex IV",
        "Complex IV",
        "Complex I",
    ]
    assert rows["Complex IV [#ffe700]"] == ("Complex IV", "Complex IV (1)", "#ffe700")
    assert rows["Complex IV [#1f77b4]"] == ("Complex IV", "Complex IV (2)", "#1f77b4")
    assert not any("[#" in text or "[#" in name for name, text, _ in rows.values())


def test_the_web_request_draws_its_rule_captions_as_written(tmp_path):
    # The Web app's request holds its rule list, in the file form: a caption
    # typed with its own color ID is drawn as the app shows it.
    path = tmp_path / "specific-colors.tsv"
    path.write_text("CDS\tgene\t^COX1$\t#ffe700\tComplex IV [#ffe700]\n", encoding="utf-8")
    document = build_session_document(_saved_caption_request(tmp_path, ColorOptions()))
    payload = document.to_dict()["drawings"][0]["renderRequest"]
    payload["diagramOptions"]["colors"]["colorTableFile"] = {
        "resourceId": "colors-color-table-file",
        "representation": "file",
    }
    with materialize_session(document, output_directory=tmp_path / "web") as materialized:
        response = render_canonical_web_request(
            payload,
            resource_paths={**materialized.resource_paths, "colors-color-table-file": path},
            output_directory=tmp_path / "web-output",
        )
    svg = tmp_path / "web.svg"
    svg.write_text(response["results"][0]["content"], encoding="utf-8")
    assert legend_rows(svg)["Complex IV [#ffe700]"] == (
        "Complex IV [#ffe700]",
        "Complex IV [#ffe700]",
        "#ffe700",
    )


@pytest.mark.parametrize("form", ["file", "table"])
def test_a_saved_session_is_read_without_the_color_ids_of_its_own_rules(form, tmp_path):
    rows = [
        ["CDS", "gene", "^COX1$", "#ffe700", "Complex IV [#ffe700]"],
        ["CDS", "gene", "^COX2$", "#1f77b4", "Complex IV [#1f77b4]"],
        ["CDS", "gene", "^ND1$", "#2ca02c", "Complex I [#ffe700]"],
    ]
    path = tmp_path / "rules.tsv"
    table(rows).to_csv(path, sep="\t", index=False, header=False)
    colors = (
        ColorOptions(color_table_file=str(path)) if form == "file" else ColorOptions(color_table=table(rows))
    )
    document = build_session_document(_saved_caption_request(tmp_path, colors))
    with materialize_session(document, output_directory=tmp_path / "replay") as materialized:
        result = render_request(session_to_request(materialized))
        assert result.request.options.colors.color_table.caption.tolist() == [
            "Complex IV",
            "Complex IV",
            "Complex I [#ffe700]",
        ]
        rows_drawn = legend_rows(result.output_paths[0])
    assert rows_drawn["Complex IV [#ffe700]"][:2] == ("Complex IV", "Complex IV (1)")
    assert rows_drawn["Complex IV [#1f77b4]"][:2] == ("Complex IV", "Complex IV (2)")
    assert rows_drawn["Complex I [#ffe700]"][:2] == ("Complex I [#ffe700]", "Complex I [#ffe700]")


@pytest.mark.parametrize(
    ("caption", "pairwise"),
    [
        ("Pairwise match identity", None),
        ("Average identity", "average"),
        ("Collinear", "orientation"),
        ("Inverted", "orientation"),
    ],
)
def test_rule_captioned_like_a_comparison_legend_row_keeps_its_own_row(caption: str, pairwise: str | None) -> None:
    """OV-646 (D-03 option A): the comparison gradient replaced the rule row of its label.

    The gradient takes a source key (``<label> (2)``) when a rule already holds
    its label, as Circular's conservation rows do; both rows are drawn and
    numbered by the naming rule.
    """

    from gbdraw.legend.table import PairwiseIdentityGradient, PairwiseIdentityLegend, legend_row_names, prepare_legend_table
    from gbdraw.legend.sources import ComparisonLegendSource

    numeric = SimpleNamespace(
        show_gc=False, show_skew=False, dinucleotide="GC", stroke_color="none", stroke_width=0,
        high_fill_color="#ffffff", low_fill_color="#000000",
    )
    features = SimpleNamespace(
        color_table=pd.DataFrame(
            [("CDS", "locus_tag", "FL1", "#ff0000", caption)],
            columns=["feature_type", "qualifier", "value", "color", "caption"],
        ),
        default_colors=pd.DataFrame([("CDS", "#54bcf8"), ("default", "#d3d3d3")], columns=["feature_type", "color"]),
        block_stroke_color="none",
        block_stroke_width=0,
    )
    legend = {
        None: PairwiseIdentityLegend(),
        "average": PairwiseIdentityLegend(label="Average identity"),
        "orientation": PairwiseIdentityLegend(gradients=(
            PairwiseIdentityGradient("Collinear", "#eeeeee", "#112233"),
            PairwiseIdentityGradient("Inverted", "#ffeeee", "#445566"),
        )),
    }[pairwise]
    table = prepare_legend_table(
        numeric, numeric, features, ["CDS"],
        comparison=ComparisonLegendSource.of(SimpleNamespace(min_color="#ffffff", max_color="#000000", identity=0)),
        used_color_rules={(caption, "#ff0000")},
        show_gc=False, show_skew=False, show_depth=False,
        pairwise_identity_legend=legend,
    )

    assert table[caption]["type"] == "solid" and table[caption]["fill"] == "#ff0000"
    gradient_key = f"{caption} (2)"
    assert table[gradient_key]["type"] == "gradient"
    names = legend_row_names(table)
    assert (names[caption].name, names[gradient_key].name) == (caption, caption)
    assert (names[caption].drawn, names[gradient_key].drawn) == (f"{caption} (1)", f"{caption} (2)")
