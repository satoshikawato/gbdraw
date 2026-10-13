"""One naming rule for every legend row (PD-OI-042 scenario revision 2, PD-OI-095).

A row's key, its name, and the name the legend draws are three facts
(gbdraw/legend/table.py): the key is an identity and keeps its strings, the
name is never rewritten, and the drawn name is derived from the rows drawn.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from gbdraw.analysis.depth_tracks import sync_depth_track_legend_entries
from gbdraw.annotations.legend import sync_annotation_legend_entries
from gbdraw.annotations.models import HatchStyle, RegionAnnotationStyle
from gbdraw.legend.circular_layout import build_circular_legend_layout
from gbdraw.legend.table import (
    LegendTable,
    drawn_legend_names,
    legend_row_names,
    legend_swatch,
)

LEGEND_NAMES = json.loads(
    (Path(__file__).parent / "fixtures" / "legend_names.json").read_text(encoding="utf-8")
)["cases"]


def solid(fill: str, stroke: str = "none", width: float = 0) -> dict[str, object]:
    return {"type": "solid", "fill": fill, "stroke": stroke, "width": width}


@pytest.mark.parametrize("case", LEGEND_NAMES, ids=[case["name"] for case in LEGEND_NAMES])
def test_drawn_names_answer_the_shared_fixture(case):
    """The Web app numbers what it displays with the same rule (tests/web/legend-row-names.test.mjs)."""
    rows = [(row["name"], legend_swatch(row["swatch"])) for row in case["rows"]]
    assert drawn_legend_names(rows) == case["drawn"]


def test_sources_of_one_label_keep_their_keys_and_share_or_number_their_rows():
    table = LegendTable()
    table.add("CDS", solid("#47b8f8", "gray", 1))
    for fill in ["#808080", "gray", "#336699"]:
        table.add_source("GC content", solid(fill))

    # The second source shares the first row; its key stays taken and unused.
    assert list(table) == ["CDS", "GC content", "GC content (3)"]
    assert table.shared_keys == ("GC content (2)",)
    assert table.names() == {
        "CDS": ("CDS", "CDS"),
        "GC content": ("GC content", "GC content (1)"),
        "GC content (3)": ("GC content", "GC content (2)"),
    }

    table.discard("GC content (3)")
    assert table.names()["GC content"] == ("GC content", "GC content")


def test_a_plain_mapping_is_read_as_rows_named_by_their_keys():
    names = legend_row_names({"CDS": solid("#47b8f8"), "GC content (2)": solid("#808080")})
    assert names == {"CDS": ("CDS", "CDS"), "GC content (2)": ("GC content (2)", "GC content (2)")}


def _depth_track(index: int, fill: str, label: str = "Coverage") -> SimpleNamespace:
    return SimpleNamespace(
        track_index=index,
        label=label,
        config=SimpleNamespace(fill_color=fill, stroke_color="none", stroke_width=0),
    )


@pytest.mark.parametrize(
    "fills,keys,drawn",
    [
        (["#111111", "#111111"], ["Coverage"], ["Coverage"]),
        (["#111111", "#222222"], ["Coverage", "Coverage (2)"], ["Coverage (1)", "Coverage (2)"]),
        (
            ["#111111", "#111111", "#222222"],
            ["Coverage", "Coverage (3)"],
            ["Coverage (1)", "Coverage (2)"],
        ),
    ],
)
def test_depth_tracks_with_one_label(fills, keys, drawn):
    table = sync_depth_track_legend_entries(
        LegendTable(), [_depth_track(index, fill) for index, fill in enumerate(fills)]
    )
    assert list(table) == keys
    assert [names.drawn for names in table.names().values()] == drawn
    assert {names.name for names in table.names().values()} == {"Coverage"}


def test_annotation_rows_follow_the_rule_and_keep_their_keys():
    hatch = HatchStyle(angle=45, spacing=6, color="#334155", width=1, cross=False)

    def annotation(style: RegionAnnotationStyle) -> SimpleNamespace:
        return SimpleNamespace(  # the fields the legend reads
            set_id="islands",
            mark="highlight",
            legend_label="Island",
            style=style,
            style_is_annotation_override=True,
        )

    styles = [
        RegionAnnotationStyle(fill="#94a3b8", stroke="#334155", stroke_width=1),
        RegionAnnotationStyle(fill="#94a3b8", stroke="#334155", stroke_width=1, fill_opacity=0.4),
        RegionAnnotationStyle(fill="#94a3b8", stroke="#334155", stroke_width=1, hatch=hatch),
    ]
    bundle = SimpleNamespace(annotations=tuple(annotation(style) for style in styles))
    slot = SimpleNamespace(enabled=True, renderer="annotations", params={"set_id": "islands"})
    table = LegendTable()
    table.add("Island", solid("#ff0000", "gray", 1))
    table = sync_annotation_legend_entries(table, bundle, [slot])

    # The legend swatch draws no opacity, so the second source shares the first's row.
    assert list(table) == ["Island", "Island (2)", "Island (4)"]
    assert table.shared_keys == ("Island (3)",)
    assert [names.drawn for names in table.names().values()] == [
        "Island (1)",
        "Island (2)",
        "Island (3)",
    ]
    assert table["Island (4)"]["annotation_hatch"] == hatch


def test_the_layout_measures_the_drawn_name():
    def width(table: LegendTable) -> float:
        return build_circular_legend_layout(
            table,
            legend_position="right",
            canvas_width=1000,
            font_family="Liberation Sans",
            font_size=16,
            dpi=96,
            color_rect_size=16,
        ).width

    one = LegendTable()
    one.add("Packaging", solid("#ffe700", "gray", 1))
    two = LegendTable()
    two.add("Packaging", solid("#ffe700", "gray", 1), key="Packaging [#ffe700]")
    two.add("Packaging", solid("#f000ff", "gray", 1), key="Packaging [#f000ff]")
    named_by_key = LegendTable({"Packaging [#ffe700]": solid("#ffe700", "gray", 1)})

    assert width(one) < width(two) < width(named_by_key)
