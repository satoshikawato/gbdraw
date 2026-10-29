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
from gbdraw.legend.sources import depth_legend_sources, legend_slot_sources
from gbdraw.tracks.linear import LinearTrackSlot
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

    # A name that returns to one row is drawn as written: explicit slots take a default track's row back.
    table.add_default_track("Depth", solid("#111111"))
    table.add_source("Depth", solid("#222222"))
    assert table.names()["Depth (2)"] == ("Depth", "Depth (2)")
    table.discard_default_track("Depth")
    assert "Depth" not in table
    assert table.names()["Depth (2)"] == ("Depth", "Depth")


def test_a_default_track_row_another_source_shares_stays_for_that_source():
    """OV-588: a color rule named like a default track, with its swatch, keeps its row under explicit slots."""
    table = LegendTable()
    table.add("GC content", solid("#a1a1a1"))
    table.add_default_track("GC content", solid("#a1a1a1"))
    table.add_default_track("GC skew", solid("#333333"))
    for key in ("GC content", "GC skew", "Depth"):
        table.discard_default_track(key)
    assert list(table) == ["GC content"]
    # A row no default track drew is not a default track's row.
    table.discard_default_track("GC content")
    assert list(table) == ["GC content"]


def test_numbers_follow_the_order_the_legend_draws_solid_rows_then_gradient_rows():
    """A gradient row is drawn after every solid row (render/groups/*/legend.py), wherever the table holds it."""
    table = LegendTable()
    table.add("CDS", solid("#47b8f8", "gray", 1))
    table.add("Pairwise match identity", {"type": "gradient", "min_color": "#ffe7e7", "max_color": "#ff7272"})
    table.add_source("Pairwise match identity", solid("#88aaff", "#404040", 1.5))

    # The table keeps its order and its keys; the solid row, drawn first, takes the first number.
    assert list(table) == ["CDS", "Pairwise match identity", "Pairwise match identity (2)"]
    assert table.names() == {
        "CDS": ("CDS", "CDS"),
        "Pairwise match identity": ("Pairwise match identity", "Pairwise match identity (2)"),
        "Pairwise match identity (2)": ("Pairwise match identity", "Pairwise match identity (1)"),
    }


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
        LegendTable(), depth_legend_sources(_depth_track(index, fill) for index, fill in enumerate(fills))
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
    slot = LinearTrackSlot(id="islands", renderer="annotations", params={"set_id": "islands"})
    table = LegendTable()
    table.add("Island", solid("#ff0000", "gray", 1))
    table = sync_annotation_legend_entries(table, bundle, legend_slot_sources([slot]))

    # The legend swatch draws no opacity, so the second source shares the first's row.
    assert list(table) == ["Island", "Island (2)", "Island (4)"]
    assert table.shared_keys == ("Island (3)",)
    assert [names.drawn for names in table.names().values()] == [
        "Island (1)",
        "Island (2)",
        "Island (3)",
    ]
    assert table["Island (4)"]["annotation_hatch"] == hatch


def test_only_a_legend_label_draws_an_annotation_color_in_the_legend():
    """O-8: the Web applies an annotation item's Fill and Line color live only
    while neither the item nor its set has a Legend label (the resolved
    ``legend_label`` is the item's or else the set's). Without one, no Legend
    row reads the colors; with one, each distinct color is a row of its own."""
    slot = LinearTrackSlot(id="islands", renderer="annotations", params={"set_id": "islands"})

    def rows(*annotations: tuple[str, str | None]) -> list[str]:
        bundle = SimpleNamespace(annotations=tuple(
            SimpleNamespace(  # the fields the legend reads
                set_id="islands",
                mark="highlight",
                legend_label=legend_label,
                style=RegionAnnotationStyle(fill=fill, stroke="#334155", stroke_width=1),
                style_is_annotation_override=True,
            )
            for fill, legend_label in annotations
        ))
        return list(sync_annotation_legend_entries(LegendTable(), bundle, legend_slot_sources([slot])))

    assert rows(("#ff0000", None)) == rows(("#00ff00", None)) == []
    assert len(rows(("#ff0000", "Island"), ("#00ff00", "Island"))) == 2
    assert len(rows(("#ff0000", "Island"), ("#ff0000", "Island"))) == 1


def test_the_layout_measures_the_drawn_name():
    def width(table: LegendTable) -> float:
        return build_circular_legend_layout(
            table,
            legend_position="right",
            canvas_width=1000,
            font_family="Liberation Sans",
            font_size=16,
            color_rect_size=16,
        ).width

    one = LegendTable()
    one.add("Packaging", solid("#ffe700", "gray", 1))
    two = LegendTable()
    two.add("Packaging", solid("#ffe700", "gray", 1), key="Packaging [#ffe700]")
    two.add("Packaging", solid("#f000ff", "gray", 1), key="Packaging [#f000ff]")
    named_by_key = LegendTable({"Packaging [#ffe700]": solid("#ffe700", "gray", 1)})

    assert width(one) < width(two) < width(named_by_key)
