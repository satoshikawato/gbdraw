"""Similarity-ring colour schemes (PD-OI-110, ``tracks.circular-conservation-ring-colors`` revision 1).

One setting picks the colour of every ring without a colour of its own: ``same`` (the configured default
gradient), ``ten`` (ten colours, repeating from ring 11) or ``per_five`` (one colour per five rings, repeating from
ring 31). A ring's position is its place among the drawn rings, outermost first.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import get_args

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.conservation import (
    CONSERVATION_COLOR_SCHEMES,
    ConservationTrack,
    conservation_track_gradient_colors,
    resolve_conservation_track_colors,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_circular_diagram_from_record, assemble_circular_diagram_from_records
from gbdraw.config.models.objects import ConservationColorScheme
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.tracks import CircularTrackSlot

TEN = ("#6e91b7", "#f4a251", "#8fc4c0", "#e67577", "#c6bebb", "#77b26f", "#be92b2", "#f0d369", "#ae8e7c", "#ffafb7")
HUES = ("#8b9cc1", "#c88890", "#62ad9a", "#bf89ab", "#7caa7e", "#bc9464")
SCHEME_PATH = "objects.conservation.color_scheme"


def _track(index: int, color: str | None = None) -> ConservationTrack:
    return ConservationTrack(
        source_index=index - 1,
        track_index=index,
        track_label=f"ring {index}",
        track_color=color,
        reference_side="subject",
        hits=pd.DataFrame({"track_color": [color or ""]}),
    )


def _slot(track_index: int, *, enabled: bool = True) -> CircularTrackSlot:
    return CircularTrackSlot(
        id=f"conservation_{track_index}",
        renderer="sequence_conservation",
        enabled=enabled,
        params={"track_index": track_index},
    )


def test_scheme_colours_are_the_receipt_cycles() -> None:
    assert CONSERVATION_COLOR_SCHEMES["same"] == ()
    assert CONSERVATION_COLOR_SCHEMES["ten"] == TEN
    assert CONSERVATION_COLOR_SCHEMES["per_five"] == tuple(hue for hue in HUES for _ in range(5))


@pytest.mark.parametrize(
    ("scheme", "ring", "expected"),
    [
        ("ten", 1, TEN[0]),
        ("ten", 10, TEN[9]),
        ("ten", 11, TEN[0]),
        ("per_five", 5, HUES[0]),
        ("per_five", 6, HUES[1]),
        ("per_five", 30, HUES[5]),
        ("per_five", 31, HUES[0]),
        ("same", 3, None),
    ],
)
def test_ring_without_a_colour_takes_the_scheme_colour_of_its_drawn_position(
    scheme: str, ring: int, expected: str | None
) -> None:
    tracks = tuple(_track(index) for index in range(1, 32))
    slots = [_slot(index) for index in range(1, 32)]

    resolved = resolve_conservation_track_colors(scheme, slots, tracks)

    assert resolved[ring - 1].track_color == expected


def test_position_is_the_drawn_order_and_an_own_colour_wins() -> None:
    tracks = (_track(1), _track(2, "#123456"), _track(3), _track(4))
    # Drawn outermost first: ring 3, ring 2, ring 1; ring 4 has no slot and is not drawn.
    slots = [CircularTrackSlot(id="features", renderer="features"), _slot(3), _slot(2), _slot(1)]

    resolved = resolve_conservation_track_colors("ten", slots, tracks)

    assert [track.track_index for track in resolved] == [3, 2, 1]
    assert [track.track_color for track in resolved] == [TEN[0], "#123456", TEN[2]]
    # The per-hit colour attribute reads the same colour.
    assert [track.hits["track_color"].tolist() for track in resolved] == [[TEN[0]], ["#123456"], [TEN[2]]]


def test_a_disabled_slot_draws_no_ring() -> None:
    resolved = resolve_conservation_track_colors("ten", [_slot(1, enabled=False), _slot(2)], (_track(1), _track(2)))

    assert [(track.track_index, track.track_color) for track in resolved] == [(2, TEN[0])]


def test_a_ring_coloured_like_the_default_max_keeps_its_own_gradient() -> None:
    """A ring's own colour draws white-tinted-to-colour, even when it equals the configured max (dev behaviour;
    saved work keeps its colours)."""
    from gbdraw.core.color import tint_color

    cfg = apply_config_overrides(None, None).objects.conservation
    default = conservation_track_gradient_colors(None, default_min_color=cfg.min_color, default_max_color=cfg.max_color)

    assert default == (cfg.min_color, cfg.max_color)
    assert conservation_track_gradient_colors(
        "#8B9CC1", default_min_color=cfg.min_color, default_max_color=cfg.max_color
    ) == (tint_color("#8b9cc1"), "#8b9cc1")


def test_config_scheme_defaults_to_same_and_rejects_other_values() -> None:
    assert apply_config_overrides(None, None).objects.conservation.color_scheme == "same"
    assert apply_config_overrides(None, {SCHEME_PATH: "per_five"}).objects.conservation.color_scheme == "per_five"
    with pytest.raises(ValidationError) as excinfo:
        apply_config_overrides(None, {SCHEME_PATH: "rainbow"})
    assert excinfo.value.diagnostic == {"code": "INPUT_INVALID", "configPath": SCHEME_PATH}


def _record(record_id: str = "rec1", length: int = 120) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, name=record_id)
    record.annotations["molecule_type"] = "DNA"
    return record


def _frame(subject: str = "rec1") -> pd.DataFrame:
    row = ("query1", subject, 95.0, 20, 0, 0, 1, 20, 5, 24, 1e-20, 100.0)
    return pd.DataFrame([row], columns=COMPARISON_COLUMNS)


def _ring_colours(svg: str) -> list[str | None]:
    return [
        element.get("data-track-color")
        for element in ET.fromstring(svg).iter()
        if element.tag.rsplit("}", 1)[-1] == "g" and element.get("data-gbdraw-slot-renderer") == "sequence_conservation"
    ]


def test_assembly_colours_rings_by_the_configured_scheme() -> None:
    canvas = assemble_circular_diagram_from_record(
        _record(),
        cfg=apply_config_overrides(None, {SCHEME_PATH: "ten"}),
        conservation_dataframes=[_frame(), _frame(), _frame()],
        conservation_reference="subject",
        conservation_colors=["", "#123456", ""],
        legend="none",
    )

    assert _ring_colours(canvas.tostring()) == [TEN[0], "#123456", TEN[2]]


def test_multi_record_canvas_colours_every_record_alike() -> None:
    canvas = assemble_circular_diagram_from_records(
        [_record("rec1"), _record("rec2")],
        cfg=apply_config_overrides(None, {SCHEME_PATH: "ten"}),
        conservation_dataframes=[_frame("rec1"), _frame("rec2")],
        conservation_reference="subject",
        legend="none",
    )

    assert _ring_colours(canvas.tostring()) == [TEN[0], TEN[1], TEN[0], TEN[1]]


def _write_inputs(tmp_path, count: int = 3):
    from Bio import SeqIO

    gbk = tmp_path / "rec1.gb"
    SeqIO.write([_record()], str(gbk), "genbank")
    blasts = []
    for index in range(count):
        path = tmp_path / f"ring{index + 1}.tsv"
        _frame().to_csv(path, sep="\t", header=False, index=False)
        blasts.append(str(path))
    return str(gbk), blasts


def _run_cli(tmp_path, monkeypatch, *args: str) -> str:
    import sys

    from gbdraw import cli

    output = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["gbdraw", "circular", *args, "-o", str(output), "-f", "svg"])
    cli.main()
    return output.with_suffix(".svg").read_text()


def test_cli_scheme_and_auto_ring_colours(tmp_path, monkeypatch) -> None:
    gbk, blasts = _write_inputs(tmp_path)

    svg = _run_cli(
        tmp_path, monkeypatch,
        "--gbk", gbk, "--conservation_blast", *blasts, "--conservation_reference", "subject",
        "--conservation_colors", "auto", "#123456", "AUTO", "--conservation_color_scheme", "ten",
    )

    assert _ring_colours(svg) == [TEN[0], "#123456", TEN[2]]


def test_scheme_names_have_one_owner() -> None:
    from gbdraw.circular import _get_args

    assert tuple(CONSERVATION_COLOR_SCHEMES) == get_args(ConservationColorScheme)
    for scheme in get_args(ConservationColorScheme):
        assert _get_args(["--gbk", "dummy.gb", "--conservation_color_scheme", scheme]).conservation_color_scheme == scheme


def test_cli_parses_the_scheme_choices() -> None:
    from gbdraw.circular import _get_args

    assert _get_args(["--gbk", "dummy.gb"]).conservation_color_scheme is None
    assert _get_args(["--gbk", "dummy.gb", "--conservation_color_scheme", "per_five"]).conservation_color_scheme == "per_five"
    with pytest.raises(SystemExit):
        _get_args(["--gbk", "dummy.gb", "--conservation_color_scheme", "rainbow"])


def test_cli_table_blank_colour_cell_takes_the_scheme_colour(tmp_path, monkeypatch) -> None:
    gbk, blasts = _write_inputs(tmp_path, 2)
    table = tmp_path / "rings.tsv"
    table.write_text(f"blast\tlabel\tcolor\n{blasts[0]}\tA\t\n{blasts[1]}\tB\tred\n")

    svg = _run_cli(
        tmp_path, monkeypatch,
        "--gbk", gbk, "--conservation_table", str(table), "--conservation_reference", "subject",
        "--conservation_color_scheme", "per_five",
    )

    assert _ring_colours(svg) == [HUES[0], "#ff0000"]


def test_python_api_scheme_and_rings_without_a_colour() -> None:
    import gbdraw

    rings = gbdraw.ComparisonRingOptions(
        tracks=[
            gbdraw.ComparisonRingTrackOptions(source=_frame()),
            gbdraw.ComparisonRingTrackOptions(source=_frame(), color="#123456"),
            gbdraw.ComparisonRingTrackOptions(source=_frame()),
        ],
        reference="subject",
        color_scheme="ten",
    )
    diagram = gbdraw.draw_circular(_record(), options=gbdraw.CircularOptions(comparison_rings=rings, legend="none"))

    assert _ring_colours(diagram.to_svg()) == [TEN[0], "#123456", TEN[2]]


@pytest.mark.parametrize(
    ("color_scheme", "expected"), [(None, [TEN[0], TEN[1]]), ("same", [None, None]), ("per_five", [HUES[0], HUES[0]])]
)
def test_python_api_scheme_none_keeps_the_config_scheme_and_a_value_overrides_it(
    color_scheme: str | None, expected: list[str | None]
) -> None:
    import gbdraw

    rings = gbdraw.ComparisonRingOptions(
        tracks=[gbdraw.ComparisonRingTrackOptions(source=_frame()), gbdraw.ComparisonRingTrackOptions(source=_frame())],
        reference="subject",
        color_scheme=color_scheme,
    )
    options = gbdraw.CircularOptions(comparison_rings=rings, legend="none", config_overrides={SCHEME_PATH: "ten"})
    diagram = gbdraw.draw_circular(_record(), options=options)

    assert [color or None for color in _ring_colours(diagram.to_svg())] == expected


def test_python_api_rejects_an_unknown_scheme() -> None:
    import gbdraw

    rings = gbdraw.ComparisonRingOptions(tracks=[gbdraw.ComparisonRingTrackOptions(source=_frame())], color_scheme="x")
    with pytest.raises(ValidationError):
        gbdraw.draw_circular(_record(), options=gbdraw.CircularOptions(comparison_rings=rings))


_SVG = {"svg": "http://www.w3.org/2000/svg"}


def _ring_block(svg: str) -> ET.Element:
    block = ET.fromstring(svg).find(".//svg:g[@id='conservation_identity_legend']", _SVG)
    assert block is not None
    return block


def _block_rows(block: ET.Element) -> list[tuple[str | None, str | None, str | None, str | None]]:
    """(name, source index, track index, swatch fill) of each ring row, top to bottom."""

    rows = []
    for row in block.findall("svg:g[@data-legend-key]", _SVG):
        swatch = row.find("svg:path", _SVG)
        rows.append((
            row.findtext("svg:text", namespaces=_SVG),
            row.get("data-source-index"),
            row.get("data-track-index"),
            None if swatch is None else swatch.get("fill"),
        ))
    return rows


def _block_texts(block: ET.Element) -> list[str | None]:
    return [text.text for text in block.iter("{http://www.w3.org/2000/svg}text")]


def _band_stops(block: ET.Element) -> list[str | None]:
    return [stop.get("stop-color") for stop in block.iter("{http://www.w3.org/2000/svg}stop")]


def _draw_legend(scheme: str, colors: list[str] | None = None, *, slots: list[CircularTrackSlot] | None = None) -> str:
    return assemble_circular_diagram_from_record(
        _record(),
        cfg=apply_config_overrides(None, {SCHEME_PATH: scheme}),
        conservation_dataframes=[_frame(), _frame(), _frame()],
        conservation_reference="subject",
        conservation_labels=["A", "B", "C"],
        conservation_colors=colors,
        circular_track_slots=slots,
        identity=90,
        legend="right",
    ).tostring()


def test_one_colour_legend_lists_the_drawn_rings_without_swatches() -> None:
    block = _ring_block(_draw_legend("same"))

    assert _block_rows(block) == [("A", "0", "1", None), ("B", "1", "2", None), ("C", "2", "3", None)]
    assert _block_texts(block) == ["outer", "A", "B", "C", "inner", "identity", "90%", "100%"]
    assert _band_stops(block) == ["#f0f1f5", "#8b9cc1"]


def test_captions_and_scale_labels_are_muted_and_ring_names_keep_the_legend_text_colour() -> None:
    block = _ring_block(_draw_legend("ten", ["", "#123456", ""]))
    font_sizes = {text.get("font-size") for text in block.iter("{http://www.w3.org/2000/svg}text")}
    names = [row.find("svg:text", _SVG) for row in block.findall("svg:g[@data-legend-key]", _SVG)]
    muted = [text for text in block.iter("{http://www.w3.org/2000/svg}text") if text not in names]

    assert [text.text for text in muted] == ["outer", "inner", "identity", "90%", "100%"]
    assert {text.get("fill") for text in muted} == {"#5a6273"}
    assert [text.get("fill") for text in names] == [None, None, None]
    assert len(font_sizes) == 1

@pytest.mark.parametrize(
    ("colors", "band_min"), [(["#8B9CC1", "", ""], "tint"), (["", "#8B9CC1", ""], "#f0f1f5")]
)
def test_a_ring_coloured_like_the_default_keeps_the_one_colour_legend(colors: list[str], band_min: str) -> None:
    """A ring without a colour counts as the configured max; the band is the first drawn ring's gradient."""
    from gbdraw.core.color import tint_color

    block = _ring_block(_draw_legend("same", colors))

    assert [row[3] for row in _block_rows(block)] == [None, None, None]
    assert _band_stops(block) == [tint_color("#8b9cc1") if band_min == "tint" else band_min, "#8b9cc1"]


def test_mixed_colour_legend_shows_swatches_and_the_grey_band() -> None:
    from gbdraw.core.color import tint_color

    block = _ring_block(_draw_legend("ten", ["", "#123456", ""]))

    assert [row[3] for row in _block_rows(block)] == [TEN[0], "#123456", TEN[2]]
    assert _band_stops(block) == [tint_color("#9c9c9c"), "#9c9c9c"]


def test_legend_lists_rings_in_drawn_order() -> None:
    slots = [
        CircularTrackSlot(id="features", renderer="features"),
        *[
            CircularTrackSlot(id=f"c{index}", renderer="sequence_conservation", params={"track_index": index})
            for index in (3, 1, 2)
        ],
    ]
    block = _ring_block(_draw_legend("ten", slots=slots))

    assert _block_rows(block) == [("C", "2", "3", TEN[0]), ("A", "0", "1", TEN[1]), ("B", "1", "2", TEN[2])]


def test_multi_record_canvas_draws_the_same_ring_block() -> None:
    def block_of(svg: str) -> list[object]:
        block = _ring_block(svg)
        return [_block_rows(block), _block_texts(block), _band_stops(block)]

    single = _draw_legend("ten", ["", "#123456", ""])
    multi = assemble_circular_diagram_from_records(
        [_record("rec1"), _record("rec1b")],
        cfg=apply_config_overrides(None, {SCHEME_PATH: "ten"}),
        conservation_dataframes=[_frame(), _frame(), _frame()],
        conservation_reference="subject",
        conservation_labels=["A", "B", "C"],
        conservation_colors=["", "#123456", ""],
        identity=90,
        legend="right",
    ).tostring()

    assert block_of(multi) == block_of(single)
