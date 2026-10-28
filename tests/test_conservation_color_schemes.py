"""Similarity-ring colour schemes (PD-OI-110, ``tracks.circular-conservation-ring-colors`` revision 1).

One setting picks the colour of every ring without a colour of its own: ``same`` (the configured default
gradient), ``ten`` (ten colours, repeating from ring 11) or ``per_five`` (one colour per five rings, repeating from
ring 31). A ring's position is its place among the drawn rings, outermost first.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

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


def test_a_ring_coloured_like_the_default_max_draws_the_default_gradient() -> None:
    cfg = apply_config_overrides(None, None).objects.conservation
    default = conservation_track_gradient_colors(None, default_min_color=cfg.min_color, default_max_color=cfg.max_color)

    assert conservation_track_gradient_colors(
        "#8B9CC1", default_min_color=cfg.min_color, default_max_color=cfg.max_color
    ) == default == (cfg.min_color, cfg.max_color)


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
