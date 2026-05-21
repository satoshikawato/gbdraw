from __future__ import annotations

from gbdraw.svg.circular_ticks import (
    generate_circular_tick_labels,
    resolve_circular_tick_label_geometry,
    set_tick_label_anchor_value,
)


def _tick_label_kwargs() -> dict:
    return {
        "center_radius_px": 100.0,
        "total_len": 4_641_652,
        "size": "large",
        "font_size": 14.0,
        "font_family": "Arial",
        "track_type": "tuckin",
        "strandedness": True,
        "dpi": 96,
        "label_side": "outside",
        "tick_side": "inside",
        "tick_length_px": 10.0,
        "tick_width": 2.0,
    }


def test_circular_tick_label_anchor_stays_middle_and_baseline_uses_edge_values_by_angle() -> None:
    assert set_tick_label_anchor_value(4_641_652, 1_000_000) == ("middle", "text-after-edge")
    assert set_tick_label_anchor_value(4_641_652, 2_500_000) == ("middle", "text-before-edge")
    assert set_tick_label_anchor_value(4_641_652, 3_000_000) == ("middle", "text-before-edge")
    assert set_tick_label_anchor_value(4_641_652, 4_500_000) == ("middle", "text-after-edge")


def test_circular_tick_label_geometry_uses_middle_anchor() -> None:
    base_kwargs = _tick_label_kwargs()

    right_side = resolve_circular_tick_label_geometry(
        **base_kwargs,
        tick=1_000_000,
        label_text="1.0 Mbp",
    )
    left_side = resolve_circular_tick_label_geometry(
        **base_kwargs,
        tick=3_000_000,
        label_text="3.0 Mbp",
    )
    top_side = resolve_circular_tick_label_geometry(
        **base_kwargs,
        tick=4_500_000,
        label_text="4.5 Mbp",
    )

    assert right_side.text_anchor == "middle"
    assert left_side.text_anchor == "middle"
    assert top_side.text_anchor == "middle"


def test_circular_tick_label_textpath_uses_middle_anchor() -> None:
    elements = generate_circular_tick_labels(
        100.0,
        4_641_652,
        "large",
        [1_000_000, 3_000_000, 4_500_000],
        "none",
        "black",
        14.0,
        "normal",
        "Arial",
        "tuckin",
        True,
        96,
        label_side="outside",
        tick_side="inside",
        tick_length_px=10.0,
        tick_width=2.0,
    )
    text_elements = [element.tostring() for element in elements if element.elementname == "text"]

    assert 'text-anchor="middle"' in text_elements[0]
    assert 'startOffset="50%"' in text_elements[0]
    assert 'dominant-baseline="text-after-edge"' in text_elements[0]
    assert 'text-anchor="middle"' in text_elements[1]
    assert 'startOffset="50%"' in text_elements[1]
    assert 'dominant-baseline="text-before-edge"' in text_elements[1]
    assert 'text-anchor="middle"' in text_elements[2]
    assert 'startOffset="50%"' in text_elements[2]
    assert 'dominant-baseline="text-after-edge"' in text_elements[2]
