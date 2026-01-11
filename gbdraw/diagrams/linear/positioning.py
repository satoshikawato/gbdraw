#!/usr/bin/env python
# coding: utf-8

"""Positioning helpers for linear diagram assembly."""

from __future__ import annotations

from typing import Tuple

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]


def calculate_record_offsets(count: int, record: SeqRecord, canvas_config: LinearCanvasConfigurator) -> Tuple[float, float]:
    """
    Calculates the vertical and horizontal offsets for placing a genomic record in the canvas.
    """
    offset: float = canvas_config.vertical_offset + (
        canvas_config.cds_height
        + canvas_config.vertical_padding
        + canvas_config.gc_padding
        + canvas_config.skew_padding
        + canvas_config.comparison_height
        + canvas_config.vertical_padding
    ) * (count - 1)
    offset_x: float = (
        (canvas_config.alignment_width * ((canvas_config.longest_genome - len(record.seq)) / canvas_config.longest_genome) / 2)
        if canvas_config.align_center
        else 0
    )
    return offset, offset_x


def position_record_group(record_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """Translates the record group to its designated position on the canvas."""
    record_group.translate(offset_x + canvas_config.horizontal_offset, offset)
    return record_group


def position_gc_content_group(
    gc_content_group: Group, offset_y: float, offset_x: float, canvas_config: LinearCanvasConfigurator
) -> Group:
    """Positions the GC content group on the canvas."""
    gc_content_group.translate(
        offset_x + canvas_config.horizontal_offset,
        offset_y + canvas_config.cds_padding + canvas_config.vertical_padding,
    )
    return gc_content_group


def position_gc_skew_group(gc_skew_group: Group, offset_y: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """Positions the GC skew group on the canvas."""
    y_offset = offset_y + canvas_config.cds_padding + canvas_config.vertical_padding
    if canvas_config.show_gc:
        y_offset += canvas_config.gc_height

    gc_skew_group.translate(offset_x + canvas_config.horizontal_offset, y_offset)
    return gc_skew_group


def position_record_definition_group(
    record_definition_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator
) -> Group:
    """Places the record definition group in the correct position on the canvas."""
    record_definition_group.translate(canvas_config.horizontal_offset - offset_x, offset)
    return record_definition_group


def position_length_bar_group(total_height: float, vertical_offset: float, vertical_padding: float) -> float:
    """Calculates the vertical position for the length bar on the canvas."""
    return total_height - (vertical_offset - vertical_padding)


def position_comparison_group(comparison_count: int, canvas_config: LinearCanvasConfigurator) -> float:
    """Calculates the vertical position for placing a comparison group on the canvas."""
    return (
        (
            canvas_config.vertical_offset
            + 0.5 * canvas_config.cds_height
            + canvas_config.vertical_padding
            + 0.9 * canvas_config.gc_padding
            + 0.9 * canvas_config.skew_padding
        )
        + (
            (
                canvas_config.comparison_height
                + canvas_config.vertical_padding
                + canvas_config.cds_height
                + canvas_config.vertical_padding
                + canvas_config.gc_padding
                + canvas_config.skew_padding
            )
            * (comparison_count - 1)
        )
    )


__all__ = [
    "calculate_record_offsets",
    "position_record_group",
    "position_gc_content_group",
    "position_gc_skew_group",
    "position_record_definition_group",
    "position_length_bar_group",
    "position_comparison_group",
]


