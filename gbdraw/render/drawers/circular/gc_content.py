#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path

from ....svg.circular_tracks import (  # type: ignore[reportMissingImports]
    generate_circular_gc_content_path_desc,
)


class GcContentDrawer:
    """
    Draws GC content track on a circular canvas.
    """

    def __init__(self, gc_config) -> None:
        self.gc_path_high_fill_color: str = gc_config.high_fill_color
        self.gc_path_low_fill_color: str = gc_config.low_fill_color
        self.gc_path_fill_color: str = gc_config.fill_color
        self.gc_path_stroke_color: str = gc_config.stroke_color
        self.gc_path_stroke_width: float = gc_config.stroke_width
        self.gc_path_fill_opacity: float = gc_config.fill_opacity

    def draw(
        self,
        radius: float,
        group: Group,
        gc_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        dinucleotide: str,
    ) -> Group:
        gc_path_desc: str = generate_circular_gc_content_path_desc(
            radius, record_len, gc_df, track_width, norm_factor, dinucleotide
        )
        gc_path: Path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd",
        )
        group.add(gc_path)
        return group


__all__ = ["GcContentDrawer"]


