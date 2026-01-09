#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path

from ....svg.linear_tracks import calculate_gc_content_path_desc


class GcContentDrawer:
    """
    Draws GC content track on a linear canvas.
    """

    def __init__(self, gc_config):
        self.gc_path_high_fill_color: str = gc_config.high_fill_color
        self.gc_path_low_fill_color: str = gc_config.low_fill_color
        self.gc_path_fill_color: str = gc_config.fill_color
        self.gc_path_stroke_color: str = gc_config.stroke_color
        self.gc_path_stroke_width: float = gc_config.stroke_width
        self.gc_path_fill_opacity: float = gc_config.fill_opacity

    def draw(
        self,
        group: Group,
        gc_df: DataFrame,
        record_len: int,
        alignment_width: float,
        genome_size_normalization_factor: float,
        track_height: float,
        start_x: float,
        start_y: float,
        dinucleotide: str,
    ) -> Group:
        gc_path_desc: str = calculate_gc_content_path_desc(
            start_x,
            start_y,
            gc_df,
            record_len,
            alignment_width,
            genome_size_normalization_factor,
            track_height,
            dinucleotide,
        )
        gc_path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd",
        )
        group.add(gc_path)
        return group


__all__ = ["GcContentDrawer"]


