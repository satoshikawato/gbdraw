#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.masking import ClipPath

from ...configurators import GcSkewConfigurator
from ...svg.linear_tracks import calculate_gc_skew_path_desc


class SkewDrawer:
    """
    Draws GC skew track on a linear canvas.
    """

    def __init__(self, skew_config: GcSkewConfigurator) -> None:
        self.skew_high_fill_color: str = skew_config.high_fill_color
        self.skew_low_fill_color: str = skew_config.low_fill_color
        self.skew_stroke_color: str = skew_config.stroke_color
        self.skew_stroke_width: float = skew_config.stroke_width
        self.skew_fill_opacity: float = skew_config.fill_opacity

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
        skew_desc: str = calculate_gc_skew_path_desc(
            start_x, start_y, gc_df, record_len, alignment_width, genome_size_normalization_factor, track_height
        )

        clip_id = f"clipper_line_{abs(hash(skew_desc))}"
        clip_path = ClipPath(id=clip_id)
        clip_path.add(
            Path(
                d=(
                    f"M{start_x},{start_y} "
                    f"L{alignment_width * genome_size_normalization_factor},{start_y} "
                    f"L{alignment_width * genome_size_normalization_factor},{start_y + track_height} "
                    f"L{start_x},{start_y + track_height} z"
                ),
                fill="white",
                stroke="none",
            )
        )

        skew_high = Path(
            d=skew_desc,
            fill=self.skew_high_fill_color,
            stroke=self.skew_stroke_color,
            stroke_width=self.skew_stroke_width,
            fill_opacity=self.skew_fill_opacity,
            fill_rule="evenodd",
        )

        skew_low = Path(
            d=skew_desc,
            fill=self.skew_low_fill_color,
            stroke=self.skew_stroke_color,
            stroke_width=self.skew_stroke_width,
            fill_opacity=self.skew_fill_opacity,
            clip_path=f"url(#{clip_id})",
            clip_rule="nonzero",
            fill_rule="evenodd",
        )

        group.add(clip_path)
        group.add(skew_high)
        group.add(skew_low)

        return group


__all__ = ["SkewDrawer"]


