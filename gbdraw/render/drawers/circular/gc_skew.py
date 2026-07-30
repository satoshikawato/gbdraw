#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.masking import ClipPath

from ....configurators import GcSkewConfigurator
from ....svg.circular_tracks import generate_circle_path_desc, generate_circular_gc_skew_path_desc
from ....svg.ids import stable_svg_id


class SkewDrawer:
    """
    Draws GC skew track on a circular canvas.
    """

    def __init__(self, skew_config: GcSkewConfigurator) -> None:
        self.skew_high_fill_color: str = skew_config.high_fill_color
        self.skew_low_fill_color: str = skew_config.low_fill_color
        self.skew_stroke_color: str = skew_config.stroke_color
        self.skew_stroke_width: float = skew_config.stroke_width
        self.skew_fill_opacity: float = skew_config.fill_opacity

    def draw(
        self,
        radius: float,
        group: Group,
        gc_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        dinucleotide: str,
        record_identifier: str | None = None,
        group_identifier: str | None = None,
    ) -> Group:
        skew_desc: str = generate_circular_gc_skew_path_desc(
            radius, gc_df, record_len, track_width, norm_factor, dinucleotide
        )
        circle_desc: str = generate_circle_path_desc(radius, norm_factor)
        clip_id = stable_svg_id(
            "clipper_circle",
            "circular-skew-clip",
            group_identifier,
            record_identifier,
            record_len,
            track_width,
            norm_factor,
            dinucleotide,
            skew_desc,
        )
        circle_path: ClipPath = ClipPath(id=clip_id)
        circle_path.add(Path(d=circle_desc, fill="white", stroke="none"))
        skew_high: Path = Path(
            d=skew_desc,
            fill=self.skew_high_fill_color,
            stroke=self.skew_stroke_color,
            fill_opacity=self.skew_fill_opacity,
            fill_rule="evenodd",
        )
        skew_low: Path = Path(
            d=skew_desc,
            fill=self.skew_low_fill_color,
            stroke=self.skew_stroke_color,
            fill_opacity=self.skew_fill_opacity,
            clip_path=f"url(#{clip_id})",
            clip_rule="nonzero",
            fill_rule="evenodd",
        )
        if (
            self.skew_stroke_width > 0
            or self.skew_stroke_color.strip().lower() != "none"
        ):
            for path in (skew_high, skew_low):
                path["stroke-width"] = self.skew_stroke_width
        group.add(circle_path)
        group.add(skew_high)
        group.add(skew_low)
        return group


__all__ = ["SkewDrawer"]


