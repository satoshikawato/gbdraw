#!/usr/bin/env python
# coding: utf-8

import hashlib

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.masking import ClipPath

from ....configurators import GcSkewConfigurator
from ....svg.circular_tracks import generate_circle_path_desc, generate_circular_gc_skew_path_desc


class SkewDrawer:
    """
    Draws GC skew track on a circular canvas.
    """

    def __init__(self, skew_config: GcSkewConfigurator) -> None:
        self.skew_high_fill_color: str = skew_config.high_fill_color
        self.skew_low_fill_color: str = skew_config.low_fill_color
        self.skew_stroke_color: str = skew_config.stroke_color
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
    ) -> Group:
        skew_desc: str = generate_circular_gc_skew_path_desc(
            radius, gc_df, record_len, track_width, norm_factor, dinucleotide
        )
        circle_desc: str = generate_circle_path_desc(radius, norm_factor)
        clip_id_source = (
            f"{record_identifier or ''}|{record_len}|{track_width:g}|{norm_factor:g}|{dinucleotide}|{skew_desc}"
        )
        clip_id = f"clipper_circle_{hashlib.md5(clip_id_source.encode('utf-8')).hexdigest()[:12]}"
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
        group.add(circle_path)
        group.add(skew_high)
        group.add(skew_low)
        return group


__all__ = ["SkewDrawer"]


