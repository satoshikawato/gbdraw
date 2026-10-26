#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from gbdraw.svg.elements import ClipPath, Group, Path

from ....configurators import GcSkewConfigurator
from ....layout.circular import RadialBand
from ....layout.circular_ink import CircularInk, CircularInkProfile, CircularInkRing, stroke_pad_px
from ....svg.circular_tracks import (
    circular_gc_skew_line,
    generate_circle_path_desc,
    generate_circular_gc_skew_path_desc,
)
from ....svg.ids import stable_svg_id


from gbdraw.layout.record_coordinates import RecordDisplayTransform

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

    def ink(
        self,
        radius: float,
        gc_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        dinucleotide: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> tuple[CircularInk, ...]:
        """What ``draw`` paints with the same arguments: the areas between the baseline circle and the skew line,
        with half the stroke."""

        line = circular_gc_skew_line(radius, gc_df, record_len, track_width, norm_factor, dinucleotide, record_transform)
        if line is None:
            return ()
        base = float(radius) * float(norm_factor)
        pad = stroke_pad_px(self.skew_stroke_color, self.skew_stroke_width)
        area = CircularInkProfile.between(line.positions, line.radii, record_length=record_len, base_px=base, pad_px=pad)
        return (area, CircularInkRing(RadialBand(base - pad, base + pad))) if pad > 0 else (area,)

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
        record_transform: RecordDisplayTransform | None = None,
    ) -> Group:
        skew_desc: str = generate_circular_gc_skew_path_desc(
            radius, gc_df, record_len, track_width, norm_factor, dinucleotide,
            record_transform=record_transform,
        )
        if not skew_desc:
            return group
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


