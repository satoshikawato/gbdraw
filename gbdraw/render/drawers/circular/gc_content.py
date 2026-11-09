#!/usr/bin/env python
# coding: utf-8


import math

from pandas import DataFrame
from gbdraw.svg.elements import Circle, Group, Path

from ....auto_sizes import depth_axis_tick_font_size_px
from ....analysis.gc import gc_content_percent_df
from ....layout.circular import RadialBand
from ....layout.circular_depth_axis import (
    DEPTH_AXIS_SMALL_TICK_SIZE_PX,
    DEPTH_AXIS_STROKE_WIDTH_PX,
    DEPTH_AXIS_TICK_SIZE_PX,
    CircularRangeLabel,
    ScalarAxisTick,
    circular_range_label,
    scalar_axis_ink,
)
from ....layout.scalar_axis import (
    percent_axis_labels,
    scalar_axis_small_tick_values,
    scaled_scalar_fraction,
)
from ....layout.circular_ink import CircularInk, CircularInkProfile, CircularInkRing, paint_shows, stroke_pad_px
from ....svg.circular_arcs import circle_arc
from ....svg.circular_tracks import (
    circle_record_turns,
    circle_unit_arc,
    circle_unit_turns,
    circular_gc_content_line,
    generate_circular_gc_content_path_desc,
    circular_scalar_area_line,
    generate_circular_scalar_area_path_desc,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform
from .scalar_axis import circular_range_label_elements, circular_scalar_axis_group


def _draws_scalar_axes(record_transform: RecordDisplayTransform) -> bool:
    """Whether the record draws its track's axis: always alone, once on a circle several records share."""

    return record_transform.circle is None or record_transform.circle.draws_scalar_axes


class GcContentDrawer:
    """
    Draws GC content track on a circular canvas.
    """

    def __init__(self, gc_config) -> None:
        self._gc_config = gc_config
        self.gc_path_high_fill_color: str = gc_config.high_fill_color
        self.gc_path_low_fill_color: str = gc_config.low_fill_color
        self.gc_path_fill_color: str = gc_config.fill_color
        self.gc_path_stroke_color: str = gc_config.stroke_color
        self.gc_path_stroke_width: float = gc_config.stroke_width
        self.gc_path_fill_opacity: float = gc_config.fill_opacity
        self.mode: str = str(getattr(gc_config, "mode", "deviation"))
        min_percent = getattr(gc_config, "min_percent", None)
        max_percent = getattr(gc_config, "max_percent", None)
        self.min_percent: float = float(0.0 if min_percent is None else min_percent)
        self.max_percent: float = float(100.0 if max_percent is None else max_percent)
        self.show_axis: bool = bool(getattr(gc_config, "show_axis", True))
        self.show_ticks: bool = bool(getattr(gc_config, "show_ticks", True))
        self.large_tick_interval: float | None = getattr(gc_config, "large_tick_interval", 20.0)
        self.small_tick_interval: float | None = getattr(gc_config, "small_tick_interval", None)
        self.tick_font_size: float | None = getattr(gc_config, "tick_font_size", None)
        self.font_family: str = gc_config.font_family
        self.percent_background_color: str = str(gc_config.percent_background_color)
        self.percent_background_opacity: float = float(getattr(gc_config, "percent_background_opacity", 1.0))
        self.percent_border_color: str = str(getattr(gc_config, "percent_border_color", "#4b5563"))
        self.percent_border_width: float = float(getattr(gc_config, "percent_border_width", 0.8))
        self.axis_stroke: str = "#4b5563"
        self.axis_stroke_width: float = DEPTH_AXIS_STROKE_WIDTH_PX
        self.axis_tick_size: float = DEPTH_AXIS_TICK_SIZE_PX
        self.axis_small_tick_size: float = DEPTH_AXIS_SMALL_TICK_SIZE_PX

    def _paints_background(self) -> bool:
        return paint_shows(self.percent_background_color, self.percent_background_opacity)

    def _paints_border(self) -> bool:
        return self.percent_border_width > 0 and paint_shows(self.percent_border_color)

    def paints_percent_band_edges(self) -> bool:
        """Whether the track paints its band's edges: in percent mode, its background ring or its border circles
        show. The tick labels next to the band keep their margin from it (OV-3116)."""

        return self.mode == "percent" and (self._paints_background() or self._paints_border())

    def _tick_font_size(self, track_width: float) -> float:
        return depth_axis_tick_font_size_px(self._gc_config, track_width)

    @staticmethod
    def _percent_radii(radius: float, track_width: float, norm_factor: float) -> tuple[float, float]:
        baseline_radius = max(
            0.0,
            (float(radius) * float(norm_factor)) - (0.5 * float(track_width)),
        )
        return baseline_radius, baseline_radius + float(track_width)

    def _add_percent_container(
        self,
        group: Group,
        radius: float,
        track_width: float,
        norm_factor: float,
        group_id: str,
        unit_arc: tuple[float, float] | None = None,
    ) -> Group:
        baseline_radius, outer_radius = self._percent_radii(radius, track_width, norm_factor)
        middle_radius = baseline_radius + (0.5 * float(track_width))
        if unit_arc is not None:
            # A record that shares its circle: the background is the record's own arc of the ring.
            group.add(
                Path(
                    id=f"{group_id}_percent_background",
                    d=circle_arc(middle_radius, *unit_arc),
                    fill="none",
                    stroke=self.percent_background_color,
                    stroke_width=float(track_width),
                    stroke_opacity=self.percent_background_opacity,
                )
            )
            return group
        group.add(
            Circle(
                id=f"{group_id}_percent_background",
                center=(0, 0),
                r=middle_radius,
                fill="none",
                stroke=self.percent_background_color,
                stroke_width=float(track_width),
                stroke_opacity=self.percent_background_opacity,
            )
        )
        return group

    def _add_percent_border(
        self,
        group: Group,
        baseline_radius: float,
        outer_radius: float,
        group_id: str,
        unit_arc: tuple[float, float] | None = None,
    ) -> Group:
        if self.percent_border_width <= 0:
            return group
        border_group = Group(id=f"{group_id}_percent_border")
        for border_radius in (baseline_radius, outer_radius):
            if unit_arc is not None:
                border_group.add(
                    Path(
                        d=circle_arc(border_radius, *unit_arc),
                        fill="none",
                        stroke=self.percent_border_color,
                        stroke_width=self.percent_border_width,
                    )
                )
                continue
            border_group.add(
                Circle(
                    center=(0, 0),
                    r=border_radius,
                    fill="none",
                    stroke=self.percent_border_color,
                    stroke_width=self.percent_border_width,
                )
            )
        group.add(border_group)
        return group

    def _percent_axis_ticks(
        self, baseline_radius: float, track_width: float, record_transform: RecordDisplayTransform | None
    ) -> tuple[list[ScalarAxisTick], CircularRangeLabel | None]:
        """The percent axis ticks as drawn, small ticks first, then large ticks with their numbers; and the range
        label this record draws: on a circle several records share, the label lies on the first record that holds
        it (``on_shared_circle``), never over a gap between records, and only that record draws it, after its own
        band, so no later record paints over it."""

        font_size = self._tick_font_size(track_width)
        axis_labels = percent_axis_labels(
            self.min_percent,
            self.max_percent,
            label_style=self._gc_config.axis_labels,
            show_ticks=self.show_ticks,
            large_tick_interval=self.large_tick_interval,
            length_px=track_width,
            label_line_px=font_size,
        )
        large_ticks = axis_labels.ticks
        ticks = [
            ScalarAxisTick(
                baseline_radius + float(track_width) * scaled_scalar_fraction(value, self.min_percent, self.max_percent),
                self.axis_small_tick_size,
                None,
            )
            for value in scalar_axis_small_tick_values(
                self.min_percent,
                self.max_percent,
                show_ticks=self.show_ticks,
                small_tick_interval=self.small_tick_interval,
                large_tick_values=[tick.value for tick in large_ticks],
            )
        ]
        ticks.extend(
            ScalarAxisTick(baseline_radius + tick.offset_px, self.axis_tick_size, tick.label) for tick in large_ticks
        )
        range_label = circular_range_label(
            axis_labels.range_label, baseline_radius + float(track_width), font_size=font_size
        )
        place = None if record_transform is None else record_transform.circle
        if range_label is not None and place is not None:
            record_turns = circle_record_turns(place)
            range_label = range_label.on_shared_circle(record_turns, self.font_family)
            own_start, own_end = record_turns[place.position]
            if range_label is not None and not own_start <= range_label.start_rad / (2.0 * math.pi) < own_end:
                range_label = None
        return ticks, range_label

    def _add_percent_axis(
        self,
        group: Group,
        radius: float,
        track_width: float,
        norm_factor: float,
        axis_group_id: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> Group:
        if not self.show_axis:
            return group

        baseline_radius, outer_radius = self._percent_radii(radius, track_width, norm_factor)
        ticks, range_label = self._percent_axis_ticks(baseline_radius, track_width, record_transform)
        axis_group = circular_scalar_axis_group(
            axis_group_id,
            baseline_radius,
            outer_radius,
            ticks,
            range_label,
            stroke=self.axis_stroke,
            stroke_width=self.axis_stroke_width,
            font_size=self._tick_font_size(track_width),
            font_family=self.font_family,
            # A record of a shared circle draws its baseline as its own border arc.
            baseline_circle=record_transform is None or record_transform.circle is None,
        )

        group.add(axis_group)
        return group

    def ink(
        self,
        radius: float,
        gc_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        dinucleotide: str,
        record_transform: RecordDisplayTransform,
    ) -> tuple[CircularInk, ...]:
        """What ``draw`` paints with the same arguments: the area between the baseline circle and the content line
        with half its stroke; in percent mode the area, the background ring when it shows, the border circles and
        the axis, its labels measured in the figure's font as ``draw`` writes them."""

        pad = stroke_pad_px(self.gc_path_stroke_color, self.gc_path_stroke_width)
        if self.mode != "percent":
            line = circular_gc_content_line(
                radius, record_len, gc_df, track_width, norm_factor, dinucleotide, record_transform
            )
            if line is None:
                return ()
            base = float(radius) * float(norm_factor)
            area = CircularInkProfile.between(line.positions, line.radii, record_length=record_len, base_px=base, pad_px=pad)
            # The baseline closes along the record's own arc (``circle_unit_turns``), the whole circle for a record alone.
            baseline = CircularInkRing(RadialBand(base - pad, base + pad), *circle_unit_turns(record_transform))
            return (area, baseline) if pad > 0 else (area,)
        plot_df = gc_content_percent_df(
            gc_df, dinucleotide=dinucleotide, min_percent=self.min_percent, max_percent=self.max_percent
        )
        scalar_line = circular_scalar_area_line(
            radius, record_len, plot_df, track_width, norm_factor,
            close_at_record_len=True, record_transform=record_transform,
        )
        if scalar_line is None:
            return ()
        positions, outer_radii, baseline_radius = scalar_line.positions, scalar_line.radii, scalar_line.baseline_radius
        inks: list[CircularInk] = [
            CircularInkProfile.between(positions, outer_radii, record_length=record_len, base_px=baseline_radius,
                                       pad_px=pad)
        ]
        outer_radius = baseline_radius + float(track_width)
        unit_arc = circle_unit_arc(record_transform)
        # The arc of the ring a record that shares its circle paints; the whole ring else.
        turns = circle_unit_turns(record_transform)
        if float(track_width) > 0 and self._paints_background():
            inks.append(CircularInkRing(RadialBand(baseline_radius, outer_radius), *turns))
        if self._paints_border():
            half = 0.5 * self.percent_border_width
            inks.extend(CircularInkRing(RadialBand(edge - half, edge + half), *turns) for edge in (baseline_radius, outer_radius))
        if self.show_axis and _draws_scalar_axes(record_transform):
            ticks, range_label = self._percent_axis_ticks(baseline_radius, track_width, record_transform)
            inks.extend(
                scalar_axis_ink(
                    baseline_radius,
                    outer_radius,
                    ticks,
                    font_size=self._tick_font_size(track_width),
                    font_family=self.font_family,
                    range_label=range_label,
                    baseline_circle=unit_arc is None,
                )
            )
        elif self.show_axis:
            _, range_label = self._percent_axis_ticks(baseline_radius, track_width, record_transform)
            if range_label is not None:
                inks.append(range_label.ink(self.font_family))
        return tuple(inks)

    def draw(
        self,
        radius: float,
        group: Group,
        gc_df: DataFrame,
        record_len: int,
        track_width: float,
        norm_factor: float,
        dinucleotide: str,
        group_id: str,
        *,
        record_transform: RecordDisplayTransform,
    ) -> Group:
        if self.mode == "percent":
            plot_df = gc_content_percent_df(
                gc_df,
                dinucleotide=dinucleotide,
                min_percent=self.min_percent,
                max_percent=self.max_percent,
            )
            gc_path_desc: str = generate_circular_scalar_area_path_desc(
                radius,
                record_len,
                plot_df,
                track_width,
                norm_factor,
                close_at_record_len=True,
                record_transform=record_transform,
            )
        else:
            gc_path_desc = generate_circular_gc_content_path_desc(
                radius, record_len, gc_df, track_width, norm_factor, dinucleotide,
                record_transform=record_transform,
            )
        if not gc_path_desc:
            return group
        unit_arc = circle_unit_arc(record_transform)
        if self.mode == "percent":
            self._add_percent_container(group, radius, track_width, norm_factor, group_id, unit_arc)
        gc_path: Path = Path(
            d=gc_path_desc,
            fill=self.gc_path_fill_color,
            stroke=self.gc_path_stroke_color,
            fill_opacity=self.gc_path_fill_opacity,
            fill_rule="evenodd",
        )
        if (
            self.gc_path_stroke_width > 0
            or self.gc_path_stroke_color.strip().lower() != "none"
        ):
            gc_path["stroke-width"] = self.gc_path_stroke_width
        group.add(gc_path)
        if self.mode == "percent":
            baseline_radius, outer_radius = self._percent_radii(radius, track_width, norm_factor)
            self._add_percent_border(group, baseline_radius, outer_radius, group_id, unit_arc)
            if _draws_scalar_axes(record_transform):
                self._add_percent_axis(group, radius, track_width, norm_factor, f"{group_id}_axis", record_transform)
            elif self.show_axis:
                self._add_shared_circle_range_label(group, baseline_radius, track_width, f"{group_id}_axis", record_transform)
        return group

    def _add_shared_circle_range_label(
        self, group: Group, baseline_radius: float, track_width: float, axis_group_id: str,
        record_transform: RecordDisplayTransform,
    ) -> None:
        """The range label of the circle's axis, drawn by a record other than the axis's when its arc holds it."""

        _, range_label = self._percent_axis_ticks(baseline_radius, track_width, record_transform)
        if range_label is None:
            return
        label_group = Group(id=axis_group_id)
        for element in circular_range_label_elements(
            range_label, axis_group_id, stroke=self.axis_stroke,
            font_size=self._tick_font_size(track_width), font_family=self.font_family,
        ):
            label_group.add(element)
        group.add(label_group)


__all__ = ["GcContentDrawer"]


