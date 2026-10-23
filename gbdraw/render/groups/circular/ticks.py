#!/usr/bin/env python
# coding: utf-8

import math

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group, Path, Text

from ....canvas import CircularCanvasGeometry
from ....config.models import CircularRenderProfile
from ....layout.circular import RadialBand
from ....layout.circular_ink import CircularInk, CircularInkPolygon, CircularInkRing
from ....svg.circular_ticks import (
    circular_tick_labels,
    circular_tick_marks,
    get_circular_tick_intervals,
    generate_circular_tick_paths,
    generate_circular_tick_labels,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class TickGroup:
    """
    Represents a group for tick marks and labels on a circular genomic plot.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        geometry: CircularCanvasGeometry,
        *,
        profile: CircularRenderProfile,
        radius: float | None = None,
        label_side: str = "legacy",
        tick_side: str = "legacy",
        tick_length_px: float | None = None,
        track_preset: str | None = None,
        group_id: str | None = None,
        slot_id: str | None = None,
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        self.record_transform = record_transform
        self.gb_record: SeqRecord = gb_record
        self.geometry = geometry
        self.radius: float = float(radius) if radius is not None else geometry.radius
        self.tick_group = Group(id=str(group_id or "tick"))
        if slot_id:
            self.tick_group.attribs["data-gbdraw-slot-id"] = str(slot_id)
            self.tick_group.attribs["data-gbdraw-slot-renderer"] = "ticks"
        self.total_len: int = len(self.gb_record.seq)
        cfg = profile.config
        ticks_config = cfg.objects.ticks
        self.tick_width = ticks_config.tick_width
        self.stroke = ticks_config.tick_labels.stroke
        self.fill = ticks_config.tick_labels.fill
        self.font_size = ticks_config.tick_labels.font_size
        self.font_weight = ticks_config.tick_labels.font_weight
        self.manual_interval = cfg.objects.scale.interval
        self.font_family = cfg.objects.text.font_family
        self.track_type = str(track_preset or cfg.canvas.circular.track_type)
        self.length_threshold = int(cfg.labels.length_threshold.circular)
        self.separate_strands = profile.strandedness
        self.label_side = str(label_side or "legacy").strip().lower()
        self.tick_side = str(tick_side or "legacy").strip().lower()
        self.tick_length_px = tick_length_px
        self.dpi = geometry.dpi
        self.set_tick_size()
        self.add_elements_to_group()

    def set_tick_size(self) -> None:
        tick_large, tick_small = get_circular_tick_intervals(
            self.total_len,
            manual_interval=self.manual_interval,
        )
        self.tick_large: int = tick_large
        self.tick_small: int = tick_small

    def add_elements_to_group(self) -> Group:
        ticks_large = list(range(0, self.total_len, self.tick_large))
        size: str = "large"
        tick_paths_large: list[Path] = generate_circular_tick_paths(
            self.radius,
            self.total_len,
            size,
            ticks_large,
            self.tick_width,
            self.track_type,
            self.separate_strands,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            length_reference_radius_px=float(self.geometry.radius),
            record_transform=self.record_transform,
            length_threshold=self.length_threshold,
        )
        ticks_large_nonzero: list[int] = [x for x in ticks_large if x != 0]
        tick_label_paths_large: list[Text] = generate_circular_tick_labels(
            self.radius,
            self.total_len,
            size,
            ticks_large_nonzero,
            self.stroke,
            self.fill,
            self.font_size,
            self.font_weight,
            self.font_family,
            self.track_type,
            self.separate_strands,
            self.dpi,
            label_side=self.label_side,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            tick_width=float(self.tick_width),
            length_reference_radius_px=float(self.geometry.radius),
            group_identifier=str(self.tick_group.attribs.get("id", "tick")),
            record_identifier=str(self.gb_record.id),
            record_transform=self.record_transform,
            tick_interval=self.tick_large,
            length_threshold=self.length_threshold,
        )
        for tick_path_large in tick_paths_large:
            self.tick_group.add(tick_path_large)
        for tick_label_path_large in tick_label_paths_large:
            self.tick_group.add(tick_label_path_large)
        return self.tick_group

    def ink(self) -> tuple[CircularInk, ...]:
        """What the group paints: each large tick mark (butt-capped, ``tick_width`` wide) and each label's letters
        along its circle."""

        ticks_large = list(range(0, self.total_len, self.tick_large))
        half = 0.5 * float(self.tick_width)
        inks: list[CircularInk] = [
            CircularInkPolygon.segment(start, end, half)
            for start, end in circular_tick_marks(
                self.radius,
                self.total_len,
                "large",
                ticks_large,
                self.track_type,
                self.separate_strands,
                tick_side=self.tick_side,
                tick_length_px=self.tick_length_px,
                length_reference_radius_px=float(self.geometry.radius),
                record_transform=self.record_transform,
                length_threshold=self.length_threshold,
            )
        ] if half > 0 else []
        for _, tick, label in circular_tick_labels(
            self.radius,
            self.total_len,
            "large",
            [x for x in ticks_large if x != 0],
            self.font_size,
            self.font_family,
            self.track_type,
            self.separate_strands,
            self.dpi,
            label_side=self.label_side,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            tick_width=float(self.tick_width),
            length_reference_radius_px=float(self.geometry.radius),
            record_transform=self.record_transform,
            tick_interval=self.tick_large,
            length_threshold=self.length_threshold,
        ):
            # The text is centred on its tick along the circle at its path radius.
            half_turn = 0.5 * label.bbox_width_px / max(1.0, label.path_radius_px) / (2.0 * math.pi)
            turn = float(tick) / self.total_len
            inks.append(
                CircularInkRing(
                    RadialBand(label.radial_inner_px, label.radial_outer_px), turn - half_turn, turn + half_turn
                )
            )
        return tuple(inks)

    def get_group(self) -> Group:
        return self.tick_group


__all__ = ["TickGroup"]
