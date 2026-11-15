#!/usr/bin/env python
# coding: utf-8

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from typing import ParamSpec

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group, Path, Text

from ....canvas import CircularCanvasGeometry
from ....config.models import CircularRenderProfile
from ....labels.label_style import scale_label_stride
from ....layout.circular import RadialBand
from ....layout.circular_ink import CircularInk, CircularInkPolygon, CircularInkRing, stroke_meets_ink
from ....render.drawers.circular.annotations import LeaderSegment
from ....render.drawers.label_backing import BACKING_PADDING_EM, SectorBacking, sector_backing
from ....render.group_content import GroupContent
from ....svg.circular_arcs import position_angle_rad, position_turn_rad
from ....svg.ids import track_slot_svg_id
from ....svg.circular_ticks import (
    CircularTickLabelGeometry,
    circular_tick_labels,
    circular_tick_marks,
    get_circular_tick_intervals,
    generate_circular_tick_paths,
    generate_circular_tick_labels,
)
from ....svg.circular_tracks import within_record_arc


from gbdraw.layout.record_coordinates import CirclePlace, RecordDisplayTransform
from gbdraw.layout.record_transform import circular_display_length


P = ParamSpec("P")


def _closest_label_gap_px(
    labels: Sequence[tuple[str, float, CircularTickLabelGeometry]], total_len: int
) -> float | None:
    """The smallest gap between neighbouring tick labels, round the circle: along the circle at the smaller of
    their path radii, between the ends of their text (each centred on its tick); None for fewer than two."""

    if len(labels) < 2:
        return None
    arcs = sorted(
        (float(tick) / total_len - half_turn, float(tick) / total_len + half_turn, radius)
        for _text, tick, geometry in labels
        for radius in [max(1.0, geometry.path_radius_px)]
        for half_turn in [0.5 * geometry.bbox_width_px / radius / (2.0 * math.pi)]
    )
    following = [*arcs[1:], (arcs[0][0] + 1.0, arcs[0][1] + 1.0, arcs[0][2])]
    return min(
        position_turn_rad(right[0] - left[1], 1.0) * min(left[2], right[2]) for left, right in zip(arcs, following)
    )


@dataclass(frozen=True)
class TickNumber:
    """One scale number: its letters along its circle, and the backing it is drawn on when a feature label's leader
    runs behind it or within the backing's padding of it (OV-3301)."""

    letters: CircularInkRing
    backing: SectorBacking


@dataclass(frozen=True)
class TickSlotLayout:
    """A tick slot's ticks and numbers before any backing (render stage ``ticks_layout``): ``content`` holds the
    ``mark_count`` tick marks, then the numbers; ``pad_px`` is how near a leader's stroke may come to a number's
    letters before the number is backed."""

    content: GroupContent
    mark_count: int
    numbers: tuple[TickNumber, ...]
    pad_px: float


def backed_tick_numbers(layout: TickSlotLayout, feature_leaders: Sequence[LeaderSegment]) -> tuple[int, ...]:
    """The numbers (their places in ``layout.numbers``) a feature label's leader runs behind or within
    ``layout.pad_px`` of (render stage ``ticks_backed_numbers``)."""

    return tuple(
        index
        for index, number in enumerate(layout.numbers)
        if any(
            stroke_meets_ink(leader.start, leader.end, leader.half_width + layout.pad_px, number.letters)
            for leader in feature_leaders
        )
    )


def backed_tick_content(layout: TickSlotLayout, backed: tuple[int, ...]) -> GroupContent:
    """The tick slot's content (render stage ``ticks_content``): the tick marks, the backings of the ``backed``
    numbers, then the numbers, so a leader drawn before the slot passes behind a backed number."""

    if not backed:
        return layout.content
    children = layout.content.children
    backings = tuple(layout.numbers[index].backing.element() for index in backed)
    return replace(
        layout.content, children=(*children[: layout.mark_count], *backings, *children[layout.mark_count :])
    )


def build_tick_slot_layout(
    group_class: Callable[P, "TickGroup"], /, *args: P.args, **kwargs: P.kwargs
) -> TickSlotLayout:
    """A render-stage function (``ticks_layout``): build a tick slot's group with ``group_class`` (``TickGroup``)
    and keep its layout, as ``build_group_content`` keeps a group's content."""

    return group_class(*args, **kwargs).layout()


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
        slot_id: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        self.record_transform = record_transform
        self.gb_record: SeqRecord = gb_record
        self.geometry = geometry
        self.radius: float = float(radius) if radius is not None else geometry.radius
        self.tick_group = Group(id=track_slot_svg_id(slot_id, renderer="ticks"))
        self.tick_group.attribs["data-gbdraw-slot-id"] = slot_id
        self.tick_group.attribs["data-gbdraw-slot-renderer"] = "ticks"
        self.total_len: int = circular_display_length(gb_record, record_transform)
        place = None if record_transform is None else record_transform.circle
        # The positions ticks run over: the record's own bases (on a shared circle, numbered from its start, at the
        # interval and in the unit of the circle's longest record, with no tick at the record's start).
        self.tick_extent: int = len(self.gb_record.seq)
        self.scale_length: int = self.total_len if place is None else place.scale_length
        self.first_tick: int = 0
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
        self.set_tick_size()
        self.set_labelled_ticks()
        self.add_elements_to_group()

    def set_tick_size(self) -> None:
        tick_large, tick_small = get_circular_tick_intervals(
            self.scale_length,
            manual_interval=self.manual_interval,
        )
        self.tick_large: int = tick_large
        self.tick_small: int = tick_small
        if self.scale_length != self.total_len:
            self.first_tick = tick_large

    def _tick_labels(self, ticks: list[int]) -> list[tuple[str, float, CircularTickLabelGeometry]]:
        return circular_tick_labels(
            self.radius,
            self.total_len,
            "large",
            ticks,
            self.font_size,
            self.font_family,
            self.track_type,
            self.separate_strands,
            label_side=self.label_side,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            tick_width=float(self.tick_width),
            length_reference_radius_px=float(self.geometry.radius),
            record_transform=self.record_transform,
            tick_interval=self.tick_large,
            length_threshold=self.length_threshold,
        )

    def set_labelled_ticks(self) -> None:
        """Which large ticks carry a label (OV-3020): every ``label_stride``-th (``scale_label_stride``), with their
        gaps measured along the circle. Every tick keeps its mark at the interval, Auto or typed."""

        ticks = [tick for tick in range(0, self.tick_extent, self.tick_large) if tick != 0]
        # On a circle several records share, every record labels the same ticks of its own numbering: the
        # stride the circle's longest record would take at the circle's scale, so neighbours' labels keep apart.
        probe = (
            self._tick_labels(ticks)
            if self.scale_length == self.total_len
            else circular_tick_labels(
                self.radius, self.total_len, "large",
                [tick for tick in range(self.tick_large, self.scale_length, self.tick_large)],
                self.font_size, self.font_family, self.track_type, self.separate_strands,
                label_side=self.label_side, tick_side=self.tick_side, tick_length_px=self.tick_length_px,
                tick_width=float(self.tick_width), length_reference_radius_px=float(self.geometry.radius),
                # The circle's longest record, from the circle's start: its labels' text and spacing.
                record_transform=RecordDisplayTransform(self.scale_length, 1, 1, None, False).with_circle(
                    CirclePlace(0, self.total_len, self.scale_length), self.scale_length
                ),
                tick_interval=self.tick_large, length_threshold=self.length_threshold,
            )
        )
        self.label_stride: int = scale_label_stride(
            lambda stride: _closest_label_gap_px(probe[stride - 1 :: stride], self.total_len), self.font_size
        )
        self.labelled_ticks: list[int] = self._labels_within_record(ticks[self.label_stride - 1 :: self.label_stride])

    def _labels_within_record(self, ticks: list[int]) -> list[int]:
        """The ticks whose label lies within its record's arc on a circle several records share (OV-3113): a label
        that would run past the record's start or end, over the gap and onto its neighbour, is left out; its tick
        keeps its mark. A record alone keeps every label."""

        transform = self.record_transform
        labels = [] if transform is None or transform.circle is None else self._tick_labels(ticks)
        if not labels:  # A record alone, or no labels drawn at all.
            return ticks
        kept = []
        for (_, position, label), tick in zip(labels, ticks, strict=True):
            if within_record_arc(*self._label_turns(position, label), transform):
                kept.append(tick)
        return kept

    def _large_ticks(self) -> list[int]:
        """The positions of the large tick marks: every interval over the record's own bases (on a shared circle,
        from one interval past the record's start)."""

        return list(range(self.first_tick, self.tick_extent, self.tick_large))

    def _label_turns(self, tick: float, label: CircularTickLabelGeometry) -> tuple[float, float]:
        """The turns of the circle from the top that a label's letters span, centred on its tick (its display
        place) at its path radius; not wrapped past 1."""

        half_turn = 0.5 * label.bbox_width_px / max(1.0, label.path_radius_px) / (2.0 * math.pi)
        turn = float(tick) / self.total_len
        return turn - half_turn, turn + half_turn

    def add_elements_to_group(self) -> Group:
        size: str = "large"
        tick_paths_large: list[Path] = generate_circular_tick_paths(
            self.radius,
            self.total_len,
            size,
            self._large_ticks(),
            self.tick_width,
            self.track_type,
            self.separate_strands,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            length_reference_radius_px=float(self.geometry.radius),
            record_transform=self.record_transform,
            length_threshold=self.length_threshold,
        )
        tick_label_paths_large: list[Text] = generate_circular_tick_labels(
            self.radius,
            self.total_len,
            size,
            self.labelled_ticks,
            self.stroke,
            self.fill,
            self.font_size,
            self.font_weight,
            self.font_family,
            self.track_type,
            self.separate_strands,
            label_side=self.label_side,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            tick_width=float(self.tick_width),
            length_reference_radius_px=float(self.geometry.radius),
            group_identifier=str(self.tick_group.attribs["id"]),
            record_identifier=str(self.gb_record.id),
            record_transform=self.record_transform,
            tick_interval=self.tick_large,
            length_threshold=self.length_threshold,
        )
        for tick_path_large in tick_paths_large:
            self.tick_group.add(tick_path_large)
        # The numbers' backings go between the marks and the numbers (``backed_tick_content``).
        self.mark_count: int = len(tick_paths_large)
        for tick_label_path_large in tick_label_paths_large:
            self.tick_group.add(tick_label_path_large)
        return self.tick_group

    def _tick_marks(self) -> list[tuple[tuple[float, float], tuple[float, float]]]:
        return circular_tick_marks(
            self.radius,
            self.total_len,
            "large",
            self._large_ticks(),
            self.track_type,
            self.separate_strands,
            tick_side=self.tick_side,
            tick_length_px=self.tick_length_px,
            length_reference_radius_px=float(self.geometry.radius),
            record_transform=self.record_transform,
            length_threshold=self.length_threshold,
        )

    def _label_rings(self) -> list[CircularInkRing]:
        """Each label's letters along its circle: centred on its tick at its path radius, as the turns of the circle
        from the top (a ring across the top runs from its start past 1 to its end)."""

        rings: list[CircularInkRing] = []
        for _, tick, label in self._tick_labels(self.labelled_ticks):
            start, end = self._label_turns(tick, label)
            band = RadialBand(label.radial_inner_px, label.radial_outer_px)
            rings.append(CircularInkRing(band, start % 1.0, end % 1.0))
        return rings

    def _numbers(self) -> tuple[TickNumber, ...]:
        """Each number's letters and its backing (OV-3301): as a region label's (``sector_backing``,
        CIR-RGN-FRONT), stopped before the tick marks and their stroke on the side they lie."""

        rings = self._label_rings()
        if not rings:
            return ()
        radii = [math.hypot(*point) for mark in self._tick_marks() for point in mark]
        half = 0.5 * float(self.tick_width)
        numbers: list[TickNumber] = []
        for ring in rings:
            outside = not radii or ring.band.inner_px >= max(radii)
            edge = (max(radii) + half if radii else 0.0) if outside else min(radii) - half
            low = position_angle_rad(ring.start_turn, 1.0)
            high = low + position_turn_rad((ring.end_turn - ring.start_turn) % 1.0, 1.0)
            backing = sector_backing(ring.band.inner_px, ring.band.outer_px, low, high, float(self.font_size), edge=edge)
            numbers.append(TickNumber(ring, backing))
        return tuple(numbers)

    def layout(self) -> TickSlotLayout:
        """The group's content, its numbers and their backings, before any number is backed."""

        return TickSlotLayout(
            GroupContent.of(self.tick_group),
            self.mark_count,
            self._numbers(),
            BACKING_PADDING_EM * float(self.font_size),
        )

    def ink(self) -> tuple[CircularInk, ...]:
        """What the group paints: each large tick mark (butt-capped, ``tick_width`` wide) and each label's letters
        along its circle."""

        half = 0.5 * float(self.tick_width)
        inks: list[CircularInk] = (
            [CircularInkPolygon.segment(start, end, half) for start, end in self._tick_marks()] if half > 0 else []
        )
        inks.extend(self._label_rings())
        return tuple(inks)

    def get_group(self) -> Group:
        return self.tick_group


__all__ = [
    "TickGroup",
    "TickNumber",
    "TickSlotLayout",
    "backed_tick_content",
    "backed_tick_numbers",
    "build_tick_slot_layout",
]
