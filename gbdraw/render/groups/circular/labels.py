#!/usr/bin/env python
# coding: utf-8

"""External (non-embedded) label layer for circular diagrams."""

from __future__ import annotations

from typing import Literal, Sequence

from gbdraw.svg.elements import Group, Line

from ....canvas import CircularCanvasGeometry
from ....labels.circular_types import CircularLabel
from ....labels.label_ink import FeatureLabelInk, circular_label_ink
from ....labels.label_style import label_leader_width_px
from ....layout.circular import CircularRecordRenderContext
from ....svg.ids import stable_svg_id
from ...drawers.circular.labels import LabelDrawer
from ...drawers.label_backing import feature_label_box, label_halo
from ...label_binding import bind_label_part

Point = tuple[float, float]


class LabelsGroup:
    """Draws the leader lines (``label_leaders``) or the text (``label_text``) of the non-embedded labels.

    The text layer is in front of every leader: per ``FeatureLabelStyle``, the
    boxes, the halos and the underlines go in it before all of its letters,
    never between a text and the leaders of ``label_leaders``.
    """

    def __init__(
        self,
        record_length: int,
        geometry: CircularCanvasGeometry,
        render_context: CircularRecordRenderContext,
        *,
        precalculated_labels: Sequence[CircularLabel] | None = None,
        phase: Literal["leaders", "text"],
    ) -> None:
        self.record_length = int(record_length)
        self.geometry = geometry
        self.render_context = render_context
        self.precalculated_labels: Sequence[CircularLabel] | None = precalculated_labels
        self.phase = phase
        self.show_labels = render_context.profile.labels_enabled
        self.labels_group: Group = self.setup_labels_group()

    def _part(self, start: Point, end: Point, label: CircularLabel) -> Line:
        """One straight part of a leader or of its underline: one style for both (LP-D95)."""

        part = Line(
            start=start,
            end=end,
            stroke=self.render_context.profile.config.labels.stroke_color.label_stroke_color,
            stroke_width=label_leader_width_px(label),
            stroke_linecap="round",
        )
        bind_label_part(part, label["feature_id"])
        return part

    def setup_labels_group(self) -> Group:
        group = Group(id="label_leaders" if self.phase == "leaders" else "label_text")
        if not self.show_labels:
            return group

        # The record assembly places every drawn label (and checks their fit) once.
        labels = [label for label in self.precalculated_labels or [] if not label.get("is_embedded")]
        style = self.render_context.profile.label_style
        drawn: list[tuple[CircularLabel, FeatureLabelInk | None, Point, tuple[Point, Point, Point] | None]] = []
        for label in labels:
            contact: Point = (label.get("leader_start_x", label["start_x"]), label.get("leader_start_y", label["start_y"]))
            ink = circular_label_ink(label, self.record_length, style) if style.draws_from_ink else None
            line = ink.underline if ink is not None else None
            # The leader ends on its underline, which continues it (LP-D88), or stops above the letters.
            drawn.append((label, ink, line[0] if line is not None else contact, line))

        if self.phase == "leaders":
            for label, _ink, contact, _line in drawn:
                bend = (label["middle_x"], label["middle_y"])
                group.add(self._part(bend, contact, label))
                feature = (
                    label.get("feature_anchor_x", label["feature_middle_x"]),
                    label.get("feature_anchor_y", label["feature_middle_y"]),
                )
                group.add(self._part(bend, feature, label))
            return group

        for label, ink, _contact, _line in drawn:
            if ink is not None and ink.draws_box(style):
                group.add(feature_label_box(ink, label["feature_id"]))
        drawer = LabelDrawer(profile=self.render_context.profile)
        feature_anchor_radius = float(self.geometry.radius)
        letters = Group()
        for label, _ink, _contact, _line in drawn:
            drawer.draw(label, letters, self.record_length, feature_anchor_radius, self.geometry.track_ratio)
        if style.protection == "halo":
            for (label, _ink, _contact, _line), text in zip(drawn, letters.elements):
                group.add(label_halo(text, stable_svg_id("label_text", label["stable_id"], label["input_order"])))
        for label, _ink, _contact, line in drawn:
            if line is not None and line[1] != line[2]:
                group.add(self._part(line[1], line[2], label))
        group.elements.extend(letters.elements)
        return group

    def get_group(self) -> Group:
        return self.labels_group


__all__ = ["LabelsGroup"]
