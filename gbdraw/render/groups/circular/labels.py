#!/usr/bin/env python
# coding: utf-8

"""External (non-embedded) label layer for circular diagrams."""

from __future__ import annotations

from typing import Literal, Sequence

from gbdraw.svg.elements import Group, Line

from ....canvas import CircularCanvasGeometry
from ....labels.circular_types import CircularLabel
from ....layout.circular import CircularRecordRenderContext
from ...drawers.circular.labels import LabelDrawer, circular_label_ink
from ...drawers.label_backing import backing_rects, feature_label_backing_mode
from ...label_binding import bind_label_part


class LabelsGroup:
    """Draws non-embedded labels + leader lines as a separate SVG group."""

    def __init__(
        self,
        record_length: int,
        geometry: CircularCanvasGeometry,
        render_context: CircularRecordRenderContext,
        *,
        precalculated_labels: Sequence[CircularLabel] | None = None,
        phase: Literal["all", "leaders", "text"] = "all",
    ) -> None:
        self.record_length = int(record_length)
        self.geometry = geometry
        self.render_context = render_context
        self.precalculated_labels: Sequence[CircularLabel] | None = precalculated_labels
        self.phase = phase if phase in {"all", "leaders", "text"} else "all"
        self.show_labels = render_context.profile.labels_enabled
        self.labels_group: Group = self.setup_labels_group()

    def setup_labels_group(self) -> Group:
        group_id = {
            "all": "labels",
            "leaders": "label_leaders",
            "text": "label_text",
        }[self.phase]
        group = Group(id=group_id)
        if not self.show_labels:
            return group

        record_length = self.record_length
        feature_anchor_radius = float(self.geometry.radius)
        # The record assembly places every drawn label (and checks their fit) once.
        label_list = self.precalculated_labels or []

        drawer = LabelDrawer(profile=self.render_context.profile)
        if self.phase in {"all", "leaders"}:
            # Only the leaders read the stroke, so a stroke change keeps the text.
            cfg = self.render_context.profile.config
            stroke_width = cfg.labels.stroke_width.for_length_param(self.geometry.length_param)
            stroke_color = cfg.labels.stroke_color.label_stroke_color
        backing_mode = feature_label_backing_mode()
        if self.phase in {"all", "text"} and backing_mode != "off":
            # The text layer is in front of every leader: backings go under all of its letters (LP-D85).
            inks = [
                circular_label_ink(label, drawer, record_length)
                for label in label_list
                if not label.get("is_embedded")
            ]
            for rect in backing_rects(inks, backing_mode):
                group.add(rect)
        for label in label_list:
            if label.get("is_embedded"):
                continue
            leader_start_x = label.get("leader_start_x", label["start_x"])
            leader_start_y = label.get("leader_start_y", label["start_y"])
            if self.phase in {"all", "leaders"}:
                line_path = Line(
                    start=(label["middle_x"], label["middle_y"]),
                    end=(leader_start_x, leader_start_y),
                    stroke=stroke_color,
                    stroke_width=stroke_width,
                    stroke_linecap="round",
                )
                bind_label_part(line_path, label["feature_id"])
                group.add(line_path)
                line_path2 = Line(
                    start=(label["middle_x"], label["middle_y"]),
                    end=(
                        label.get("feature_anchor_x", label["feature_middle_x"]),
                        label.get("feature_anchor_y", label["feature_middle_y"]),
                    ),
                    stroke=stroke_color,
                    stroke_width=stroke_width,
                    stroke_linecap="round",
                )
                bind_label_part(line_path2, label["feature_id"])
                group.add(line_path2)
            if self.phase in {"all", "text"}:
                group = drawer.draw(label, group, record_length, feature_anchor_radius, self.geometry.track_ratio)

        return group

    def get_group(self) -> Group:
        return self.labels_group


__all__ = ["LabelsGroup"]
