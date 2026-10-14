#!/usr/bin/env python
# coding: utf-8

"""External (non-embedded) label layer for circular diagrams."""

from __future__ import annotations

from typing import Literal, Optional

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group, Line

from ....canvas import CircularCanvasGeometry
from ....labels.circular_types import CircularLabel
from ....layout.circular import CircularRecordRenderContext
from ...drawers.circular.labels import LabelDrawer
from ...label_binding import bind_label_part


class LabelsGroup:
    """Draws non-embedded labels + leader lines as a separate SVG group."""

    def __init__(
        self,
        gb_record: SeqRecord,
        geometry: CircularCanvasGeometry,
        render_context: CircularRecordRenderContext,
        *,
        outer_arena: tuple[float, float] | None = None,
        precalculated_labels: Optional[list[CircularLabel]] = None,
        feature_track_ratio_factor_override: float | None = None,
        feature_anchor_radius_px: float | None = None,
        phase: Literal["all", "leaders", "text"] = "all",
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.geometry = geometry
        self.render_context = render_context
        self.outer_arena = outer_arena
        self.precalculated_labels: Optional[list[CircularLabel]] = precalculated_labels
        self.feature_track_ratio_factor_override = feature_track_ratio_factor_override
        self.feature_anchor_radius_px = feature_anchor_radius_px
        self.phase = phase if phase in {"all", "leaders", "text"} else "all"
        cfg = render_context.profile.config

        self.show_labels = render_context.profile.labels_enabled

        self.label_stroke_width = cfg.labels.stroke_width.for_length_param(geometry.length_param)
        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color

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

        record_length: int = len(self.gb_record.seq)
        feature_anchor_radius = (
            float(self.feature_anchor_radius_px)
            if self.feature_anchor_radius_px is not None
            else float(self.geometry.radius)
        )
        # The record assembly places every drawn label (and checks their fit) once.
        label_list = self.precalculated_labels or []

        drawer = LabelDrawer(profile=self.render_context.profile)
        for label in label_list:
            if label.get("is_embedded"):
                continue
            leader_start_x = label.get("leader_start_x", label["start_x"])
            leader_start_y = label.get("leader_start_y", label["start_y"])
            if self.phase in {"all", "leaders"}:
                line_path = Line(
                    start=(label["middle_x"], label["middle_y"]),
                    end=(leader_start_x, leader_start_y),
                    stroke=self.label_stroke_color,
                    stroke_width=self.label_stroke_width,
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
                    stroke=self.label_stroke_color,
                    stroke_width=self.label_stroke_width,
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
