#!/usr/bin/env python
# coding: utf-8

"""External (non-embedded) label layer for circular diagrams."""

from __future__ import annotations

from typing import Literal, Optional

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group, Line, Rect

from ....canvas import CircularCanvasConfigurator
from ....features.factory import FeatureBuildResult
from ....labels.circular import prepare_label_list
from ....layout.circular import CircularRecordRenderContext
from ...drawers.circular.labels import LabelDrawer
from ...label_binding import bind_label_part


class LabelsGroup:
    """Draws non-embedded labels + leader lines as a separate SVG group."""

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        feature_layers: FeatureBuildResult,
        render_context: CircularRecordRenderContext,
        *,
        outer_arena: tuple[float, float] | None = None,
        precalculated_labels: Optional[list[dict]] = None,
        feature_track_ratio_factor_override: float | None = None,
        feature_anchor_radius_px: float | None = None,
        phase: Literal["all", "leaders", "text"] = "all",
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.feature_layers = feature_layers
        self.render_context = render_context
        self.outer_arena = outer_arena
        self.precalculated_labels: Optional[list[dict]] = precalculated_labels
        self.feature_track_ratio_factor_override = feature_track_ratio_factor_override
        self.feature_anchor_radius_px = feature_anchor_radius_px
        self.phase = phase if phase in {"all", "leaders", "text"} else "all"
        cfg = render_context.profile.config

        self.show_labels = render_context.profile.labels_enabled

        self.label_stroke_width = cfg.labels.stroke_width.for_length_param(self.canvas_config.length_param)
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

        feature_dict = self.feature_layers.foreground_features

        record_length: int = len(self.gb_record.seq)
        feature_anchor_radius = (
            float(self.feature_anchor_radius_px)
            if self.feature_anchor_radius_px is not None
            else float(self.canvas_config.radius)
        )
        if self.precalculated_labels is not None:
            label_list = self.precalculated_labels
        else:
            label_list = prepare_label_list(
                feature_dict,
                record_length,
                feature_anchor_radius,
                self.canvas_config.track_ratio,
                self.render_context.profile,
                outer_arena=self.outer_arena,
                feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                feature_layout=self.render_context.feature_layout,
                radial_layout=self.render_context.radial_layout,
                track_preset=self.render_context.track_preset,
                feature_lane_direction=str(
                    self.render_context.feature_lane_direction
                ),
            )

        drawer = LabelDrawer(profile=self.render_context.profile)
        if any(label.get("proto_bg") is not None or label.get("proto_ul") or label.get("proto_gap") is not None for label in label_list):
            return self._study_group(group, label_list, drawer, record_length, feature_anchor_radius)
        for label in label_list:
            if label.get("is_embedded"):
                continue
            leader_start_x = label.get("leader_start_x", label["start_x"])
            leader_start_y = label.get("leader_start_y", label["start_y"])
            if self.phase in {"all", "leaders"} and label.get("leader_points") is not None:
                # study prototype: polyline leaders (and extra strokes such as brackets)
                polylines = [label["leader_points"], *label.get("leader_extra", ())]
                for polyline in polylines:
                    for seg_start, seg_end in zip(polyline, polyline[1:]):
                        segment = Line(
                            start=(float(seg_start[0]), float(seg_start[1])),
                            end=(float(seg_end[0]), float(seg_end[1])),
                            stroke=self.label_stroke_color,
                            stroke_width=self.label_stroke_width,
                            stroke_linecap="round",
                        )
                        bind_label_part(segment, label["feature_id"])
                        group.add(segment)
            elif self.phase in {"all", "leaders"}:
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
                group = drawer.draw(label, group, record_length, feature_anchor_radius, self.canvas_config.track_ratio)

        return group

    def _study_group(self, group: Group, label_list: list[dict], drawer: LabelDrawer, record_length: int, anchor_radius: float) -> Group:
        """Study prototype K: all leaders, then the filled label boxes, then underlines and all text (z-order)."""
        shown = [label for label in label_list if not label.get("is_embedded")]

        def line(a, b, label):
            segment = Line(start=(float(a[0]), float(a[1])), end=(float(b[0]), float(b[1])), stroke=self.label_stroke_color,
                           stroke_width=self.label_stroke_width, stroke_linecap="round")
            bind_label_part(segment, label["feature_id"])
            return segment

        def pieces(a, b, others, pad):
            """The parts of segment a-b outside the boxes `others` grown by pad (leaders interrupted behind labels)."""
            cuts = []
            for x0, y0, x1, y1 in others:
                x0, y0, x1, y1 = x0 - pad, y0 - pad, x1 + pad, y1 + pad
                t0, t1 = 0.0, 1.0
                for d, lo, hi, st in ((b[0] - a[0], x0, x1, a[0]), (b[1] - a[1], y0, y1, a[1])):
                    if abs(d) < 1e-12:
                        if st <= lo or st >= hi:
                            t0, t1 = 1.0, 0.0
                        continue
                    ta, tb = sorted(((lo - st) / d, (hi - st) / d))
                    t0, t1 = max(t0, ta), min(t1, tb)
                if t0 < t1:
                    cuts.append((t0, t1))
            out, t = [], 0.0
            for c0, c1 in sorted(cuts):
                if c0 > t:
                    out.append((t, c0))
                t = max(t, c1)
            if t < 1.0:
                out.append((t, 1.0))
            return [((a[0] + (b[0] - a[0]) * u0, a[1] + (b[1] - a[1]) * u0), (a[0] + (b[0] - a[0]) * u1, a[1] + (b[1] - a[1]) * u1)) for u0, u1 in out if u1 - u0 > 1e-6]

        if self.phase in {"all", "leaders"}:
            for label in shown:
                if label.get("leader_points") is None:
                    continue
                gap = label.get("proto_gap")
                others = [other["aabb_local"] for other in shown if other is not label] if gap is not None else []
                for polyline in [label["leader_points"], *label.get("leader_extra", ())]:
                    for a, b in zip(polyline, polyline[1:]):
                        for u, v in (pieces(a, b, others, float(gap)) if gap is not None else [(a, b)]):
                            group.add(line(u, v, label))
        if self.phase in {"all", "text"}:
            for label in shown:
                if label.get("proto_bg") is not None:
                    x0, y0, x1, y1 = label["aabb_local"]
                    pad = float(label.get("proto_pad", 0.0))
                    box = Rect(insert=(x0 - pad, y0 - pad), size=(x1 - x0 + 2 * pad, y1 - y0 + 2 * pad), fill="#ffffff",
                               fill_opacity=float(label["proto_bg"]), stroke="none")
                    bind_label_part(box, label["feature_id"])
                    group.add(box)
            for label in shown:
                if label.get("proto_ul") and label.get("leader_points"):
                    x0, y0, x1, y1 = label["aabb_local"]
                    end = label["leader_points"][-1]
                    uy = y1 + 0.5
                    group.add(line((x0, uy), (x1, uy), label))
                    if abs(float(end[1]) - uy) > 0.3:
                        group.add(line(end, (end[0], uy), label))
            for label in label_list:
                if not label.get("is_embedded"):
                    group = drawer.draw(label, group, record_length, anchor_radius, self.canvas_config.track_ratio)
        return group

    def get_group(self) -> Group:
        return self.labels_group


__all__ = ["LabelsGroup"]
