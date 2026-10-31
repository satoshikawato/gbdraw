#!/usr/bin/env python
# coding: utf-8

from typing import Tuple

from gbdraw.svg.elements import Path, Text, TextPath
from gbdraw.svg.circular_arcs import circle_arc, position_angle

from ....config.models import CircularRenderProfile
from ....auto_sizes import determine_length_parameter
from ....layout.circular import calculate_feature_position_factors_circular
from ....layout.common import calculate_cds_ratio
from ....svg.ids import stable_svg_id
from ....svg.text_path import generate_text_path
from ....labels.circular_leaders import TextAnchor, TextBaseline, horizontal_label_alignment, rim_text_placement
from ....labels.circular_types import CircularLabel
from ...label_binding import bind_label_part


class LabelDrawer:
    def __init__(self, *, profile: CircularRenderProfile) -> None:
        cfg = profile.config
        self._profile = profile
        self._cfg = cfg
        self.strandedness = profile.strandedness

    def set_feature_label_anchor_value(
        self,
        total_len: int,
        tick: float,
        start_x: float,
        is_inner: bool = False,
    ) -> Tuple[TextAnchor, TextBaseline]:
        return horizontal_label_alignment(total_len, tick, start_x, is_inner=is_inner)

    def embed_label(
        self,
        group,
        label: CircularLabel,
        radius,
        record_length,
        track_ratio,
        feature_track_ratio_factor_override: float | None = None,
    ):
        length_param = determine_length_parameter(record_length, self._cfg.labels.length_threshold.circular)
        track_ratio_factor = (
            float(feature_track_ratio_factor_override)
            if feature_track_ratio_factor_override is not None
            else float(self._cfg.canvas.circular.track_ratio_factors[length_param][0])
        )
        cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
        # Get track_id from label for overlap resolution
        track_id = label.get("track_id", 0)
        explicit_center_radius = label.get("feature_center_radius_px")
        if explicit_center_radius is None:
            factors: list[float] = calculate_feature_position_factors_circular(
                record_length,
                label["strand"],
                track_ratio,
                cds_ratio,
                offset,
                self.track_type,
                self.strandedness,
                track_id,
            )
            label_radius = float(radius) * float(factors[1])
        else:
            label_radius = float(explicit_center_radius)
        angle = 360.0 * (label["middle"] / record_length)
        # A quarter of the label's line box: the same offset for every label of one face and size (OV-1743).
        center_offset = float(label["height_px"]) / 4

        # The text path runs in reading order: clockwise from the label start, or anticlockwise from its end
        # in the lower half so the text stays upright.
        span = 360.0 * (float(label["end"]) - float(label["start"])) / record_length
        if 90 <= angle < 270:
            label_axis_path_desc = circle_arc(
                label_radius + center_offset, position_angle(label["end"], record_length), -span
            )
        else:
            label_axis_path_desc = circle_arc(
                label_radius - center_offset, position_angle(label["start"], record_length), span
            )

        group_identifier = str(getattr(group, "attribs", {}).get("id", "labels"))
        label_axis_path = Path(
            id=stable_svg_id(
                "circular_feature_label_path",
                group_identifier,
                label.get("stable_id"),
                label.get("input_order"),
                len(getattr(group, "elements", ()) or ()),
                record_length,
                track_id,
                label["label_text"],
                label_axis_path_desc,
            ),
            d=label_axis_path_desc,
            stroke="none",
            fill="none",
        )
        text_path = Text("")
        text_path.add(
            TextPath(
                label_axis_path,
                text=label["label_text"],
                startOffset="50%",
                method="align",
                text_anchor="middle",
                font_size=self.font_size,
                font_style=label["font_style"],
                font_weight="normal",
                font_family=self.font_family,
                dominant_baseline="auto",
            )
        )
        group.add(label_axis_path)
        group.add(text_path)
        return group

    def add_label_on_the_rim(self, group, label: CircularLabel, radius, record_length):
        (x, y), anchor_value, baseline_value, rotation_deg = rim_text_placement(label, record_length)
        label_path = generate_text_path(
            label["label_text"],
            x,
            y,
            interval=0,
            font_size=self.font_size,
            font_weight="normal",
            font=self.font_family,
            dominant_baseline=baseline_value,
            text_anchor=anchor_value,
            font_style=label["font_style"],
        )
        if label.get("placement") == "radial":
            label_path.rotate(rotation_deg, center=(x, y))
        group.add(label_path)
        return group

    def draw(
        self,
        label: CircularLabel,
        group,
        record_length,
        radius,
        track_ratio,
        feature_track_ratio_factor_override: float | None = None,
        track_preset: str | None = None,
    ):
        cfg = self._cfg
        length_threshold = cfg.labels.length_threshold.circular
        length_param = determine_length_parameter(record_length, length_threshold)
        fallback_font_size = cfg.labels.font_size.for_length_param(length_param)
        self.font_size = label.get("font_size", fallback_font_size)
        self.font_family = label.get("font_family", cfg.objects.text.font_family)
        self.track_type = str(track_preset or label.get("track_preset") or cfg.canvas.circular.track_type)
        self.strandedness = self._profile.strandedness
        if label["is_embedded"] is True:
            group = self.embed_label(
                group,
                label,
                radius,
                record_length,
                track_ratio,
                feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            )
        else:
            group = self.add_label_on_the_rim(group, label, radius, record_length)
        # Each branch appends its text last. Bind that label to the same logical
        # rendered feature as its geometry, irrespective of multipart positions.
        text = group.elements[-1]
        bind_label_part(text, label["feature_id"], complete=True)
        return group


__all__ = ["LabelDrawer"]


