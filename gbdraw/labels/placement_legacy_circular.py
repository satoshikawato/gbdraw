#!/usr/bin/env python
# coding: utf-8

"""
Legacy circular label placement utilities.

This module contains the circular/arc-based label placement logic that used to live in
`gbdraw.labels.placement_legacy`.

Note: These functions are preserved for backward compatibility with older imports.
New code should prefer `gbdraw.labels.placement_circular`.
"""

import math

from .filtering import get_label_text, preprocess_label_filtering  # type: ignore[reportMissingImports]
from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ..features.coordinates import get_strand  # type: ignore[reportMissingImports]
from ..geometry import calculate_angle_degrees  # type: ignore[reportMissingImports]
from ..layout.circular import (  # type: ignore[reportMissingImports]
    calculate_feature_position_factors_circular,
)
from ..layout.common import calculate_cds_ratio  # type: ignore[reportMissingImports]
from ..text import calculate_bbox_dimensions  # type: ignore[reportMissingImports]


def sort_labels(labels):
    return sorted(labels, key=lambda x: x["middle"])


def y_overlap(label1, label2, total_len, minimum_margin):
    # Adjusted to consider absolute values for y coordinates
    label1_start_y = label1["start_y"]
    label2_start_y = label2["start_y"]
    label1_angle = (360.0 * (label1["middle"] / total_len)) % 360
    label2_angle = (360.0 * (label2["middle"] / total_len)) % 360

    # Helper to calculate max/min Y based on angle
    def get_y_bounds(lbl, angle, start_y):
        h = lbl["height_px"]
        if lbl.get("is_inner", False) is False:
            if 0 <= angle < 10:
                return start_y, start_y - h - 0.5 * minimum_margin
            elif 10 <= angle < 170:
                return (
                    start_y + 0.5 * h + 0.5 * minimum_margin,
                    start_y - 0.5 * h - 0.5 * minimum_margin,
                )
            elif 170 <= angle < 190:
                return start_y + h + 0.5 * minimum_margin, start_y - 0.5 * minimum_margin
            elif 190 <= angle < 350:
                return (
                    start_y + 0.5 * h + 0.5 * minimum_margin,
                    start_y - 0.5 * h - 0.5 * minimum_margin,
                )
            else:
                return start_y, start_y - h - 0.5 * minimum_margin
        else:
            if 0 <= angle < 10:
                return start_y + h + 0.5 * minimum_margin, start_y
            elif 10 <= angle < 170:
                return (
                    start_y + 0.5 * h + 0.5 * minimum_margin,
                    start_y - 0.5 * h - 0.5 * minimum_margin,
                )
            elif 170 <= angle < 190:
                return (
                    start_y + 0.5 * h + 0.5 * minimum_margin,
                    start_y - 0.5 * h - 0.5 * minimum_margin,
                )
            elif 190 <= angle < 350:
                return (
                    start_y + 0.5 * h + 0.5 * minimum_margin,
                    start_y - 0.5 * h - 0.5 * minimum_margin,
                )
            else:
                return start_y, start_y - h - 0.5 * minimum_margin

    max_y1, min_y1 = get_y_bounds(label1, label1_angle, label1_start_y)
    max_y2, min_y2 = get_y_bounds(label2, label2_angle, label2_start_y)

    if min_y1 < min_y2:
        return max_y1 >= min_y2
    else:
        return max_y2 >= min_y1


def x_overlap(label1, label2, minimum_margin=1.0):
    def get_x_bounds(lbl):
        w = lbl["width_px"]
        sx = lbl["start_x"]
        if lbl.get("is_inner", False) is False:
            if sx > 0:
                return sx + w + 0.5 * minimum_margin, sx - 0.5 * minimum_margin
            else:
                return sx + 0.5 * minimum_margin, sx - w - 0.5 * minimum_margin
        else:
            if sx > 0:
                return sx + 0.5 * minimum_margin, sx - w - 0.5 * minimum_margin
            else:
                return sx + w + 0.5 * minimum_margin, sx - 0.5 * minimum_margin

    max_x1, min_x1 = get_x_bounds(label1)
    max_x2, min_x2 = get_x_bounds(label2)

    if min_x1 < min_x2:
        return max_x1 >= min_x2 or max_x1 >= max_x2
    else:
        return max_x2 >= min_x1 or max_x2 >= max_x1


def _calculate_coordinates_ellipse(center_x: float, center_y: float, x_radius: float, y_radius: float, angle_degrees: float):
    """
    Legacy helper: coordinates on an ellipse (x_radius/y_radius) for an angle in degrees.

    `gbdraw.geometry.calculate_coordinates()` is circle-only (single radius) and has a different
    signature; the legacy placement code expects ellipse radii.
    """
    angle_radians = math.radians(angle_degrees)
    x = center_x + x_radius * math.cos(angle_radians)
    y = center_y + y_radius * math.sin(angle_radians)
    return x, y


def place_labels_on_arc_fc(
    labels: list[dict],
    center_x: float,
    center_y: float,
    x_radius: float,
    y_radius: float,
    start_angle: float,
    end_angle: float,
    total_length: int,
) -> list[dict]:
    def check_overlap(label1, label2, total_length, margin):
        return y_overlap(label1, label2, total_length, margin) and x_overlap(
            label1, label2, minimum_margin=2
        )

    if not labels:
        return []

    rearranged_labels = []
    labels = sort_labels(labels)
    current_angle = -75
    increment = 0.1

    for i, label in enumerate(labels):
        if i == 0:
            label["start_x"], label["start_y"] = _calculate_coordinates_ellipse(
                center_x, center_y, x_radius, y_radius, current_angle
            )
            rearranged_labels.append(label)
        else:
            new_angle = current_angle + increment
            if new_angle < -75:
                new_angle = -75
            elif -75 <= new_angle < 85:
                if label["middle"] > (total_length / 2) or i >= len(labels) * (2 / 3):
                    new_angle = 85
                else:
                    new_agle = new_angle
            elif 85 < new_angle < 90:
                new_angle = 90
            else:
                new_angle = new_angle
            if new_angle > 269:
                new_angle = 269
            label["start_x"], label["start_y"] = _calculate_coordinates_ellipse(
                center_x, center_y, x_radius, y_radius, new_angle
            )

            while rearranged_labels and check_overlap(
                label, rearranged_labels[-1], total_length, 0.1
            ):
                new_angle += 0.01
                label["start_x"], label["start_y"] = _calculate_coordinates_ellipse(
                    center_x, center_y, x_radius, y_radius, new_angle
                )

            rearranged_labels.append(label)
            current_angle = new_angle

    return rearranged_labels


def improved_label_placement_fc(
    labels,
    center_x,
    center_y,
    x_radius,
    y_radius,
    feature_radius,
    total_length,
    start_angle,
    end_angle,
    y_margin=0.1,
    max_iterations=10000,
):
    def move_label(label, angle):
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def calculate_angle_of_three_points(x1, y1, x2, y2, x3, y3):
        v1 = (x1 - x2, y1 - y2)
        v2 = (x3 - x2, y3 - y2)
        dot_product = v1[0] * v2[0] + v1[1] * v2[1]
        mag1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2)
        mag2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2)
        try:
            cos_angle = dot_product / (mag1 * mag2)
            angle = math.acos(max(-1, min(1, cos_angle)))
            return math.degrees(angle)
        except ZeroDivisionError:
            return 0

    def check_overlap(label1, label2, total_length):
        return y_overlap(label1, label2, total_length, y_margin) and x_overlap(
            label1, label2, minimum_margin=1
        )

    labels = sort_labels(labels)
    for iteration in range(max_iterations):
        changes_made = False
        for i, label in enumerate(reversed(labels)):
            reverse_i = len(labels) - 1 - i
            normalize = False
            current_angle = calculate_angle_degrees(
                center_x,
                center_y,
                label["start_x"],
                label["start_y"],
                label["middle"],
                start_angle,
                end_angle,
                total_length,
                x_radius,
                y_radius,
                normalize=normalize,
            )

            current_score = calculate_angle_of_three_points(
                label["feature_middle_x"],
                label["feature_middle_y"],
                0,
                0,
                label["start_x"],
                label["start_y"],
            )

            overlaps_prev = False
            overlaps_next = False

            if i == 0:
                overlaps_prev = check_overlap(label, labels[reverse_i - 1], total_length)
                overlaps_next = check_overlap(label, labels[0], total_length)
            elif 0 < i < len(labels) - 1:
                overlaps_prev = check_overlap(labels[reverse_i - 1], label, total_length)
                overlaps_next = check_overlap(label, labels[reverse_i + 1], total_length)
            elif i == len(labels) - 1:
                overlaps_prev = check_overlap(label, labels[-1], total_length)
                overlaps_next = check_overlap(label, labels[reverse_i + 1], total_length)

            if overlaps_prev and overlaps_next:
                continue

            if overlaps_prev:
                direction = 1
            elif overlaps_next:
                direction = -1
            else:
                test_angle_plus = current_angle + 1
                test_x_plus, test_y_plus = move_label(label, test_angle_plus)
                score_plus = calculate_angle_of_three_points(
                    label["feature_middle_x"], label["feature_middle_y"], 0, 0, test_x_plus, test_y_plus
                )

                test_angle_minus = current_angle - 1
                test_x_minus, test_y_minus = move_label(label, test_angle_minus)
                score_minus = calculate_angle_of_three_points(
                    label["feature_middle_x"], label["feature_middle_y"], 0, 0, test_x_minus, test_y_minus
                )
                direction = 1 if abs(score_plus) < abs(score_minus) else -1

            creates_new_overlap = False
            while True:
                new_angle = current_angle + direction * 0.05
                new_x, new_y = move_label(label, new_angle)
                new_score = calculate_angle_of_three_points(
                    label["feature_middle_x"], label["feature_middle_y"], 0, 0, new_x, new_y
                )

                label_copy = label.copy()
                label_copy["start_x"], label_copy["start_y"] = new_x, new_y
                if i == 0:
                    if overlaps_prev:
                        creates_new_overlap = check_overlap(label_copy, labels[0], total_length)
                    elif overlaps_next:
                        creates_new_overlap = check_overlap(label_copy, labels[reverse_i - 1], total_length)
                    else:
                        creates_new_overlap = check_overlap(label_copy, labels[reverse_i - 1], total_length) or check_overlap(
                            label_copy, labels[0], total_length
                        )
                elif 0 < i < len(labels) - 1:
                    prev_label = labels[reverse_i - 1]
                    next_label = labels[reverse_i + 1]
                    if overlaps_prev:
                        creates_new_overlap = check_overlap(label_copy, next_label, total_length)
                    elif overlaps_next:
                        creates_new_overlap = check_overlap(label_copy, prev_label, total_length)
                    else:
                        creates_new_overlap = check_overlap(label_copy, prev_label, total_length) or check_overlap(
                            label_copy, next_label, total_length
                        )
                elif i == len(labels) - 1:
                    if overlaps_prev:
                        creates_new_overlap = check_overlap(label_copy, labels[reverse_i + 1], total_length)
                    elif overlaps_next:
                        creates_new_overlap = check_overlap(label_copy, labels[-1], total_length)
                    else:
                        creates_new_overlap = check_overlap(label_copy, labels[-1], total_length) or check_overlap(
                            label_copy, labels[reverse_i + 1], total_length
                        )

                if (abs(new_score) <= abs(current_score)) and not creates_new_overlap:
                    label["start_x"], label["start_y"] = new_x, new_y
                    current_angle = new_angle
                    current_score = new_score
                    changes_made = True
                else:
                    break

            if not changes_made:
                break

    return labels


def rearrange_labels_fc(
    labels,
    feature_radius,
    total_length,
    genome_len,
    config_dict,
    strands,
    is_outer,
    cfg: GbdrawConfig | None = None,
):
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    track_type = cfg.canvas.circular.track_type

    if is_outer:
        offset_config = cfg.labels.unified_adjustment.outer_labels
        x_radius_factor = cfg.labels.arc_x_radius_factor[track_type][strands][genome_len]
        y_radius_factor = cfg.labels.arc_y_radius_factor[track_type][strands][genome_len]
        default_center_x = cfg.labels.arc_center_x[track_type][genome_len]
    else:
        offset_config = cfg.labels.unified_adjustment.inner_labels
        x_radius_factor = cfg.labels.inner_arc_x_radius_factor[track_type][strands][genome_len]
        y_radius_factor = cfg.labels.inner_arc_y_radius_factor[track_type][strands][genome_len]
        default_center_x = cfg.labels.inner_arc_center_x[track_type][genome_len]

    x_radius_offset = offset_config.x_radius_offset
    y_radius_offset = offset_config.y_radius_offset

    if is_outer:
        x_radius = feature_radius * x_radius_factor * x_radius_offset
        y_radius = feature_radius * y_radius_factor * y_radius_offset
    else:
        x_radius = feature_radius * x_radius_factor * (2 - x_radius_offset)
        y_radius = feature_radius * y_radius_factor * (2 - y_radius_offset)

    center_y = 0
    center_x = default_center_x
    start_angle = 0
    end_angle = 360

    labels = sorted(labels, key=lambda x: x["middle"])
    labels = place_labels_on_arc_fc(labels, center_x, center_y, x_radius, y_radius, start_angle, end_angle, total_length)
    labels = improved_label_placement_fc(
        labels, center_x, center_y, x_radius, y_radius, feature_radius, total_length, start_angle, end_angle
    )

    return labels


def prepare_label_list(feature_dict, total_length, radius, track_ratio, config_dict, cfg: GbdrawConfig | None = None):
    embedded_labels = []
    outer_labels = []
    inner_labels = []

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    label_filtering = cfg.labels.filtering.as_dict()
    label_filtering = preprocess_label_filtering(label_filtering)
    length_threshold = cfg.labels.length_threshold.circular
    length_param = determine_length_parameter(total_length, length_threshold)
    track_type = cfg.canvas.circular.track_type
    strandedness = cfg.canvas.strandedness

    strands = "separate" if strandedness else "single"
    allow_inner_labels = cfg.canvas.circular.allow_inner_labels
    radius_factor = cfg.labels.radius_factor[track_type][strands][length_param]
    inner_radius_factor = cfg.labels.inner_radius_factor[track_type][strands][length_param]

    font_family = cfg.objects.text.font_family
    font_size: float = cfg.labels.font_size.for_length_param(length_param)
    interval = cfg.canvas.dpi

    track_ratio_factor = cfg.canvas.circular.track_ratio_factors[length_param][0]

    for feature_object in feature_dict.values():
        feature_label_text = get_label_text(feature_object, label_filtering)
        if feature_label_text == "":
            continue

        label_entry = dict()
        longest_segment_length = 0
        is_embedded = False
        coordinate_strand: str = "undefined"
        feature_location_list = feature_object.location
        list_of_coordinates = feature_object.coordinates
        # `FeatureObject.strand` is derived at creation time; keep placement pure.

        feature_location_count = 0
        longest_segment_start = 0
        longest_segment_end = 0
        longeset_segment_middle = 0

        for coordinate in list_of_coordinates:
            if feature_location_list[feature_location_count].kind == "line":
                feature_location_count += 1
                continue
            else:
                coordinate_start = int(coordinate.start)
                coordinate_end = int(coordinate.end)
                coordinate_strand = get_strand(coordinate.strand)
                interval_length = abs(int(coordinate_end - coordinate_start) + 1)
                interval_middle = int(coordinate_end + coordinate_start) / 2
                feature_location_count += 1
                if interval_length > longest_segment_length:
                    longest_segment_start = coordinate_start
                    longest_segment_end = coordinate_end
                    longeset_segment_middle = interval_middle
                    longest_segment_length = interval_length

        cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
        factors: list[float] = calculate_feature_position_factors_circular(
            total_length, coordinate_strand, track_ratio, cds_ratio, offset, track_type, strandedness
        )

        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(
            feature_label_text, font_family, font_size, interval
        )

        label_middle = longeset_segment_middle
        label_as_feature_length = total_length * (1.1 * bbox_width_px) / (2 * math.pi * radius)
        label_start = label_middle - (label_as_feature_length / 2)
        label_end = label_middle + (label_as_feature_length / 2)
        feature_middle_x: float = (radius * factors[1]) * math.cos(
            math.radians(360.0 * (label_middle / total_length) - 90)
        )
        feature_middle_y: float = (radius * factors[1]) * math.sin(
            math.radians(360.0 * (label_middle / total_length) - 90)
        )

        if feature_object.strand == "positive" or allow_inner_labels is False:
            middle_x: float = (radius_factor * radius) * math.cos(
                math.radians(360.0 * (label_middle / total_length) - 90)
            )
            middle_y: float = (radius_factor * radius) * math.sin(
                math.radians(360.0 * (label_middle / total_length) - 90)
            )
        else:
            middle_x = (inner_radius_factor * radius) * math.cos(
                math.radians(360.0 * (label_middle / total_length) - 90)
            )
            middle_y = (inner_radius_factor * radius) * math.sin(
                math.radians(360.0 * (label_middle / total_length) - 90)
            )

        is_embedded = label_start > longest_segment_start and label_end < longest_segment_end

        label_entry["label_text"] = feature_label_text
        label_entry["middle"] = label_middle
        label_entry["start"] = label_start
        label_entry["end"] = label_end
        label_entry["middle_x"] = middle_x
        label_entry["middle_y"] = middle_y
        label_entry["feature_middle_x"] = feature_middle_x
        label_entry["feature_middle_y"] = feature_middle_y
        label_entry["width_px"] = bbox_width_px
        label_entry["height_px"] = bbox_height_px
        label_entry["strand"] = coordinate_strand
        label_entry["is_embedded"] = is_embedded
        label_entry["font_size"] = font_size
        label_entry["font_family"] = font_family

        if is_embedded:
            embedded_labels.append(label_entry)
        else:
            if feature_object.strand == "positive" or allow_inner_labels is False:
                label_entry["is_inner"] = False
                outer_labels.append(label_entry)
            else:
                label_entry["is_inner"] = True
                inner_labels.append(label_entry)

    outer_labels_rearranged = rearrange_labels_fc(
        outer_labels,
        radius,
        total_length,
        length_param,
        config_dict,
        strands,
        is_outer=True,
        cfg=cfg,
    )
    inner_labels_rearranged = []
    if allow_inner_labels:
        inner_labels_rearranged = rearrange_labels_fc(
            inner_labels,
            radius,
            total_length,
            length_param,
            config_dict,
            strands,
            is_outer=False,
            cfg=cfg,
        )

    label_list_fc = embedded_labels + outer_labels_rearranged + inner_labels_rearranged
    return label_list_fc


__all__ = [
    "improved_label_placement_fc",
    "place_labels_on_arc_fc",
    "prepare_label_list",
    "rearrange_labels_fc",
    "sort_labels",
    "x_overlap",
    "y_overlap",
]


