#!/usr/bin/env python
# coding: utf-8

"""
Linear label placement utilities.

This module contains the linear track-based label placement logic that used to live in
`gbdraw.labels.placement`.
"""

import math
from collections import defaultdict

from .filtering import get_label_text  # type: ignore[reportMissingImports]
from .policy import normalize_label_rendering
from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..features.coordinates import get_strand  # type: ignore[reportMissingImports]
from ..core.text import calculate_bbox_dimensions  # type: ignore[reportMissingImports]
from ..core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ..layout.linear_coords import normalize_position_to_linear_track  # type: ignore[reportMissingImports]
from ..layout.linear import calculate_feature_position_factors_linear  # type: ignore[reportMissingImports]
from ..layout.spatial import Aabb, AabbIndex, Interval, IntervalIndex


def check_label_overlap(label1, label2):
    """Check if two labels overlap horizontally"""
    return not (label1["end"] < label2["start"] or label2["end"] < label1["start"])


def find_lowest_available_track(track_dict, label):
    """Find the lowest track (closest to track_1) where the label can be placed without overlap"""
    track_num = 1
    while True:
        track_id = f"track_{track_num}"
        # Check if this track has overlaps
        has_overlap = False
        if track_id in track_dict:
            for existing_label in track_dict[track_id]:
                if check_label_overlap(label, existing_label):
                    has_overlap = True
                    break

        if not has_overlap:
            return track_num
        track_num += 1


def _linear_label_interval(label: dict) -> Interval:
    return Interval(float(label["start"]), float(label["end"]))


def _linear_label_aabb(label: dict) -> Aabb:
    left, right, top, bottom = calculate_label_bounds(label)
    return Aabb(left, top, right, bottom)


def _external_label_bucket_size(alignment_width: float) -> float:
    return max(16.0, float(alignment_width) / 256.0)


def _above_feature_label_bucket_size(labels: list[dict], min_gap_px: float) -> float:
    max_extent = max(
        (
            max(float(label.get("width_px", 0.0)), float(label.get("height_px", 0.0)))
            for label in labels
        ),
        default=16.0,
    )
    return max(16.0, max_extent + float(min_gap_px))


def _find_lowest_available_track_indexed(
    track_dict: dict,
    track_indexes: dict[str, IntervalIndex],
    label_by_id: dict[int, dict],
    label: dict,
    bucket_size: float,
) -> int:
    """Find the lowest non-overlapping external label track using index candidates."""
    label_interval = _linear_label_interval(label)
    track_num = 1
    while True:
        track_id = f"track_{track_num}"
        if track_id not in track_dict or not track_dict[track_id]:
            return track_num

        index = track_indexes.get(track_id)
        if index is None:
            candidates = track_dict[track_id]
        else:
            candidates = [label_by_id[label_id] for label_id in index.query(label_interval)]

        if not any(check_label_overlap(label, existing_label) for existing_label in candidates):
            return track_num
        track_num += 1


def _insert_external_label_index(
    track_indexes: dict[str, IntervalIndex],
    track_id: str,
    label_id: int,
    label: dict,
    bucket_size: float,
) -> None:
    index = track_indexes.setdefault(track_id, IntervalIndex(bucket_size=bucket_size))
    index.insert(label_id, _linear_label_interval(label))


def _anchor_x_values(width_px: float, text_anchor: str) -> tuple[float, float]:
    if text_anchor == "start":
        return 0.0, width_px
    if text_anchor == "end":
        return -width_px, 0.0
    return -width_px / 2.0, width_px / 2.0


def _rotated_corner_offsets(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
) -> list[tuple[float, float]]:
    half_height = height_px / 2.0
    x_values = _anchor_x_values(width_px, text_anchor)
    y_values = (-half_height, half_height)
    radians = math.radians(rotation_deg)
    sin_theta = math.sin(radians)
    cos_theta = math.cos(radians)
    return [
        (
            (x * cos_theta) - (y * sin_theta),
            (x * sin_theta) + (y * cos_theta),
        )
        for x in x_values
        for y in y_values
    ]


def _rotated_y_bounds_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
) -> tuple[float, float]:
    """Return min/max y-offset after rotation around the text anchor point."""
    _, _, y_min, y_max = _rotated_bounds_from_anchor(width_px, height_px, rotation_deg, text_anchor)
    return y_min, y_max


def _rotated_bounds_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
) -> tuple[float, float, float, float]:
    """Return min/max x/y offsets after rotation around the text anchor point."""
    offsets = _rotated_corner_offsets(width_px, height_px, rotation_deg, text_anchor)
    x_offsets = [point[0] for point in offsets]
    y_offsets = [point[1] for point in offsets]
    return min(x_offsets), max(x_offsets), min(y_offsets), max(y_offsets)


def _rotated_extreme_y_point_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
    *,
    choose_max: bool,
) -> tuple[float, float]:
    """Return the rotated corner offset at the top/bottom edge used as a leader contact."""
    points = _rotated_corner_offsets(width_px, height_px, rotation_deg, text_anchor)
    return max(points, key=lambda point: point[1]) if choose_max else min(points, key=lambda point: point[1])


def calculate_label_bounds(label: dict) -> tuple[float, float, float, float]:
    """Return absolute left/right/top/bottom coordinates for a label after rotation."""
    x_min_offset, x_max_offset, y_min_offset, y_max_offset = _rotated_bounds_from_anchor(
        float(label["width_px"]),
        float(label["height_px"]),
        float(label.get("rotation_deg", 0.0)),
        str(label.get("text_anchor", "middle")),
    )
    middle_x = float(label["middle_x"])
    middle_y = float(label["middle_y"])
    return (
        middle_x + x_min_offset,
        middle_x + x_max_offset,
        middle_y + y_min_offset,
        middle_y + y_max_offset,
    )


def calculate_label_y_bounds(label: dict) -> tuple[float, float]:
    """Return absolute top/bottom y coordinates for a label after rotation."""
    _, _, top_y, bottom_y = calculate_label_bounds(label)
    return top_y, bottom_y


def _label_bounds_overlap(label1: dict, label2: dict, min_gap_px: float) -> bool:
    left1, right1, top1, bottom1 = calculate_label_bounds(label1)
    left2, right2, top2, bottom2 = calculate_label_bounds(label2)
    return not (
        right1 + min_gap_px <= left2
        or right2 + min_gap_px <= left1
        or bottom1 + min_gap_px <= top2
        or bottom2 + min_gap_px <= top1
    )


def _update_linear_label_leader_end(label: dict) -> None:
    contact_x_offset = float(label.get("label_contact_x_offset", 0.0))
    contact_y_offset = float(label.get("label_contact_y_offset", 0.0))
    label["leader_end_x"] = float(label["middle_x"]) + contact_x_offset
    label["leader_end_y"] = float(label["middle_y"]) + contact_y_offset


def _place_linear_label_without_overlap(
    label: dict,
    placed: list[dict],
    min_gap_px: float,
    direction: float,
    placed_index: AabbIndex | None = None,
    placed_by_id: dict[int, dict] | None = None,
) -> bool:
    """Move one linear label away from the axis until its bbox clears placed labels."""
    original_y = float(label["middle_y"])
    for _ in range(len(placed) + 1):
        required_shift = 0.0
        _, _, candidate_top, candidate_bottom = calculate_label_bounds(label)

        if placed_index is None or placed_by_id is None:
            candidates = placed
        else:
            candidates = [
                placed_by_id[label_id]
                for label_id in placed_index.query(_linear_label_aabb(label), padding=min_gap_px)
            ]

        for other in candidates:
            if not _label_bounds_overlap(label, other, min_gap_px):
                continue

            _, _, placed_top, placed_bottom = calculate_label_bounds(other)
            if direction < 0.0:
                shift = (placed_top - float(min_gap_px)) - candidate_bottom
                required_shift = min(required_shift, shift)
            else:
                shift = (placed_bottom + float(min_gap_px)) - candidate_top
                required_shift = max(required_shift, shift)

        if required_shift == 0.0:
            return abs(float(label["middle_y"]) - original_y) > 0.001

        label["middle_y"] = float(label["middle_y"]) + required_shift

    return abs(float(label["middle_y"]) - original_y) > 0.001


def _resolve_above_feature_label_overlaps(labels: list[dict], min_gap_px: float) -> None:
    """Stack rotated above-feature labels away from the axis when their bboxes collide."""
    if not labels:
        return
    if len(labels) > 300:
        return

    for place_above in (True, False):
        side_labels = [label for label in labels if bool(label.get("above_feature_place_above", True)) is place_above]
        if not side_labels:
            continue
        placed: list[dict] = []
        placed_by_id: dict[int, dict] = {}
        placed_index = AabbIndex(bucket_size=_above_feature_label_bucket_size(side_labels, min_gap_px))
        direction = -1.0 if place_above else 1.0
        for label in sorted(side_labels, key=lambda item: float(item.get("feature_anchor_x", item["middle_x"]))):
            if _place_linear_label_without_overlap(
                label,
                placed,
                min_gap_px,
                direction,
                placed_index=placed_index,
                placed_by_id=placed_by_id,
            ):
                label["leader_line"] = True
                _update_linear_label_leader_end(label)
            placed_id = len(placed_by_id)
            placed_by_id[placed_id] = label
            placed_index.insert(placed_id, _linear_label_aabb(label))
            placed.append(label)


def prepare_label_list_linear(
    feature_dict,
    genome_length,
    alignment_width,
    genome_size_normalization_factor,
    cds_height,
    strandedness,
    track_layout,
    track_axis_gap,
    config_dict,
    cfg: GbdrawConfig | None = None,
):
    """
    Prepares a list of labels for linear genome visualization with proper track organization.
    """
    embedded_labels = []
    external_labels = []
    track_dict = defaultdict(list)
    external_track_indexes: dict[str, IntervalIndex] = {}
    external_label_by_id: dict[int, dict] = {}
    feature_track_positions = {}  # Store feature track positions
    track_layout_normalized = str(track_layout).strip().lower()
    axis_gap_factor = (
        (float(track_axis_gap) / float(cds_height))
        if (track_axis_gap is not None and float(cds_height) > 0.0)
        else None
    )

    # Get configuration values
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    length_threshold = cfg.labels.length_threshold.circular
    length_param = determine_length_parameter(genome_length, length_threshold)
    font_family = cfg.objects.text.font_family
    font_size = cfg.labels.font_size.linear.for_length_param(length_param)
    linear_label_cfg = cfg.labels.linear
    label_rendering = normalize_label_rendering(cfg.labels.rendering)
    label_spacing_px = float(cfg.labels.spacing.linear)
    external_bucket_size = _external_label_bucket_size(alignment_width)
    force_above_feature = linear_label_cfg.placement == "above_feature"
    if force_above_feature and label_rendering != "auto":
        raise ValueError(
            "label_rendering embedded_only|external_only cannot be used with label_placement above_feature"
        )
    base_rotation_deg = linear_label_cfg.rotation
    interval = cfg.canvas.dpi
    label_filtering = cfg.labels.filtering.as_dict()
    # First pass: Calculate feature track positions
    for feature_id, feature_object in feature_dict.items():
        if len(feature_object.coordinates) == 0:
            continue

        # Get feature track info
        coordinate = feature_object.coordinates[0]  # Use first coordinate for strand
        strand = get_strand(coordinate.strand)
        feature_track_id = feature_object.feature_track_id

        # Calculate track position using the same logic as for features
        factors = calculate_feature_position_factors_linear(
            strand,
            feature_track_id,
            strandedness,
            track_layout=track_layout,
            axis_gap_factor=axis_gap_factor,
        )

        track_y_position = cds_height * factors[1]  # Use middle factor
        track_top_y = cds_height * factors[0]
        track_bottom_y = cds_height * factors[2]

        # Store track position for this feature
        feature_track_positions[feature_id] = {
            "middle_y": track_y_position,
            "top_y": track_top_y,
            "bottom_y": track_bottom_y,
        }

    if feature_track_positions:
        top_feature_y_limit = min(pos["top_y"] for pos in feature_track_positions.values())
        bottom_feature_y_limit = max(pos["bottom_y"] for pos in feature_track_positions.values())
    else:
        top_feature_y_limit = 0.0
        bottom_feature_y_limit = 0.0

    # Second pass: Process labels
    max_bbox_height = 0
    for feature_id, feature_object in reversed(list(feature_dict.items())):
        feature_label_text = get_label_text(feature_object, label_filtering)
        feature_track_id = feature_object.feature_track_id
        if not feature_label_text:
            continue

        # Calculate label dimensions and positions
        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
        if bbox_height_px > max_bbox_height:
            max_bbox_height = bbox_height_px
        # Find the longest segment and its middle point
        longest_segment_length = 0
        coordinate_strand = None
        factors = None
        longest_segment_start = 0
        longest_segment_end = 0

        for coordinate in feature_object.coordinates:
            coordinate_strand = get_strand(coordinate.strand)
            factors = calculate_feature_position_factors_linear(
                coordinate_strand,
                feature_track_id,
                strandedness,
                track_layout=track_layout,
                axis_gap_factor=axis_gap_factor,
            )
            start = int(coordinate.start)
            end = int(coordinate.end)
            segment_length = abs(end - start + 1)
            if segment_length > longest_segment_length:
                longest_segment_start = start
                longest_segment_end = end
                longest_segment_length = segment_length

        # Calculate normalized positions
        normalized_start = normalize_position_to_linear_track(
            longest_segment_start, genome_length, alignment_width, genome_size_normalization_factor
        )
        normalized_end = normalize_position_to_linear_track(
            longest_segment_end, genome_length, alignment_width, genome_size_normalization_factor
        )
        longest_segment_length_in_pixels = abs(normalized_end - normalized_start) + 1

        if force_above_feature:
            # Above-feature mode tilts labels in the opposite direction of the user-provided angle.
            # In separate-strands mode, negative strand labels are mirrored and placed below features.
            is_negative_separate = strandedness and coordinate_strand == "negative"
            label_rotation_deg = base_rotation_deg if is_negative_separate else -base_rotation_deg
            label_text_anchor = "start" if label_rotation_deg != 0.0 else "middle"
        else:
            # Keep auto mode horizontal regardless of configured rotation.
            label_rotation_deg = 0.0
            label_text_anchor = "middle"

        feature_anchor_x = (normalized_start + normalized_end) / 2.0
        label_anchor_x = feature_anchor_x
        label_contact_x_offset = 0.0
        label_contact_y_offset = 0.0
        label_leader_start_y = 0.0
        if force_above_feature:
            place_above_feature = not (strandedness and coordinate_strand == "negative")
            if label_rotation_deg == 0.0:
                y_min_offset, y_max_offset = _rotated_y_bounds_from_anchor(
                    bbox_width_px, bbox_height_px, label_rotation_deg, label_text_anchor
                )
                label_contact_y_offset = y_max_offset if place_above_feature else y_min_offset
            else:
                label_contact_x_offset, label_contact_y_offset = _rotated_extreme_y_point_from_anchor(
                    bbox_width_px,
                    bbox_height_px,
                    label_rotation_deg,
                    label_text_anchor,
                    choose_max=place_above_feature,
                )
            label_anchor_x = feature_anchor_x - label_contact_x_offset
        bbox_start = label_anchor_x - (bbox_width_px / 2)
        bbox_end = label_anchor_x + (bbox_width_px / 2)

        # Get actual feature track position
        feature_position = feature_track_positions.get(
            feature_id,
            {
                "middle_y": cds_height * factors[1],
                "top_y": cds_height * factors[0],
                "bottom_y": cds_height * factors[2],
            },
        )
        feature_y = feature_position["middle_y"]
        feature_top_y = feature_position["top_y"]
        feature_bottom_y = feature_position["bottom_y"]

        # Create base label entry
        label_entry = {
            "label_text": feature_label_text,
            "middle": label_anchor_x,
            "start": bbox_start,
            "end": bbox_end,
            "middle_x": label_anchor_x,
            "width_px": bbox_width_px,
            "height_px": bbox_height_px,
            "strand": coordinate_strand,
            "feature_middle_y": feature_y,  # Use actual feature position
            "feature_top_y": feature_top_y,
            "feature_bottom_y": feature_bottom_y,
            "feature_start_x": normalized_start,
            "feature_end_x": normalized_end,
            "feature_anchor_x": feature_anchor_x,
            "font_size": font_size,
            "font_family": font_family,
            "rotation_deg": label_rotation_deg,
            "text_anchor": label_text_anchor,
        }

        # Determine if label should be embedded
        if force_above_feature:
            y_min_offset, y_max_offset = _rotated_y_bounds_from_anchor(
                bbox_width_px, bbox_height_px, label_rotation_deg, label_text_anchor
            )
            label_vertical_gap = max(1.0, bbox_height_px * 0.05)
            is_negative_separate = strandedness and coordinate_strand == "negative"
            if is_negative_separate:
                # Keep rotated label top edge below the feature bottom.
                label_y = feature_bottom_y + label_vertical_gap - y_min_offset
                label_leader_start_y = feature_bottom_y
            else:
                # Keep rotated label bottom edge above the feature top.
                label_y = feature_top_y - label_vertical_gap - y_max_offset
                label_leader_start_y = feature_top_y
            label_entry.update(
                {
                    "middle_y": label_y,
                    "is_embedded": True,
                    "track_id": "track_0",
                    "above_feature_place_above": not is_negative_separate,
                    "label_contact_x_offset": label_contact_x_offset,
                    "label_contact_y_offset": label_contact_y_offset,
                    "leader_line": False,
                    "leader_start_x": feature_anchor_x,
                    "leader_start_y": label_leader_start_y,
                    "leader_end_x": feature_anchor_x,
                    "leader_end_y": label_y + label_contact_y_offset,
                }
            )
            track_dict["track_0"].append(label_entry)
        elif label_rendering != "external_only" and bbox_width_px < longest_segment_length_in_pixels:
            label_entry.update({"middle_y": feature_y, "is_embedded": True, "track_id": "track_0"})
            track_dict["track_0"].append(label_entry)
        elif label_rendering != "embedded_only":
            # For external labels: find lowest possible track
            label_entry.update({"middle_y": 0, "is_embedded": False})

            # Find lowest track where label can be placed without overlaps
            best_track = _find_lowest_available_track_indexed(
                track_dict,
                external_track_indexes,
                external_label_by_id,
                label_entry,
                external_bucket_size,
            )
            label_entry["track_id"] = f"track_{best_track}"
            track_dict[f"track_{best_track}"].append(label_entry)
            external_label_id = len(external_label_by_id)
            external_label_by_id[external_label_id] = label_entry
            _insert_external_label_index(
                external_track_indexes,
                f"track_{best_track}",
                external_label_id,
                label_entry,
                external_bucket_size,
            )

    # Process embedded labels
    if "track_0" in track_dict:
        if force_above_feature:
            above_feature_collision_gap = max(1.0, max_bbox_height * 0.05)
            _resolve_above_feature_label_overlaps(track_dict["track_0"], above_feature_collision_gap)
        for label in track_dict["track_0"]:
            embedded_labels.append(label)

    track_height = max_bbox_height + label_spacing_px
    for track_id in sorted(track_dict.keys()):
        if track_id == "track_0":
            continue

        track_labels = track_dict[track_id]
        track_num = int(track_id.split("_")[1])

        for label in track_labels:
            # Compact vertical positioning for external labels
            if track_layout_normalized in {"below", "tuckin"}:
                label["middle_y"] = bottom_feature_y_limit + (track_height * track_num)
            else:
                label["middle_y"] = top_feature_y_limit - (track_height * track_num)
            external_labels.append(label)

    return embedded_labels + external_labels


__all__ = [
    "calculate_label_bounds",
    "calculate_label_y_bounds",
    "check_label_overlap",
    "find_lowest_available_track",
    "prepare_label_list_linear",
]


