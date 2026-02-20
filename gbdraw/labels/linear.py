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
from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..features.coordinates import get_strand  # type: ignore[reportMissingImports]
from ..core.text import calculate_bbox_dimensions  # type: ignore[reportMissingImports]
from ..core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ..layout.linear_coords import normalize_position_to_linear_track  # type: ignore[reportMissingImports]
from ..layout.linear import calculate_feature_position_factors_linear  # type: ignore[reportMissingImports]


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


def _rotated_y_bounds_from_anchor(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: str,
) -> tuple[float, float]:
    """Return min/max y-offset after rotation around the text anchor point."""
    half_height = height_px / 2.0
    if text_anchor == "start":
        x_values = (0.0, width_px)
    elif text_anchor == "end":
        x_values = (-width_px, 0.0)
    else:
        x_values = (-width_px / 2.0, width_px / 2.0)
    y_values = (-half_height, half_height)
    radians = math.radians(rotation_deg)
    sin_theta = math.sin(radians)
    cos_theta = math.cos(radians)
    y_offsets = [(x * sin_theta) + (y * cos_theta) for x in x_values for y in y_values]
    return min(y_offsets), max(y_offsets)


def calculate_label_y_bounds(label: dict) -> tuple[float, float]:
    """Return absolute top/bottom y coordinates for a label after rotation."""
    y_min_offset, y_max_offset = _rotated_y_bounds_from_anchor(
        float(label["width_px"]),
        float(label["height_px"]),
        float(label.get("rotation_deg", 0.0)),
        str(label.get("text_anchor", "middle")),
    )
    middle_y = float(label["middle_y"])
    return middle_y + y_min_offset, middle_y + y_max_offset


def prepare_label_list_linear(
    feature_dict,
    genome_length,
    alignment_width,
    genome_size_normalization_factor,
    cds_height,
    strandedness,
    track_layout,
    config_dict,
    cfg: GbdrawConfig | None = None,
):
    """
    Prepares a list of labels for linear genome visualization with proper track organization.
    """
    embedded_labels = []
    external_labels = []
    track_dict = defaultdict(list)
    feature_track_positions = {}  # Store feature track positions

    # Get configuration values
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    length_threshold = cfg.labels.length_threshold.circular
    length_param = determine_length_parameter(genome_length, length_threshold)
    font_family = cfg.objects.text.font_family
    font_size = cfg.labels.font_size.linear.for_length_param(length_param)
    linear_label_cfg = cfg.labels.linear
    force_above_feature = linear_label_cfg.placement == "above_feature"
    base_rotation_deg = linear_label_cfg.rotation
    interval = cfg.canvas.dpi
    label_filtering = cfg.labels.filtering.as_dict()
    # Find the maximum feature track ID (positive tracks grow upwards, same as labels)
    max_feature_track = 0
    for feature_object in feature_dict.values():
        # Only consider positive (upward) tracks
        if feature_object.feature_track_id > max_feature_track:
            max_feature_track = feature_object.feature_track_id

    # Get the y-coordinate of the top of the highest feature track
    top_factors = calculate_feature_position_factors_linear(
        strand="positive",
        track_id=max_feature_track,
        separate_strands=strandedness,
        track_layout=track_layout,
    )
    # top_factors[0] is the 'top' y-factor. This will be a negative number (i.e., high on the canvas).
    top_feature_y_limit = cds_height * top_factors[0]

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

        label_anchor_x = (normalized_start + normalized_end) / 2.0
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
            else:
                # Keep rotated label bottom edge above the feature top.
                label_y = feature_top_y - label_vertical_gap - y_max_offset
            label_entry.update(
                {
                    "middle_y": label_y,
                    "is_embedded": True,
                    "track_id": "track_0",
                }
            )
            track_dict["track_0"].append(label_entry)
        elif bbox_width_px < longest_segment_length_in_pixels:
            label_entry.update({"middle_y": feature_y, "is_embedded": True, "track_id": "track_0"})
            track_dict["track_0"].append(label_entry)
        else:
            # For external labels: find lowest possible track
            label_entry.update({"middle_y": 0, "is_embedded": False})

            # Find lowest track where label can be placed without overlaps
            best_track = find_lowest_available_track(track_dict, label_entry)
            label_entry["track_id"] = f"track_{best_track}"
            track_dict[f"track_{best_track}"].append(label_entry)

    # Process embedded labels
    if "track_0" in track_dict:
        for label in track_dict["track_0"]:
            embedded_labels.append(label)

    track_height = max_bbox_height * 1.1  # Reduced height for separate strands mode
    for track_id in sorted(track_dict.keys()):
        if track_id == "track_0":
            continue

        track_labels = track_dict[track_id]
        track_num = int(track_id.split("_")[1])

        for label in track_labels:
            # Compact vertical positioning for external labels
            label["middle_y"] = top_feature_y_limit - (track_height * track_num)
            external_labels.append(label)

    return embedded_labels + external_labels


__all__ = [
    "calculate_label_y_bounds",
    "check_label_overlap",
    "find_lowest_available_track",
    "prepare_label_list_linear",
]


