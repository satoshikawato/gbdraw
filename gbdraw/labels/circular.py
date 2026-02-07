#!/usr/bin/env python
# coding: utf-8

"""
Circular label placement utilities.

This module contains the circular/arc-based label placement logic that used to live in
`gbdraw.labels.placement`.
"""

import math

from .filtering import get_label_text, preprocess_label_filtering
from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..features.coordinates import get_strand
from ..core.text import calculate_bbox_dimensions
from ..core.sequence import determine_length_parameter
from ..layout.common import calculate_cds_ratio
from ..layout.circular import calculate_feature_position_factors_circular


def y_overlap(label1, label2, total_len, minimum_margin):
    # Adjusted to consider absolute values for y coordinates
    label1_start_y = label1["start_y"]
    label2_start_y = label2["start_y"]
    label1_angle = (360.0 * (label1["middle"] / total_len)) % 360
    label2_angle = (360.0 * (label2["middle"] / total_len)) % 360
    if label1["is_inner"] is False:
        if 0 <= label1_angle < 10:  # baseline_value = "text-after-edge"
            max_y1 = label1_start_y
            min_y1 = label1_start_y - 1.0 * label1["height_px"] - 0.5 * minimum_margin
        elif 10 <= label1_angle < 170:  # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        elif 170 <= label1_angle < 190:  # baseline_value = "hanging"
            max_y1 = label1_start_y + 1.0 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * minimum_margin
        elif 190 <= label1_angle < 350:  # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        else:  # baseline_value = "text-after-edge"
            max_y1 = label1_start_y
            min_y1 = label1_start_y - 1.0 * label1["height_px"] - 0.5 * minimum_margin
    else:
        if 0 <= label1_angle < 10:  # baseline_value = "hanging"
            max_y1 = label1_start_y + 1.0 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y
        elif 10 <= label1_angle < 170:  # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        elif 170 <= label1_angle < 190:  # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        elif 190 <= label1_angle < 350:  # baseline_value = "middle"
            max_y1 = label1_start_y + 0.5 * label1["height_px"] + 0.5 * minimum_margin
            min_y1 = label1_start_y - 0.5 * label1["height_px"] - 0.5 * minimum_margin
        else:  # baseline_value = "hanging"
            max_y1 = label1_start_y
            min_y1 = label1_start_y - 1.0 * label1["height_px"] - 0.5 * minimum_margin
    if label2["is_inner"] is False:
        if 0 <= label2_angle < 10:  # baseline_value = "text-after-edge"
            max_y2 = label2_start_y
            min_y2 = label2_start_y - 1.0 * label2["height_px"] - 0.5 * minimum_margin
        elif 10 <= label2_angle < 170:  # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        elif 170 <= label2_angle < 190:  # baseline_value = "hanging"
            max_y2 = label2_start_y + 1.0 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * minimum_margin
        elif 190 <= label2_angle < 350:  # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        else:  # baseline_value = "text-after-edge"
            max_y2 = label2_start_y
            min_y2 = label2_start_y - 1.0 * label2["height_px"] - 0.5 * minimum_margin
    else:
        if 0 <= label2_angle < 10:  # baseline_value = "hanging"
            max_y2 = label2_start_y + 1.0 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y
        elif 10 <= label2_angle < 170:  # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        elif 170 <= label2_angle < 190:  # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        elif 190 <= label2_angle < 350:  # baseline_value = "middle"
            max_y2 = label2_start_y + 0.5 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y - 0.5 * label2["height_px"] - 0.5 * minimum_margin
        else:  # baseline_value = "hanging"
            max_y2 = label2_start_y + 1.0 * label2["height_px"] + 0.5 * minimum_margin
            min_y2 = label2_start_y
    if min_y1 < min_y2:
        return max_y1 >= min_y2
    else:
        return max_y2 >= min_y1


def x_overlap(label1, label2, minimum_margin=1.0):
    # Adjusted to directly return the evaluated condition

    if label1["is_inner"] is False:
        if label1["start_x"] > 0:
            max_x1 = label1["start_x"] + label1["width_px"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - 0.5 * minimum_margin
        else:
            max_x1 = label1["start_x"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - label1["width_px"] - 0.5 * minimum_margin
    else:
        if label1["start_x"] > 0:
            max_x1 = label1["start_x"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - label1["width_px"] - 0.5 * minimum_margin
        else:
            max_x1 = label1["start_x"] + label1["width_px"] + 0.5 * minimum_margin
            min_x1 = label1["start_x"] - 0.5 * minimum_margin
    if label2["is_inner"] is False:
        if label2["start_x"] > 0:
            max_x2 = label2["start_x"] + label2["width_px"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - 0.5 * minimum_margin
        else:
            max_x2 = label2["start_x"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - label2["width_px"] - 0.5 * minimum_margin
    else:
        if label2["start_x"] > 0:
            max_x2 = label2["start_x"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - label2["width_px"] - 0.5 * minimum_margin
        else:
            max_x2 = label2["start_x"] + label2["width_px"] + 0.5 * minimum_margin
            min_x2 = label2["start_x"] - 0.5 * minimum_margin
    if min_x1 < min_x2:
        if max_x1 >= min_x2:
            return True
        # when nested also true
        elif max_x1 >= max_x2:
            return True
        else:
            return False
    else:
        if max_x2 >= min_x1:
            return True
        # when nested also true
        elif max_x2 >= max_x1:
            return True
        else:
            return False


def calculate_angle_degrees(
    center_x, center_y, x, y, middle, start_angle, end_angle, total_length, x_radius, y_radius, normalize
):
    """Calculate the angle in degrees from a point (x, y) relative to the center of an ellipse."""
    # Normalize coordinates to a unit circle
    x_normalized = (x - center_x) / x_radius
    y_normalized = (y - center_y) / y_radius

    angle_radians = math.atan2(y_normalized, x_normalized)
    angle_degrees = math.degrees(angle_radians)
    return angle_degrees


def calculate_coordinates(center_x, center_y, x_radius, y_radius, angle_degrees, middle, total_length):
    angle_radians = math.radians(angle_degrees)
    y = center_y + y_radius * math.sin(angle_radians)
    x = center_x + x_radius * math.cos(angle_radians)

    return x, y


def calculate_angle_for_y(center_y, y_radius, y):
    """Calculate the angle in degrees for a given y-coordinate on an ellipse."""
    if center_y - y_radius <= y <= center_y + y_radius:
        # Calculate the arcsine of the normalized y-coordinate
        angle_radians = math.asin((y - center_y) / y_radius)
        angle_degrees = math.degrees(angle_radians)
        return angle_degrees
    else:
        return None  # Indicates the y-coordinate is outside the ellipse's bounds


def angle_from_middle(middle: float, total_length: int) -> float:
    """Convert genome midpoint to circular angle (degrees)."""
    return (360.0 * (middle / total_length) - 90.0) % 360.0


def angle_from_middle_unwrapped(middle: float, total_length: int) -> float:
    """Convert genome midpoint to monotonic circular angle (degrees, unwrapped)."""
    return 360.0 * (middle / total_length) - 90.0


def angle_delta_deg(a: float, b: float) -> float:
    """Return the signed smallest angle delta from b to a in [-180, 180)."""
    return ((a - b + 180.0) % 360.0) - 180.0


def normalize_angle_near_reference(angle_deg: float, reference_unwrapped: float) -> float:
    """Project an angle to the nearest 360-degree branch around a reference angle."""
    normalized = angle_deg % 360.0
    return normalized + 360.0 * round((reference_unwrapped - normalized) / 360.0)


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

    rearranged_labels: list[dict] = []
    labels = sort_labels(labels)
    angle_step = 0.25
    max_shift = 25.0
    max_steps = int(max_shift / angle_step)
    min_order_gap = 0.05
    prev_angle_unwrapped = -float("inf")

    for label in labels:
        target_angle_unwrapped = angle_from_middle_unwrapped(label["middle"], total_length)
        label["target_angle_unwrapped"] = target_angle_unwrapped
        label["target_angle"] = target_angle_unwrapped % 360.0
        lower_bound = prev_angle_unwrapped + min_order_gap
        best_x = None
        best_y = None
        best_overlaps = None
        best_angle_unwrapped = None

        for step in range(max_steps + 1):
            if step == 0:
                candidate_offsets = [0.0]
            else:
                offset = step * angle_step
                candidate_offsets = [offset, -offset]

            for offset in candidate_offsets:
                candidate_angle_unwrapped = target_angle_unwrapped + offset
                if candidate_angle_unwrapped < lower_bound:
                    continue
                candidate_x, candidate_y = calculate_coordinates(
                    center_x, center_y, x_radius, y_radius, candidate_angle_unwrapped, label["middle"], total_length
                )
                candidate = label.copy()
                candidate["start_x"] = candidate_x
                candidate["start_y"] = candidate_y

                overlap_count = 0
                for existing in rearranged_labels:
                    if check_overlap(candidate, existing, total_length, 0.1):
                        overlap_count += 1

                if overlap_count == 0:
                    label["start_x"] = candidate_x
                    label["start_y"] = candidate_y
                    best_overlaps = 0
                    best_angle_unwrapped = candidate_angle_unwrapped
                    break

                if best_overlaps is None or overlap_count < best_overlaps:
                    best_x = candidate_x
                    best_y = candidate_y
                    best_overlaps = overlap_count
                    best_angle_unwrapped = candidate_angle_unwrapped

            if best_overlaps == 0:
                break

        if best_overlaps != 0:
            if best_angle_unwrapped is None:
                # Fallback: preserve ordering even when local placement is crowded.
                best_angle_unwrapped = lower_bound
                best_x, best_y = calculate_coordinates(
                    center_x, center_y, x_radius, y_radius, best_angle_unwrapped, label["middle"], total_length
                )
            label["start_x"] = best_x if best_x is not None else 0.0
            label["start_y"] = best_y if best_y is not None else 0.0

        label["angle_unwrapped"] = best_angle_unwrapped if best_angle_unwrapped is not None else lower_bound
        prev_angle_unwrapped = label["angle_unwrapped"]
        rearranged_labels.append(label)

    return rearranged_labels


def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def sort_labels(labels):
    return sorted(labels, key=lambda x: x["middle"])


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
    y_margin=0.3,
    max_angle_shift_deg=25.0,
    max_iterations=60,
):
    def move_label(label, angle):
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def check_overlap(label1, label2, total_length):
        return y_overlap(label1, label2, total_length, y_margin) and x_overlap(
            label1, label2, minimum_margin=1.5
        )

    def overlap_count(label_idx, candidate_label=None):
        current = candidate_label if candidate_label is not None else labels[label_idx]
        n_labels = len(labels)
        if n_labels <= 1:
            return 0
        # In angular-sorted order, overlaps are local. Limit checks to nearby labels
        # to avoid quadratic work on dense records.
        neighbor_span = 6
        overlaps = 0
        checked: set[int] = set()
        for delta in range(1, neighbor_span + 1):
            prev_idx = (label_idx - delta) % n_labels
            next_idx = (label_idx + delta) % n_labels
            checked.add(prev_idx)
            checked.add(next_idx)
        for other_idx in checked:
            if other_idx == label_idx:
                continue
            if check_overlap(current, labels[other_idx], total_length):
                overlaps += 1
        return overlaps

    def update_label_coordinates(label: dict, angle_unwrapped: float) -> None:
        label["angle_unwrapped"] = angle_unwrapped
        label["start_x"], label["start_y"] = move_label(label, angle_unwrapped)

    def local_density_weight(label_idx: int) -> float:
        dense_threshold_deg = 8.0
        n_labels = len(labels)
        target = labels[label_idx]["target_angle_unwrapped"]
        close_neighbors = 0
        if label_idx > 0:
            prev_target = labels[label_idx - 1]["target_angle_unwrapped"]
            if (target - prev_target) <= dense_threshold_deg:
                close_neighbors += 1
        if label_idx < n_labels - 1:
            next_target = labels[label_idx + 1]["target_angle_unwrapped"]
            if (next_target - target) <= dense_threshold_deg:
                close_neighbors += 1
        return 1.0 + close_neighbors

    def is_dense_center(label_idx: int) -> bool:
        if label_idx <= 0 or label_idx >= len(labels) - 1:
            return False
        dense_threshold_deg = 8.0
        left_gap = labels[label_idx]["target_angle_unwrapped"] - labels[label_idx - 1]["target_angle_unwrapped"]
        right_gap = labels[label_idx + 1]["target_angle_unwrapped"] - labels[label_idx]["target_angle_unwrapped"]
        return left_gap <= dense_threshold_deg and right_gap <= dense_threshold_deg

    def middle_proximity_penalty(label_idx: int, candidate_angle: float) -> float:
        """
        Penalize configurations where edge labels in a dense triplet are closer to
        their features than the center label.
        """
        target = labels[label_idx]["target_angle_unwrapped"]
        edge_delta = abs(candidate_angle - target)
        penalty = 0.0

        center_right = label_idx + 1
        if center_right < len(labels) - 1 and is_dense_center(center_right):
            center_delta = abs(
                labels[center_right]["angle_unwrapped"] - labels[center_right]["target_angle_unwrapped"]
            )
            if edge_delta < center_delta:
                penalty += center_delta - edge_delta

        center_left = label_idx - 1
        if center_left > 0 and is_dense_center(center_left):
            center_delta = abs(labels[center_left]["angle_unwrapped"] - labels[center_left]["target_angle_unwrapped"])
            if edge_delta < center_delta:
                penalty += center_delta - edge_delta

        return penalty

    def enforce_order(min_gap: float) -> bool:
        changed = False
        for i in range(1, len(labels)):
            minimum_allowed = labels[i - 1]["angle_unwrapped"] + min_gap
            if labels[i]["angle_unwrapped"] < minimum_allowed:
                update_label_coordinates(labels[i], minimum_allowed)
                changed = True
        return changed

    labels = sort_labels(labels)
    angle_step = 0.25
    max_steps = int(max_angle_shift_deg / angle_step)
    min_order_gap = 0.05

    for label in labels:
        target_angle_unwrapped = label.get("target_angle_unwrapped")
        if target_angle_unwrapped is None:
            target_angle_unwrapped = angle_from_middle_unwrapped(label["middle"], total_length)
            label["target_angle_unwrapped"] = target_angle_unwrapped
        label["target_angle"] = target_angle_unwrapped % 360.0

        current_angle_mod = calculate_angle_degrees(
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
            normalize=False,
        )
        angle_unwrapped = normalize_angle_near_reference(current_angle_mod, target_angle_unwrapped)
        update_label_coordinates(label, angle_unwrapped)

    enforce_order(min_order_gap)

    for _ in range(max_iterations):
        changes_made = False
        if enforce_order(min_order_gap):
            changes_made = True

        iteration_order = sorted(range(len(labels)), key=local_density_weight)
        for label_idx in iteration_order:
            label = labels[label_idx]
            current_angle = label["angle_unwrapped"]
            target_angle = label["target_angle_unwrapped"]
            current_overlaps = overlap_count(label_idx)
            current_target_delta = abs(current_angle - target_angle)
            current_middle_penalty = middle_proximity_penalty(label_idx, current_angle)

            if current_overlaps == 0 and current_middle_penalty == 0:
                continue

            lower_bound = target_angle - max_angle_shift_deg
            upper_bound = target_angle + max_angle_shift_deg
            if label_idx > 0:
                lower_bound = max(lower_bound, labels[label_idx - 1]["angle_unwrapped"] + min_order_gap)
            if label_idx < len(labels) - 1:
                upper_bound = min(upper_bound, labels[label_idx + 1]["angle_unwrapped"] - min_order_gap)
            if lower_bound > upper_bound:
                lower_bound = labels[label_idx - 1]["angle_unwrapped"] + min_order_gap if label_idx > 0 else lower_bound
                upper_bound = labels[label_idx + 1]["angle_unwrapped"] - min_order_gap if label_idx < len(labels) - 1 else upper_bound
                if lower_bound > upper_bound:
                    continue

            anchor_weight = local_density_weight(label_idx)
            best_angle = current_angle
            best_x = label["start_x"]
            best_y = label["start_y"]
            best_overlaps = current_overlaps
            best_target_delta = current_target_delta
            best_move = 0.0
            best_middle_penalty = current_middle_penalty

            search_complete = False
            for step in range(max_steps + 1):
                if step == 0:
                    offsets = [0.0]
                else:
                    offset = step * angle_step
                    offsets = [offset, -offset]

                for offset in offsets:
                    candidate_angle = target_angle + offset
                    if not (lower_bound <= candidate_angle <= upper_bound):
                        continue
                    candidate_x, candidate_y = move_label(label, candidate_angle)
                    label_copy = label.copy()
                    label_copy["start_x"] = candidate_x
                    label_copy["start_y"] = candidate_y

                    candidate_overlaps = overlap_count(label_idx, label_copy)
                    candidate_target_delta = abs(candidate_angle - target_angle)
                    candidate_move = abs(candidate_angle - current_angle)
                    candidate_middle_penalty = middle_proximity_penalty(label_idx, candidate_angle)
                    proximity_penalty_factor = 4.0
                    current_score = (
                        best_overlaps,
                        anchor_weight * best_target_delta + proximity_penalty_factor * best_middle_penalty,
                        best_move,
                    )
                    candidate_score = (
                        candidate_overlaps,
                        anchor_weight * candidate_target_delta + proximity_penalty_factor * candidate_middle_penalty,
                        candidate_move,
                    )

                    if candidate_score < current_score:
                        best_angle = candidate_angle
                        best_x = candidate_x
                        best_y = candidate_y
                        best_overlaps = candidate_overlaps
                        best_target_delta = candidate_target_delta
                        best_move = candidate_move
                        best_middle_penalty = candidate_middle_penalty

                    if best_overlaps == 0 and best_target_delta == 0 and best_middle_penalty == 0:
                        search_complete = True
                        break
                if search_complete:
                    break

            improved = False
            if best_overlaps < current_overlaps:
                improved = True
            elif best_overlaps == current_overlaps and best_middle_penalty < current_middle_penalty:
                improved = True

            if improved:
                label["start_x"] = best_x
                label["start_y"] = best_y
                label["angle_unwrapped"] = best_angle
                changes_made = True
        if not changes_made:
            break  # No changes were made in this iteration, so we can stop

    enforce_order(min_order_gap)
    return labels


def rearrange_labels_fc(
    labels,
    feature_radius,
    total_length,
    genome_len,
    config_dict,
    strands,
    is_outer,
    arena_outer_radius: float | None = None,
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

    # Optional "arena" override: scale the ellipse so it fits the desired outer radius (px).
    if arena_outer_radius is not None:
        arena_outer_radius = float(arena_outer_radius)
        if arena_outer_radius > 0:
            current_max = max(abs(x_radius), abs(y_radius))
            if current_max == 0:
                x_radius = arena_outer_radius
                y_radius = arena_outer_radius
            else:
                scale = arena_outer_radius / current_max
                x_radius *= scale
                y_radius *= scale

    center_y = 0
    center_x = default_center_x
    start_angle = 0
    end_angle = 360

    labels = sorted(labels, key=lambda x: x["middle"])
    # Initial placement of labels
    labels = place_labels_on_arc_fc(labels, center_x, center_y, x_radius, y_radius, start_angle, end_angle, total_length)

    # Apply improved label placement
    labels = improved_label_placement_fc(
        labels, center_x, center_y, x_radius, y_radius, feature_radius, total_length, start_angle, end_angle
    )

    return labels


def prepare_label_list(
    feature_dict,
    total_length,
    radius,
    track_ratio,
    config_dict,
    cfg: GbdrawConfig | None = None,
    *,
    outer_arena: tuple[float, float] | None = None,
):
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    embedded_labels = []
    outer_labels = []
    inner_labels = []
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
        else:
            label_entry = dict()
            longest_segment_length = 0
            is_embedded = False
            label_middle = 0
            coordinate_strand: str = "undefined"
            feature_location_list = feature_object.location
            list_of_coordinates = feature_object.coordinates
            # `FeatureObject.strand` is derived at creation time; keep placement pure.

            feature_location_count = 0
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
            # Get track_id for overlap resolution
            track_id = getattr(feature_object, 'feature_track_id', 0)
            factors: list[float] = calculate_feature_position_factors_circular(
                total_length, coordinate_strand, track_ratio, cds_ratio, offset, track_type, strandedness, track_id
            )
            # Store track_id in label entry for embedded label drawing
            label_entry["track_id"] = track_id
            bbox_width_px, bbox_height_px = calculate_bbox_dimensions(feature_label_text, font_family, font_size, interval)
            label_middle = longeset_segment_middle
            label_as_feature_length = total_length * (1.1 * bbox_width_px) / (2 * math.pi * radius)
            label_start = label_middle - (label_as_feature_length / 2)
            label_end = label_middle + (label_as_feature_length / 2)
            feature_middle_x: float = (radius * factors[1]) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
            feature_middle_y: float = (radius * factors[1]) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            is_outer_label = (feature_object.strand == "positive") or (allow_inner_labels is False)
            if is_outer_label:
                if outer_arena is not None:
                    # Use the inner edge of the arena as the "anchor" radius for leader lines.
                    anchor_radius = float(outer_arena[0])
                    middle_x = anchor_radius * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                    middle_y = anchor_radius * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
                else:
                    middle_x = (radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                    middle_y = (radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            else:
                middle_x = (inner_radius_factor * radius) * math.cos(math.radians(360.0 * (label_middle / total_length) - 90))
                middle_y = (inner_radius_factor * radius) * math.sin(math.radians(360.0 * (label_middle / total_length) - 90))
            if label_start > longest_segment_start and label_end < longest_segment_end:
                is_embedded = True
            else:
                is_embedded = False
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
                if label_entry["middle"] > (total_length / 2):
                    if feature_object.strand == "positive" or allow_inner_labels is False:
                        label_entry["is_inner"] = False
                        outer_labels.append(label_entry)
                    else:
                        label_entry["is_inner"] = True
                        inner_labels.append(label_entry)

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
        arena_outer_radius=(float(outer_arena[1]) if outer_arena is not None else None),
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
    "calculate_angle_degrees",
    "calculate_angle_for_y",
    "calculate_coordinates",
    "angle_from_middle",
    "angle_from_middle_unwrapped",
    "angle_delta_deg",
    "normalize_angle_near_reference",
    "euclidean_distance",
    "improved_label_placement_fc",
    "place_labels_on_arc_fc",
    "prepare_label_list",
    "rearrange_labels_fc",
    "sort_labels",
    "x_overlap",
    "y_overlap",
]


