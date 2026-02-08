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

MIN_BBOX_GAP_RATIO = 0.05
HEAVY_CLUSTER_RELAX_MAX_LABELS = 70
FULL_SCAN_LABEL_LIMIT = 70


def minimum_bbox_gap_px(label1: dict, label2: dict, base_margin_px: float = 0.0) -> float:
    """Return a minimum spacing margin in pixels derived from label bbox size."""
    width_scale = max(float(label1.get("width_px", 0.0)), float(label2.get("width_px", 0.0)))
    height_scale = max(float(label1.get("height_px", 0.0)), float(label2.get("height_px", 0.0)))
    ratio_gap = MIN_BBOX_GAP_RATIO * max(width_scale, height_scale)
    return max(float(base_margin_px), ratio_gap)


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
        effective_margin = minimum_bbox_gap_px(label1, label2, base_margin_px=max(margin, 2.0))
        return y_overlap(label1, label2, total_length, effective_margin) and x_overlap(
            label1, label2, minimum_margin=effective_margin
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
    labels = sort_labels(labels)
    if not labels:
        return []

    # Preserve original sorted order; optimization may rotate labels to avoid
    # placing a dense cluster across the 0/360-degree seam.
    for idx, label in enumerate(labels):
        label["_opt_original_idx"] = idx

    if len(labels) > 2:
        raw_targets = [angle_from_middle_unwrapped(label["middle"], total_length) for label in labels]
        gap_values = [raw_targets[i + 1] - raw_targets[i] for i in range(len(labels) - 1)]
        gap_values.append(raw_targets[0] + 360.0 - raw_targets[-1])
        largest_gap_idx = max(range(len(gap_values)), key=gap_values.__getitem__)
        if largest_gap_idx < len(labels) - 1:
            rotate_by = largest_gap_idx + 1
            labels = labels[rotate_by:] + labels[:rotate_by]

    # Reassign targets as a monotonic unwrapped sequence in optimization order.
    prev_target_unwrapped: float | None = None
    for label in labels:
        raw_target = angle_from_middle_unwrapped(label["middle"], total_length)
        if prev_target_unwrapped is None:
            target_angle_unwrapped = raw_target
        else:
            target_angle_unwrapped = normalize_angle_near_reference(raw_target, prev_target_unwrapped)
            while target_angle_unwrapped <= prev_target_unwrapped:
                target_angle_unwrapped += 360.0
        label["target_angle_unwrapped"] = target_angle_unwrapped
        label["target_angle"] = target_angle_unwrapped % 360.0
        prev_target_unwrapped = target_angle_unwrapped

    pair_margins: list[list[float]] = []

    def move_label(label, angle):
        new_x = center_x + x_radius * math.cos(math.radians(angle))
        new_y = center_y + y_radius * math.sin(math.radians(angle))
        return new_x, new_y

    def check_overlap(label1, label2, total_length, effective_margin: float):
        return y_overlap(label1, label2, total_length, effective_margin) and x_overlap(
            label1, label2, minimum_margin=effective_margin
        )

    def pair_margin(idx1: int, idx2: int) -> float:
        if idx1 < idx2:
            return pair_margins[idx1][idx2]
        return pair_margins[idx2][idx1]

    def overlap_count(label_idx, candidate_label=None, *, full_scan: bool = False):
        current = candidate_label if candidate_label is not None else labels[label_idx]
        n_labels = len(labels)
        if n_labels <= 1:
            return 0
        checked_indices: set[int] = set()
        if full_scan:
            checked_indices = set(range(n_labels))
            checked_indices.discard(label_idx)
        else:
            # In angular-sorted order, most overlaps are local. Start with nearby labels.
            neighbor_span = 6
            for delta in range(1, neighbor_span + 1):
                prev_idx = (label_idx - delta) % n_labels
                next_idx = (label_idx + delta) % n_labels
                checked_indices.add(prev_idx)
                checked_indices.add(next_idx)
        overlaps = 0
        for other_idx in checked_indices:
            if other_idx == label_idx:
                continue
            if check_overlap(current, labels[other_idx], total_length, pair_margin(label_idx, other_idx)):
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

    def count_remaining_overlaps() -> int:
        overlap_count_all = 0
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if check_overlap(labels[i], labels[j], total_length, pair_margin(i, j)):
                    overlap_count_all += 1
        return overlap_count_all

    def count_overlaps_for_label_list(label_list: list[dict]) -> int:
        overlap_count_all = 0
        for i in range(len(label_list)):
            for j in range(i + 1, len(label_list)):
                if check_overlap(label_list[i], label_list[j], total_length, pair_margin(i, j)):
                    overlap_count_all += 1
        return overlap_count_all

    def find_overlap_components() -> list[list[int]]:
        n_labels = len(labels)
        adjacency: dict[int, set[int]] = {i: set() for i in range(n_labels)}
        for i in range(n_labels):
            for j in range(i + 1, n_labels):
                if check_overlap(labels[i], labels[j], total_length, pair_margin(i, j)):
                    adjacency[i].add(j)
                    adjacency[j].add(i)

        components: list[list[int]] = []
        visited: set[int] = set()
        for idx in range(n_labels):
            if idx in visited or not adjacency[idx]:
                continue
            stack = [idx]
            component: list[int] = []
            visited.add(idx)
            while stack:
                current = stack.pop()
                component.append(current)
                for neigh in adjacency[current]:
                    if neigh not in visited:
                        visited.add(neigh)
                        stack.append(neigh)
            if component:
                components.append(sorted(component))
        return components

    def merge_close_components(components: list[list[int]], max_gap: int = 2) -> list[list[int]]:
        if not components:
            return []
        sorted_components = sorted((sorted(comp) for comp in components), key=lambda comp: comp[0])
        merged: list[list[int]] = [sorted_components[0]]
        for comp in sorted_components[1:]:
            if comp[0] - merged[-1][-1] <= max_gap:
                merged[-1] = sorted(set(merged[-1] + comp))
            else:
                merged.append(comp)
        return merged

    def relax_overlapping_clusters(shift_limit_deg: float) -> bool:
        """
        When local single-label optimization stalls, nudge clusters as a group.
        This allows neighboring labels to move together while keeping the
        cluster center label relatively close to its feature.
        """
        if len(labels) > HEAVY_CLUSTER_RELAX_MAX_LABELS:
            return False

        components = merge_close_components(find_overlap_components(), max_gap=2)
        if not components:
            return False

        cluster_step_limit = min(int(shift_limit_deg / angle_step), 40)
        n_labels = len(labels)
        changed = False

        for component in components:
            component_center = component[len(component) // 2]
            cluster_padding = min(6, max(2, len(component) + 1))
            cluster_start = max(0, component[0] - cluster_padding)
            cluster_end = min(n_labels - 1, component[-1] + cluster_padding)
            cluster_indices = list(range(cluster_start, cluster_end + 1))

            current_score = (
                count_remaining_overlaps(),
                abs(labels[component_center]["angle_unwrapped"] - labels[component_center]["target_angle_unwrapped"]),
                sum(abs(labels[idx]["angle_unwrapped"] - labels[idx]["target_angle_unwrapped"]) for idx in cluster_indices),
            )
            best_score = current_score
            best_ordered: list[float] | None = None

            center_candidates = [component_center]
            if len(component) >= 2:
                center_candidates.append(component[0])
                center_candidates.append(component[-1])
            center_candidates = sorted(set(center_candidates))
            center_shift_limit_steps = min(int(shift_limit_deg / angle_step), 48)
            center_shift_steps = range(-center_shift_limit_steps, center_shift_limit_steps + 1, 4)

            for center_idx in center_candidates:
                for step in range(0, cluster_step_limit + 1):
                    spread = step * angle_step
                    for shift_step in center_shift_steps:
                        center_shift = shift_step * angle_step
                        ordered = [float(labels[i]["angle_unwrapped"]) for i in range(n_labels)]

                        for idx in cluster_indices:
                            rel = idx - center_idx
                            distance_weight = 1.0 + 0.25 * max(0, abs(rel) - 1)
                            target = float(labels[idx]["target_angle_unwrapped"])
                            lower_bound = target - shift_limit_deg
                            upper_bound = target + shift_limit_deg
                            candidate = target + center_shift + (rel * spread * distance_weight)
                            if candidate < lower_bound:
                                candidate = lower_bound
                            elif candidate > upper_bound:
                                candidate = upper_bound
                            ordered[idx] = candidate

                        # Keep global ordering constraints valid.
                        for i in range(1, n_labels):
                            minimum_allowed = ordered[i - 1] + min_order_gap
                            if ordered[i] < minimum_allowed:
                                ordered[i] = minimum_allowed
                        for i in range(n_labels - 2, -1, -1):
                            maximum_allowed = ordered[i + 1] - min_order_gap
                            if ordered[i] > maximum_allowed:
                                ordered[i] = maximum_allowed

                        # Respect per-label shift budget.
                        violates_budget = False
                        for idx in range(n_labels):
                            target = float(labels[idx]["target_angle_unwrapped"])
                            if abs(ordered[idx] - target) > shift_limit_deg + 1e-6:
                                violates_budget = True
                                break
                        if violates_budget:
                            continue

                        temp_labels = [label.copy() for label in labels]
                        for idx in range(n_labels):
                            temp_labels[idx]["angle_unwrapped"] = ordered[idx]
                            temp_labels[idx]["start_x"], temp_labels[idx]["start_y"] = move_label(
                                temp_labels[idx], ordered[idx]
                            )

                        candidate_score = (
                            count_overlaps_for_label_list(temp_labels),
                            abs(
                                temp_labels[component_center]["angle_unwrapped"]
                                - temp_labels[component_center]["target_angle_unwrapped"]
                            ),
                            sum(
                                abs(
                                    temp_labels[idx]["angle_unwrapped"]
                                    - temp_labels[idx]["target_angle_unwrapped"]
                                )
                                for idx in cluster_indices
                            ),
                        )
                        if candidate_score < best_score:
                            best_score = candidate_score
                            best_ordered = ordered
                        if best_score[0] == 0:
                            break
                    if best_score[0] == 0:
                        break
                if best_score[0] == 0:
                    break

            if best_ordered is not None and best_score < current_score:
                for idx in range(n_labels):
                    update_label_coordinates(labels[idx], best_ordered[idx])
                changed = True

        return changed

    def optimize_with_shift_limit(shift_limit_deg: float, *, full_scan_checks: bool) -> None:
        max_steps = int(shift_limit_deg / angle_step)
        for _ in range(max_iterations):
            changes_made = False
            if enforce_order(min_order_gap):
                changes_made = True

            iteration_order = sorted(range(len(labels)), key=local_density_weight)
            for label_idx in iteration_order:
                label = labels[label_idx]
                current_angle = label["angle_unwrapped"]
                target_angle = label["target_angle_unwrapped"]
                current_overlaps = overlap_count(label_idx, full_scan=full_scan_checks)
                current_target_delta = abs(current_angle - target_angle)
                current_middle_penalty = middle_proximity_penalty(label_idx, current_angle)

                if current_overlaps == 0 and current_middle_penalty == 0:
                    continue

                lower_bound = target_angle - shift_limit_deg
                upper_bound = target_angle + shift_limit_deg
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

                        candidate_overlaps = overlap_count(label_idx, label_copy, full_scan=full_scan_checks)
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
                break
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            for _ in range(3):
                if not relax_overlapping_clusters(shift_limit_deg):
                    break
                if count_remaining_overlaps() == 0:
                    break

    def separate_remaining_adjacent_overlaps(shift_limit_deg: float) -> bool:
        changed = False
        n_labels = len(labels)
        if n_labels < 2:
            return False

        for _ in range(160):
            pass_changed = False
            for i in range(n_labels - 1):
                if not check_overlap(labels[i], labels[i + 1], total_length, pair_margin(i, i + 1)):
                    continue

                pair_changed = False
                for _ in range(80):
                    can_shift_left = True
                    for k in range(0, i + 1):
                        target = float(labels[k]["target_angle_unwrapped"])
                        if labels[k]["angle_unwrapped"] - angle_step < target - shift_limit_deg - 1e-6:
                            can_shift_left = False
                            break

                    can_shift_right = True
                    for k in range(i + 1, n_labels):
                        target = float(labels[k]["target_angle_unwrapped"])
                        if labels[k]["angle_unwrapped"] + angle_step > target + shift_limit_deg + 1e-6:
                            can_shift_right = False
                            break

                    if not can_shift_left and not can_shift_right:
                        break

                    if can_shift_left:
                        for k in range(0, i + 1):
                            update_label_coordinates(labels[k], labels[k]["angle_unwrapped"] - angle_step)
                    if can_shift_right:
                        for k in range(i + 1, n_labels):
                            update_label_coordinates(labels[k], labels[k]["angle_unwrapped"] + angle_step)

                    enforce_order(min_order_gap)
                    if not check_overlap(labels[i], labels[i + 1], total_length, pair_margin(i, i + 1)):
                        pair_changed = True
                        break

                if pair_changed:
                    pass_changed = True
                    changed = True

            if not pass_changed:
                break

        return changed

    n_labels = len(labels)
    pair_margins = [[0.0 for _ in range(n_labels)] for _ in range(n_labels)]
    base_margin = max(y_margin, 1.5)
    for i in range(n_labels):
        for j in range(i + 1, n_labels):
            pair_margins[i][j] = minimum_bbox_gap_px(labels[i], labels[j], base_margin_px=base_margin)

    angle_step = 0.25
    min_order_gap = 0.05

    def snapshot_angles() -> list[float]:
        return [float(label["angle_unwrapped"]) for label in labels]

    def restore_angles(angle_snapshot: list[float]) -> None:
        for idx, angle in enumerate(angle_snapshot):
            update_label_coordinates(labels[idx], angle)

    def reset_angles_to_targets() -> None:
        for label in labels:
            update_label_coordinates(label, float(label["target_angle_unwrapped"]))
        enforce_order(min_order_gap)

    def placement_score() -> tuple[int, float]:
        return (
            count_remaining_overlaps(),
            sum(abs(float(label["angle_unwrapped"]) - float(label["target_angle_unwrapped"])) for label in labels),
        )

    for label in labels:
        target_angle_unwrapped = float(label["target_angle_unwrapped"])

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

    shift_schedule = [max_angle_shift_deg]
    if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS:
        for extra_shift in (35.0, 45.0, 60.0):
            if extra_shift > max_angle_shift_deg:
                shift_schedule.append(extra_shift)

    for phase_idx, shift_limit in enumerate(shift_schedule):
        run_full_scan = (
            len(labels) <= FULL_SCAN_LABEL_LIMIT
            and len(shift_schedule) > 1
            and phase_idx == (len(shift_schedule) - 1)
        )
        optimize_with_shift_limit(shift_limit, full_scan_checks=run_full_scan)
        if count_remaining_overlaps() == 0:
            break

    if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
        separate_remaining_adjacent_overlaps(max(shift_schedule))

    if count_remaining_overlaps() > 0:
        baseline_angles = snapshot_angles()
        baseline_score = placement_score()
        optimize_with_shift_limit(
            max(shift_schedule),
            full_scan_checks=(len(labels) <= FULL_SCAN_LABEL_LIMIT),
        )
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            separate_remaining_adjacent_overlaps(max(shift_schedule))
        continued_score = placement_score()
        if continued_score > baseline_score:
            restore_angles(baseline_angles)
            continued_score = baseline_score

    if count_remaining_overlaps() > 0:
        baseline_angles = snapshot_angles()
        baseline_score = placement_score()
        reset_angles_to_targets()
        optimize_with_shift_limit(
            max(shift_schedule),
            full_scan_checks=(len(labels) <= FULL_SCAN_LABEL_LIMIT),
        )
        if len(labels) <= HEAVY_CLUSTER_RELAX_MAX_LABELS and count_remaining_overlaps() > 0:
            separate_remaining_adjacent_overlaps(max(shift_schedule))
        retry_score = placement_score()
        if retry_score > baseline_score:
            restore_angles(baseline_angles)

    enforce_order(min_order_gap)
    labels.sort(key=lambda x: x["_opt_original_idx"])
    for label in labels:
        label.pop("_opt_original_idx", None)
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
    "minimum_bbox_gap_px",
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


