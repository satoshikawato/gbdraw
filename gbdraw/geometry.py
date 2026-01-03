#!/usr/bin/env python
# coding: utf-8

import math


def deg_to_rad(deg):
    return deg * math.pi / 180.0


def calculate_coordinates(center_x, center_y, radius, angle_degrees):
    angle_radians = deg_to_rad(angle_degrees)
    x = center_x + radius * math.cos(angle_radians)
    y = center_y + radius * math.sin(angle_radians)
    return x, y


def calculate_angle_degrees(
    center_x,
    center_y,
    x,
    y,
    middle,
    start_angle,
    end_angle,
    total_length,
    x_radius,
    y_radius,
    normalize,
):
    """Calculate the angle in degrees from a point (x, y) relative to the center of an ellipse."""
    x_normalized = (x - center_x) / x_radius
    y_normalized = (y - center_y) / y_radius

    angle_radians = math.atan2(y_normalized, x_normalized)
    angle_degrees = math.degrees(angle_radians)
    return angle_degrees


def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


__all__ = [
    "calculate_angle_degrees",
    "calculate_coordinates",
    "deg_to_rad",
    "euclidean_distance",
]


