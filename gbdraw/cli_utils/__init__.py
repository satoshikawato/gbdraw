#!/usr/bin/env python
# coding: utf-8

"""CLI module with shared argument definitions and utilities."""

from .common import (
    add_analysis_args,
    add_color_args,
    add_feature_args,
    add_input_args,
    add_label_args,
    calculate_window_step,
    handle_output_formats,
    setup_logging,
    validate_input_args,
    validate_label_args,
)

__all__ = [
    "add_analysis_args",
    "add_color_args",
    "add_feature_args",
    "add_input_args",
    "add_label_args",
    "calculate_window_step",
    "handle_output_formats",
    "setup_logging",
    "validate_input_args",
    "validate_label_args",
]
