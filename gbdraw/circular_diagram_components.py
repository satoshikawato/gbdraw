#!/usr/bin/env python
# coding: utf-8

"""Backwards-compatibility shim.

The implementation has been moved to `gbdraw.diagrams.circular`.
This module re-exports the public API for backwards compatibility.
"""

from gbdraw.diagrams.circular.assemble import (
    add_record_on_circular_canvas,
    assemble_circular_diagram,
    plot_circular_diagram,
)
from gbdraw.diagrams.circular.builders import (
    add_axis_group_on_canvas,
    add_gc_content_group_on_canvas,
    add_gc_skew_group_on_canvas,
    add_labels_group_on_canvas,
    add_legend_group_on_canvas,
    add_record_definition_group_on_canvas,
    add_record_group_on_canvas,
    add_tick_group_on_canvas,
)
from gbdraw.diagrams.circular.positioning import (
    center_group_on_canvas,
    place_legend_on_canvas,
)

__all__ = [
    # Main assembly functions
    "add_record_on_circular_canvas",
    "assemble_circular_diagram",
    "plot_circular_diagram",
    # Builder functions
    "add_axis_group_on_canvas",
    "add_gc_content_group_on_canvas",
    "add_gc_skew_group_on_canvas",
    "add_labels_group_on_canvas",
    "add_legend_group_on_canvas",
    "add_record_definition_group_on_canvas",
    "add_record_group_on_canvas",
    "add_tick_group_on_canvas",
    # Positioning functions
    "center_group_on_canvas",
    "place_legend_on_canvas",
]
