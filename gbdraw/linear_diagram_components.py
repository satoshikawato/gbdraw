#!/usr/bin/env python
# coding: utf-8

"""Linear diagram assembly utilities (compatibility façade).

Historically this module contained all linear diagram composition logic.
To improve cohesion and make the codebase easier to use as a library, implementations
were moved into `gbdraw.diagrams.linear.*`.

Public functions remain available from this module via re-exports.
"""

from .diagrams.linear.assemble import (  # type: ignore[reportMissingImports]
    assemble_linear_diagram,
    plot_linear_diagram,
)
from .diagrams.linear.builders import (  # type: ignore[reportMissingImports]
    add_comparison_on_linear_canvas,
    add_gc_content_group,
    add_gc_skew_group,
    add_legends_on_linear_canvas,
    add_length_bar_on_linear_canvas,
    add_record_definition_group,
    add_record_group,
)
from .diagrams.linear.positioning import (  # type: ignore[reportMissingImports]
    calculate_record_offsets,
    position_comparison_group,
    position_gc_content_group,
    position_gc_skew_group,
    position_length_bar_group,
    position_record_definition_group,
    position_record_group,
)
from .diagrams.linear.precalc import (  # type: ignore[reportMissingImports]
    _precalculate_definition_widths,
    _precalculate_label_dimensions,
)

__all__ = [
    "_precalculate_definition_widths",
    "_precalculate_label_dimensions",
    "add_comparison_on_linear_canvas",
    "add_gc_content_group",
    "add_gc_skew_group",
    "add_legends_on_linear_canvas",
    "add_length_bar_on_linear_canvas",
    "add_record_definition_group",
    "add_record_group",
    "assemble_linear_diagram",
    "calculate_record_offsets",
    "plot_linear_diagram",
    "position_comparison_group",
    "position_gc_content_group",
    "position_gc_skew_group",
    "position_length_bar_group",
    "position_record_definition_group",
    "position_record_group",
]


