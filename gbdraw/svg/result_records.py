"""The records a Result keeps of what Python drew (generated into the Web's ``services/result-records.generated.js``).

The Web's Result executor writes them when it edits a Result; Python writes the Legend ones when it draws a request's
Legend row edits (``gbdraw.legend.table.apply_legend_row_edits``), as the executor would have written them, so the Web
reads both Results alike. ``data-gbdraw-base-<attribute>`` holds Python's value of a changed attribute (empty: Python
drew none); ``LEGEND_ORDER_RECORD`` on a feature Legend group holds Python's row order (its keys, as JSON) once the
rows are reordered; ``LEGEND_MERGED_ROW`` on a hidden row holds the key of the earlier row of its name and swatch that
the Legend draws for both (D-28).
"""

from __future__ import annotations

RESULT_BASE_ATTRIBUTE_PREFIX = "data-gbdraw-base-"
# The attributes whose drawn value a Result records, in the Web's order.
RESULT_BASE_ATTRIBUTES = ("fill", "stroke", "stroke-width", "display", "data-legend-key")
LEGEND_ORDER_RECORD = f"{RESULT_BASE_ATTRIBUTE_PREFIX}legend-order"
LEGEND_MERGED_ROW = "data-legend-merged"


def result_base_attribute(attribute: str) -> str:
    """The record of ``attribute``: ``data-gbdraw-base-<attribute>``."""

    return f"{RESULT_BASE_ATTRIBUTE_PREFIX}{attribute}"


__all__ = [
    "LEGEND_MERGED_ROW",
    "LEGEND_ORDER_RECORD",
    "RESULT_BASE_ATTRIBUTES",
    "RESULT_BASE_ATTRIBUTE_PREFIX",
    "result_base_attribute",
]
