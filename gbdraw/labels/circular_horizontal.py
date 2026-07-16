"""Compatibility adapter for the established horizontal circular label strategy."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from .circular_types import CircularLabelLayoutResult


def keep_horizontal_layout(labels: Sequence[Mapping[str, object]]) -> CircularLabelLayoutResult:
    """Wrap legacy horizontal dictionaries in the common placement result contract."""
    return CircularLabelLayoutResult(labels=tuple(labels), content_bounds=None)


__all__ = ["keep_horizontal_layout"]
