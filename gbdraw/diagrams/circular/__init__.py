"""Circular diagram assembly helpers (internal)."""

from __future__ import annotations

from typing import Any


def __getattr__(name: str) -> Any:
    if name in {"assemble_circular_diagram", "plot_circular_diagram"}:
        from .assemble import assemble_circular_diagram, plot_circular_diagram

        return {
            "assemble_circular_diagram": assemble_circular_diagram,
            "plot_circular_diagram": plot_circular_diagram,
        }[name]
    raise AttributeError(name)

__all__ = [
    "assemble_circular_diagram",
    "plot_circular_diagram",
]
