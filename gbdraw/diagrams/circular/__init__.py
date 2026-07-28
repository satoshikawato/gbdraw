"""Circular diagram assembly helpers (internal)."""

from __future__ import annotations

from typing import Any


def __getattr__(name: str) -> Any:
    if name == "assemble_circular_diagram":
        from .assemble import assemble_circular_diagram

        return assemble_circular_diagram
    raise AttributeError(name)

__all__ = ["assemble_circular_diagram"]
