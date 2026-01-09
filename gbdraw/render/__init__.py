"""Rendering backends (internal).

This package contains:
- drawers/: Low-level SVG element builders
- groups/: High-level SVG group assemblers
- export.py: Output format handling (SVG, PNG, PDF)
"""

from .export import save_figure, parse_formats

__all__ = ["save_figure", "parse_formats"]
