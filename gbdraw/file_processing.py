#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.file_processing` contained multiple responsibilities:
  - genome I/O (GenBank, GFF3+FASTA)
  - color tables / palettes
  - comparison (BLAST) loading
  - export (SVG + optional CairoSVG conversion)
  - config loading

As part of responsibility separation, implementations were moved to:
  - `gbdraw.io.*`
  - `gbdraw.config.*`
  - `gbdraw.render.*`

New code should prefer importing from those modules directly.
"""

from .io.genome import (  # noqa: F401
    load_gbks,
    load_gff_fasta,
    merge_gff_fasta_records,
    scan_features_recursive,
    filter_features_by_type,
)
from .io.comparisons import load_comparisons  # noqa: F401
from .io.colors import (  # noqa: F401
    load_default_colors,
    read_color_table,
    resolve_color_to_hex,
)
from .render.export import parse_formats, save_figure  # noqa: F401
from .config.toml import load_config_toml  # noqa: F401

__all__ = [
    # genome
    "load_gbks",
    "load_gff_fasta",
    "merge_gff_fasta_records",
    "scan_features_recursive",
    "filter_features_by_type",
    # comparisons
    "load_comparisons",
    # colors
    "load_default_colors",
    "read_color_table",
    "resolve_color_to_hex",
    # export
    "parse_formats",
    "save_figure",
    # config
    "load_config_toml",
]


