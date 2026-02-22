"""Option bundles for the public API.

These dataclasses provide a lighter-weight entry point than the long list of
keyword arguments in assemble_* helpers. They are optional and additive.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.tracks import TrackSpec  # type: ignore[reportMissingImports]


@dataclass(frozen=True)
class ColorOptions:
    """Color table and palette inputs."""

    color_table: DataFrame | None = None
    color_table_file: str | None = None
    default_colors: DataFrame | None = None
    default_colors_palette: str = "default"
    default_colors_file: str | None = None


@dataclass(frozen=True)
class TrackOptions:
    """Track layout options (currently circular-only)."""

    track_specs: Sequence[str | TrackSpec] | None = None


@dataclass(frozen=True)
class OutputOptions:
    """High-level output/layout options used by the assemblers."""

    output_prefix: str = "out"
    legend: str = "right"


@dataclass(frozen=True)
class DiagramOptions:
    """Bundled options for diagram assembly helpers."""

    config: GbdrawConfig | dict | None = None
    config_overrides: Mapping[str, object] | None = None
    colors: ColorOptions | None = None
    tracks: TrackOptions | None = None
    output: OutputOptions | None = None
    selected_features_set: Sequence[str] | None = None
    feature_shapes: Mapping[str, str] | None = None
    dinucleotide: str = "GC"
    window: int | None = None
    step: int | None = None
    species: str | None = None
    strain: str | None = None
    blast_files: Sequence[str] | None = None
    evalue: float = 1e-5
    bitscore: float = 50.0
    identity: float = 70.0


__all__ = [
    "ColorOptions",
    "DiagramOptions",
    "OutputOptions",
    "TrackOptions",
]
