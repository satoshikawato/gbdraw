#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from svgwrite import Drawing

from ..config.models import LinearRenderProfile
from ..core.sequence import determine_length_parameter

# The configured arrow coefficients are calibrated to these non-stranded heights.
_ARROW_LENGTH_REFERENCE_CDS_HEIGHT = {
    "short": 40.0,
    "long": 10.0,
}


@dataclass(frozen=True)
class _LinearPlotTrack:
    key: str
    top_extent: float
    bottom_extent: float
    gap_after: float = 0.0


def _calculate_linear_plot_track_layout(tracks: list[_LinearPlotTrack]) -> tuple[dict[str, float], dict[str, float], float]:
    """Return anchor offsets, legacy padding segments, and visual stack span."""
    if not tracks:
        return {}, {}, 0.0

    anchor_offsets: dict[str, float] = {}
    current_anchor = 0.0
    previous_track: _LinearPlotTrack | None = None
    for track in tracks:
        if previous_track is not None:
            current_anchor += previous_track.bottom_extent + previous_track.gap_after + track.top_extent
        anchor_offsets[track.key] = current_anchor
        previous_track = track

    first_track = tracks[0]
    last_track = tracks[-1]
    stack_span = (
        first_track.top_extent
        + anchor_offsets[last_track.key]
        + last_track.bottom_extent
        + last_track.gap_after
    )

    padding_segments: dict[str, float] = {}
    consumed_span = 0.0
    for index, track in enumerate(tracks):
        if index < len(tracks) - 1:
            next_track = tracks[index + 1]
            segment = anchor_offsets[next_track.key] - anchor_offsets[track.key]
        else:
            segment = stack_span - consumed_span
        padding_segments[track.key] = segment
        consumed_span += segment

    return anchor_offsets, padding_segments, stack_span


def _calculate_linear_plot_track_visual_bottom(
    tracks: list[_LinearPlotTrack],
    anchor_offsets: dict[str, float],
) -> float:
    """Return the drawn stack bottom relative to the shared plot-track anchor."""
    if not tracks:
        return 0.0
    return max(
        anchor_offsets[track.key] + track.bottom_extent
        for track in tracks
    )


class LinearCanvasConfigurator:
    """
    Configures the settings for a linear canvas used for genomic data visualization.

    Attributes:
    output_prefix (str): Prefix for the output file.
    fig_width (int): Width of the figure.
    vertical_offset (float): Vertical offset for alignment.
    horizontal_offset (float): Horizontal offset for alignment.
    vertical_padding (float): Vertical padding between elements.
    configured_track_spacing (float): Default spacing between linear data tracks.
    comparison_height (float): Height for comparison tracks.
    definition_gap (float): Minimum gap between definition text and record axis.
    default_cds_height (float): Default height for coding sequences.
    default_gc_height (float): Default height for GC content.
    show_gc (bool): Flag to display GC content.
    strandedness (bool): Flag to display strandedness.
    num_of_entries (int): Number of entries to visualize.
    align_center (bool): Flag to align content to the center.
    longest_genome (int): Length of the longest genome in the dataset.

    Methods:
    set_gc_height_and_gc_padding(): Sets plot track heights, offsets, and spacing.
    set_cds_height_and_cds_padding(): Sets CDS height and padding.
    set_arrow_length(): Sets the arrow length for representation.
    calculate_dimensions(): Calculates dimensions for the canvas.
    create_svg_canvas(): Creates and returns an SVG canvas for drawing.
    """

    def __init__(
        self,
        num_of_entries: int,
        longest_genome: int,
        profile: LinearRenderProfile,
        legend: str,
        output_prefix="out",
        has_comparisons: bool = False,
        depth_track_count: int = 1,
        depth_track_heights: list[float | None] | tuple[float | None, ...] | None = None,
    ):
        """
        Initializes the linear canvas configurator with given settings.

        Args:
        num_of_entries (int): Number of entries to visualize.
        longest_genome (int): Length of the longest genome in the dataset.
        profile (LinearRenderProfile): Resolved linear render settings.
        output_prefix (str, optional): Prefix for the output file. Defaults to 'out'.
        show_gc (bool, optional): Flag to display GC content. Defaults to False.
        strandedness (bool, optional): Flag to display strandedness. Defaults to False.
        align_center (bool, optional): Flag to align content to the center. Defaults to False.
        has_comparisons (bool, optional): Flag indicating if BLAST comparisons are present. Defaults to False.
        """
        cfg = profile.config
        self.output_prefix: str = output_prefix
        self._profile = profile
        self._cfg = cfg
        self.longest_genome: int = longest_genome
        self.fig_width: int = cfg.canvas.linear.width
        self.original_vertical_offset: float = cfg.canvas.linear.vertical_offset
        self.vertical_offset: float = cfg.canvas.linear.vertical_offset
        self.horizontal_offset: float = cfg.canvas.linear.horizontal_offset
        self.vertical_padding: float = cfg.canvas.linear.vertical_padding
        self.configured_track_spacing: float = cfg.canvas.linear.track_spacing
        self.has_comparisons: bool = has_comparisons
        # Store the configured comparison_height for BLAST ribbons
        self.configured_comparison_height: float = cfg.canvas.linear.comparison_height
        # For initial dimension calculation, use comparison_height only when BLAST is present
        # The actual inter-record spacing will be recalculated in assemble_linear_diagram
        # based on feature track heights
        self.comparison_height: float = cfg.canvas.linear.comparison_height if has_comparisons else 0
        self.definition_gap: float = cfg.canvas.linear.definition_gap
        self.length_threshold = cfg.labels.length_threshold.linear
        self.length_param = determine_length_parameter(self.longest_genome, self.length_threshold)
        self.default_cds_height: float = getattr(cfg.canvas.linear.default_cds_height, self.length_param)
        self.default_gc_height: float = cfg.canvas.linear.default_gc_height
        self.default_depth_height: float = cfg.canvas.linear.depth_height
        self.configured_depth_padding: float = cfg.canvas.linear.depth_padding
        self.dpi: int = cfg.canvas.dpi
        self.depth_track_count: int = max(1, int(depth_track_count))
        self.configured_depth_track_heights = self._normalize_depth_track_heights(depth_track_heights)
        self.track_axis_gap: float | None = cfg.canvas.linear.track_axis_gap
        self.align_center: bool = cfg.canvas.linear.align_center
        self.keep_definition_left_aligned: bool = cfg.canvas.linear.keep_definition_left_aligned
        self.normalize_length: bool = cfg.canvas.linear.normalize_length
        self.legend_position = legend
        self.num_of_entries: int = num_of_entries
        self.plot_tracks_top_extent = 0.0
        self.plot_tracks_bottom_extent = 0.0

        self.calculate_dimensions()
        self.set_arrow_length()

    @property
    def profile(self) -> LinearRenderProfile:
        return self._profile

    @property
    def show_gc(self) -> bool:
        return self._profile.show_gc

    @property
    def show_skew(self) -> bool:
        return self._profile.show_skew

    @property
    def show_depth(self) -> bool:
        return self._profile.show_depth

    @property
    def strandedness(self) -> bool:
        return self._profile.strandedness

    @property
    def resolve_overlaps(self) -> bool:
        return self._profile.resolve_overlaps

    @property
    def track_layout(self) -> str:
        return self._profile.track_layout

    def _normalize_depth_track_heights(
        self,
        depth_track_heights: list[float | None] | tuple[float | None, ...] | None,
    ) -> list[float | None]:
        if not depth_track_heights:
            return []
        normalized: list[float | None] = []
        for raw_height in depth_track_heights:
            if raw_height is None:
                normalized.append(None)
                continue
            value = float(raw_height)
            normalized.append(value if value > 0 else None)
        return normalized

    def _depth_track_height_at(self, depth_index: int) -> float:
        if 0 <= depth_index < len(self.configured_depth_track_heights):
            configured_height = self.configured_depth_track_heights[depth_index]
            if configured_height is not None:
                return float(configured_height)
        return float(self.default_depth_height)

    def set_gc_height_and_gc_padding(self) -> None:
        """
        Sets linear plot track heights, anchor offsets, and legacy padding segments.
        """

        active_depth_track_count = self.depth_track_count if self.show_depth else 0
        self.depth_track_heights: list[float] = [
            self._depth_track_height_at(depth_index)
            for depth_index in range(active_depth_track_count)
        ]
        self.depth_height: float = self.depth_track_heights[0] if self.depth_track_heights else 0.0
        self.gc_height: float = self.default_gc_height if self.show_gc else 0.0
        self.skew_height: float = self.default_gc_height if self.show_skew else 0.0
        gc_content_mode = str(getattr(self._cfg.objects.gc_content, "mode", "deviation"))

        tracks: list[_LinearPlotTrack] = []
        for depth_index in range(active_depth_track_count):
            tracks.append(
                _LinearPlotTrack(
                    key="depth" if active_depth_track_count == 1 else f"depth_{depth_index + 1}",
                    top_extent=0.0,
                    bottom_extent=self.depth_track_heights[depth_index],
                    gap_after=self.configured_depth_padding,
                )
            )
        if self.show_gc:
            if gc_content_mode == "percent":
                gc_top_extent = 0.0
                gc_bottom_extent = self.gc_height
            else:
                gc_top_extent = 0.5 * self.gc_height
                gc_bottom_extent = 0.5 * self.gc_height
            tracks.append(
                _LinearPlotTrack(
                    key="gc_content",
                    top_extent=gc_top_extent,
                    bottom_extent=gc_bottom_extent,
                    gap_after=self.configured_track_spacing if self.show_skew else 0.0,
                )
            )
        if self.show_skew:
            tracks.append(
                _LinearPlotTrack(
                    key="gc_skew",
                    top_extent=0.5 * self.skew_height,
                    bottom_extent=0.5 * self.skew_height,
                )
            )

        track_offsets, padding_segments, stack_span = _calculate_linear_plot_track_layout(tracks)
        self.plot_track_offsets = track_offsets
        self.plot_tracks_height = stack_span
        self.plot_tracks_top_extent = 0.0
        self.plot_tracks_bottom_extent = stack_span
        self.plot_tracks_visual_bottom = _calculate_linear_plot_track_visual_bottom(
            tracks,
            track_offsets,
        )
        if active_depth_track_count <= 1:
            self.depth_track_offsets = [track_offsets.get("depth", 0.0)] if self.show_depth else []
        else:
            self.depth_track_offsets = [
                track_offsets.get(f"depth_{index + 1}", 0.0)
                for index in range(active_depth_track_count)
            ]
        self.depth_track_offset = self.depth_track_offsets[0] if self.depth_track_offsets else 0.0
        self.gc_content_track_offset = track_offsets.get("gc_content", 0.0)
        self.gc_skew_track_offset = track_offsets.get("gc_skew", 0.0)
        self.depth_padding: float = (
            sum(
                padding_segments.get(
                    "depth" if active_depth_track_count == 1 else f"depth_{index + 1}",
                    0.0,
                )
                for index in range(active_depth_track_count)
            )
            if active_depth_track_count > 0
            else 0.0
        )
        self.gc_padding: float = padding_segments.get("gc_content", 0.0)
        self.skew_padding: float = padding_segments.get("gc_skew", 0.0)

    def set_cds_height_and_cds_padding(self) -> None:
        """
        Sets the height and padding for the coding sequences (CDS) track based on configuration settings.
        This method adjusts the cds_height and cds_padding attributes.
        """

        if self.strandedness:
            self.cds_height: float = 1 * self.default_cds_height
            self.cds_padding: float = 0.6 * self.cds_height + 5
        else:
            self.cds_height: float = 0.5 * self.default_cds_height
            self.cds_padding: float = 0.6 * self.cds_height + 5

    def set_arrow_length(self) -> None:
        """
        Sets the length of the arrow used in the representation based on the longest genome.
        This method adjusts the arrow_length attribute.
        """

        self.arrow_length_param = getattr(self._cfg.canvas.linear.arrow_length_parameter, self.length_param)
        base_arrow_length = self.arrow_length_param * self.longest_genome
        baseline_cds_height = _ARROW_LENGTH_REFERENCE_CDS_HEIGHT[self.length_param]
        if baseline_cds_height <= 0.0:
            self.arrow_length = base_arrow_length
            return
        current_non_stranded_cds_height = 0.5 * float(self.default_cds_height)
        self.arrow_length = base_arrow_length * (current_non_stranded_cds_height / baseline_cds_height)

    def calculate_dimensions(self) -> None:
        """
        Calculates the dimensions for the linear canvas including the total width and height,
        considering all the elements and padding. This method updates total_width and total_height attributes.
        """

        self.set_gc_height_and_gc_padding()
        self.set_cds_height_and_cds_padding()
        self.add_margin: float | Literal[0] = (
            2 * self.cds_height if ((self.show_gc or self.show_depth) and not self.strandedness) else 0
        )
        # Keep the record axis width fixed to the configured figure width from the start.
        # Horizontal offsets reposition the plotted record; they do not shorten its scale.
        self.alignment_width: float = self.fig_width
        self.total_width = int(self.horizontal_offset + self.fig_width)
        self.total_height = int(
            2 * self.vertical_offset
            + (self.cds_height + self.plot_tracks_height)
            + (
                self.vertical_padding
                + self.comparison_height
                + self.vertical_padding
                + self.cds_height
                + self.plot_tracks_height
            )
            * (self.num_of_entries - 1)
        )

    def configure_plot_width(self, definition_reserve_width: float) -> None:
        """Resolve the legend-independent horizontal plot-space origin."""

        reserve = max(0.0, float(definition_reserve_width))
        self.alignment_width = float(self.fig_width)
        self.horizontal_offset = reserve + float(self.definition_gap)
        self.total_width = self.horizontal_offset + self.alignment_width

    def create_svg_canvas(self) -> Drawing:
        """
        Creates and returns an SVG canvas for the linear representation based on the configurator's settings.

        Returns:
        Drawing: An SVG Drawing object representing the linear canvas.
        """

        return Drawing(
            filename=self.output_prefix + ".svg",
            size=(str(self.total_width) + "px", str(self.total_height) + "px"),
            viewBox=("0 0 " + str(self.total_width) + " " + str(self.total_height)),
            debug=False,
        )


__all__ = ["LinearCanvasConfigurator"]
