#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

from svgwrite import Drawing

from ..config.models import GbdrawConfig
from ..config.toml import load_config_toml
from ..core.sequence import determine_length_parameter


@lru_cache(maxsize=1)
def _get_default_linear_non_stranded_cds_heights() -> dict[str, float]:
    """Return baseline CDS heights from the unmodified default config."""
    default_config = load_config_toml("gbdraw.data", "config.toml")
    default_cfg = GbdrawConfig.from_dict(default_config)
    return {
        "short": 0.5 * float(default_cfg.canvas.linear.default_cds_height.short),
        "long": 0.5 * float(default_cfg.canvas.linear.default_cds_height.long),
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


class LinearCanvasConfigurator:
    """
    Configures the settings for a linear canvas used for genomic data visualization.

    Attributes:
    output_prefix (str): Prefix for the output file.
    fig_width (int): Width of the figure.
    vertical_offset (float): Vertical offset for alignment.
    horizontal_offset (float): Horizontal offset for alignment.
    vertical_padding (float): Vertical padding between elements.
    comparison_height (float): Height for comparison tracks.
    canvas_padding (float): Padding around the canvas.
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
        config_dict: dict,
        legend: str,
        output_prefix="out",
        cfg: GbdrawConfig | None = None,
        has_comparisons: bool = False,
    ):
        """
        Initializes the linear canvas configurator with given settings.

        Args:
        num_of_entries (int): Number of entries to visualize.
        longest_genome (int): Length of the longest genome in the dataset.
        config_dict (dict): Configuration dictionary with canvas settings.
        output_prefix (str, optional): Prefix for the output file. Defaults to 'out'.
        show_gc (bool, optional): Flag to display GC content. Defaults to False.
        strandedness (bool, optional): Flag to display strandedness. Defaults to False.
        align_center (bool, optional): Flag to align content to the center. Defaults to False.
        has_comparisons (bool, optional): Flag indicating if BLAST comparisons are present. Defaults to False.
        """

        self.output_prefix: str = output_prefix
        self.config_dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.longest_genome: int = longest_genome
        self.fig_width: int = cfg.canvas.linear.width
        self.original_vertical_offset: float = cfg.canvas.linear.vertical_offset
        self.vertical_offset: float = cfg.canvas.linear.vertical_offset
        self.horizontal_offset: float = cfg.canvas.linear.horizontal_offset
        self.original_horizontal_offset = self.horizontal_offset
        self.vertical_padding: float = cfg.canvas.linear.vertical_padding
        self.has_comparisons: bool = has_comparisons
        # Store the configured comparison_height for BLAST ribbons
        self.configured_comparison_height: float = cfg.canvas.linear.comparison_height
        # For initial dimension calculation, use comparison_height only when BLAST is present
        # The actual inter-record spacing will be recalculated in assemble_linear_diagram
        # based on feature track heights
        self.comparison_height: float = cfg.canvas.linear.comparison_height if has_comparisons else 0
        self.canvas_padding: float = cfg.canvas.linear.canvas_padding
        self.length_threshold = cfg.labels.length_threshold.linear
        self.length_param = determine_length_parameter(self.longest_genome, self.length_threshold)
        self.default_cds_height: float = getattr(cfg.canvas.linear.default_cds_height, self.length_param)
        self.baseline_non_stranded_cds_heights = _get_default_linear_non_stranded_cds_heights()
        self.default_gc_height: float = cfg.canvas.linear.default_gc_height
        self.default_depth_height: float = cfg.canvas.linear.depth_height
        self.configured_depth_padding: float = cfg.canvas.linear.depth_padding
        self.dpi: int = cfg.canvas.dpi
        self.show_gc: bool = cfg.canvas.show_gc
        self.show_skew: bool = cfg.canvas.show_skew
        self.show_depth: bool = cfg.canvas.show_depth
        self.strandedness: bool = cfg.canvas.strandedness
        self.resolve_overlaps: bool = cfg.canvas.resolve_overlaps
        self.track_layout: str = cfg.canvas.linear.track_layout
        self.track_axis_gap: float | None = cfg.canvas.linear.track_axis_gap
        self.ruler_on_axis: bool = cfg.canvas.linear.ruler_on_axis
        self.align_center: bool = cfg.canvas.linear.align_center
        self.keep_definition_left_aligned: bool = cfg.canvas.linear.keep_definition_left_aligned
        self.normalize_length: bool = cfg.canvas.linear.normalize_length
        _label_setting = cfg.canvas.show_labels
        if isinstance(_label_setting, str):
            self.show_labels = _label_setting in ["all", "first"]
        else:
            self.show_labels = bool(_label_setting)
        self.legend_position = legend
        self.num_of_entries: int = num_of_entries
        self.total_height = 0

        self.calculate_dimensions()
        self.set_arrow_length()

    def set_gc_height_and_gc_padding(self) -> None:
        """
        Sets linear plot track heights, anchor offsets, and legacy padding segments.
        """

        self.depth_height: float = self.default_depth_height if self.show_depth else 0.0
        self.gc_height: float = self.default_gc_height if self.show_gc else 0.0
        self.skew_height: float = self.default_gc_height if self.show_skew else 0.0

        tracks: list[_LinearPlotTrack] = []
        if self.show_depth:
            tracks.append(
                _LinearPlotTrack(
                    key="depth",
                    top_extent=0.0,
                    bottom_extent=self.depth_height,
                    gap_after=self.configured_depth_padding,
                )
            )
        if self.show_gc:
            tracks.append(
                _LinearPlotTrack(
                    key="gc_content",
                    top_extent=0.5 * self.gc_height,
                    bottom_extent=0.5 * self.gc_height,
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
        self.depth_track_offset = track_offsets.get("depth", 0.0)
        self.gc_content_track_offset = track_offsets.get("gc_content", 0.0)
        self.gc_skew_track_offset = track_offsets.get("gc_skew", 0.0)
        self.depth_padding: float = padding_segments.get("depth", 0.0)
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
        baseline_cds_height = float(self.baseline_non_stranded_cds_heights.get(self.length_param, 0.0))
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
        self.total_width = int(self.fig_width + 2 * self.canvas_padding)
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

    def recalculate_canvas_dimensions(self, legend_group, max_definition_width):
        """
        Calculates final canvas dimensions and legend offsets, ensuring the legend fits within the canvas.
        """
        # Keep the alignment width fixed to the configured figure width so the longest record fits.
        # Horizontal offsets and legend widths expand the total canvas width around that baseline.
        self.alignment_width = self.fig_width

        def calculate_optimal_legend_y():
            genome_area_top = self.vertical_offset
            genome_area_bottom = self.total_height - self.vertical_offset - self.vertical_padding
            genome_area_center_y = genome_area_top + (genome_area_bottom - genome_area_top) / 2
            legend_y = genome_area_center_y - (legend_group.legend_height / 2)

            if legend_y < 0 or (legend_y + legend_group.legend_height) > self.total_height:
                legend_y = max((self.total_height - legend_group.legend_height) / 2, self.vertical_padding)

            return legend_y

        padding = self.canvas_padding
        legend_width = legend_group.legend_width

        if self.legend_position == "right":
            self.horizontal_offset = 2 * padding + max_definition_width
            self.total_width = (
                self.horizontal_offset
                + self.alignment_width
                + 2 * padding
                + legend_width
                + 1 * padding
            )
            self.legend_offset_x = self.horizontal_offset + self.alignment_width + 2 * padding
            self.legend_offset_y = calculate_optimal_legend_y()

        elif self.legend_position == "left":
            self.horizontal_offset = 1 * padding + legend_width + 2 * padding + max_definition_width
            self.total_width = self.horizontal_offset + self.alignment_width + 2 * padding
            self.legend_offset_x = padding
            self.legend_offset_y = calculate_optimal_legend_y()

        elif self.legend_position in ["top", "bottom"]:
            self.horizontal_offset = 2 * padding + max_definition_width
            self.total_width = self.horizontal_offset + self.alignment_width + 2 * padding
            self.legend_offset_x = (self.total_width - legend_group.legend_width) / 2
            if self.legend_position == "top":
                self.legend_offset_y = self.original_vertical_offset + 2 * self.vertical_padding
            elif self.legend_position == "bottom":
                self.legend_offset_y = self.total_height - self.original_vertical_offset - legend_group.legend_height

        else:
            self.horizontal_offset = 2 * padding + max_definition_width
            self.total_width = self.horizontal_offset + self.alignment_width + 2 * padding
            self.legend_offset_x = 0
            self.legend_offset_y = 0

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


