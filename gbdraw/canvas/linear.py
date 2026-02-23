#!/usr/bin/env python
# coding: utf-8

from typing import Literal

from svgwrite import Drawing

from ..config.models import GbdrawConfig
from ..core.sequence import determine_length_parameter


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
    set_gc_height_and_gc_padding(): Sets GC height and padding.
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
        self.default_gc_height: float = cfg.canvas.linear.default_gc_height
        self.dpi: int = cfg.canvas.dpi
        self.show_gc: bool = cfg.canvas.show_gc
        self.show_skew: bool = cfg.canvas.show_skew
        self.strandedness: bool = cfg.canvas.strandedness
        self.resolve_overlaps: bool = cfg.canvas.resolve_overlaps
        self.track_layout: str = cfg.canvas.linear.track_layout
        self.track_axis_gap: float | None = cfg.canvas.linear.track_axis_gap
        self.ruler_on_axis: bool = cfg.canvas.linear.ruler_on_axis
        self.align_center: bool = cfg.canvas.linear.align_center
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
        Sets the height and padding for the GC content track based on configuration settings.
        This method adjusts the gc_height and gc_padding attributes.
        """

        if self.show_gc:
            self.gc_height: float = self.default_gc_height
            if self.show_skew:
                self.gc_padding: float = self.gc_height
            else:
                self.gc_padding: float = self.gc_height
        else:
            self.gc_height: float = 0
            self.gc_padding: float = 0
        if self.show_skew:
            self.skew_height: float = self.default_gc_height
            if self.show_gc:
                self.skew_padding: float = self.skew_height
            else:
                self.skew_padding: float = self.skew_height
        else:
            self.skew_height: float = 0
            self.skew_padding: float = 0

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
        self.arrow_length: float = self.arrow_length_param * self.longest_genome

    def calculate_dimensions(self) -> None:
        """
        Calculates the dimensions for the linear canvas including the total width and height,
        considering all the elements and padding. This method updates total_width and total_height attributes.
        """

        self.set_gc_height_and_gc_padding()
        self.set_cds_height_and_cds_padding()
        self.add_margin: float | Literal[0] = 2 * self.cds_height if (self.show_gc and not self.strandedness) else 0
        self.alignment_width: float = self.fig_width - self.horizontal_offset
        self.total_width = int(self.fig_width + 2 * self.canvas_padding)
        self.total_height = int(
            2 * self.vertical_offset
            + (self.cds_height + self.gc_padding)
            + (
                self.vertical_padding
                + self.comparison_height
                + self.vertical_padding
                + self.cds_height
                + self.gc_padding
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


