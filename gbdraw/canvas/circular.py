#!/usr/bin/env python
# coding: utf-8

from typing import Literal

from svgwrite import Drawing  # type: ignore[reportMissingImports]

from ..config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ..core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]


def resolve_circular_side_legend_geometry(
    *,
    canvas_height: float,
    legend_width: float,
    color_rect_size: float,
) -> tuple[float, float, float]:
    """Return (inner_gap, edge_margin, reserved_width) for circular side legends."""
    side_inner_gap_px = 0.05 * float(legend_width)
    side_edge_margin_px = max(0.0, (0.05 * float(canvas_height)) - (0.5 * float(color_rect_size)))
    side_reserved_width_px = float(legend_width) + side_inner_gap_px + side_edge_margin_px
    return side_inner_gap_px, side_edge_margin_px, side_reserved_width_px


class CircularCanvasConfigurator:
    """
    Configures the settings for a circular canvas used for genomic data visualization.

    Attributes:
    output_prefix (str): Prefix for the output file.
    total_width (int): Total width of the canvas.
    total_height (int): Total height of the canvas.
    radius (float): Radius of the circular canvas.
    track_ratio (float): Ratio to determine track width.
    show_gc (bool): Flag to display GC content.
    show_skew (bool): Flag to display GC skew.
    strandedness (bool): Flag to display strandedness.

    Methods:
    calculate_dimensions(): Calculates dimensions for the canvas.
    create_svg_canvas(): Creates and returns an SVG canvas for drawing.
    get_track_ids(): Determines the track IDs for visualization.
    """

    def __init__(
        self,
        output_prefix: str,
        config_dict: dict,
        legend: str,
        gb_record,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the circular canvas configurator with given settings.

        Args:
        output_prefix (str): Prefix for the output file.
        config_dict (dict): Configuration dictionary with canvas settings.
        show_gc (bool, optional): Flag to display GC content. Defaults to True.
        strandedness (bool, optional): Flag to display strandedness. Defaults to True.
        show_skew (bool, optional): Flag to display GC skew. Defaults to True.
        """

        self.output_prefix: str = output_prefix
        self.config_dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg

        raw_show_labels = cfg.canvas.show_labels
        if isinstance(raw_show_labels, str):
            self.show_labels = raw_show_labels != "none"
        else:
            self.show_labels = bool(raw_show_labels)

        if self.show_labels:
            label_setting = "with_labels"
        else:
            label_setting = "without_labels"
        self.default_width: int = (
            cfg.canvas.circular.width.with_labels
            if label_setting == "with_labels"
            else cfg.canvas.circular.width.without_labels
        )
        self.default_height: int = cfg.canvas.circular.height
        self.radius: float = cfg.canvas.circular.radius
        self.track_ratio: float = cfg.canvas.circular.track_ratio
        self.show_gc: bool = cfg.canvas.show_gc
        self.show_skew: bool = cfg.canvas.show_skew
        self.strandedness: bool = cfg.canvas.strandedness
        self.dpi: int = cfg.canvas.dpi
        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        # track_width (legacy) used to be config-driven; now derived from radius/track_ratio factors.
        self.track_ratio_factors = cfg.canvas.circular.track_ratio_factors[self.length_param]
        self.legend_position: str = legend

        self.calculate_dimensions()
        self.get_track_ids()

    def calculate_dimensions(self) -> None:
        """
        Calculates the dimensions and offsets for the circular canvas based on the configuration.
        """

        self.total_height = self.default_height
        self.offset_y: float = self.total_height * 0.5
        if self.legend_position == "left":
            self.total_width = self.default_width * 1.2
            self.offset_x: float = self.default_width * 0.6
        elif self.legend_position == "right":
            self.total_width = self.default_width * 1.2
            self.offset_x: float = self.default_width * 0.5
        else:
            self.total_width = self.default_width
            self.offset_x: float = self.default_width * 0.5

        # Create linear canvas

    def recalculate_canvas_dimensions(self, legend_config):
        legend_local_top = -0.5 * float(legend_config.color_rect_size)
        legend_edge_margin = 16.0
        legend_content_gap = 12.0
        top_bottom_reserved_height = (
            float(legend_config.legend_height) + legend_edge_margin + legend_content_gap
        )
        side_inner_gap, side_edge_margin, side_reserved_width = resolve_circular_side_legend_geometry(
            canvas_height=float(self.total_height),
            legend_width=float(legend_config.legend_width),
            color_rect_size=float(legend_config.color_rect_size),
        )

        if self.legend_position == "right":
            self.total_width = self.default_width + side_reserved_width
            self.legend_offset_x = self.default_width + side_inner_gap
            self.legend_offset_y = (self.total_height - legend_config.legend_height) / 2
        elif self.legend_position == "left":
            self.total_width = self.default_width + side_reserved_width
            self.legend_offset_x = side_edge_margin
            self.offset_x: float = (self.default_width * 0.5) + side_reserved_width
            self.legend_offset_y = (self.total_height - legend_config.legend_height) / 2
        elif self.legend_position == "top":
            self.total_height = self.default_height + top_bottom_reserved_height
            self.offset_y = (self.default_height * 0.5) + top_bottom_reserved_height
            self.legend_offset_x = (self.total_width - legend_config.legend_width) / 2
            self.legend_offset_y = legend_edge_margin - legend_local_top
        elif self.legend_position == "bottom":
            self.total_height = self.default_height + top_bottom_reserved_height
            self.offset_y = self.default_height * 0.5
            self.legend_offset_x = (self.total_width - legend_config.legend_width) / 2
            self.legend_offset_y = self.default_height + legend_content_gap - legend_local_top
        elif self.legend_position == "upper_left":
            self.legend_offset_x: float = 0.025 * self.total_width
            self.legend_offset_y: float = 0.05 * self.total_height
        elif self.legend_position == "upper_right":
            self.legend_offset_x: float = 0.85 * self.total_width
            self.legend_offset_y: float = 0.05 * self.total_height
        elif self.legend_position == "lower_left":
            self.legend_offset_x: float = 0.025 * self.total_width
            self.legend_offset_y: float = 0.78 * self.total_height
        elif self.legend_position == "lower_right":
            self.legend_offset_x: float = 0.875 * self.total_width
            self.legend_offset_y: float = 0.75 * self.total_height
        elif self.legend_position == "none":
            self.legend_offset_x: float = 0
            self.legend_offset_y: float = 0
        else:
            self.legend_offset_x: float = 0
            self.legend_offset_y: float = 0

    def create_svg_canvas(self) -> Drawing:
        """
        Creates and returns an SVG canvas based on the configurator's settings.

        Returns:
        Drawing: An SVG Drawing object.
        """

        return Drawing(
            filename=self.output_prefix + ".svg",
            size=(str(self.total_width) + "px", str(self.total_height) + "px"),
            viewBox=("0 0 " + str(self.total_width) + " " + str(self.total_height)),
            debug=False,
        )

    def get_track_ids(self) -> None:
        """
        Determines and assigns track IDs for the visualization based on the configurator's settings.
        """

        self.track_ids: dict = {}
        gc_track_id: Literal[2] | None = 2 if self.show_gc or not self.show_skew else None
        skew_track_id: Literal[3, 2] | None = (3 if self.show_gc else 2) if self.show_skew else None

        if gc_track_id is not None:
            self.track_ids["gc_track"] = gc_track_id
        if skew_track_id is not None:
            self.track_ids["skew_track"] = skew_track_id

__all__ = ["CircularCanvasConfigurator", "resolve_circular_side_legend_geometry"]


