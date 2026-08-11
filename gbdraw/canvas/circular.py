#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import Literal

from svgwrite import Drawing  # type: ignore[reportMissingImports]

from ..config.models import CircularRenderProfile  # type: ignore[reportMissingImports]
from ..core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]


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
        profile: CircularRenderProfile,
        legend: str,
        gb_record,
    ) -> None:
        """
        Initializes the circular canvas configurator with given settings.

        Args:
        output_prefix (str): Prefix for the output file.
        profile (CircularRenderProfile): Resolved circular render settings.
        show_gc (bool, optional): Flag to display GC content. Defaults to True.
        strandedness (bool, optional): Flag to display strandedness. Defaults to True.
        show_skew (bool, optional): Flag to display GC skew. Defaults to True.
        """
        cfg = profile.config
        self.output_prefix: str = output_prefix
        self._profile = profile
        self._cfg = cfg

        if profile.labels_enabled:
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
        self.dpi: int = cfg.canvas.dpi
        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        # track_width (legacy) used to be config-driven; now derived from radius/track_ratio factors.
        self.track_ratio_factors = cfg.canvas.circular.track_ratio_factors[self.length_param]
        self.legend_position: str = legend

        self.calculate_dimensions()
        self.get_track_ids()

    @property
    def profile(self) -> CircularRenderProfile:
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

    def calculate_dimensions(self) -> None:
        """
        Calculates the dimensions and offsets for the circular canvas based on the configuration.
        """

        self.total_width = self.default_width
        self.total_height = self.default_height
        self.offset_x: float = self.default_width * 0.5
        self.offset_y: float = self.default_height * 0.5
        # Decorations are positioned only by the final composition plan.  Keep
        # these attributes as neutral compatibility state for older callers.
        self.legend_offset_x: float = 0.0
        self.legend_offset_y: float = 0.0

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
        depth_track_id: Literal[2] | None = 2 if self.show_depth else None
        if self.show_depth:
            gc_track_id: Literal[3] | None = 3 if self.show_gc else None
            skew_track_id: Literal[4, 3] | None = (4 if self.show_gc else 3) if self.show_skew else None
        else:
            gc_track_id: Literal[2] | None = 2 if self.show_gc or not self.show_skew else None
            skew_track_id: Literal[3, 2] | None = (3 if self.show_gc else 2) if self.show_skew else None

        if depth_track_id is not None:
            self.track_ids["depth_track"] = depth_track_id
        if gc_track_id is not None:
            self.track_ids["gc_track"] = gc_track_id
        if skew_track_id is not None:
            self.track_ids["skew_track"] = skew_track_id

__all__ = ["CircularCanvasConfigurator"]
