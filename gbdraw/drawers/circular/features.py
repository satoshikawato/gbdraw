#!/usr/bin/env python
# coding: utf-8

import math
from typing import Optional, Union, List, Dict

from svgwrite.container import Group
from svgwrite.path import Path

from ...feature_objects import FeatureObject
from ...layout.common import calculate_cds_ratio
from ...configurators import FeatureDrawingConfigurator
from ...svg.circular_features import (
    generate_circular_arrowhead_path,
    generate_circular_intron_path,
    generate_circular_rectangle_path,
)


class FeatureDrawer:
    """
    Draws genomic features (blocks + intron lines) on a circular canvas.
    """

    def __init__(self, feature_config: FeatureDrawingConfigurator) -> None:
        self.default_feature_color: str = feature_config.block_fill_color
        self.default_stroke_color: str = feature_config.block_stroke_color
        self.default_stroke_width: float = feature_config.block_stroke_width
        self.intron_stroke_color: str = feature_config.line_stroke_color
        self.intron_stroke_width: float = feature_config.line_stroke_width

    def draw_path(
        self,
        path_data: str,
        group: Group,
        fill_color: str,
        stroke_color_specified: Optional[str] = None,
        stroke_width_specified: Optional[float] = None,
    ) -> None:
        stroke_color: str = (
            stroke_color_specified if stroke_color_specified is not None else self.default_stroke_color
        )
        stroke_width: float = (
            stroke_width_specified if stroke_width_specified is not None else self.default_stroke_width
        )
        path: Path = Path(
            d=path_data,
            fill=fill_color,
            stroke=stroke_color,
            stroke_width=stroke_width,
            stroke_linejoin="round",
            stroke_linecap="round",
            stroke_miterlimit=4,
        )
        group.add(path)

    def draw(
        self,
        feature_object: FeatureObject,
        group: Group,
        total_length: int,
        radius: float,
        track_ratio: float,
        track_ratio_factor,
        track_type: str,
        strandedness: bool,
        length_param,
    ) -> Group:
        """
        Draw a feature on the circular canvas.
        
        Args:
            feature_object: The feature to draw
            group: SVG group to add the feature to
            total_length: Total genome length
            radius: Base radius of the canvas
            track_ratio: Track ratio from config
            track_ratio_factor: Track ratio factor for length parameter
            track_type: "tuckin", "middle", or "spreadout"
            strandedness: Whether strands are separated
            length_param: Length parameter ("short" or "long")
        
        Returns:
            Updated SVG group with the feature added
        """
        cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
        
        # Get the track_id from the feature for overlap resolution
        track_id = getattr(feature_object, 'feature_track_id', 0)
        
        gene_paths = FeaturePathGenerator(
            radius, total_length, track_ratio, cds_ratio, offset, track_type, strandedness, track_id
        ).generate_circular_gene_path(feature_object)
        
        for gene_path in gene_paths:
            if not gene_path or not gene_path[0]:
                continue
            path_type: str = gene_path[0]
            path_data = gene_path[1]
            if path_type == "block":
                self.draw_path(path_data, group, fill_color=feature_object.color)
            elif path_type == "line":
                self.draw_path(
                    path_data,
                    group,
                    fill_color="none",
                    stroke_color_specified=self.intron_stroke_color,
                    stroke_width_specified=self.intron_stroke_width,
                )
            elif path_type == "label":
                for element in path_data:
                    group.add(element)
        return group


class FeaturePathGenerator:
    """
    Generates SVG path data for genomic features on a circular canvas.
    """

    def __init__(
        self,
        radius: float,
        total_length: int,
        track_ratio: float,
        cds_ratio: float,
        offset: float,
        track_type: str,
        strandedness: bool,
        track_id: int = 0,
    ) -> None:
        """
        Initialize the path generator.
        
        Args:
            radius: Base radius of the circular canvas
            total_length: Total genome length
            track_ratio: Track ratio from config
            cds_ratio: Calculated CDS ratio
            offset: Base offset
            track_type: "tuckin", "middle", or "spreadout"
            strandedness: Whether strands are separated
            track_id: Track number for overlap resolution (0 = default track)
        """
        self.radius: float = radius
        self.total_length: int = total_length
        self.track_ratio: float = track_ratio
        self.cds_ratio: float = cds_ratio
        self.offset: float = offset
        self.track_type = track_type
        self.strandedness = strandedness
        self.track_id = track_id
        self.set_arrow_length()

    def set_arrow_length(self) -> None:
        MIN_ARROW_LENGTH = 30
        MAX_ARROW_LENGTH = 700
        PARAM_A = 3
        PARAM_B = 5
        self.arrow_length: float = (
            MIN_ARROW_LENGTH
            + (MAX_ARROW_LENGTH - MIN_ARROW_LENGTH)
            * 1
            / (1 + math.exp(-PARAM_A * (math.log10(self.total_length) - PARAM_B)))
        )

    def generate_circular_gene_path(self, feature_object: FeatureObject):
        """
        Generate SVG path data for a feature.
        
        Args:
            feature_object: The feature to generate paths for
        
        Returns:
            List of [path_type, path_data] pairs
        """
        coords = feature_object.location
        coordinates_paths: List[List[str]] = []
        for coord in coords:
            coord_dict: Dict[str, Union[str, int]] = {
                "coord_type": str(coord.kind),
                "coord_strand": coord.strand,
                "coord_start": coord.start,
                "coord_end": coord.end,
            }
            coord_type: str = str(coord_dict["coord_type"])
            if coord_type == "line":
                coord_path: List[str] = generate_circular_intron_path(
                    self.radius,
                    coord_dict,
                    self.total_length,
                    self.track_ratio,
                    self.cds_ratio,
                    self.offset,
                    self.track_type,
                    self.strandedness,
                    self.track_id,
                )
            elif coord_type == "block":
                if coord.is_last and feature_object.is_directional is True:
                    coord_path = generate_circular_arrowhead_path(
                        self.radius,
                        coord_dict,
                        self.total_length,
                        self.arrow_length,
                        self.track_ratio,
                        self.cds_ratio,
                        self.offset,
                        self.track_type,
                        self.strandedness,
                        self.track_id,
                    )
                else:
                    coord_path = generate_circular_rectangle_path(
                        self.radius,
                        coord_dict,
                        self.total_length,
                        self.track_ratio,
                        self.cds_ratio,
                        self.offset,
                        self.track_type,
                        self.strandedness,
                        self.track_id,
                    )
            else:
                coord_path = []
            coordinates_paths.append(coord_path)
        return coordinates_paths


__all__ = ["FeatureDrawer", "FeaturePathGenerator"]
