#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union, List, Dict

from svgwrite.container import Group
from svgwrite.path import Path

from ....features.objects import FeatureObject
from ....features.ids import compute_feature_object_hash
from ....layout.common import calculate_cds_ratio
from ....configurators import FeatureDrawingConfigurator
from ....svg.circular_features import (
    generate_circular_arrowhead_path,
    generate_circular_arrowhead_path_with_radii,
    generate_circular_intron_path,
    generate_circular_intron_path_with_radii,
    generate_circular_rectangle_path,
    generate_circular_rectangle_path_with_radii,
)
from ....svg.arrows import calculate_circular_arrow_length

if TYPE_CHECKING:
    from ....diagrams.circular.radial_layout import CircularFeatureLayout


class FeatureDrawer:
    """
    Draws genomic features (blocks + intron lines) on a circular canvas.
    """

    def __init__(
        self,
        feature_config: FeatureDrawingConfigurator,
        feature_layout: CircularFeatureLayout | None = None,
    ) -> None:
        self.default_feature_color: str = feature_config.block_fill_color
        self.default_stroke_color: str = feature_config.block_stroke_color
        self.default_stroke_width: float = feature_config.block_stroke_width
        self.intron_stroke_color: str = feature_config.line_stroke_color
        self.intron_stroke_width: float = feature_config.line_stroke_width
        self.feature_layout = feature_layout

    @staticmethod
    def get_feature_data_id(feature_object: FeatureObject) -> Optional[str]:
        return compute_feature_object_hash(feature_object)

    def draw_path(
        self,
        path_data: str,
        group: Group,
        fill_color: str,
        stroke_color_specified: Optional[str] = None,
        stroke_width_specified: Optional[float] = None,
        feature_data_id: Optional[str] = None,
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
            debug=False,
        )
        if feature_data_id:
            path.attribs["data-gbdraw-feature-id"] = feature_data_id
            path.attribs["id"] = feature_data_id
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
            radius,
            total_length,
            track_ratio,
            cds_ratio,
            offset,
            track_type,
            strandedness,
            track_id,
            feature_layout=self.feature_layout,
        ).generate_circular_gene_path(feature_object)

        # Get feature identifier for instant preview support
        feature_data_id = self.get_feature_data_id(feature_object)

        # Keep intron lines behind blocks for the same feature by drawing lines first.
        ordered_types: tuple[str, ...] = ("line", "block", "label")
        for expected_type in ordered_types:
            for gene_path in gene_paths:
                if not gene_path or not gene_path[0]:
                    continue
                path_type: str = gene_path[0]
                if path_type != expected_type:
                    continue
                path_data = gene_path[1]
                if path_type == "block":
                    self.draw_path(path_data, group, fill_color=feature_object.color, feature_data_id=feature_data_id)
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
        feature_layout: CircularFeatureLayout | None = None,
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
        self.feature_layout = feature_layout
        self.set_arrow_length()

    def set_arrow_length(self) -> None:
        self.arrow_length = calculate_circular_arrow_length(self.total_length)

    def _coalesce_origin_spanning_block(self, feature_object: FeatureObject) -> Optional[Dict[str, Union[str, int]]]:
        """
        Collapse two-block origin-spanning features into one block.

        This avoids visual splitting at the start/end boundary in circular mode
        (e.g. D-loop represented as join(end..genome,1..start)).
        """
        block_coords = [coord for coord in feature_object.location if coord.kind == "block"]
        if len(block_coords) != 2:
            return None

        starts = [int(coord.start) for coord in block_coords]
        ends = [int(coord.end) for coord in block_coords]

        if min(starts) > 1 or max(ends) < self.total_length:
            return None

        left_blocks = [coord for coord in block_coords if int(coord.start) <= 1]
        right_blocks = [coord for coord in block_coords if int(coord.end) >= self.total_length]

        if len(left_blocks) != 1 or len(right_blocks) != 1 or left_blocks[0] is right_blocks[0]:
            return None

        merged_start = max(starts)
        merged_end = min(ends)
        if merged_start <= merged_end:
            return None

        return {
            "coord_type": "block",
            "coord_strand": block_coords[0].strand,
            "coord_start": merged_start,
            "coord_end": merged_end,
        }

    def generate_circular_gene_path(self, feature_object: FeatureObject):
        """
        Generate SVG path data for a feature.

        Args:
            feature_object: The feature to generate paths for

        Returns:
            List of [path_type, path_data] pairs
        """
        lane = None
        if self.feature_layout is not None:
            lane = self.feature_layout.lane_for_track_id(int(getattr(feature_object, "feature_track_id", 0)))
        merged_coord = self._coalesce_origin_spanning_block(feature_object)
        if merged_coord is not None:
            merged_strand = str(merged_coord["coord_strand"])
            merged_coord_for_draw = merged_coord
            if merged_strand not in {"positive", "negative"}:
                merged_coord_for_draw = dict(merged_coord)
                merged_coord_for_draw["coord_strand"] = "positive"
            if feature_object.is_directional and merged_strand in {"positive", "negative"}:
                if lane is None:
                    merged_path = generate_circular_arrowhead_path(
                        self.radius,
                        merged_coord_for_draw,
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
                    merged_path = generate_circular_arrowhead_path_with_radii(
                        merged_coord_for_draw,
                        self.total_length,
                        self.arrow_length,
                        lane.inner_px,
                        lane.center_px,
                        lane.outer_px,
                    )
            else:
                # Fallback to rectangle for undefined strand to avoid arrow path errors.
                if lane is None:
                    merged_path = generate_circular_rectangle_path(
                        self.radius,
                        merged_coord_for_draw,
                        self.total_length,
                        self.track_ratio,
                        self.cds_ratio,
                        self.offset,
                        self.track_type,
                        self.strandedness,
                        self.track_id,
                    )
                else:
                    merged_path = generate_circular_rectangle_path_with_radii(
                        merged_coord_for_draw,
                        self.total_length,
                        lane.inner_px,
                        lane.center_px,
                        lane.outer_px,
                    )
            return [merged_path]

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
                if lane is None:
                    coord_path = generate_circular_intron_path(
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
                    coord_path = generate_circular_intron_path_with_radii(
                        coord_dict,
                        self.total_length,
                        lane.center_px,
                    )
            elif coord_type == "block":
                if coord.is_last and feature_object.is_directional is True:
                    if lane is None:
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
                        coord_path = generate_circular_arrowhead_path_with_radii(
                            coord_dict,
                            self.total_length,
                            self.arrow_length,
                            lane.inner_px,
                            lane.center_px,
                            lane.outer_px,
                        )
                else:
                    if lane is None:
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
                        coord_path = generate_circular_rectangle_path_with_radii(
                            coord_dict,
                            self.total_length,
                            lane.inner_px,
                            lane.center_px,
                            lane.outer_px,
                        )
            else:
                coord_path = []
            coordinates_paths.append(coord_path)
        return coordinates_paths


__all__ = ["FeatureDrawer", "FeaturePathGenerator"]
