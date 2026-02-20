#!/usr/bin/env python
# coding: utf-8

import hashlib
from typing import Optional

from svgwrite.container import Group
from svgwrite.path import Path

from ....configurators import FeatureDrawingConfigurator
from ....svg.linear_features import (
    create_arrowhead_path_linear,
    create_intron_path_linear,
    create_rectangle_path_linear,
)


class FeatureDrawer:
    """
    Draws genomic features on a linear canvas.
    """

    def __init__(self, feature_config: FeatureDrawingConfigurator) -> None:
        self.default_feature_color: str = feature_config.block_fill_color
        self.default_stroke_color: str = feature_config.block_stroke_color
        self.default_stroke_width: float = feature_config.block_stroke_width
        self.intron_stroke_color: str = feature_config.line_stroke_color
        self.intron_stroke_width: float = feature_config.line_stroke_width

    @staticmethod
    def get_feature_data_id(feature_object) -> Optional[str]:
        """
        Generate a stable identifier for SVG id attribute.
        Hashes feature properties (record id, type, start, end, strand) to create a unique ID.
        Must use original coordinates (not processed location) to match colors.py compute_feature_hash().
        """
        # Use original coordinates (SimpleLocation from BioPython) for consistent hashing
        coords = getattr(feature_object, 'coordinates', None)
        if not coords:
            return None
        coord = coords[0]  # First SimpleLocation
        feat_type = getattr(feature_object, 'feature_type', '') or ''
        record_id = getattr(feature_object, 'record_id', None)
        # Use integer strand directly from SimpleLocation
        if record_id is not None:
            key = f"{record_id}:{feat_type}:{int(coord.start)}:{int(coord.end)}:{coord.strand}"
        else:
            key = f"{feat_type}:{int(coord.start)}:{int(coord.end)}:{coord.strand}"
        return "f" + hashlib.md5(key.encode()).hexdigest()[:8]

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
        path = Path(
            d=path_data,
            fill=fill_color,
            stroke=stroke_color,
            stroke_width=stroke_width,
            stroke_linejoin="round",
            stroke_linecap="round",
        )
        if feature_data_id:
            # Use feature's internal ID directly (e.g., "gene_000000001")
            path.attribs['id'] = feature_data_id
        group.add(path)

    def draw(
        self,
        feature_object,
        group: Group,
        genome_length: int,
        cds_height: float,
        alignment_width: float,
        normalization_factor: float,
        feature_strand: str,
        separate_strands: bool,
        arrow_length: float,
        track_layout: str = "middle",
        track_axis_gap: float | None = None,
    ) -> Group:
        path_generator = FeaturePathGenerator(
            genome_length=genome_length,
            alignment_width=alignment_width,
            cds_height=cds_height,
            genome_size_normalization_factor=normalization_factor,
            feature_strand=feature_strand,
            separate_strands=separate_strands,
            arrow_length=arrow_length,
            track_layout=track_layout,
            track_axis_gap=track_axis_gap,
        )

        gene_paths = path_generator.generate_linear_gene_path(feature_object)

        # Get feature identifier for instant preview support
        feature_data_id = self.get_feature_data_id(feature_object)

        for gene_path in gene_paths:
            path_type, path_data = gene_path[0], gene_path[1]
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
        return group


class FeaturePathGenerator:
    """
    Generates SVG path data for genomic features on a linear canvas.
    """

    def __init__(
        self,
        genome_length: int,
        alignment_width: float,
        cds_height: float,
        genome_size_normalization_factor: float,
        feature_strand: str,
        separate_strands: bool,
        arrow_length: float,
        track_layout: str = "middle",
        track_axis_gap: float | None = None,
    ) -> None:
        self.genome_length = genome_length
        self.alignment_width = alignment_width
        self.cds_height = cds_height
        self.genome_size_normalization_factor = genome_size_normalization_factor
        self.feature_strand = feature_strand
        self.separate_strands = separate_strands
        self.arrow_length = arrow_length
        self.track_layout = track_layout
        self.track_axis_gap = track_axis_gap

    def generate_linear_gene_path(self, gene_object):
        feature_track_id = gene_object.feature_track_id
        coords = gene_object.location
        coordinates_paths = []

        for coord in coords:
            coord_dict = {
                "feat_type": coord.kind,
                "feat_strand": coord.strand,
                "feat_start": coord.start,
                "feat_end": coord.end,
            }
            feat_type = coord_dict["feat_type"]

            if feat_type == "line":
                coord_path = create_intron_path_linear(
                    coord_dict=coord_dict,
                    genome_length=self.genome_length,
                    alignment_width=self.alignment_width,
                    genome_size_normalization_factor=self.genome_size_normalization_factor,
                    cds_height=self.cds_height,
                    feature_strand=self.feature_strand,
                    separate_strands=self.separate_strands,
                    feature_track_id=feature_track_id,
                    track_layout=self.track_layout,
                    track_axis_gap=self.track_axis_gap,
                )
            elif feat_type == "block":
                if coord.is_last and gene_object.is_directional:
                    coord_path = create_arrowhead_path_linear(
                        coord_dict=coord_dict,
                        arrow_length=self.arrow_length,
                        cds_height=self.cds_height,
                        feature_strand=self.feature_strand,
                        genome_length=self.genome_length,
                        alignment_width=self.alignment_width,
                        genome_size_normalization_factor=self.genome_size_normalization_factor,
                        separate_strands=self.separate_strands,
                        feature_track_id=feature_track_id,
                        track_layout=self.track_layout,
                        track_axis_gap=self.track_axis_gap,
                    )
                else:
                    coord_path = create_rectangle_path_linear(
                        coord_dict=coord_dict,
                        genome_length=self.genome_length,
                        alignment_width=self.alignment_width,
                        genome_size_normalization_factor=self.genome_size_normalization_factor,
                        cds_height=self.cds_height,
                        feature_strand=self.feature_strand,
                        separate_strands=self.separate_strands,
                        feature_track_id=feature_track_id,
                        track_layout=self.track_layout,
                        track_axis_gap=self.track_axis_gap,
                    )
            else:
                continue

            coordinates_paths.append(coord_path)

        return coordinates_paths


__all__ = ["FeatureDrawer", "FeaturePathGenerator"]


