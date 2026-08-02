#!/usr/bin/env python
# coding: utf-8

from typing import Optional

from svgwrite.container import Group
from svgwrite.path import Path

from ....configurators import FeatureDrawingConfigurator
from ....features.ids import compute_feature_object_hash, make_linear_rendered_feature_id
from ....layout.linear import LinearFeatureLaneGeometry
from ....svg.ids import instance_svg_id
from ....svg.linear_features import (
    create_arrow_path_linear,
    create_arrow_shaft_path_linear,
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
        arrow_geometry = getattr(feature_config, "arrow_geometry", None)
        self.head_length_ratio = getattr(
            arrow_geometry,
            "head_length_ratio",
            "auto",
        )
        self.shaft_width_ratio: float = float(
            getattr(arrow_geometry, "shaft_width_ratio", 1.0)
        )

    @staticmethod
    def get_feature_data_id(feature_object) -> Optional[str]:
        return compute_feature_object_hash(feature_object)

    def draw_path(
        self,
        path_data: str,
        group: Group,
        fill_color: str,
        stroke_color_specified: Optional[str] = None,
        stroke_width_specified: Optional[float] = None,
        feature_data_id: Optional[str] = None,
        dom_element_id: Optional[str] = None,
        rendered_feature_id: Optional[str] = None,
        stable_feature_id: Optional[str] = None,
        record_id: Optional[str] = None,
        record_index: int | None = None,
        feature_part: Optional[str] = None,
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
            debug=False,
        )
        if feature_data_id:
            path.attribs["data-gbdraw-feature-id"] = feature_data_id
            path.attribs["id"] = dom_element_id or feature_data_id
            if rendered_feature_id:
                path.attribs["data-gbdraw-rendered-feature-id"] = (
                    rendered_feature_id
                )
            if stable_feature_id:
                path.attribs["data-gbdraw-stable-feature-id"] = stable_feature_id
            if record_id:
                path.attribs["data-gbdraw-record-id"] = record_id
            if record_index is not None:
                path.attribs["data-gbdraw-record-index"] = str(int(record_index))
            if feature_part:
                path.attribs["data-gbdraw-feature-part"] = feature_part
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
        feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
        record_index: int = 0,
        record_count: int = 1,
        feature_instance_id: str | None = None,
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
            feature_lane_geometry=feature_lane_geometry,
            head_length_ratio=self.head_length_ratio,
            shaft_width_ratio=self.shaft_width_ratio,
        )

        gene_paths = path_generator.generate_linear_gene_path(feature_object)

        # Get feature identifier for instant preview support
        stable_feature_id = self.get_feature_data_id(feature_object)
        feature_data_id = make_linear_rendered_feature_id(
            record_index=record_index,
            stable_feature_id=stable_feature_id,
            record_count=record_count,
        )
        feature_dom_id = (
            instance_svg_id(feature_data_id, feature_instance_id)
            if feature_data_id and feature_instance_id is not None
            else feature_data_id
        )
        block_count = sum(1 for path_type, *_rest in gene_paths if path_type == "block")
        block_index = 0
        line_index = 0

        for gene_path in gene_paths:
            path_type, path_data = gene_path[0], gene_path[1]
            if path_type == "block":
                block_index += 1
                dom_element_id = None
                if feature_dom_id:
                    dom_element_id = (
                        feature_dom_id
                        if block_count <= 1
                        else f"{feature_dom_id}__part{block_index}"
                    )
                self.draw_path(
                    path_data,
                    group,
                    fill_color=feature_object.color,
                    feature_data_id=feature_data_id,
                    dom_element_id=dom_element_id,
                    rendered_feature_id=(
                        feature_dom_id
                        if feature_instance_id is not None
                        else None
                    ),
                    stable_feature_id=stable_feature_id,
                    record_id=getattr(feature_object, "record_id", None),
                    record_index=record_index,
                    feature_part="block",
                )
            elif path_type == "line":
                line_index += 1
                dom_element_id = f"{feature_dom_id}__line{line_index}" if feature_dom_id else None
                self.draw_path(
                    path_data,
                    group,
                    fill_color="none",
                    stroke_color_specified=self.intron_stroke_color,
                    stroke_width_specified=self.intron_stroke_width,
                    feature_data_id=feature_data_id,
                    dom_element_id=dom_element_id,
                    rendered_feature_id=(
                        feature_dom_id
                        if feature_instance_id is not None
                        else None
                    ),
                    stable_feature_id=stable_feature_id,
                    record_id=getattr(feature_object, "record_id", None),
                    record_index=record_index,
                    feature_part="connector",
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
        feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
        head_length_ratio: str | float = "auto",
        shaft_width_ratio: float = 1.0,
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
        self.feature_lane_geometry = feature_lane_geometry
        self.head_length_ratio = head_length_ratio
        self.shaft_width_ratio = shaft_width_ratio

    def generate_linear_gene_path(self, gene_object):
        feature_track_id = gene_object.feature_track_id
        coords = gene_object.location
        coordinates_paths = []
        feature_y_positions = None
        if self.feature_lane_geometry is not None:
            feature_y_positions = self.feature_lane_geometry.lane_for(
                strand=self.feature_strand,
                track_id=feature_track_id,
                separate_strands=self.separate_strands,
            ).positions

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
                    feature_y_positions=feature_y_positions,
                )
            elif feat_type == "block":
                glyph_kind = getattr(
                    gene_object,
                    "glyph_kind",
                    "arrow" if gene_object.is_directional else "rectangle",
                )
                has_defined_strand = coord.strand in {"positive", "negative"}
                if coord.is_last and glyph_kind == "arrow" and has_defined_strand:
                    coord_path = create_arrow_path_linear(
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
                        feature_y_positions=feature_y_positions,
                        head_length_ratio=self.head_length_ratio,
                        shaft_width_ratio=self.shaft_width_ratio,
                    )
                elif (
                    glyph_kind == "arrow"
                    and has_defined_strand
                    and self.shaft_width_ratio < 1.0
                ):
                    coord_path = create_arrow_shaft_path_linear(
                        coord_dict=coord_dict,
                        genome_length=self.genome_length,
                        alignment_width=self.alignment_width,
                        genome_size_normalization_factor=self.genome_size_normalization_factor,
                        cds_height=self.cds_height,
                        feature_strand=self.feature_strand,
                        separate_strands=self.separate_strands,
                        feature_track_id=feature_track_id,
                        shaft_width_ratio=self.shaft_width_ratio,
                        track_layout=self.track_layout,
                        track_axis_gap=self.track_axis_gap,
                        feature_y_positions=feature_y_positions,
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
                        feature_y_positions=feature_y_positions,
                    )
            else:
                continue

            coordinates_paths.append(coord_path)

        return coordinates_paths


__all__ = ["FeatureDrawer", "FeaturePathGenerator"]

