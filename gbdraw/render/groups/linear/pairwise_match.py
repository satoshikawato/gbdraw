#!/usr/bin/env python
# coding: utf-8

import math
from typing import Tuple, Dict

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path

from ....canvas import LinearCanvasConfigurator
from ....layout.linear_coords import normalize_position_linear
from ....core.color import (
    DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
    DEFAULT_COLLINEAR_ORIENTATION_COLORS,
    interpolate_color,
)
from ....features.ids import make_linear_rendered_feature_id
from ....layout.linear_multi_record import LinearRecordPlacement


def _row_value(row: object, name: str, default: object = "") -> object:
    value = getattr(row, name, default)
    if value is None:
        return default
    return value


def _is_missing_value(value: object) -> bool:
    if value is None:
        return True
    try:
        return bool(value != value)
    except Exception:
        return False


def _attribute_text(value: object) -> str:
    if _is_missing_value(value):
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.12g}"
    return str(value).strip()


def _required_attribute_text(value: object) -> str:
    text = _attribute_text(value)
    return text if text else " "


def _row_float(row: object, name: str, default: float = 0.0) -> float:
    try:
        return float(_row_value(row, name, default))
    except (TypeError, ValueError):
        return float(default)


def _row_span(row: object, start_name: str, end_name: str) -> float:
    return abs(_row_float(row, end_name) - _row_float(row, start_name))


def _match_draw_order_key(row: object) -> tuple[int, float, float, float, float, float, str]:
    block_id = str(_row_value(row, "collinearity_block_id", "") or "").strip()
    orientation = str(_row_value(row, "collinearity_orientation", "") or "").strip().lower()
    orientation_rank = 0 if orientation == "plus" else 2 if orientation == "minus" else 1
    query_span = _row_span(row, "qstart", "qend")
    subject_span = _row_span(row, "sstart", "send")
    area = query_span * subject_span
    total_span = query_span + subject_span
    return (
        orientation_rank,
        -area,
        -total_span,
        _row_float(row, "identity"),
        _row_float(row, "qstart"),
        _row_float(row, "sstart"),
        block_id,
    )


def _normalize_collinearity_color_mode(value: object) -> str:
    return str(value or "").strip().lower().replace("-", "_")


def _match_kind(row: object) -> str:
    collinearity_block_id = _attribute_text(_row_value(row, "collinearity_block_id", ""))
    if collinearity_block_id:
        return "collinear"
    orthogroup_id = _attribute_text(_row_value(row, "orthogroup_id", ""))
    if orthogroup_id:
        return "orthogroup"
    return "pairwise"


class PairWiseMatchGroup:
    """
    Handles the visualization of pairwise matches (comparisons) in a linear layout.

    Attributes:
        canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
        sequence_length_dict (Dict[str, int]): Dictionary of sequence lengths.
        comparison_df (DataFrame): DataFrame containing comparison data.
        comparison_height (float): Height of the comparison track.
        comparison_count (int): Counter for the number of comparisons.
        config_dict (dict): Configuration dictionary with styling parameters.
    """

    def __init__(
        self,
        canvas_config: LinearCanvasConfigurator,
        sequence_length_dict: dict,
        comparison_df: DataFrame,
        actual_comparison_height: float,
        comparison_count: int,
        blast_config,
        records,
        record_offsets_x: dict[int, float] | None = None,
        query_record_index: int | None = None,
        subject_record_index: int | None = None,
        query_placement: LinearRecordPlacement | None = None,
        subject_placement: LinearRecordPlacement | None = None,
        query_y: float | None = None,
        subject_y: float | None = None,
    ) -> None:
        """
        Initializes the PairWiseMatchGroup with necessary data and configurations.

        Args:
            canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
            sequence_length_dict (Dict[str, int]): Dictionary of sequence lengths.
            comparison_df (DataFrame): DataFrame containing comparison data.
            comparison_height (float): Height of the comparison track.
            comparison_count (int): Counter for the number of comparisons.
            config_dict (dict): Configuration dictionary with styling parameters.
        """
        self.canvas_config: LinearCanvasConfigurator = canvas_config
        self.sequence_length_dict: Dict[str, int] = sequence_length_dict
        self.comparison_df: DataFrame = comparison_df
        self.comparison_height: float = actual_comparison_height
        self.match_fill_color: str = blast_config.fill_color
        self.min_identity: float = blast_config.identity
        self.match_min_color: str = blast_config.min_color
        self.match_max_color: str = blast_config.max_color
        self.match_fill_opacity: float = blast_config.fill_opacity
        self.match_stroke_color: str = blast_config.stroke_color
        self.match_stroke_width: float = blast_config.stroke_width
        self.collinearity_orientation_colors: dict[str, str] = {
            **DEFAULT_COLLINEAR_ORIENTATION_COLORS,
            **(getattr(blast_config, "collinearity_orientation_colors", {}) or {}),
        }
        self.collinearity_orientation_min_colors: dict[str, str] = {
            **DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
            **(getattr(blast_config, "collinearity_orientation_min_colors", {}) or {}),
        }
        self.match_style: str = str(getattr(blast_config, "match_style", "ribbon")).lower()
        self.curve_tension: float = float(getattr(blast_config, "curve_tension", 0.5))
        self.comparison_count: int = comparison_count
        self.records = records
        self.query_record_index = (
            int(query_record_index)
            if query_record_index is not None
            else self.comparison_count - 1
        )
        self.subject_record_index = (
            int(subject_record_index)
            if subject_record_index is not None
            else self.comparison_count
        )
        self.query_placement = query_placement
        self.subject_placement = subject_placement
        self.query_y = 0 if query_y is None else float(query_y)
        self.subject_y = self.comparison_height if subject_y is None else float(subject_y)
        self.record_offsets_x = record_offsets_x or {}
        self.track_id: str = "comparison" + str(self.comparison_count)
        self._pairwise_match_counter = 0
        self.calculate_query_subject_offsets()
        self.match_group = Group(id=self.track_id, debug=False)
        if self.query_placement is not None:
            self.match_group.attribs["data-query-record-index"] = self.query_record_index
            self.match_group.attribs["data-subject-record-index"] = self.subject_record_index
            self.match_group.attribs["data-query-row"] = self.query_placement.row
        if self.subject_placement is not None:
            self.match_group.attribs["data-subject-row"] = self.subject_placement.row
        self.add_elements_to_group()

    def calculate_query_subject_offsets(self) -> Tuple[float, float]:
        """
        Calculates the horizontal offsets for the query and subject sequences.

        This method determines the starting x-coordinates for the query and subject in the alignment,
        ensuring they are centered if specified in the configuration.

        Args:
            query_id (str): Identifier for the query sequence.
            subject_id (str): Identifier for the subject sequence.

        Returns:
            Tuple[float, float]: The x-coordinate offsets for the query and subject sequences.
        """
        if self.query_placement is not None and self.subject_placement is not None:
            self.query_offset_x = 0.0
            self.subject_offset_x = 0.0
        elif self.canvas_config.normalize_length:
            self.query_offset_x = 0
            self.subject_offset_x = 0
        else:
            if self.canvas_config.align_center:
                qlen = len(self.records[self.query_record_index].seq)
                slen = len(self.records[self.subject_record_index].seq)
                self.query_offset_x = (self.canvas_config.longest_genome - qlen) / 2
                self.subject_offset_x = (self.canvas_config.longest_genome - slen) / 2
            else:
                self.query_offset_x = self.subject_offset_x = 0
        self.query_alignment_offset_x = float(self.record_offsets_x.get(self.query_record_index, 0.0))
        self.subject_alignment_offset_x = float(self.record_offsets_x.get(self.subject_record_index, 0.0))

    def _record_id_for_index(self, index: int) -> str:
        try:
            record = self.records[index]
        except (AttributeError, IndexError, TypeError):
            return ""
        return _attribute_text(
            getattr(record, "id", None)
            or getattr(record, "name", None)
            or getattr(record, "description", None)
            or ""
        )

    def _next_pairwise_match_id(self, match_index: int | None = None) -> str:
        if match_index is None:
            self._pairwise_match_counter = int(getattr(self, "_pairwise_match_counter", 0)) + 1
            match_index = self._pairwise_match_counter
        track_id = _attribute_text(getattr(self, "track_id", "")) or f"comparison{self.comparison_count}"
        return f"{track_id}_match{int(match_index)}"

    def _rendered_feature_svg_id_values(self, value: object, *, record_index: int) -> str:
        stable_ids = [
            part.strip()
            for part in _attribute_text(value).split(";")
            if part.strip()
        ]
        rendered_ids = [
            make_linear_rendered_feature_id(
                record_index=record_index,
                stable_feature_id=stable_id,
                record_count=len(self.records),
            )
            for stable_id in stable_ids
        ]
        return ";".join(rendered_id for rendered_id in rendered_ids if rendered_id)

    def generate_linear_match_path(self, row: DataFrame, match_index: int | None = None) -> Path:
        """
        Generates an SVG path for a pairwise match based on the provided row data.

        Args:
            row (DataFrame): A row from the DataFrame containing match information.

        Returns:
            Path: An SVG path element representing the pairwise match.
        """
        query_start: float
        query_end: float
        subject_start: float
        subject_end: float
        query_start_x: float
        query_start_y: float
        query_end_x: float
        query_end_y: float
        subject_start_x: float
        subject_start_y: float
        subject_end_x: float
        subject_end_y: float
        identity_percent = float(row.identity)
        factor = (identity_percent - self.min_identity) / (100 - self.min_identity)
        default_gradient_color = interpolate_color(self.match_min_color, self.match_max_color, factor)
        dynamic_fill_color = self.resolve_match_fill_color(row, factor, default_gradient_color)
        query_start, query_end, subject_start, subject_end = self.calculate_offsets(row)
        query_start_x, query_start_y, query_end_x, query_end_y = self.normalize_positions(
            query_start, query_end, getattr(self, "query_y", 0), is_query=True
        )
        subject_start_x, subject_start_y, subject_end_x, subject_end_y = self.normalize_positions(
            subject_start,
            subject_end,
            getattr(self, "subject_y", self.comparison_height),
            is_query=False,
        )

        match_style = str(getattr(self, "match_style", "ribbon")).lower()
        if match_style == "curve":
            match_path_desc: str = self.construct_curved_ribbon_path_description(
                query_start_x,
                query_start_y,
                query_end_x,
                query_end_y,
                subject_start_x,
                subject_start_y,
                subject_end_x,
                subject_end_y,
            )
        else:
            match_style = "ribbon"
            match_path_desc = self.construct_path_description(
                query_start_x,
                query_start_y,
                query_end_x,
                query_end_y,
                subject_start_x,
                subject_start_y,
                subject_end_x,
                subject_end_y,
            )
        path = Path(
            d=match_path_desc,
            fill=dynamic_fill_color,
            fill_opacity=self.match_fill_opacity,
            stroke=self.match_stroke_color,
            stroke_width=self.match_stroke_width,
            debug=False,
        )
        path.attribs["data-pairwise-match-style"] = match_style
        path.attribs["data-identity-factor"] = f"{factor:.6g}"
        self.add_required_metadata_attributes(path, row, match_index=match_index)
        self.add_optional_metadata_attributes(path, row)
        return path

    def add_required_metadata_attributes(
        self,
        path: Path,
        row: DataFrame,
        match_index: int | None = None,
    ) -> None:
        query_record_index = int(
            getattr(self, "query_record_index", self.comparison_count - 1)
        )
        subject_record_index = int(
            getattr(self, "subject_record_index", self.comparison_count)
        )
        query_record_id = (
            _attribute_text(_row_value(row, "query", ""))
            or _attribute_text(_row_value(row, "qseqid", ""))
            or self._record_id_for_index(query_record_index)
        )
        subject_record_id = (
            _attribute_text(_row_value(row, "subject", ""))
            or _attribute_text(_row_value(row, "sseqid", ""))
            or self._record_id_for_index(subject_record_index)
        )
        match_id = self._next_pairwise_match_id(match_index)
        required_attributes = {
            "data-gbdraw-match-id": match_id,
            "data-gbdraw-pairwise-match-id": match_id,
            "data-match-kind": _match_kind(row),
            "data-query-record-index": query_record_index,
            "data-subject-record-index": subject_record_index,
            "data-query-record-id": query_record_id,
            "data-subject-record-id": subject_record_id,
            "data-qstart": _row_value(row, "qstart", ""),
            "data-qend": _row_value(row, "qend", ""),
            "data-sstart": _row_value(row, "sstart", ""),
            "data-send": _row_value(row, "send", ""),
            "data-identity": _row_value(row, "identity", ""),
            "data-alignment-length": _row_value(row, "alignment_length", ""),
            "data-evalue": _row_value(row, "evalue", ""),
            "data-bitscore": _row_value(row, "bitscore", ""),
            "data-mismatches": _row_value(row, "mismatches", ""),
            "data-gap-opens": _row_value(row, "gap_opens", ""),
        }
        for attribute, value in required_attributes.items():
            path.attribs[attribute] = _required_attribute_text(value)

    def resolve_match_fill_color(
        self,
        row: DataFrame,
        factor: float,
        default_gradient_color: str,
    ) -> str:
        collinearity_block_id = str(_row_value(row, "collinearity_block_id", "") or "").strip()
        if not collinearity_block_id:
            return default_gradient_color

        color_mode = _normalize_collinearity_color_mode(
            _row_value(row, "collinearity_color_mode", "")
        )
        if color_mode == "average_identity":
            return default_gradient_color

        orientation = str(_row_value(row, "collinearity_orientation", "") or "").strip().lower()
        orientation_colors = getattr(
            self,
            "collinearity_orientation_colors",
            DEFAULT_COLLINEAR_ORIENTATION_COLORS,
        )
        orientation_color = orientation_colors.get(orientation)
        if not orientation_color:
            return default_gradient_color

        if color_mode == "orientation":
            return orientation_color
        if color_mode == "orientation_identity":
            orientation_min_colors = getattr(
                self,
                "collinearity_orientation_min_colors",
                DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
            )
            min_color = orientation_min_colors.get(orientation, self.match_min_color)
            return interpolate_color(min_color, orientation_color, factor)

        return default_gradient_color

    def add_optional_metadata_attributes(self, path: Path, row: DataFrame) -> None:
        has_collinearity_metadata = bool(_attribute_text(_row_value(row, "collinearity_block_id", "")))
        has_orthogroup_metadata = bool(_attribute_text(_row_value(row, "orthogroup_id", "")))
        if not (has_collinearity_metadata or has_orthogroup_metadata):
            return

        metadata_columns = {
            "collinearity_block_id": "data-collinearity-block-id",
            "collinearity_block_kind": "data-collinearity-block-kind",
            "collinearity_orientation": "data-collinearity-orientation",
            "collinearity_block_score": "data-collinearity-block-score",
            "collinearity_block_evalue": "data-collinearity-block-evalue",
            "collinearity_anchor_index": "data-collinearity-anchor-index",
            "collinearity_anchor_count": "data-collinearity-anchor-count",
            "collinearity_color_mode": "data-collinearity-color-mode",
            "group_kind": "data-group-kind",
            "group_scope": "data-group-scope",
            "collinear_group_scope": "data-collinear-group-scope",
            "orthogroup_id": "data-orthogroup-id",
            "rbh_orthogroup_id": "data-rbh-orthogroup-id",
            "ortholog_path_id": "data-ortholog-path-id",
            "edge_kind": "data-edge-kind",
            "render_role": "data-render-role",
            "query_orthogroup_representative": "data-query-orthogroup-representative",
            "subject_orthogroup_representative": "data-subject-orthogroup-representative",
            "query_orthogroup_member_count": "data-query-orthogroup-member-count",
            "subject_orthogroup_member_count": "data-subject-orthogroup-member-count",
            "query_orthogroup_role": "data-query-orthogroup-role",
            "subject_orthogroup_role": "data-subject-orthogroup-role",
            "query_orthogroup_confidence": "data-query-orthogroup-confidence",
            "subject_orthogroup_confidence": "data-subject-orthogroup-confidence",
            "query_orthogroup_assignment_reason": "data-query-orthogroup-assignment-reason",
            "subject_orthogroup_assignment_reason": "data-subject-orthogroup-assignment-reason",
            "query_protein_id": "data-query-protein-id",
            "subject_protein_id": "data-subject-protein-id",
            "query_feature_svg_id": "data-query-feature-svg-id",
            "subject_feature_svg_id": "data-subject-feature-svg-id",
            "query_unit_id": "data-query-unit-id",
            "subject_unit_id": "data-subject-unit-id",
            "query_locus_id": "data-query-locus-id",
            "subject_locus_id": "data-subject-locus-id",
            "query_display_name": "data-query-display-name",
            "subject_display_name": "data-subject-display-name",
        }
        query_record_index = int(
            getattr(self, "query_record_index", self.comparison_count - 1)
        )
        subject_record_index = int(
            getattr(self, "subject_record_index", self.comparison_count)
        )
        for column, attribute in metadata_columns.items():
            text = _attribute_text(_row_value(row, column, ""))
            if text:
                if column == "query_feature_svg_id":
                    path.attribs["data-query-stable-feature-svg-id"] = text
                    path.attribs[attribute] = self._rendered_feature_svg_id_values(
                        text,
                        record_index=query_record_index,
                    )
                elif column == "subject_feature_svg_id":
                    path.attribs["data-subject-stable-feature-svg-id"] = text
                    path.attribs[attribute] = self._rendered_feature_svg_id_values(
                        text,
                        record_index=subject_record_index,
                    )
                else:
                    path.attribs[attribute] = text

    def calculate_offsets(self, row: DataFrame) -> tuple[float, float, float, float]:
        """
        Calculates the start and end positions for a pairwise match.

        Args:
            row (DataFrame): A row from the DataFrame containing match information.

        Returns:
            Tuple[float, float, float, float]: Start and end positions for the query and subject matches.
        """
        query_start: float = row.qstart + self.query_offset_x
        query_end: float = row.qend + self.query_offset_x
        subject_start: float = row.sstart + self.subject_offset_x
        subject_end: float = row.send + self.subject_offset_x
        return query_start, query_end, subject_start, subject_end

    def normalize_positions(
        self, start: float, end: float, y_position: float, is_query: bool
    ) -> tuple[float, float, float, float]:
        """
        Normalizes the start and end positions for display on the linear canvas.
        """
        placement = (
            getattr(self, "query_placement", None)
            if is_query
            else getattr(self, "subject_placement", None)
        )
        if placement is not None:
            start_x = placement.x_for_position(start)
            end_x = placement.x_for_position(end)
        elif self.canvas_config.normalize_length:
            # Normalize based on the length of each record itself
            if is_query:
                query_record_index = int(
                    getattr(self, "query_record_index", self.comparison_count - 1)
                )
                genome_length = len(self.records[query_record_index].seq)
            else:
                subject_record_index = int(
                    getattr(self, "subject_record_index", self.comparison_count)
                )
                genome_length = len(self.records[subject_record_index].seq)

            start_x = self.canvas_config.alignment_width * (start / genome_length)
            end_x = self.canvas_config.alignment_width * (end / genome_length)
        else:
            # Normalize based on the longest genome as before
            start_x = normalize_position_linear(start, self.canvas_config.longest_genome, self.canvas_config.alignment_width)
            end_x = normalize_position_linear(end, self.canvas_config.longest_genome, self.canvas_config.alignment_width)

        if placement is None:
            alignment_offset = self.query_alignment_offset_x if is_query else self.subject_alignment_offset_x
            start_x += alignment_offset
            end_x += alignment_offset
        return start_x, y_position, end_x, y_position

    def construct_path_description(
        self,
        query_start_x: float,
        query_start_y: float,
        query_end_x: float,
        query_end_y: float,
        subject_start_x: float,
        subject_start_y: float,
        subject_end_x: float,
        subject_end_y: float,
    ) -> str:
        """
        Constructs the path description for an SVG path element representing a pairwise match.

        Args:
            query_start_x (float): x-coordinate of the query start.
            query_start_y (float): y-coordinate of the query start.
            query_end_x (float): x-coordinate of the query end.
            query_end_y (float): y-coordinate of the query end.
            subject_start_x (float): x-coordinate of the subject start.
            subject_start_y (float): y-coordinate of the subject start.
            subject_end_x (float): x-coordinate of the subject end.
            subject_end_y (float): y-coordinate of the subject end.

        Returns:
            str: SVG path description string.
        """
        return f"M {query_start_x},{query_start_y}L{query_end_x},{query_end_y} L{subject_end_x},{subject_end_y}L{subject_start_x},{subject_start_y} z"

    def construct_curved_ribbon_path_description(
        self,
        query_start_x: float,
        query_start_y: float,
        query_end_x: float,
        query_end_y: float,
        subject_start_x: float,
        subject_start_y: float,
        subject_end_x: float,
        subject_end_y: float,
    ) -> str:
        """
        Constructs a closed curved ribbon path preserving both aligned intervals.
        """
        tension = float(getattr(self, "curve_tension", 0.5))
        dy = subject_start_y - query_start_y
        c1y = query_start_y + dy * tension
        c2y = subject_start_y - dy * tension
        return (
            f"M {query_start_x},{query_start_y}"
            f"C{query_start_x},{c1y} {subject_start_x},{c2y} {subject_start_x},{subject_start_y} "
            f"L{subject_end_x},{subject_end_y}"
            f"C{subject_end_x},{c2y} {query_end_x},{c1y} {query_end_x},{query_end_y} z"
        )

    def add_elements_to_group(self) -> Group:
        """
        Adds all the pairwise match paths to the SVG group.

        Iterates through the pairwise matches and creates SVG path elements for each,
        adding them to the group.

        Returns:
            Group: The SVG group with all match paths added.
        """
        rows = sorted(self.comparison_df.itertuples(), key=_match_draw_order_key)
        for index, row in enumerate(rows, start=1):
            match_path: Path = self.generate_linear_match_path(row, match_index=index)
            self.match_group.add(match_path)
        return self.match_group

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing all the pairwise match visualizations.

        Returns:
            Group: The SVG group with pairwise match elements.
        """
        return self.match_group


__all__ = ["PairWiseMatchGroup"]


