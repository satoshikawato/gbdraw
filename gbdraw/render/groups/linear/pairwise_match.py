#!/usr/bin/env python
# coding: utf-8

from typing import Tuple, Dict

from pandas import DataFrame
from svgwrite.container import Group
from svgwrite.path import Path

from ....canvas import LinearCanvasConfigurator
from ....layout.linear_coords import normalize_position_linear
from ....core.color import interpolate_color


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
        self.comparison_count: int = comparison_count
        self.records = records
        self.track_id: str = "comparison" + str(self.comparison_count)
        self.calculate_query_subject_offsets()
        self.match_group = Group(id=self.track_id)
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
        if self.canvas_config.normalize_length:
            self.query_offset_x = 0
            self.subject_offset_x = 0
        else:
            if self.canvas_config.align_center:
                qlen = len(self.records[self.comparison_count - 1].seq)
                slen = len(self.records[self.comparison_count].seq)
                self.query_offset_x = (self.canvas_config.longest_genome - qlen) / 2
                self.subject_offset_x = (self.canvas_config.longest_genome - slen) / 2
            else:
                self.query_offset_x = self.subject_offset_x = 0

    def generate_linear_match_path(self, row: DataFrame) -> Path:
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
        dynamic_fill_color = interpolate_color(self.match_min_color, self.match_max_color, factor)
        query_start, query_end, subject_start, subject_end = self.calculate_offsets(row)
        query_start_x, query_start_y, query_end_x, query_end_y = self.normalize_positions(
            query_start, query_end, 0, is_query=True
        )
        subject_start_x, subject_start_y, subject_end_x, subject_end_y = self.normalize_positions(
            subject_start, subject_end, self.comparison_height, is_query=False
        )

        match_path_desc: str = self.construct_path_description(
            query_start_x,
            query_start_y,
            query_end_x,
            query_end_y,
            subject_start_x,
            subject_start_y,
            subject_end_x,
            subject_end_y,
        )
        return Path(
            d=match_path_desc,
            fill=dynamic_fill_color,
            fill_opacity=self.match_fill_opacity,
            stroke=self.match_stroke_color,
            stroke_width=self.match_stroke_width,
        )

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
        if self.canvas_config.normalize_length:
            # Normalize based on the length of each record itself
            if is_query:
                genome_length = len(self.records[self.comparison_count - 1].seq)
            else:
                genome_length = len(self.records[self.comparison_count].seq)

            start_x = self.canvas_config.alignment_width * (start / genome_length)
            end_x = self.canvas_config.alignment_width * (end / genome_length)
        else:
            # Normalize based on the longest genome as before
            start_x = normalize_position_linear(start, self.canvas_config.longest_genome, self.canvas_config.alignment_width)
            end_x = normalize_position_linear(end, self.canvas_config.longest_genome, self.canvas_config.alignment_width)

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

    def add_elements_to_group(self) -> Group:
        """
        Adds all the pairwise match paths to the SVG group.

        Iterates through the pairwise matches and creates SVG path elements for each,
        adding them to the group.

        Returns:
            Group: The SVG group with all match paths added.
        """
        for row in self.comparison_df.itertuples():
            match_path: Path = self.generate_linear_match_path(row)
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


