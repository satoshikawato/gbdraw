#!/usr/bin/env python
# coding: utf-8

from collections.abc import Sequence

from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group

from ....analysis.skew import skew_df
from ...drawers.linear.gc_content import GcContentDrawer
from ...drawers.linear.wrap_rows import draw_track_rows
from ....configurators import GcContentConfigurator


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class GcContentGroup:
    """
    Manages the visualization of GC content for a genomic sequence in a linear layout.

    Attributes:
        gb_record (SeqRecord): GenBank record containing genomic data.
        track_height (float): Height of the GC content track.
        sequence_width (float): The record's drawn axis width, from its placement.
        gc_config (GcContentConfigurator): Configuration for GC content visualization.
        start_x (float): Starting x-coordinate for the GC content visualization.
        start_y (float): Starting y-coordinate for the GC content visualization.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        track_height: float,
        gc_config: GcContentConfigurator,
        *,
        sequence_width: float,
        start_x: float = 0,
        start_y: float = 0,
        gc_df: DataFrame | None = None,
        group_id: str = "gc_content",
        record_transform: RecordDisplayTransform | None = None,
        row_shifts: Sequence[float] = (),
    ) -> None:
        """
        Initializes the GcContentGroup with the given parameters and configurations.

        Args:
            gb_record (SeqRecord): GenBank record containing genomic data.
            track_height (float): Height of the GC content track.
            sequence_width (float): The record's drawn axis width, from its placement.
            gc_config (GcContentConfigurator): Configuration for GC content visualization.
            start_x (float): Starting x-coordinate for the GC content visualization.
            start_y (float): Starting y-coordinate for the GC content visualization.
        """
        self.record_transform = record_transform
        self.row_shifts = row_shifts
        self.group_id = group_id
        self.gc_group: Group = Group(id=self.group_id)
        self.start_x: float = start_x
        self.start_y: float = start_y
        self.gc_config: GcContentConfigurator = gc_config
        self.gb_record: SeqRecord = gb_record
        self.window: int = self.gc_config.window
        self.step: int = self.gc_config.step
        self.dinucleotide: str = self.gc_config.dinucleotide
        self.track_height: float = track_height
        self.sequence_width: float = float(sequence_width)
        self.record_len: int = len(self.gb_record.seq)
        self.generate_gc_df(gc_df)
        self.add_elements_to_group()

    def generate_gc_df(self, gc_df: DataFrame | None = None) -> None:
        """
        Generates a DataFrame for GC content based on the genomic record.

        This method uses the `skew_df` function to calculate the GC content of the genomic record.
        The resulting DataFrame is used for visualizing the GC content.
        """
        if gc_df is not None:
            self.gc_df = gc_df
            return
        self.gc_df = skew_df(self.gb_record, self.window, self.step, self.dinucleotide)

    def add_elements_to_group(self) -> None:
        """
        Adds the GC content visualization elements to the SVG group.

        This method calls the GcContentDrawer to draw the GC content based on the calculated
        DataFrame and adds the resulting SVG elements to the group.
        """
        drawer = GcContentDrawer(self.gc_config)
        self.gc_group = draw_track_rows(
            self.gc_group,
            self.record_transform,
            self.record_len,
            self.sequence_width,
            self.row_shifts,
            lambda group, length, width, row: drawer.draw(
                group,
                self.gc_df,
                length,
                width,
                self.track_height,
                self.start_x,
                self.start_y,
                self.dinucleotide,
                self.group_id,
                record_transform=self.record_transform,
                row=row,
            ),
        )

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the GC content visualization.

        Returns:
            Group: The SVG group with GC content visualization elements.
        """
        return self.gc_group


__all__ = ["GcContentGroup"]


