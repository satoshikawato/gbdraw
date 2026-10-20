#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group

from ....analysis.skew import skew_df
from ...drawers.linear.gc_skew import SkewDrawer
from ....configurators import GcSkewConfigurator


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class GcSkewGroup:
    """
    Manages the visualization of GC skew for a genomic sequence in a linear layout.
    (This is the new independent class)
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        track_height: float,
        skew_config: GcSkewConfigurator,
        *,
        sequence_width: float,
        start_x: float = 0,
        start_y: float = 0,
        gc_df: DataFrame | None = None,
        group_id: str = "gc_skew",
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        """
        Initializes the GcSkewGroup with the given parameters and configurations.
        """
        self.record_transform = record_transform
        self.skew_group: Group = Group(id=group_id)
        self.start_x: float = start_x
        self.start_y: float = start_y
        self.skew_config: GcSkewConfigurator = skew_config
        self.group_id = group_id
        self.gb_record: SeqRecord = gb_record
        self.window: int = self.skew_config.window
        self.step: int = self.skew_config.step
        self.dinucleotide: str = self.skew_config.dinucleotide
        self.track_height: float = track_height
        # The record's drawn axis width, from its placement.
        self.sequence_width: float = float(sequence_width)
        self.record_len: int = len(self.gb_record.seq)
        self.generate_gc_df(gc_df)
        self.add_elements_to_group()

    def generate_gc_df(self, gc_df: DataFrame | None = None) -> None:
        if gc_df is not None:
            self.skew_df = gc_df
            return
        self.skew_df = skew_df(self.gb_record, self.window, self.step, self.dinucleotide)

    def add_elements_to_group(self) -> None:
        """
        Adds the GC skew visualization elements to the SVG group.
        """
        self.skew_group = SkewDrawer(self.skew_config).draw(
            self.skew_group,
            self.skew_df,
            self.record_len,
            self.sequence_width,
            self.track_height,
            self.start_x,
            self.start_y,
            self.dinucleotide,
            self.group_id,
            record_transform=self.record_transform,
        )

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the GC skew visualization.
        """
        return self.skew_group


__all__ = ["GcSkewGroup"]


