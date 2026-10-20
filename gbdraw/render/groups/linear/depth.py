#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from gbdraw.svg.elements import Group

from ....analysis.depth import depth_df as build_depth_df
from ....configurators import DepthConfigurator
from ...drawers.linear.depth import DepthDrawer


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class DepthGroup:
    """Manages a linear depth coverage group."""

    def __init__(
        self,
        gb_record: SeqRecord,
        track_height: float,
        depth_config: DepthConfigurator,
        *,
        sequence_width: float,
        depth_table: DataFrame | None = None,
        start_x: float = 0,
        start_y: float = 0,
        depth_df: DataFrame | None = None,
        group_id: str = "depth",
        axis_group_id: str = "depth_axis",
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        self.record_transform = record_transform
        self.depth_group = Group(id=group_id)
        self.start_x = float(start_x)
        self.start_y = float(start_y)
        self.depth_config = depth_config
        self.gb_record = gb_record
        self.track_height = float(track_height)
        # The record's drawn axis width, from its placement.
        self.sequence_width = float(sequence_width)
        self.record_len = len(self.gb_record.seq)
        self.axis_group_id = str(axis_group_id)
        if depth_df is not None:
            self.depth_df = depth_df
        elif depth_table is not None:
            self.depth_df = build_depth_df(
                self.gb_record,
                depth_table,
                self.depth_config.window,
                self.depth_config.step,
                normalize=self.depth_config.normalize,
                min_depth=self.depth_config.min_depth,
                max_depth=self.depth_config.max_depth,
            )
        else:
            self.depth_df = DataFrame(columns=["position", "depth", "depth_normalized"])
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        self.depth_group = DepthDrawer(self.depth_config).draw(
            self.depth_group,
            self.depth_df,
            self.record_len,
            self.sequence_width,
            1.0,
            self.track_height,
            self.start_x,
            self.start_y,
            self.axis_group_id,
            record_transform=self.record_transform,
        )

    def get_group(self) -> Group:
        return self.depth_group


__all__ = ["DepthGroup"]
