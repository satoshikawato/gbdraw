#!/usr/bin/env python
# coding: utf-8

from collections.abc import Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from gbdraw.svg.elements import Group

from ....configurators import DepthConfigurator
from ...drawers.linear.depth import DepthDrawer
from ...drawers.linear.wrap_rows import draw_track_rows


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
        start_x: float = 0,
        start_y: float = 0,
        depth_df: DataFrame | None = None,
        group_id: str,
        axis_group_id: str,
        record_transform: RecordDisplayTransform | None = None,
        row_shifts: Sequence[float] = (),
    ) -> None:
        self.record_transform = record_transform
        self.row_shifts = row_shifts
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
        else:
            self.depth_df = DataFrame(columns=["position", "depth"])
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        drawer = DepthDrawer(self.depth_config)
        self.depth_group = draw_track_rows(
            self.depth_group,
            self.record_transform,
            self.record_len,
            self.sequence_width,
            self.row_shifts,
            lambda group, length, width, row: drawer.draw(
                group,
                self.depth_df,
                length,
                width,
                self.track_height,
                self.start_x,
                self.start_y,
                self.axis_group_id,
                record_transform=self.record_transform,
                row=row,
            ),
        )

    def get_group(self) -> Group:
        return self.depth_group


__all__ = ["DepthGroup"]
