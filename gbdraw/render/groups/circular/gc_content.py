#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from gbdraw.svg.elements import Group

from ...drawers.circular.gc_content import GcContentDrawer
from ....configurators import GcContentConfigurator


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class GcContentGroup:
    """
    This class is responsible for creating a group for GC content visualization on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        gc_df: DataFrame,
        radius: float,
        track_width: float,
        gc_config: GcContentConfigurator,
        *,
        norm_factor_override: float,
        group_id: str | None = None,
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        self.record_transform = record_transform
        self.group_id = group_id or "gc_content"
        self.gc_group: Group = Group(id=self.group_id)
        self.radius: float = radius
        self.gc_config: GcContentConfigurator = gc_config
        self.gb_record: SeqRecord = gb_record
        self.record_len: int = len(self.gb_record.seq)
        self.gc_df: DataFrame = gc_df
        self.track_width: float = track_width
        self.norm_factor = float(norm_factor_override)
        self.dinucleotide: str = self.gc_config.dinucleotide
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        self.gc_group = GcContentDrawer(self.gc_config).draw(
            self.radius,
            self.gc_group,
            self.gc_df,
            self.record_len,
            self.track_width,
            self.norm_factor,
            self.dinucleotide,
            self.group_id,
            record_transform=self.record_transform,
        )

    def get_group(self) -> Group:
        return self.gc_group


__all__ = ["GcContentGroup"]


