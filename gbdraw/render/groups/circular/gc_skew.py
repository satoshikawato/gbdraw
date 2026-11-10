#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from gbdraw.svg.elements import Group

from ...drawers.circular.gc_skew import SkewDrawer
from ....configurators import GcSkewConfigurator


from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import circular_display_length, unset_record_transform

class GcSkewGroup:
    """
    Represents a group for GC skew visualization on a circular genomic plot.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        gc_df: DataFrame,
        radius: float,
        track_width: float,
        skew_config: GcSkewConfigurator,
        *,
        norm_factor_override: float,
        group_id: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        # The transform the record is drawn by: whether its lines wrap round the origin or hold to 12 o'clock.
        self.record_transform = record_transform if record_transform is not None else unset_record_transform(gb_record)
        self.gb_record: SeqRecord = gb_record
        self.gc_df: DataFrame = gc_df
        self.radius: float = radius
        self.track_width: float = track_width
        self.skew_config: GcSkewConfigurator = skew_config
        self.record_len: int = circular_display_length(self.gb_record, self.record_transform)
        self.skew_group: Group = Group(id=group_id)
        self.norm_factor = float(norm_factor_override)
        self.dinucleotide: str = self.skew_config.dinucleotide
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        self.skew_group = SkewDrawer(self.skew_config).draw(
            self.radius,
            self.skew_group,
            self.gc_df,
            self.record_len,
            self.track_width,
            self.norm_factor,
            self.dinucleotide,
            record_identifier=self.gb_record.id,
            group_identifier=str(self.skew_group.attribs["id"]),
            record_transform=self.record_transform,
        )

    def get_group(self) -> Group:
        return self.skew_group


__all__ = ["GcSkewGroup"]


