#!/usr/bin/env python
# coding: utf-8

from typing import Dict

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ....core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.circular.gc_skew import SkewDrawer  # type: ignore[reportMissingImports]
from ....configurators import GcSkewConfigurator  # type: ignore[reportMissingImports]


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
        config_dict: Dict,
        track_id: str,
        norm_factor_override: float | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.gc_df: DataFrame = gc_df
        self.radius: float = radius
        self.track_width: float = track_width
        self.skew_config: GcSkewConfigurator = skew_config
        self.record_len: int = len(self.gb_record.seq)
        self.skew_group = Group(id="skew")
        self.config_dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.track_type = cfg.canvas.circular.track_type
        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.norm_factor = (
            float(norm_factor_override)
            if norm_factor_override is not None
            else cfg.canvas.circular.track_dict[self.length_param][self.track_type][str(track_id)]
        )
        self.dinucleotide: str = self.skew_config.dinucleotide
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        self.skew_group: Group = SkewDrawer(self.skew_config).draw(
            self.radius,
            self.skew_group,
            self.gc_df,
            self.record_len,
            self.track_width,
            self.norm_factor,
            self.dinucleotide,
        )

    def get_group(self) -> Group:
        return self.skew_group


__all__ = ["GcSkewGroup"]


