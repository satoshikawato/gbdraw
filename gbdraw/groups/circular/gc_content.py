#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.circular.gc_content import GcContentDrawer  # type: ignore[reportMissingImports]
from ...configurators import GcContentConfigurator  # type: ignore[reportMissingImports]


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
        config_dict: dict,
        track_id: str,
        norm_factor_override: float | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gc_group = Group(id="gc_content")
        self.radius: float = radius
        self.gc_config: GcContentConfigurator = gc_config
        self.gb_record: SeqRecord = gb_record
        self.record_len: int = len(self.gb_record.seq)
        self.config_dict: dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.gc_df: DataFrame = gc_df
        self.track_width: float = track_width
        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.track_type = cfg.canvas.circular.track_type
        self.norm_factor = (
            float(norm_factor_override)
            if norm_factor_override is not None
            else cfg.canvas.circular.track_dict[self.length_param][self.track_type][str(track_id)]
        )
        self.dinucleotide: str = self.gc_config.dinucleotide
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        self.gc_group: Group = GcContentDrawer(self.gc_config).draw(
            self.radius,
            self.gc_group,
            self.gc_df,
            self.record_len,
            self.track_width,
            self.norm_factor,
            self.dinucleotide,
        )

    def get_group(self) -> Group:
        return self.gc_group


__all__ = ["GcContentGroup"]


