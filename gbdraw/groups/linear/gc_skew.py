#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group

from ...analysis.skew import skew_df
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.linear.gc_skew import SkewDrawer
from ...configurators import GcSkewConfigurator


class GcSkewGroup:
    """
    Manages the visualization of GC skew for a genomic sequence in a linear layout.
    (This is the new independent class)
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        longest_record_len: int,
        track_height: float,
        alignment_width: float,
        skew_config: GcSkewConfigurator,
        config_dict: dict,
        start_x: float = 0,
        start_y: float = 0,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the GcSkewGroup with the given parameters and configurations.
        """
        self.skew_group = Group(id="gc_skew")
        self.start_x: float = start_x
        self.start_y: float = start_y
        self.longest_record_len: int = longest_record_len
        self.skew_config: GcSkewConfigurator = skew_config
        self.gb_record: SeqRecord = gb_record
        self.window: int = self.skew_config.window
        self.step: int = self.skew_config.step
        self.dinucleotide: str = self.skew_config.dinucleotide
        self.track_height: float = track_height
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.bool_normalize_length = cfg.canvas.linear.normalize_length
        self.alignment_width: float = alignment_width
        self.generate_gc_df()
        self.normalize_length()
        self.add_elements_to_group()

    def normalize_length(self) -> None:
        """
        Normalizes the length of the genomic record relative to the longest record.
        """
        self.record_len: int = len(self.gb_record.seq)
        if self.bool_normalize_length:
            self.genome_size_normalization_factor: float = 1.0
        else:
            self.genome_size_normalization_factor: float = self.record_len / self.longest_record_len

    def generate_gc_df(self) -> None:
        self.skew_df: DataFrame = skew_df(self.gb_record, self.window, self.step, self.dinucleotide)

    def add_elements_to_group(self) -> None:
        """
        Adds the GC skew visualization elements to the SVG group.
        """
        self.skew_group: Group = SkewDrawer(self.skew_config).draw(
            self.skew_group,
            self.skew_df,
            self.record_len,
            self.alignment_width,
            self.genome_size_normalization_factor,
            self.track_height,
            self.start_x,
            self.start_y,
            self.dinucleotide,
        )

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the GC skew visualization.
        """
        return self.skew_group


__all__ = ["GcSkewGroup"]


