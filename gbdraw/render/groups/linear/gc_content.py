#!/usr/bin/env python
# coding: utf-8

from pandas import DataFrame  # type: ignore[reportMissingImports]
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ....analysis.skew import skew_df  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.linear.gc_content import GcContentDrawer  # type: ignore[reportMissingImports]
from ....configurators import GcContentConfigurator  # type: ignore[reportMissingImports]


class GcContentGroup:
    """
    Manages the visualization of GC content for a genomic sequence in a linear layout.

    Attributes:
        gb_record (SeqRecord): GenBank record containing genomic data.
        longest_record_len (int): Length of the longest record for normalization.
        track_height (float): Height of the GC content track.
        alignment_width (float): Width of the alignment area.
        gc_config (GcContentConfigurator): Configuration for GC content visualization.
        config_dict (dict): Configuration dictionary with styling parameters.
        start_x (float): Starting x-coordinate for the GC content visualization.
        start_y (float): Starting y-coordinate for the GC content visualization.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        longest_record_len: int,
        track_height: float,
        alignment_width: float,
        gc_config: GcContentConfigurator,
        config_dict: dict,
        start_x: float = 0,
        start_y: float = 0,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the GcContentGroup with the given parameters and configurations.

        Args:
            gb_record (SeqRecord): GenBank record containing genomic data.
            longest_record_len (int): Length of the longest record for normalization.
            track_height (float): Height of the GC content track.
            alignment_width (float): Width of the alignment area.
            gc_config (GcContentConfigurator): Configuration for GC content visualization.
            config_dict (dict): Configuration dictionary with styling parameters.
            start_x (float): Starting x-coordinate for the GC content visualization.
            start_y (float): Starting y-coordinate for the GC content visualization.
        """
        self.gc_group = Group(id="gc_content")
        self.start_x: float = start_x
        self.start_y: float = start_y
        self.longest_record_len: int = longest_record_len
        self.gc_config: GcContentConfigurator = gc_config
        self.gb_record: SeqRecord = gb_record
        self.window: int = self.gc_config.window
        self.step: int = self.gc_config.step
        self.dinucleotide: str = self.gc_config.dinucleotide
        self.track_height: float = track_height
        self.alignment_width: float = alignment_width
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.bool_normalize_length = cfg.canvas.linear.normalize_length
        self.normalize_length()
        self.generate_gc_df()
        self.add_elements_to_group()

    def normalize_length(self) -> None:
        """
        Normalizes the length of the genomic record relative to the longest record.

        This method calculates the normalization factor based on the length of the current
        genomic record and the length of the longest record in the dataset. This factor is used
        to scale the visualization appropriately.
        """
        self.record_len: int = len(self.gb_record.seq)
        if self.bool_normalize_length:
            self.genome_size_normalization_factor: float = 1.0
        else:
            self.genome_size_normalization_factor: float = self.record_len / self.longest_record_len

    def generate_gc_df(self) -> None:
        """
        Generates a DataFrame for GC content based on the genomic record.

        This method uses the `skew_df` function to calculate the GC content of the genomic record.
        The resulting DataFrame is used for visualizing the GC content.
        """
        self.gc_df: DataFrame = skew_df(self.gb_record, self.window, self.step, self.dinucleotide)

    def add_elements_to_group(self) -> None:
        """
        Adds the GC content visualization elements to the SVG group.

        This method calls the GcContentDrawer to draw the GC content based on the calculated
        DataFrame and adds the resulting SVG elements to the group.
        """
        self.gc_group: Group = GcContentDrawer(self.gc_config).draw(
            self.gc_group,
            self.gc_df,
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
        Retrieves the SVG group containing the GC content visualization.

        Returns:
            Group: The SVG group with GC content visualization elements.
        """
        return self.gc_group


__all__ = ["GcContentGroup"]


