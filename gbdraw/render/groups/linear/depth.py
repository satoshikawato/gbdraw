#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ....analysis.depth import depth_df as build_depth_df  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....configurators import DepthConfigurator  # type: ignore[reportMissingImports]
from ...drawers.linear.depth import DepthDrawer


class DepthGroup:
    """Manages a linear depth coverage group."""

    def __init__(
        self,
        gb_record: SeqRecord,
        longest_record_len: int,
        track_height: float,
        alignment_width: float,
        depth_config: DepthConfigurator,
        config_dict: dict,
        depth_table: DataFrame | None = None,
        start_x: float = 0,
        start_y: float = 0,
        cfg: GbdrawConfig | None = None,
        depth_df: DataFrame | None = None,
    ) -> None:
        self.depth_group = Group(id="depth")
        self.start_x = float(start_x)
        self.start_y = float(start_y)
        self.longest_record_len = int(longest_record_len)
        self.depth_config = depth_config
        self.gb_record = gb_record
        self.track_height = float(track_height)
        self.alignment_width = float(alignment_width)
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.bool_normalize_length = cfg.canvas.linear.normalize_length
        self.record_len = len(self.gb_record.seq)
        self.genome_size_normalization_factor = (
            1.0 if self.bool_normalize_length else self.record_len / self.longest_record_len
        )
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
            self.alignment_width,
            self.genome_size_normalization_factor,
            self.track_height,
            self.start_x,
            self.start_y,
        )

    def get_group(self) -> Group:
        return self.depth_group


__all__ = ["DepthGroup"]
