#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....configurators import DepthConfigurator  # type: ignore[reportMissingImports]
from ....core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ...drawers.circular.depth import DepthDrawer


class DepthGroup:
    """Creates a depth coverage group for a circular canvas."""

    def __init__(
        self,
        gb_record: SeqRecord,
        depth_df: DataFrame,
        radius: float,
        track_width: float,
        depth_config: DepthConfigurator,
        config_dict: dict,
        track_id: str | int,
        norm_factor_override: float | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.depth_group = Group(id="depth")
        self.radius = float(radius)
        self.depth_config = depth_config
        self.gb_record = gb_record
        self.record_len = len(self.gb_record.seq)
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.depth_df = depth_df
        self.track_width = float(track_width)
        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.track_type = cfg.canvas.circular.track_type
        self.norm_factor = (
            float(norm_factor_override)
            if norm_factor_override is not None
            else cfg.canvas.circular.track_dict[self.length_param][self.track_type][str(track_id)]
        )
        self.add_elements_to_group()

    def add_elements_to_group(self) -> None:
        self.depth_group = DepthDrawer(self.depth_config).draw(
            self.radius,
            self.depth_group,
            self.depth_df,
            self.record_len,
            self.track_width,
            self.norm_factor,
        )

    def get_group(self) -> Group:
        return self.depth_group


__all__ = ["DepthGroup"]
