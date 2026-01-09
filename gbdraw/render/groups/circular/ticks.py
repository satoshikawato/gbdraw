#!/usr/bin/env python
# coding: utf-8

from typing import Optional, Literal

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.path import Path  # type: ignore[reportMissingImports]
from svgwrite.text import Text  # type: ignore[reportMissingImports]

from ....canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....svg.circular_ticks import (  # type: ignore[reportMissingImports]
    generate_circular_tick_paths,
    generate_circular_tick_labels,
)


class TickGroup:
    """
    Represents a group for tick marks and labels on a circular genomic plot.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        config_dict: dict,
        radius: float | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.radius: float = float(radius) if radius is not None else self.canvas_config.radius
        self.tick_group = Group(id="tick")
        self.total_len: int = len(self.gb_record.seq)
        self.config_dict: dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        ticks_config = cfg.objects.ticks
        self.tick_width = ticks_config.tick_width
        self.stroke = ticks_config.tick_labels.stroke
        self.fill = ticks_config.tick_labels.fill
        self.font_size = ticks_config.tick_labels.font_size
        self.font_weight = ticks_config.tick_labels.font_weight
        self.manual_interval = cfg.objects.scale.interval
        self.font_family = cfg.objects.text.font_family
        self.track_type = cfg.canvas.circular.track_type
        self.separate_strands = cfg.canvas.strandedness
        self.dpi = self.canvas_config.dpi
        self.set_tick_size()
        self.add_elements_to_group()

    def set_tick_size(self) -> None:
        if self.manual_interval is not None and self.manual_interval > 0:
            self.tick_large = self.manual_interval
            self.tick_small = self.manual_interval // 10
        else:
            if self.total_len <= 30000:
                tick_large, tick_small = 1000, 100
            elif 30000 < self.total_len <= 50000:
                tick_large, tick_small = 5000, 1000
            elif 50000 < self.total_len <= 150000:
                tick_large, tick_small = 10000, 1000
            elif 150000 < self.total_len <= 1000000:
                tick_large, tick_small = 50000, 10000
            elif 1000000 < self.total_len <= 10000000:
                tick_large, tick_small = 500000, 100000
            else:
                tick_large, tick_small = 1000000, 200000
            self.tick_large: Literal[10000, 50000, 200000, 1000000] = tick_large
            self.tick_small: Literal[1000, 10000, 50000, 200000] = tick_small

    def add_elements_to_group(self) -> Group:
        ticks_large = list(range(0, self.total_len, self.tick_large))
        size: str = "large"
        tick_paths_large: list[Path] = generate_circular_tick_paths(
            self.radius,
            self.total_len,
            size,
            ticks_large,
            self.tick_width,
            self.track_type,
            self.separate_strands,
        )
        ticks_large_nonzero: list[int] = [x for x in ticks_large if x != 0]
        tick_label_paths_large: list[Text] = generate_circular_tick_labels(
            self.radius,
            self.total_len,
            size,
            ticks_large_nonzero,
            self.stroke,
            self.fill,
            self.font_size,
            self.font_weight,
            self.font_family,
            self.track_type,
            self.separate_strands,
            self.dpi,
        )
        for tick_path_large in tick_paths_large:
            self.tick_group.add(tick_path_large)
        for tick_label_path_large in tick_label_paths_large:
            self.tick_group.add(tick_label_path_large)
        return self.tick_group

    def get_group(self) -> Group:
        return self.tick_group


__all__ = ["TickGroup"]


