from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from .canvas import CanvasConfig
from .labels import LabelsConfig
from .objects import ObjectsConfig


@dataclass(frozen=True)
class GbdrawConfig:
    canvas: CanvasConfig
    labels: LabelsConfig
    objects: ObjectsConfig

    @classmethod
    def from_dict(cls, config_dict: Mapping[str, Any]) -> "GbdrawConfig":
        return cls(
            canvas=CanvasConfig.from_dict(config_dict["canvas"]),
            labels=LabelsConfig.from_dict(config_dict["labels"]),
            objects=ObjectsConfig.from_dict(config_dict["objects"]),
        )


__all__ = ["GbdrawConfig"]


