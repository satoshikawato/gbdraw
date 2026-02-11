"""Typed config models (internal).

These dataclasses provide a structured view over the nested `config_dict` used across gbdraw.
They are introduced incrementally; callers can still pass raw dicts, but core modules can
opt into these models for clearer typing and separation.
"""

from .canvas import CanvasConfig, CircularCanvasConfig, LinearCanvasConfig
from .labels import LabelsConfig, LabelsFilteringConfig, LabelsLengthThresholdConfig, LabelsLinearConfig
from .objects import ObjectsConfig
from .root import GbdrawConfig

__all__ = [
    "CanvasConfig",
    "CircularCanvasConfig",
    "LinearCanvasConfig",
    "LabelsConfig",
    "LabelsFilteringConfig",
    "LabelsLengthThresholdConfig",
    "LabelsLinearConfig",
    "ObjectsConfig",
    "GbdrawConfig",
]


