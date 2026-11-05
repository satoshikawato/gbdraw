"""The drawing configurators one diagram is drawn with, built once by the API for either mode.

The canvas configurator is mode-specific and travels beside this bundle, in
the mode's assembler call.
"""

from __future__ import annotations

from dataclasses import dataclass

from .depth import DepthConfigurator
from .features import FeatureDrawingConfigurator
from .gc import GcContentConfigurator, GcSkewConfigurator
from .legend import LegendDrawingConfigurator


@dataclass(frozen=True, slots=True)
class DiagramConfigurators:
    """The feature, GC content, GC skew, legend and depth configurators of one diagram.

    ``depth`` is ``None`` when the diagram draws no depth track.
    """

    feature: FeatureDrawingConfigurator
    gc: GcContentConfigurator
    skew: GcSkewConfigurator
    legend: LegendDrawingConfigurator
    depth: DepthConfigurator | None
