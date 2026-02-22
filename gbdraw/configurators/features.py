#!/usr/bin/env python
# coding: utf-8

from typing import Dict, List, Mapping, Optional

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.features.shapes import (
    normalize_feature_shape_overrides,
    resolve_directional_feature_types,
)


class FeatureDrawingConfigurator:
    """
    Configurator for drawing genomic features.

    This class provides configurations for how genomic features should be visualized, including
    color mapping and selection of features to display.

    Attributes:
        color_table (Optional[DataFrame]): Table defining color codes for different features.
        default_colors (Optional[DataFrame]): Default color settings for features.
        selected_features_set (str): Identifier for a set of features to be visualized.
    """

    def __init__(
        self,
        color_table: Optional[DataFrame],
        default_colors: DataFrame,
        selected_features_set: List[str],
        config_dict: Dict,
        canvas_config,
        feature_shapes: Mapping[str, str] | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the FeatureDrawingConfigurator with color settings and feature selection.

        Args:
            color_table (Optional[DataFrame]): Color table for features.
            default_colors (Optional[DataFrame]): Default colors for features.
            selected_features_set (str): Set identifier for selecting features to display.
        """

        self.color_table: Optional[DataFrame] = color_table
        self.default_colors: DataFrame = default_colors
        self.selected_features_set: List[str] = selected_features_set
        self.feature_shapes = normalize_feature_shape_overrides(feature_shapes)
        self.directional_feature_types: set[str] = resolve_directional_feature_types(
            self.feature_shapes
        )
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.block_fill_color: str = default_colors[default_colors["feature_type"] == "default"]["color"].values[0]
        self.block_stroke_color: str = cfg.objects.features.block_stroke_color
        self.block_stroke_width: float = cfg.objects.features.block_stroke_width.for_length_param(self.length_param)
        self.line_stroke_color: str = cfg.objects.features.line_stroke_color
        self.line_stroke_width: float = cfg.objects.features.line_stroke_width.for_length_param(self.length_param)
        self.qualifier_priority = cfg.labels.filtering.qualifier_priority
        self.blacklist_keywords: List[str] = cfg.labels.filtering.blacklist_keywords


__all__ = ["FeatureDrawingConfigurator"]


