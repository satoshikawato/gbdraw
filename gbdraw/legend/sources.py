"""What the legend reads from the drawing settings: its legend-only sources.

The legend table is a render stage keyed by its arguments, so it takes the few
values its rows show, not the configurators that also carry geometry, ticks,
arrows or label filters. A change to one of those then leaves the legend a hit.
Each source is built from its configurator or slot once per request by the
factory below it; ``ConservationLegendSource`` is the conservation sibling.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.comparisons.identity_scale import IdentityColorDomain
from gbdraw.core.sequence import check_feature_presence
from gbdraw.features.colors import precompute_used_color_rules
from gbdraw.tracks.parsing import SlotView, slot_skew_colors

if TYPE_CHECKING:
    from gbdraw.analysis.depth_tracks import DepthTrackData
    from gbdraw.configurators import (
        BlastMatchConfigurator,
        DepthConfigurator,
        FeatureDrawingConfigurator,
        GcContentConfigurator,
        GcSkewConfigurator,
    )
    from gbdraw.features.placement import ResolvedRecordFeatureInputs


@dataclass(frozen=True, eq=False)
class FeatureLegendSource:
    """The feature rows: colors, block stroke, and what decides which rows exist."""

    color_table: DataFrame | None
    default_colors: DataFrame
    block_stroke_color: str
    block_stroke_width: float
    selected_features: tuple[str, ...]
    feature_visibility_rules: list[dict[str, Any]] | None
    specific_color_rules: Mapping[str, Any]
    default_color_map: Mapping[str, str]
    record_features: tuple[ResolvedRecordFeatureInputs, ...]

    @classmethod
    def of(cls, config: FeatureDrawingConfigurator) -> FeatureLegendSource:
        return cls(
            color_table=config.color_table,
            default_colors=config.default_colors,
            block_stroke_color=config.block_stroke_color,
            block_stroke_width=config.block_stroke_width,
            selected_features=tuple(config.selected_features_set),
            feature_visibility_rules=config.feature_visibility_rules,
            specific_color_rules=config.specific_color_rules,
            default_color_map=config.default_color_map,
            record_features=config.record_features,
        )


def feature_legend_presence(
    records: Sequence[SeqRecord],
    features: FeatureLegendSource,
) -> tuple[list[str], set[tuple[str, str]], set[str]]:
    """The feature types drawn, the color rules features took, and the types on default colors."""

    records = list(records)
    features_present = check_feature_presence(
        records,
        list(features.selected_features),
        feature_visibility_rules=features.feature_visibility_rules,
        specific_color_rules=features.specific_color_rules,
        record_features=features.record_features,
    )
    used_color_rules, default_used_features = precompute_used_color_rules(
        records,
        features.specific_color_rules,
        features.default_color_map,
        set(features.selected_features),
        feature_visibility_rules=features.feature_visibility_rules,
        record_features=features.record_features,
    )
    return features_present, used_color_rules, default_used_features


@dataclass(frozen=True)
class NumericLegendSource:
    """The swatch of a GC content or GC skew track."""

    dinucleotide: str
    high_fill_color: str
    low_fill_color: str
    stroke_color: str
    stroke_width: float

    @classmethod
    def of(cls, config: GcContentConfigurator | GcSkewConfigurator) -> NumericLegendSource:
        return cls(
            dinucleotide=config.dinucleotide,
            high_fill_color=config.high_fill_color,
            low_fill_color=config.low_fill_color,
            stroke_color=config.stroke_color,
            stroke_width=config.stroke_width,
        )

    def for_skew_slot(self, slot: LegendSlotSource, dinucleotide: str) -> NumericLegendSource:
        """The swatch a skew slot draws: its dinucleotide and its own colors."""

        high, low = slot_skew_colors(slot)
        return replace(
            self,
            dinucleotide=str(dinucleotide).upper(),
            high_fill_color=high or self.high_fill_color,
            low_fill_color=low or self.low_fill_color,
        )


@dataclass(frozen=True)
class LegendSources:
    """The feature, GC content and GC skew sources of one drawing's legend."""

    features: FeatureLegendSource
    gc: NumericLegendSource
    skew: NumericLegendSource

    @classmethod
    def of(
        cls,
        feature_config: FeatureDrawingConfigurator,
        gc_config: GcContentConfigurator,
        skew_config: GcSkewConfigurator,
    ) -> LegendSources:
        return cls(
            FeatureLegendSource.of(feature_config),
            NumericLegendSource.of(gc_config),
            NumericLegendSource.of(skew_config),
        )


@dataclass(frozen=True)
class DepthLegendSource:
    """The swatch and name of one depth track (``track_index`` 0 for a lone track)."""

    label: str
    fill_color: str
    stroke_color: str
    stroke_width: float
    track_index: int = 0

    @classmethod
    def of_config(cls, config: DepthConfigurator) -> DepthLegendSource:
        return cls("", config.fill_color, config.stroke_color, config.stroke_width)


def depth_legend_sources(tracks: Iterable[DepthTrackData] | None) -> tuple[DepthLegendSource, ...]:
    return tuple(
        DepthLegendSource(
            track.label,
            track.config.fill_color,
            track.config.stroke_color,
            track.config.stroke_width,
            track.track_index,
        )
        for track in tracks or ()
    )


@dataclass(frozen=True)
class ComparisonLegendSource:
    """The pairwise identity gradient of a comparison: its end colors and its identity domain."""

    min_color: str
    max_color: str
    domain: IdentityColorDomain

    @classmethod
    def of(cls, config: BlastMatchConfigurator) -> ComparisonLegendSource:
        return cls(config.min_color, config.max_color, config.identity_domain)


@dataclass(frozen=True)
class LegendSlotSource:
    """What the legend reads from one enabled track slot."""

    id: str
    renderer: str
    params: Mapping[str, object]


def legend_slot_sources(slots: Iterable[SlotView] | None) -> tuple[LegendSlotSource, ...] | None:
    """The enabled slots in order; ``None`` (no slot list) stays ``None``."""

    if slots is None:
        return None
    return tuple(
        LegendSlotSource(slot.id, slot.renderer, dict(slot.params))
        for slot in slots
        if slot.enabled
    )


__all__ = [
    "ComparisonLegendSource",
    "DepthLegendSource",
    "FeatureLegendSource",
    "LegendSlotSource",
    "LegendSources",
    "NumericLegendSource",
    "depth_legend_sources",
    "feature_legend_presence",
    "legend_slot_sources",
]
