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
from typing import TYPE_CHECKING, Any, ClassVar, Literal, Protocol

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.comparisons.identity_scale import IdentityColorDomain
from gbdraw.core.sequence import check_feature_presence
from gbdraw.features.colors import TypeRuleUse, precompute_used_color_rules
from gbdraw.layout.series_style import SeriesPaint
from gbdraw.tracks.parsing import SlotView, slot_skew_fill_colors, slots_draw_features

if TYPE_CHECKING:
    from gbdraw.annotations.models import AnnotationTrackParams
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
    slots: Iterable[SlotView] | None,
) -> tuple[list[str], set[TypeRuleUse], set[str]]:
    """The feature types drawn, the color rules features of each type took, and the types on default colors.

    ``slots`` is the drawing's track stack; without a features slot
    (:func:`slots_draw_features`) no feature is drawn, so none is named.
    """

    if not slots_draw_features(slots):
        return [], set(), set()
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
class DinucleotideLegendSource:
    """The rows of a GC content or GC skew track: named ``<dinucleotide> content`` or ``<dinucleotide> skew``
    (``row_name``), in the paint the track draws with (its configurator's ``series_paint``). A paint of one fill
    shows one row; one whose fill below its baseline differs (GC skew) shows a ``(+)`` and a ``(-)`` row. Each kind
    is a type of its own (``GcContentLegendSource``, ``GcSkewLegendSource``), so a source of one kind passed where
    the other is expected fails the type check."""

    dinucleotide: str
    paint: SeriesPaint
    kind: ClassVar[Literal["content", "skew"]]

    def row_name(self, dinucleotide: str) -> str:
        """The name of the track's rows for a slot of ``dinucleotide``, before its ``legend_label``."""

        return f"{dinucleotide} {self.kind}"

    @property
    def name(self) -> str:
        """The name of the track's rows without a slot of its own."""

        return self.row_name(self.dinucleotide)

    def slot_paint(self, slot: LegendSlotSource) -> SeriesPaint:
        """The paint a slot of this track draws with."""

        return self.paint


@dataclass(frozen=True)
class GcContentLegendSource(DinucleotideLegendSource):
    """The rows of a GC content track; every slot draws in the track's paint."""

    kind: ClassVar[Literal["content", "skew"]] = "content"

    @classmethod
    def of(cls, config: GcContentConfigurator) -> GcContentLegendSource:
        return cls(config.dinucleotide, config.series_paint)


@dataclass(frozen=True)
class GcSkewLegendSource(DinucleotideLegendSource):
    """The rows of a GC skew track; a slot draws in its own fills (``slot_skew_fill_colors``)."""

    kind: ClassVar[Literal["content", "skew"]] = "skew"

    @classmethod
    def of(cls, config: GcSkewConfigurator) -> GcSkewLegendSource:
        return cls(config.dinucleotide, config.series_paint)

    def slot_paint(self, slot: LegendSlotSource) -> SeriesPaint:
        """The paint a slot of this track draws with: its own fills (``slot_skew_fill_colors``)."""

        fill, below_fill = slot_skew_fill_colors(slot, self.paint.fill, self.paint.low_fill)
        return replace(self.paint, fill=fill, below_fill=below_fill)


@dataclass(frozen=True)
class DepthLegendSource:
    """The row of a depth track, in its colour (its configurator's ``legend_paint``: one row whatever its style and
    negative colour). ``label`` and ``track_index`` name one of several logical tracks (``depth_track_legend_name``);
    a lone track's row is ``Depth`` (``lone``)."""

    label: str
    track_index: int
    paint: SeriesPaint

    @classmethod
    def of(cls, config: DepthConfigurator, label: str, track_index: int) -> DepthLegendSource:
        return cls(label, track_index, config.legend_paint)

    @classmethod
    def lone(cls, config: DepthConfigurator) -> DepthLegendSource:
        """The source of a record's only depth track, which has no label of its own."""

        return cls.of(config, "", 0)


def depth_legend_sources(tracks: Iterable[DepthTrackData] | None) -> tuple[DepthLegendSource, ...]:
    """The legend sources of the drawn depth tracks, each with its name and ring."""

    return tuple(DepthLegendSource.of(track.config, track.label, track.track_index) for track in tracks or ())


@dataclass(frozen=True)
class FeatureDensityLegendSource:
    """The row of one feature density track: its name (the query in words, or its caption), its style, and its
    colors and scale. A heatmap's row is a gradient from ``low_color`` at 0 to ``color`` at ``max_text``; a
    histogram's a swatch of ``color``."""

    track_id: str
    caption: str
    style: str
    color: str
    low_color: str
    max_value: float
    max_text: str
    stroke_color: str
    stroke_width: float


@dataclass(frozen=True)
class LegendSources:
    """The feature, GC content, GC skew and feature density sources of one drawing's legend."""

    features: FeatureLegendSource
    gc: GcContentLegendSource
    skew: GcSkewLegendSource
    feature_density: tuple[FeatureDensityLegendSource, ...] = ()

    @classmethod
    def of(
        cls,
        feature_config: FeatureDrawingConfigurator,
        gc_config: GcContentConfigurator,
        skew_config: GcSkewConfigurator,
        feature_density: Sequence[FeatureDensityLegendSource] = (),
    ) -> LegendSources:
        return cls(
            FeatureLegendSource.of(feature_config),
            GcContentLegendSource.of(gc_config),
            GcSkewLegendSource.of(skew_config),
            tuple(feature_density),
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

    @property
    def enabled(self) -> bool:
        """Always true: :func:`legend_slot_sources` keeps only the enabled slots."""

        return True


class SidedSlotView(SlotView, Protocol):
    """A normalized track slot, in either mode: its side and, for an annotation slot, its parsed params."""

    @property
    def side(self) -> str: ...

    @property
    def annotation(self) -> AnnotationTrackParams | None: ...


# The rank of each side in the Legend: the outer side, the Axis, the inner side. A slot other than an annotation
# overlay on side ``overlay`` lies on the Axis.
_SIDE_RANKS = {"outside": 0, "above": 0, "overlay": 1, "inside": 2, "below": 2}


def legend_slot_sources(slots: Iterable[SidedSlotView] | None) -> tuple[LegendSlotSource, ...] | None:
    """The enabled slots in the order the Legend lists their entries; ``None`` (no slot list) stays ``None``.

    The order is the order the rows lie in from the outer edge, in both modes: the outer side (outside, above), the
    Axis, the inner side (inside, below), each in the written order. An annotation overlay comes right after its
    anchor, and the overlays on one anchor in their written order. So the Legend reads the same order whatever the
    order of the sides in the list, as the Web stack holds it (OV-3187).
    """

    if slots is None:
        return None
    enabled = [slot for slot in slots if slot.enabled]
    place = {slot.id: index for index, slot in enumerate(enabled)}

    def row(slot: SidedSlotView) -> tuple[int, int]:
        return _SIDE_RANKS[slot.side], place[slot.id]

    def key(slot: SidedSlotView) -> tuple[int, int, int]:
        anchor_id = slot.annotation.anchor_slot if slot.annotation is not None and slot.side == "overlay" else None
        if anchor_id is None or anchor_id not in place:
            return (*row(slot), -1)
        return (*row(enabled[place[anchor_id]]), place[slot.id])

    return tuple(LegendSlotSource(slot.id, slot.renderer, dict(slot.params)) for slot in sorted(enabled, key=key))


__all__ = [
    "ComparisonLegendSource",
    "DepthLegendSource",
    "DinucleotideLegendSource",
    "FeatureDensityLegendSource",
    "FeatureLegendSource",
    "GcContentLegendSource",
    "GcSkewLegendSource",
    "LegendSlotSource",
    "LegendSources",
    "SidedSlotView",
    "depth_legend_sources",
    "feature_legend_presence",
    "legend_slot_sources",
]
