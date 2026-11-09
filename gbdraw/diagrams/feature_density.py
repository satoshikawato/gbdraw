"""What a feature density slot draws in either mode: its layer by the slot's track id, the histogram's depth
frame and settings, its color, and its legend rows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from ..analysis.feature_density import HEATMAP_LOW_COLOR, FeatureDensityLayer, FeatureDensityTrack, feature_density_layers
from ..config.models import RenderProfile
from ..configurators import DepthConfigurator
from ..io.colors import resolve_color_to_hex
from ..layout.record_coordinates import RecordDisplayTransform
from ..layout.record_transform import spans_may_cross_origin
from ..legend.sources import FeatureDensityLegendSource
from ..tracks.parsing import SlotWithParams, feature_density_slot_track_id


def slot_feature_density_layer(
    slot_id: str, slot: SlotWithParams, layers: Mapping[str, FeatureDensityLayer] | None
) -> FeatureDensityLayer | None:
    """The layer a ``feature_density`` slot draws: its track's on this record, ``None`` when the record has none."""

    return (layers or {}).get(feature_density_slot_track_id(slot_id, slot.params))


def feature_density_color(layer: FeatureDensityLayer, profile: RenderProfile) -> str:
    """The track's own color, else the depth track's fill color."""

    return str(layer.track.color or profile.config.objects.depth.fill_color)


def feature_density_histogram_source(
    layer: FeatureDensityLayer, profile: RenderProfile
) -> tuple[DataFrame, DepthConfigurator]:
    """A histogram's values as the depth drawers read them: the window values from 0 to the shared top of the
    scale, without an axis (its Legend row names the top)."""

    frame = DataFrame({"position": layer.values.samples, "depth": layer.values.values})
    config = DepthConfigurator(
        window=layer.window,
        step=layer.step,
        profile=profile,
        fill_color=feature_density_color(layer, profile),
        min_depth=0.0,
        max_depth=float(layer.scale_max),
        normalize=False,
        show_axis=False,
        show_ticks=False,
    )
    return frame, config


def feature_density_legend(
    layers: Mapping[str, FeatureDensityLayer] | None, profile: RenderProfile
) -> tuple[FeatureDensityLegendSource, ...]:
    """The legend rows of a record's feature density tracks, from their drawn layers and colors."""

    depth = profile.config.objects.depth
    return tuple(
        FeatureDensityLegendSource(
            track_id=layer.track.id,
            caption=layer.caption(),
            style=layer.track.style,
            color=resolve_color_to_hex(feature_density_color(layer, profile)),
            low_color=HEATMAP_LOW_COLOR,
            max_value=float(layer.scale_max),
            max_text=layer.scale_text(layer.scale_max),
            stroke_color=depth.stroke_color,
            stroke_width=depth.stroke_width,
        )
        for layer in (layers or {}).values()
    )


def feature_density_record_layers(
    records: Sequence[SeqRecord],
    tracks: Sequence[FeatureDensityTrack] | None,
    *,
    windows: Sequence[tuple[int, int]],
    transforms: Sequence[RecordDisplayTransform],
) -> list[dict[str, FeatureDensityLayer]]:
    """Each record's feature density layers, on the scale the records drawn together share; a record's windows
    wrap round its origin when its spans may (``spans_may_cross_origin``)."""

    if not tracks:
        return [{} for _ in records]
    return feature_density_layers(
        records,
        tracks,
        windows=windows,
        wraps=[spans_may_cross_origin(record, transform) for record, transform in zip(records, transforms)],
    )
