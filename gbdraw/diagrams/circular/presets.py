"""Preset expansion for circular track layouts.

The public ``track_type`` value is treated as a preset name at this boundary.
Downstream layout code consumes explicit slot geometry/parameters rather than
branching on the legacy preset name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...tracks.circular import CircularTrackSlot  # type: ignore[reportMissingImports]
from ...tracks.spec import CircularTrackPlacement, ScalarSpec  # type: ignore[reportMissingImports]
from ...svg.circular_ticks import (  # type: ignore[reportMissingImports]
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_ratio_bounds,
)


CircularTrackPreset = Literal["tuckin", "middle", "spreadout"]
CircularFeatureLaneDirection = Literal["inside", "split", "outside"]

_VALID_PRESETS: frozenset[str] = frozenset({"tuckin", "middle", "spreadout"})


@dataclass(frozen=True)
class CircularPresetContext:
    cfg: GbdrawConfig
    canvas_config: CircularCanvasConfigurator
    total_length: int
    strandedness: bool
    show_features: bool
    show_ticks: bool
    show_depth: bool
    show_gc: bool
    show_skew: bool
    dinucleotide: str = "GC"


@dataclass(frozen=True)
class CircularFeatureSlotDefaults:
    lane_direction: CircularFeatureLaneDirection
    placement: CircularTrackPlacement | None
    width: ScalarSpec | None


@dataclass(frozen=True)
class CircularLabelArenaDefaults:
    preset: CircularTrackPreset


def normalize_circular_track_preset(raw: str | None) -> CircularTrackPreset:
    preset = str(raw or "tuckin").strip().lower()
    if preset not in _VALID_PRESETS:
        raise ValueError(
            "Circular track preset must be one of: tuckin, middle, spreadout"
        )
    return preset  # type: ignore[return-value]


def circular_feature_lane_direction_for_preset(
    preset: str | None,
) -> CircularFeatureLaneDirection:
    normalized = normalize_circular_track_preset(preset)
    if normalized == "middle":
        return "split"
    if normalized == "spreadout":
        return "outside"
    return "inside"


def _scalar_factor(value: float) -> ScalarSpec:
    return ScalarSpec(float(value), "factor")


def _scalar_px(value: float) -> ScalarSpec:
    return ScalarSpec(float(value), "px")


def _default_feature_width_px(context: CircularPresetContext) -> float:
    length_param = str(context.canvas_config.length_param)
    return (
        float(context.canvas_config.radius)
        * float(context.canvas_config.track_ratio)
        * float(context.cfg.canvas.circular.track_ratio_factors[length_param][0])
    )


def _default_numeric_width_px(
    renderer: str,
    context: CircularPresetContext,
) -> float:
    length_param = str(context.canvas_config.length_param)
    base = float(context.canvas_config.radius) * float(context.canvas_config.track_ratio)
    if renderer == "depth":
        return base * float(context.cfg.canvas.circular.track_ratio_factors[length_param][1]) * 0.5
    if renderer == "dinucleotide_skew":
        return base * float(context.cfg.canvas.circular.track_ratio_factors[length_param][2])
    return base * float(context.cfg.canvas.circular.track_ratio_factors[length_param][1])


def _numeric_slot(
    *,
    slot_id: str,
    renderer: str,
    context: CircularPresetContext,
    params: dict[str, object] | None = None,
) -> CircularTrackSlot:
    del context
    return CircularTrackSlot(
        id=slot_id,
        renderer=renderer,
        params={"_preset_generated": True, **dict(params or {})},
    )


def circular_feature_slot_defaults_for_preset(
    preset: str,
    context: CircularPresetContext,
) -> CircularFeatureSlotDefaults:
    normalized = normalize_circular_track_preset(preset)
    lane_direction = circular_feature_lane_direction_for_preset(normalized)
    return CircularFeatureSlotDefaults(
        lane_direction=lane_direction,
        placement=CircularTrackPlacement(radius=_scalar_factor(1.0)),
        width=_scalar_px(_default_feature_width_px(context)),
    )


def _tick_slot_for_preset(
    preset: CircularTrackPreset,
    context: CircularPresetContext,
) -> CircularTrackSlot:
    base_radius = float(context.canvas_config.radius)
    tick_inner_ratio, tick_outer_ratio = get_circular_tick_path_ratio_bounds(
        int(context.total_length),
        preset,
        bool(context.strandedness),
    )
    tick_inner_ratio, tick_outer_ratio = sorted((float(tick_inner_ratio), float(tick_outer_ratio)))
    tick_side = "inside" if tick_outer_ratio <= 1.0 else "outside"
    anchor_ratio = tick_outer_ratio if tick_side == "inside" else tick_inner_ratio
    tick_width_px = max(0.0, (tick_outer_ratio - tick_inner_ratio) * base_radius)

    label_side = "outside"
    label_bounds = get_circular_tick_label_radius_bounds(
        center_radius_px=anchor_ratio * base_radius,
        total_len=int(context.total_length),
        track_type=preset,
        strandedness=bool(context.strandedness),
        font_size=float(context.cfg.objects.ticks.tick_labels.font_size),
        font_family=str(context.cfg.objects.text.font_family),
        dpi=int(context.canvas_config.dpi),
        manual_interval=context.cfg.objects.scale.interval,
        tick_width=float(context.cfg.objects.ticks.tick_width),
        tick_side=tick_side,
        tick_length_px=tick_width_px,
        length_reference_radius_px=base_radius,
    )
    if label_bounds is not None:
        label_center = (float(label_bounds[0]) + float(label_bounds[1])) / 2.0
        label_side = "inside" if label_center < (anchor_ratio * base_radius) else "outside"

    return CircularTrackSlot(
        id="ticks",
        renderer="ticks",
        placement=CircularTrackPlacement(radius=_scalar_factor(anchor_ratio)),
        width=_scalar_px(tick_width_px),
        params={
            "_preset_generated": True,
            "tick_side": tick_side,
            "label_side": label_side,
            "preset": preset,
        },
    )


def circular_track_slots_for_preset(
    preset: str,
    context: CircularPresetContext,
) -> list[CircularTrackSlot]:
    normalized = normalize_circular_track_preset(preset)
    slots: list[CircularTrackSlot] = []
    nt = str(context.dinucleotide or "GC").upper()

    if context.show_features:
        feature_defaults = circular_feature_slot_defaults_for_preset(normalized, context)
        slots.append(
            CircularTrackSlot(
                id="features",
                renderer="features",
                placement=feature_defaults.placement,
                width=feature_defaults.width,
                params={
                    "_preset_generated": True,
                    "lane_direction": feature_defaults.lane_direction,
                },
            )
        )
    if context.show_ticks:
        slots.append(_tick_slot_for_preset(normalized, context))
    if context.show_depth:
        slots.append(
            _numeric_slot(
                slot_id="depth",
                renderer="depth",
                context=context,
            )
        )
    if context.show_gc:
        slots.append(
            _numeric_slot(
                slot_id="gc_content",
                renderer="dinucleotide_content",
                context=context,
                params={"nt": nt},
            )
        )
    if context.show_skew:
        slots.append(
            _numeric_slot(
                slot_id="gc_skew",
                renderer="dinucleotide_skew",
                context=context,
                params={"nt": nt},
            )
        )
    return slots


def circular_label_arena_defaults_for_preset(
    preset: str,
    context: CircularPresetContext,
) -> CircularLabelArenaDefaults:
    del context
    return CircularLabelArenaDefaults(preset=normalize_circular_track_preset(preset))


__all__ = [
    "CircularFeatureLaneDirection",
    "CircularFeatureSlotDefaults",
    "CircularLabelArenaDefaults",
    "CircularPresetContext",
    "CircularTrackPreset",
    "circular_feature_lane_direction_for_preset",
    "circular_feature_slot_defaults_for_preset",
    "circular_label_arena_defaults_for_preset",
    "circular_track_slots_for_preset",
    "normalize_circular_track_preset",
]
