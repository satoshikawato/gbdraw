"""Preset expansion for circular track layouts.

The public ``track_type`` value is treated as a preset name at this boundary.
Downstream layout code consumes explicit slot geometry/parameters rather than
branching on the legacy preset name.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal, Sequence

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...tracks.circular import (  # type: ignore[reportMissingImports]
    CircularTrackSlot,
    NUMERIC_CIRCULAR_TRACK_RENDERERS,
)
from ...tracks.scalars import ScalarSpec  # type: ignore[reportMissingImports]
from ...svg.circular_ticks import (  # type: ignore[reportMissingImports]
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_ratio_bounds,
)


CircularTrackPreset = Literal["tuckin", "middle", "spreadout"]
CircularFeatureLaneDirection = Literal["inside", "split", "outside"]

_VALID_PRESETS: frozenset[str] = frozenset({"tuckin", "middle", "spreadout"})
_NUMERIC_PRESET_SLOT_IDS: frozenset[str] = frozenset({"depth", "gc_content", "gc_skew"})
_NON_NUMERIC_PRESET_SLOT_IDS: frozenset[str] = frozenset({"features", "ticks"})
_RENDERER_ALIASES: dict[str, str] = {
    "gc_content": "dinucleotide_content",
    "content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "skew": "dinucleotide_skew",
}


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
    radius: ScalarSpec | None
    width: ScalarSpec | None


@dataclass(frozen=True)
class CircularLabelArenaDefaults:
    preset: CircularTrackPreset


@dataclass(frozen=True)
class CircularPresetRadialPlan:
    slots: tuple[CircularTrackSlot, ...]
    preferred_anchor_slot_ids: frozenset[str]


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


def _normalized_renderer(raw: object) -> str:
    renderer = str(raw).strip().lower()
    return _RENDERER_ALIASES.get(renderer, renderer)


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
    preset: CircularTrackPreset,
    context: CircularPresetContext,
    params: dict[str, object] | None = None,
) -> CircularTrackSlot:
    track_id_key = {
        "depth": "depth_track",
        "gc_content": "gc_track",
        "gc_skew": "skew_track",
    }.get(slot_id)
    track_id = (
        context.canvas_config.track_ids.get(track_id_key)
        if track_id_key is not None
        else None
    )
    length_param = str(context.canvas_config.length_param)
    radius = None
    if track_id is not None:
        radius = _scalar_factor(
            float(context.cfg.canvas.circular.track_dict[length_param][preset][str(track_id)])
        )
    return CircularTrackSlot(
        id=slot_id,
        renderer=renderer,
        side="inside",
        radius=radius,
        width=_scalar_px(_default_numeric_width_px(renderer, context)),
        spacing=_scalar_px(max(1.0, 0.01 * float(context.canvas_config.radius))),
        params=dict(params or {}),
    )


def circular_feature_slot_defaults_for_preset(
    preset: str,
    context: CircularPresetContext,
) -> CircularFeatureSlotDefaults:
    normalized = normalize_circular_track_preset(preset)
    lane_direction = circular_feature_lane_direction_for_preset(normalized)
    return CircularFeatureSlotDefaults(
        lane_direction=lane_direction,
        radius=_scalar_factor(1.0),
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
        side=tick_side,
        radius=_scalar_factor(anchor_ratio),
        width=_scalar_px(tick_width_px),
        spacing=_scalar_px(max(1.0, 0.01 * base_radius)),
        params={
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
                side="overlay" if feature_defaults.lane_direction == "split" else feature_defaults.lane_direction,
                radius=feature_defaults.radius,
                width=feature_defaults.width,
                spacing=_scalar_px(max(1.0, 0.01 * float(context.canvas_config.radius))),
                reserve=True if feature_defaults.lane_direction == "split" else None,
                params={
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
                preset=normalized,
                context=context,
            )
        )
    if context.show_gc:
        slots.append(
            _numeric_slot(
                slot_id="gc_content",
                renderer="dinucleotide_content",
                preset=normalized,
                context=context,
                params={"nt": nt},
            )
        )
    if context.show_skew:
        slots.append(
            _numeric_slot(
                slot_id="gc_skew",
                renderer="dinucleotide_skew",
                preset=normalized,
                context=context,
                params={"nt": nt},
            )
        )
    return slots


def circular_radial_plan_for_preset(
    preset: str,
    context: CircularPresetContext,
) -> CircularPresetRadialPlan:
    slots = tuple(circular_track_slots_for_preset(preset, context))
    preferred_ids = frozenset(
        slot.id
        for slot in slots
        if slot.renderer in {"depth", "dinucleotide_content", "dinucleotide_skew"}
        and slot.side == "inside"
        and slot.radius is not None
    )
    return CircularPresetRadialPlan(
        slots=slots,
        preferred_anchor_slot_ids=preferred_ids,
    )


def _slot_requests_pure_auto(slot: CircularTrackSlot, renderer: str) -> bool:
    """Return True for the legacy explicit auto-packed numeric/depth shape."""

    return (
        renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
        and slot.radius is None
        and slot.compress is True
        and (slot.side is None or str(slot.side).strip().lower() == "inside")
    )


def _slot_uses_builtin_preset_lane(slot: CircularTrackSlot, renderer: str) -> bool:
    if renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
        return str(slot.id) in _NUMERIC_PRESET_SLOT_IDS
    if renderer in {"features", "ticks"}:
        return str(slot.id) in _NON_NUMERIC_PRESET_SLOT_IDS
    return False


def _slot_is_blank_unmatched_numeric_duplicate(slot: CircularTrackSlot, renderer: str) -> bool:
    return (
        renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
        and str(slot.id) not in _NUMERIC_PRESET_SLOT_IDS
        and slot.enabled
        and slot.side is None
        and slot.radius is None
        and slot.width is None
        and slot.spacing is None
        and slot.strict is None
        and slot.compress is None
        and slot.reserve is None
    )


def _inherited_params_for_slot(
    slot: CircularTrackSlot,
    preset_slot: CircularTrackSlot | None,
) -> dict[str, object]:
    inherited = dict(preset_slot.params or {}) if preset_slot is not None else {}
    explicit = dict(slot.params or {})

    # Programmatic callers may still use accepted aliases. Do not let an
    # inherited canonical key hide an explicit alias supplied by the caller.
    if "lanes" in explicit and "lane_direction" not in explicit:
        inherited.pop("lane_direction", None)
    if "dinucleotide" in explicit and "nt" not in explicit:
        inherited.pop("nt", None)

    inherited.update(explicit)
    return inherited


def _inherited_width_for_renderer(
    renderer: str,
    params_slot: CircularTrackSlot | None,
    context: CircularPresetContext,
) -> ScalarSpec | None:
    if renderer == "features":
        return _scalar_px(_default_feature_width_px(context))
    if renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS:
        return _scalar_px(_default_numeric_width_px(renderer, context))
    if params_slot is not None:
        return params_slot.width
    return None


def _overlay_slot_on_preset_lane(
    slot: CircularTrackSlot,
    *,
    renderer: str,
    geometry_slot: CircularTrackSlot | None,
    params_slot: CircularTrackSlot | None,
    context: CircularPresetContext,
) -> CircularTrackSlot:
    params = _inherited_params_for_slot(slot, params_slot)
    side = slot.side
    if side is None and not (
        renderer == "features"
        and ("lane_direction" in params or "lanes" in params)
    ):
        side = geometry_slot.side if geometry_slot is not None else None
    return replace(
        slot,
        renderer=renderer,
        side=side,
        radius=slot.radius if slot.radius is not None else (geometry_slot.radius if geometry_slot is not None else None),
        width=slot.width if slot.width is not None else _inherited_width_for_renderer(renderer, params_slot, context),
        spacing=slot.spacing if slot.spacing is not None else (geometry_slot.spacing if geometry_slot is not None else None),
        reserve=slot.reserve if slot.reserve is not None else (params_slot.reserve if params_slot is not None else None),
        params=params,
    )


def circular_track_slots_from_preset_order(
    slots: Sequence[CircularTrackSlot],
    preset: str,
    context: CircularPresetContext,
) -> CircularPresetRadialPlan:
    """Overlay user slot order/overrides onto record-local preset defaults.

    The returned slots are transient render-time inputs. The caller's slot
    objects remain unchanged, so preset-derived geometry is not persisted back
    into web/session/API state.
    """

    normalized_preset = normalize_circular_track_preset(preset)
    legacy_auto_mode = any(
        slot.enabled
        and _slot_requests_pure_auto(slot, _normalized_renderer(slot.renderer))
        and _slot_uses_builtin_preset_lane(slot, _normalized_renderer(slot.renderer))
        for slot in slots
    )
    if legacy_auto_mode:
        return CircularPresetRadialPlan(
            slots=tuple(
                replace(slot, renderer=_normalized_renderer(slot.renderer))
                for slot in slots
            ),
            preferred_anchor_slot_ids=frozenset(),
        )

    preset_slots = tuple(circular_track_slots_for_preset(normalized_preset, context))
    preset_by_id = {str(slot.id): slot for slot in preset_slots}
    preset_by_renderer: dict[str, CircularTrackSlot] = {}
    for preset_slot in preset_slots:
        renderer = _normalized_renderer(preset_slot.renderer)
        preset_by_renderer.setdefault(renderer, preset_slot)

    renderer_has_blank_unmatched_duplicates: dict[str, bool] = {}
    for slot in slots:
        if slot.enabled:
            renderer = _normalized_renderer(slot.renderer)
            if _slot_is_blank_unmatched_numeric_duplicate(slot, renderer):
                renderer_has_blank_unmatched_duplicates[renderer] = True

    geometry_lane_index = 0
    layout_slots: list[CircularTrackSlot] = []
    preferred_ids: set[str] = set()

    for slot in slots:
        renderer = _normalized_renderer(slot.renderer)
        geometry_slot: CircularTrackSlot | None = None
        params_slot: CircularTrackSlot | None = None
        inherited_radius = False

        has_extra_numeric_renderer = bool(renderer_has_blank_unmatched_duplicates.get(renderer, False))
        if (
            slot.enabled
            and not has_extra_numeric_renderer
            and not _slot_requests_pure_auto(slot, renderer)
            and _slot_uses_builtin_preset_lane(slot, renderer)
        ):
            if geometry_lane_index < len(preset_slots):
                geometry_slot = preset_slots[geometry_lane_index]
                geometry_lane_index += 1
            params_slot = preset_by_id.get(str(slot.id), preset_by_renderer.get(renderer))

        if geometry_slot is not None and slot.radius is None and geometry_slot.radius is not None:
            inherited_radius = True

        overlaid = _overlay_slot_on_preset_lane(
            slot,
            renderer=renderer,
            geometry_slot=geometry_slot,
            params_slot=params_slot,
            context=context,
        )
        layout_slots.append(overlaid)

        if (
            inherited_radius
            and renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
            and overlaid.radius is not None
            and str(overlaid.side or "inside").strip().lower() == "inside"
        ):
            preferred_ids.add(str(overlaid.id))

    return CircularPresetRadialPlan(
        slots=tuple(layout_slots),
        preferred_anchor_slot_ids=frozenset(preferred_ids),
    )


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
    "CircularPresetRadialPlan",
    "CircularTrackPreset",
    "circular_radial_plan_for_preset",
    "circular_feature_lane_direction_for_preset",
    "circular_feature_slot_defaults_for_preset",
    "circular_label_arena_defaults_for_preset",
    "circular_track_slots_from_preset_order",
    "circular_track_slots_for_preset",
    "normalize_circular_track_preset",
]
