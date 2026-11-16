from collections.abc import Callable, MutableMapping
from types import ModuleType
from typing import Any

import pytest


def make_feature_draw_ratio_capture(captured: list[float]):
    def capture(
        _drawer,
        _feature,
        group,
        _total_length,
        _radius,
        _track_ratio,
        track_ratio_factor,
        _track_type,
        _strandedness,
        _length_param,
        **_kwargs,
    ):
        captured.append(float(track_ratio_factor))
        return group

    return capture


def make_numeric_track_capture(
    captured: MutableMapping[str, Any],
    *,
    width_key: str | None = None,
    norm_key: str | None = None,
):
    def capture(
        canvas,
        _record,
        _dataframe,
        _canvas_config,
        _track_config,
        *,
        track_width_override=None,
        norm_factor_override=None,
        **_kwargs,
    ):
        if width_key is not None:
            captured[width_key] = track_width_override
        if norm_key is not None:
            captured[norm_key] = norm_factor_override
        return canvas

    return capture


def make_numeric_slot_capture(capture_numeric_slot, fallback_group_id: str):
    def capture(
        canvas,
        _record,
        _dataframe,
        canvas_config,
        _track_config,
        *,
        track_width_override=None,
        norm_factor_override=None,
        group_id=None,
        **_kwargs,
    ):
        capture_numeric_slot(
            str(group_id or fallback_group_id),
            canvas_config,
            track_width_override,
            norm_factor_override,
        )
        return canvas

    return capture


def make_numeric_slot_geometry_capture(captured: MutableMapping[str, Any]):
    def capture(slot_id, canvas_config, track_width_override, norm_factor_override) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

    return capture

# The render stage of each kind, as the Circular assembly names it.
STAGE_KINDS = {"gc_content_content": "gc_content", "gc_skew_content": "gc_skew", "depth_content": "depth"}


def patch_circular_series_builders(
    monkeypatch: pytest.MonkeyPatch, module: ModuleType, **fakes: Callable[..., object]
) -> None:
    """Replace ``add_numeric_series_group_on_canvas`` in ``module`` so each kind named in ``fakes`` (``gc_content``,
    ``gc_skew``, ``depth``) is drawn by its fake, called with the record's frame, the canvas config and the track's
    configurator, and the slot's fields as keywords; the other kinds draw as before."""

    real = module.add_numeric_series_group_on_canvas

    def dispatch(canvas, stage, gb_record, frame, canvas_config, config, slot):
        fake = fakes.get(STAGE_KINDS[stage])
        if fake is None:
            return real(canvas, stage, gb_record, frame, canvas_config, config, slot)
        return fake(
            canvas,
            gb_record,
            frame,
            canvas_config,
            config,
            group_key=slot.group_key,
            track_width_override=slot.track_width,
            norm_factor_override=slot.norm_factor,
            group_id=slot.group_id,
            slot_group_id=slot.slot_group_id,
            slot_id=slot.slot_id,
            slot_renderer=slot.slot_renderer,
            record_transform=slot.record_transform,
        )

    monkeypatch.setattr(module, "add_numeric_series_group_on_canvas", dispatch)
