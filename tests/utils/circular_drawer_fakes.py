from collections.abc import Callable, MutableMapping
from types import ModuleType
from typing import Any, get_args

import pytest

from gbdraw.diagrams.circular.builders import CircularSeriesSlot
from gbdraw.render.numeric_series import NumericSeriesStage


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
    """A fake ``add_numeric_series_group_on_canvas`` that records the slot's track width and anchor."""

    def capture(canvas, _stage, _record, _frame, _canvas_config, _config, slot: CircularSeriesSlot):
        if width_key is not None:
            captured[width_key] = slot.track_width
        if norm_key is not None:
            captured[norm_key] = slot.norm_factor
        return canvas

    return capture


def make_numeric_slot_capture(capture_numeric_slot, fallback_group_id: str):
    """A fake ``add_numeric_series_group_on_canvas`` that hands the slot's geometry to ``capture_numeric_slot``."""

    def capture(canvas, _stage, _record, _frame, canvas_config, _config, slot: CircularSeriesSlot):
        capture_numeric_slot(
            str(slot.group_id or fallback_group_id),
            canvas_config,
            slot.track_width,
            slot.norm_factor,
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


def patch_series_groups(
    monkeypatch: pytest.MonkeyPatch, module: ModuleType, **fakes: Callable[..., object]
) -> None:
    """Replace ``add_numeric_series_group_on_canvas`` in ``module`` so each render stage named in ``fakes``
    (``gc_content_content``, ``gc_skew_content``, ``depth_content``) is drawn by its fake, called with the same
    arguments; the other stages draw as before."""

    unknown = set(fakes) - set(get_args(NumericSeriesStage))
    assert not unknown, unknown
    real = module.add_numeric_series_group_on_canvas

    def route(canvas, stage, *args):
        return fakes.get(stage, real)(canvas, stage, *args)

    monkeypatch.setattr(module, "add_numeric_series_group_on_canvas", route)
