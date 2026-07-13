from collections.abc import MutableMapping
from typing import Any


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
        _config_dict,
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
        _config_dict,
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
