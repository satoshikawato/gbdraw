from __future__ import annotations

from gbdraw.config.models import LinearRenderProfile
from gbdraw.layout.linear import (
    LinearRecordRenderContext,
    LinearResolvedTrack,
    LinearTrackLayout,
)


def make_linear_record_render_context(
    profile: LinearRenderProfile,
    *,
    depth_available: bool = False,
) -> LinearRecordRenderContext:
    """Build the resolved feature-track context used by focused unit tests."""

    side = "overlay" if profile.track_layout == "middle" else profile.track_layout
    feature_track = LinearResolvedTrack(
        slot_index=0,
        id="features",
        renderer="features",
        side=side,
        y_offset=0.0,
        height=0.0,
        spacing_after_px=0.0,
        top_extent=0.0,
        bottom_extent=0.0,
        z=0,
        params={"nested": {"values": [1]}},
    )
    track_layout = LinearTrackLayout(
        slots=[feature_track],  # type: ignore[arg-type]
        top_extent=0.0,
        bottom_extent=0.0,
        plot_tracks_height=0.0,
        plot_tracks_visual_bottom=0.0,
        depth_track_offsets=[],  # type: ignore[arg-type]
        depth_track_heights=[],  # type: ignore[arg-type]
        gc_content_track_offset=0.0,
        gc_skew_track_offset=0.0,
    )
    return LinearRecordRenderContext(
        profile=profile,
        track_layout=track_layout,
        depth_available=depth_available,
    )


__all__ = ["make_linear_record_render_context"]
