from __future__ import annotations

import pytest

from gbdraw.diagrams.circular.radial_layout import build_circular_feature_layout


class _Feature:
    def __init__(self, track_id: int) -> None:
        self.feature_track_id = track_id


def test_middle_combined_feature_center_is_axis_radius() -> None:
    layout = build_circular_feature_layout(
        {"a": _Feature(0)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="middle",
        strandedness=False,
    )

    assert layout is not None
    assert layout.lanes_by_track_id[0].center_px == pytest.approx(390.0)


def test_spreadout_and_tuckin_feature_lanes_stay_on_expected_side() -> None:
    spreadout = build_circular_feature_layout(
        {"a": _Feature(0), "b": _Feature(1)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="spreadout",
        strandedness=False,
    )
    tuckin = build_circular_feature_layout(
        {"a": _Feature(0), "b": _Feature(1)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="tuckin",
        strandedness=False,
    )

    assert spreadout is not None
    assert tuckin is not None
    assert all(lane.inner_px > 390.0 for lane in spreadout.lanes_by_track_id.values())
    assert all(lane.outer_px < 390.0 for lane in tuckin.lanes_by_track_id.values())
