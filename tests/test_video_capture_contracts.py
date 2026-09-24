"""Timing and edit contracts for the recorded gbdraw walkthrough videos."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

CAPTURE_ROOT = Path(__file__).resolve().parents[1] / 'docs' / 'capture'
sys.path.insert(0, str(CAPTURE_ROOT))
from video.walkthrough_render import FPS, HIGHLIGHTS, Camera, TimeMap, _moments  # noqa: E402


def test_highlights_cut_only_at_marks_the_journey_records() -> None:
    source = (CAPTURE_ROOT / 'video' / 'walkthrough.py').read_text(encoding='utf-8')
    marks = set(re.findall(r'rec\.mark\("([^"]+)"\)', source))
    if 'rec.mark(f"rule-{index + 1}")' in source:
        marks |= {f'rule-{index}' for index in range(1, 5)}
    wanted = {name for first, _, last, _ in HIGHLIGHTS for name in (first, last)}
    assert wanted <= marks, sorted(wanted - marks)


def test_walkthrough_time_map_compresses_only_marked_waits() -> None:
    events = [
        {'t': 2.0, 'kind': 'speed_start', 'factor': None, 'target': 1.0},
        {'t': 10.0, 'kind': 'speed_end'},
        {'t': 12.0, 'kind': 'speed_start', 'factor': 3.0, 'target': None},
        {'t': 15.0, 'kind': 'speed_end'},
    ]
    mapping = TimeMap.build(events, 0.0, 20.0)
    assert mapping.out(2.0) == pytest.approx(2.0)
    assert mapping.out(10.0) == pytest.approx(3.0)  # 8 s wait shown in 1 s
    assert mapping.out(15.0) == pytest.approx(6.0)  # 3 s at 3x
    assert mapping.duration == pytest.approx(11.0)
    for moment in (0.5, 4.0, 11.0, 13.5, 19.0):
        assert mapping.source(mapping.out(moment)) == pytest.approx(moment)
    assert mapping.factor(mapping.out(4.0)) == pytest.approx(8.0)


def test_walkthrough_camera_eases_and_stays_inside_the_page() -> None:
    events = [
        {'t': 0.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0},
        {'t': 1.0, 'kind': 'camera', 'cx': 10, 'cy': 10, 'zoom': 2.0, 'duration': 1.0},
    ]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), 16 / 9, False)
    assert camera.at(0.5) == pytest.approx((960, 540, 1080))
    assert camera.at(3.0) == pytest.approx((480, 270, 540))  # clamped to the top-left corner
    x, y, height = camera.at(1.5)
    assert 480 < x < 960 and 270 < y < 540 and 540 < height < 1080


def test_portrait_camera_frames_the_diagram_for_whole_page_shots() -> None:
    page = {'x': 0, 'y': 0, 'width': 1920, 'height': 1080}
    diagram = {'x': 900, 'y': 150, 'width': 800, 'height': 800}
    button = {'x': 1500, 'y': 70, 'width': 60, 'height': 20}
    events = [
        {'t': 0.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0,
         'rect': page, 'subject': diagram, 'portrait': None},
        {'t': 2.0, 'kind': 'camera', 'cx': 1530, 'cy': 80, 'zoom': 2.1, 'duration': 0,
         'rect': button, 'subject': diagram, 'portrait': None},
    ]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), 936 / 1060, True)
    x, y, height = camera.at(1.0)
    assert (x, y) == pytest.approx((1300, 540)) and height == pytest.approx(1080)
    x, y, height = camera.at(3.0)
    assert height == pytest.approx(1080 / 2.6)  # small targets stop at the maximum zoom
    assert y == pytest.approx(height / 2)       # and stay inside the page


def test_highlight_moments_skim_waits_but_not_fast_forwards() -> None:
    events = [
        {'t': 1.0, 'kind': 'speed_start', 'factor': None, 'target': 1.0},
        {'t': 5.0, 'kind': 'speed_end'},
        {'t': 6.0, 'kind': 'speed_start', 'factor': 3.0, 'target': None},
        {'t': 9.0, 'kind': 'speed_end'},
    ]
    mapping = TimeMap.build(events, 0.0, 12.0)
    assert [kind for *_, kind in mapping.segments] == ['wait', 'fast']
    moments = _moments(mapping, 0.0, mapping.duration, 1.0, 2.0)
    # 1 s + 1 s wait at 2x + 1 s + 1 s fast-forward + 3 s at normal pace.
    assert len(moments) / FPS == pytest.approx(0.5 + 6.0, abs=0.05)
    assert moments == tuple(sorted(moments))
