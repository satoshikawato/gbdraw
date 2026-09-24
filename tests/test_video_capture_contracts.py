"""Input-boundary checks for the reproducible Meet gbdraw video."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'docs' / 'capture'))
from video.model import contained_file, load_assets, load_storyboard, sha256  # noqa: E402


STORYBOARD = Path(__file__).resolve().parents[1] / 'docs/videos/meet-gbdraw/storyboard.json'


def test_approved_storyboard_has_13_scenes_and_1140_frames() -> None:
    rows = load_storyboard(STORYBOARD)['scenes']
    assert len(rows) == 13
    assert rows[-1]['start_frame'] + rows[-1]['frames'] == 1140


@pytest.mark.parametrize('mutation', [
    lambda data: data['scenes'][1].update(id='intro'),
    lambda data: data['scenes'][0].update(frames=59),
    lambda data: data['scenes'][2].update(assets=['unknown.figure']),
])
def test_storyboard_rejects_invalid_timeline(tmp_path: Path, mutation) -> None:
    data = json.loads(STORYBOARD.read_text(encoding='utf-8'))
    mutation(data)
    path = tmp_path / 'storyboard.json'
    path.write_text(json.dumps(data), encoding='utf-8')
    with pytest.raises(ValueError):
        load_storyboard(path)


def test_asset_bundle_rejects_escape_and_changed_bytes(tmp_path: Path) -> None:
    source = tmp_path / 'figure.png'
    source.write_bytes(b'approved figure')
    manifest = tmp_path / 'assets.json'
    asset = {'kind': 'image', 'path': 'figure.png', 'sha256': sha256(source)}
    manifest.write_text(json.dumps({'schema_version': 1, 'assets': {'human.circular': asset}}), encoding='utf-8')
    load_assets(manifest, complete=False)
    source.write_bytes(b'changed figure')
    with pytest.raises(ValueError, match='checksum mismatch'):
        load_assets(manifest, complete=False)
    asset['path'] = '../figure.png'
    manifest.write_text(json.dumps({'schema_version': 1, 'assets': {'human.circular': asset}}), encoding='utf-8')
    with pytest.raises(ValueError, match='Unsafe relative path'):
        load_assets(manifest, complete=False)


def test_contained_file_rejects_symlink_escape(tmp_path: Path) -> None:
    outside = tmp_path.parent / 'outside-video-asset.bin'
    outside.write_bytes(b'outside')
    root = tmp_path / 'bundle'
    root.mkdir()
    (root / 'link.bin').symlink_to(outside)
    with pytest.raises(ValueError, match='escaped file'):
        contained_file(root, 'link.bin')


def test_walkthrough_time_map_compresses_only_marked_waits() -> None:
    from video.walkthrough_render import TimeMap

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
    from video.walkthrough_render import Camera, TimeMap

    events = [
        {'t': 0.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0},
        {'t': 1.0, 'kind': 'camera', 'cx': 10, 'cy': 10, 'zoom': 2.0, 'duration': 1.0},
    ]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080))
    assert camera.at(0.5) == pytest.approx((960, 540, 1.0))
    assert camera.at(3.0) == pytest.approx((480, 270, 2.0))  # clamped to the top-left corner
    x, y, zoom = camera.at(1.5)
    assert 480 < x < 960 and 270 < y < 540 and 1.0 < zoom < 2.0
