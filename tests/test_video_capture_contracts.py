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
