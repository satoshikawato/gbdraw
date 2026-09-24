"""Build, re-edit, and verify the Meet gbdraw introduction video."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

from video.model import sha256
from video.render import FONT, render_video
from video.validate import check_video

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STORYBOARD = REPO_ROOT / 'docs/videos/meet-gbdraw/storyboard.json'


def _require_empty(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f'Choose a new or empty run directory: {path}')
    path.mkdir(parents=True, exist_ok=True)


def _version(*args: str) -> str:
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout.strip()


def _provenance() -> dict:
    wheel = next((REPO_ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl'), None)
    if wheel is None:
        raise FileNotFoundError('Prepare the generated browser wheel first')
    source_paths = (
        'docs/capture/video/capture.py', 'docs/capture/video/human_edits.py',
        'docs/capture/flows/bgc_losatp.py', 'docs/capture/flows/tutorials/gui_first_circular.py',
        'docs/capture/flows/tutorials/gui_annotated_chloroplast.py',
        'docs/capture/flows/tutorials/gui_losatn.py', 'docs/capture/flows/tutorials/gui_losatp_groups.py',
        'docs/capture/flows/web_capture.py',
    )
    return {
        'git_head': _version('git', 'rev-parse', 'HEAD'),
        'git_branch': _version('git', 'branch', '--show-current'),
        'git_diff_sha256': hashlib.sha256(_version('git', 'diff', '--binary').encode()).hexdigest(),
        'capture_source_sha256': {name: sha256(REPO_ROOT / name) for name in source_paths},
        'wheel_sha256': sha256(wheel),
        'input_manifest_sha256': sha256(REPO_ROOT / 'gbdraw/web/tutorial-data/manifest.json'),
        'platform': platform.platform(), 'python': sys.version.split()[0],
        'playwright': version('playwright'),
        'ffmpeg': _version('ffmpeg', '-version').splitlines()[0],
        'font': {'name': 'Lato Bold', 'sha256': sha256(FONT)},
        'viewport': [1920, 1080],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    build = commands.add_parser('build')
    build.add_argument('--storyboard', type=Path, default=DEFAULT_STORYBOARD)
    build.add_argument('--out', type=Path, required=True)
    render = commands.add_parser('render')
    render.add_argument('--assets', type=Path, required=True)
    render.add_argument('--storyboard', type=Path, default=DEFAULT_STORYBOARD)
    render.add_argument('--out', type=Path, required=True)
    check = commands.add_parser('check')
    check.add_argument('--out', type=Path, required=True)
    check.add_argument('--visual', action='store_true')
    check.add_argument('--baseline', type=Path, default=REPO_ROOT / 'docs/videos/meet-gbdraw/reference-frames')
    args = parser.parse_args()
    if args.command == 'build':
        _require_empty(args.out)
        # Browser dependencies stay behind the build command. Render and check
        # can run from a captured bundle without starting an app or LOSAT.
        from video.capture import capture_intro_assets
        from video.human_edits import capture_human_edits
        assets = {}
        capture_intro_assets(args.out, assets)
        capture_human_edits(args.out, assets)
        manifest = args.out / 'assets.json'
        data = json.loads(manifest.read_text(encoding='utf-8'))
        data['provenance'] = _provenance()
        manifest.write_text(json.dumps(data, indent=2), encoding='utf-8')
        video = render_video(manifest, args.storyboard, args.out)
        print(video)
        print(json.dumps(check_video(args.out), indent=2))
    elif args.command == 'render':
        _require_empty(args.out)
        video = render_video(args.assets, args.storyboard, args.out)
        print(video)
        print(json.dumps(check_video(args.out), indent=2))
    else:
        print(json.dumps(check_video(args.out, visual=args.visual, baseline=args.baseline), indent=2))


if __name__ == '__main__':
    main()
