"""Record, re-render, and check the gbdraw walkthrough and highlights videos."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

from video.walkthrough_render import check_walkthrough, render_walkthrough, sha256

REPO_ROOT = Path(__file__).resolve().parents[2]


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
    return {
        'git_head': _version('git', 'rev-parse', 'HEAD'),
        'git_diff_sha256': hashlib.sha256(_version('git', 'diff', '--binary').encode()).hexdigest(),
        'wheel_sha256': sha256(wheel),
        'recorder_sha256': sha256(REPO_ROOT / 'docs/capture/video/walkthrough.py'),
        'platform': platform.platform(), 'python': sys.version.split()[0],
        'playwright': version('playwright'),
        'ffmpeg': _version('ffmpeg', '-version').splitlines()[0],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    build = commands.add_parser('build', help='Record the Web app, render both videos, and check them')
    build.add_argument('--out', type=Path, required=True)
    render = commands.add_parser('render', help='Render both videos again from an existing recording')
    render.add_argument('--recording', type=Path, required=True)
    render.add_argument('--out', type=Path, required=True)
    check = commands.add_parser('check', help='Decode and verify rendered videos')
    check.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'build':
        _require_empty(args.out)
        # Only recording needs the browser; rendering reads the saved bundle.
        from video.walkthrough import record_walkthrough
        manifest_path = record_walkthrough(args.out)
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        manifest['provenance'] = _provenance()
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        render_walkthrough(args.out, args.out)
    elif args.command == 'render':
        _require_empty(args.out)
        render_walkthrough(args.recording, args.out)
    print(json.dumps(check_walkthrough(args.out), indent=2))


if __name__ == '__main__':
    main()
