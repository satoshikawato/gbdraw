"""Record, re-render, and check the gbdraw walkthrough and Gallery tutorial videos."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path

from video.walkthrough_render import STORIES, check_walkthrough, render_walkthrough, sha256

REPO_ROOT = Path(__file__).resolve().parents[2]


def _is_empty(path: Path) -> bool:
    return not path.exists() or not any(path.iterdir())


def _require_empty(path: Path) -> None:
    if not _is_empty(path):
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
        'tutorials_sha256': sha256(REPO_ROOT / 'docs/capture/video/tutorials.py'),
        'platform': platform.platform(), 'python': sys.version.split()[0],
        'playwright': version('playwright'),
        'ffmpeg': _version('ffmpeg', '-version').splitlines()[0],
    }


def collect_final(final: Path, target: Path) -> list[Path]:
    """Copy one run's videos, subtitles, and covers to ``target``, never over a file."""

    target.mkdir(parents=True, exist_ok=True)
    copied = []
    for source in sorted(final.iterdir()):
        destination = target / source.name
        if destination.exists():
            raise FileExistsError(f'{destination} already exists')
        shutil.copyfile(source, destination)
        if sha256(destination) != sha256(source):
            raise OSError(f'{destination} does not match {source}')
        copied.append(destination)
    return copied


def build_all(names: list[str], out: Path, *, jobs: int = 1, collect: Path | None = None,
              prepare_wheel: bool = False) -> dict[str, Path]:
    """Build each journey into ``out/<journey>``, then collect the finished files.

    Every destination is checked before anything runs, so a long build never
    fails at the end on an existing file.
    """

    out = out.resolve()
    runs = {name: out / name for name in names}
    targets = {name: collect.resolve() / name for name in names} if collect else {}
    busy = [path for path in (out, *targets.values()) if not _is_empty(path)]
    if busy:
        raise FileExistsError(f'Choose new or empty directories: {", ".join(map(str, busy))}')
    if jobs < 1:
        raise ValueError('--jobs must be at least 1')
    if prepare_wheel:
        # The wheel build deletes build/, so it runs before any output exists.
        subprocess.run([sys.executable, 'tools/prepare_browser_wheel.py'], cwd=REPO_ROOT, check=True)
    logs = out / 'logs'
    logs.mkdir(parents=True)

    def build(name: str) -> tuple[str, int]:
        # One process per journey: each owns its browser and its encoder.
        with (logs / f'{name}.log').open('w', encoding='utf-8') as log:
            result = subprocess.run(
                [sys.executable, str(Path(__file__).resolve()), 'build', '--journey', name, '--out', str(runs[name])],
                cwd=REPO_ROOT, stdout=log, stderr=subprocess.STDOUT)
        print(f'{name}: {"PASS" if result.returncode == 0 else "FAIL"} ({runs[name]})', flush=True)
        return name, result.returncode

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        failed = sorted(name for name, code in pool.map(build, names) if code)
    if failed:
        raise SystemExit(f'Failed journeys: {", ".join(failed)}; see {logs}')
    for name, target in targets.items():
        collect_final(runs[name] / 'final', target)
        print(f'{name}: collected in {target}', flush=True)
    return runs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    build = commands.add_parser('build', help="Record one journey in the Web app, render its videos, and check them")
    build.add_argument('--out', type=Path, required=True)
    build.add_argument('--journey', choices=sorted(STORIES), default='meet-gbdraw',
                       help='meet-gbdraw (default) or a Gallery example id')
    render = commands.add_parser('render', help="Render a recording's videos again")
    render.add_argument('--recording', type=Path, required=True)
    render.add_argument('--out', type=Path, required=True)
    check = commands.add_parser('check', help='Decode and verify rendered videos')
    check.add_argument('--out', type=Path, required=True)
    every = commands.add_parser('build-all', help='Build several journeys, each in OUT/<journey>, and collect their videos')
    every.add_argument('--out', type=Path, required=True)
    every.add_argument('--journey', action='append', choices=sorted(STORIES),
                       help='journey to build; repeat for more (default: every journey)')
    every.add_argument('--jobs', type=int, default=1, help='journeys to build at once (default: 1)')
    every.add_argument('--collect', type=Path, help='copy each final/ folder to COLLECT/<journey>/')
    every.add_argument('--prepare-wheel', action='store_true',
                       help='run tools/prepare_browser_wheel.py first; it deletes build/')
    args = parser.parse_args()
    if args.command == 'build-all':
        build_all(list(dict.fromkeys(args.journey or STORIES)), args.out, jobs=args.jobs, collect=args.collect,
                  prepare_wheel=args.prepare_wheel)
        return
    if args.command == 'build':
        _require_empty(args.out)
        # Only recording needs the browser; rendering reads the saved bundle.
        from video.walkthrough import record_walkthrough
        manifest_path = record_walkthrough(args.out, args.journey)
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
