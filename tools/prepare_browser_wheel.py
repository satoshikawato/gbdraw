#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"


def _load_build_support_module():
    build_support_path = REPO_ROOT / "gbdraw" / "_build_support.py"
    spec = spec_from_file_location("gbdraw_build_support", build_support_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load build support module from {build_support_path}")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_browser_wheel(*, refresh_cache_bust: bool = False) -> int:
    build_support = _load_build_support_module()
    expected_name = build_support.expected_browser_wheel_name()

    with tempfile.TemporaryDirectory(prefix="gbdraw-browser-wheel-") as tmpdir:
        outdir = Path(tmpdir) / "dist"
        stash_dir = Path(tmpdir) / "stashed-web-wheels"
        outdir.mkdir(parents=True, exist_ok=True)
        stash_dir.mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        env[build_support.BROWSER_WHEEL_BUILD_ENV] = "1"
        stashed_wheels: list[tuple[Path, Path]] = []
        try:
            for existing_wheel in sorted(WEB_ROOT.glob("gbdraw-*.whl")):
                stashed_path = stash_dir / existing_wheel.name
                shutil.move(existing_wheel, stashed_path)
                stashed_wheels.append((stashed_path, existing_wheel))

            shutil.rmtree(REPO_ROOT / "build", ignore_errors=True)

            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "build",
                    "--wheel",
                    "--no-isolation",
                    "--outdir",
                    str(outdir),
                ],
                cwd=REPO_ROOT,
                env=env,
                check=True,
            )

            built_wheel = outdir / expected_name
            if not built_wheel.exists():
                available = ", ".join(path.name for path in sorted(outdir.glob("gbdraw-*.whl")))
                raise FileNotFoundError(
                    f"Expected browser wheel {expected_name} was not produced. Available wheels: {available or 'none'}"
                )
            browser_wheel_member = f"gbdraw/web/{expected_name}"
            with zipfile.ZipFile(built_wheel) as zf:
                if browser_wheel_member in set(zf.namelist()):
                    raise RuntimeError(f"Browser wheel recursively contains itself: {browser_wheel_member}")

            target_path = WEB_ROOT / expected_name
            shutil.copy2(built_wheel, target_path)
        except Exception:
            for stashed_path, original_path in stashed_wheels:
                if stashed_path.exists():
                    shutil.move(stashed_path, original_path)
            raise

    cache_bust = build_support.generate_cache_bust_token() if refresh_cache_bust else None
    build_support.update_browser_wheel_config(wheel_name=expected_name, cache_bust=cache_bust)
    build_support.validate_browser_wheel_prepared()

    print(f"Prepared browser wheel: {target_path.relative_to(REPO_ROOT)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the browser wheel into gbdraw/web while preventing recursive self-inclusion. "
            "By default this keeps the existing cache-bust token unchanged."
        )
    )
    parser.add_argument(
        "--refresh-cache-bust",
        action="store_true",
        help="Refresh GBDRAW_WHEEL_CACHE_BUST in gbdraw/web/js/config.js after preparing the wheel.",
    )
    args = parser.parse_args(argv)
    return prepare_browser_wheel(refresh_cache_bust=args.refresh_cache_bust)


if __name__ == "__main__":
    raise SystemExit(main())
