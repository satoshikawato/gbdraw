#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from runpy import run_path
from types import SimpleNamespace


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"


def _load_build_support_module():
    return SimpleNamespace(**run_path(str(REPO_ROOT / "gbdraw" / "_build_support.py")))


def add_pyodide_bytecode(wheel: Path) -> None:
    """Add the vendored Pyodide's bytecode of every module to the wheel (tools/browser_wheel_bytecode.py)."""

    node = shutil.which("node")
    if node is None:
        raise RuntimeError(
            "The browser wheel needs Node.js (`node`): tools/compile_browser_wheel_bytecode.mjs compiles its "
            "bytecode with the vendored Pyodide."
        )
    subprocess.run([node, str(REPO_ROOT / "tools" / "compile_browser_wheel_bytecode.mjs"), str(wheel)], check=True)


def prepare_browser_wheel(*, build_isolation: bool = True) -> int:
    build_support = _load_build_support_module()
    expected_name = build_support.expected_browser_wheel_name()
    build_support.refresh_open_source_notices()

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

            wheel_command = [
                sys.executable,
                "-m",
                "pip",
                "wheel",
                "--no-deps",
                "--no-cache-dir",
            ]
            if not build_isolation:
                wheel_command.append("--no-build-isolation")
            wheel_command.extend(["--wheel-dir", str(outdir), str(REPO_ROOT)])
            subprocess.run(
                wheel_command,
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
            add_pyodide_bytecode(built_wheel)

            target_path = WEB_ROOT / expected_name
            shutil.copy2(built_wheel, target_path)
        except Exception:
            for stashed_path, original_path in stashed_wheels:
                if stashed_path.exists():
                    shutil.move(stashed_path, original_path)
            raise

    build_support.update_browser_wheel_config(wheel_name=expected_name)
    build_support.validate_browser_wheel_prepared()

    print(f"Prepared browser wheel: {target_path.relative_to(REPO_ROOT)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the browser wheel into gbdraw/web while preventing recursive self-inclusion. "
            "The hosted bundle (tools/prepare_cloudflare_pages.py) names it by its content hash."
        )
    )
    parser.add_argument(
        "--no-build-isolation",
        action="store_true",
        help="Use build requirements already installed in the current environment.",
    )
    args = parser.parse_args(argv)
    return prepare_browser_wheel(build_isolation=not args.no_build_isolation)


if __name__ == "__main__":
    raise SystemExit(main())
