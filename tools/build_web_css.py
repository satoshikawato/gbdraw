#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_CSS = REPO_ROOT / "gbdraw" / "web" / "css" / "tailwind.source.css"
OUTPUT_CSS = REPO_ROOT / "gbdraw" / "web" / "vendor" / "tailwind" / "app.min.css"
CONFIG_JS = REPO_ROOT / "gbdraw" / "web" / "tailwind.config.js"
LINUX_TAILWIND = REPO_ROOT / "tools" / "tailwind" / "bin" / "tailwindcss-linux-x64"
WINDOWS_TAILWIND = REPO_ROOT / "tools" / "tailwind" / "bin" / "tailwindcss-windows-x64.exe"


def get_tailwind_binary() -> Path:
    if os.name == "nt":
        return WINDOWS_TAILWIND
    return LINUX_TAILWIND


def ensure_inputs() -> None:
    missing = [
        p
        for p in (INPUT_CSS, CONFIG_JS)
        if not p.exists()
    ]
    if missing:
        missing_text = "\n".join(f"- {p}" for p in missing)
        raise RuntimeError(f"Missing required web CSS source file(s):\n{missing_text}")


def ensure_tailwind_binary(binary: Path) -> None:
    if not binary.exists():
        raise RuntimeError(
            "Tailwind standalone CLI is missing.\n"
            f"Expected binary: {binary}\n"
            "Place the standalone executable at the path above."
        )
    if os.name != "nt" and not os.access(binary, os.X_OK):
        raise RuntimeError(
            "Tailwind standalone CLI is not executable.\n"
            f"Expected executable file: {binary}\n"
            f"Run: chmod +x {binary}"
        )


def run_tailwind(binary: Path, output_css: Path) -> None:
    cmd = [
        str(binary),
        "-i",
        str(INPUT_CSS),
        "-o",
        str(output_css),
        "--config",
        str(CONFIG_JS),
        "--minify",
    ]
    try:
        subprocess.run(cmd, check=True)
    except OSError as exc:
        raise RuntimeError(
            "Failed to execute Tailwind standalone CLI.\n"
            f"Binary: {binary}\n"
            f"Details: {exc}"
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "Tailwind CSS build failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Exit code: {exc.returncode}"
        ) from exc


def check_mode(binary: Path) -> int:
    with tempfile.TemporaryDirectory(prefix="gbdraw-tailwind-check-") as tmpdir:
        tmp_output = Path(tmpdir) / "app.min.css"
        run_tailwind(binary, tmp_output)
        if not OUTPUT_CSS.exists():
            print(
                f"[FAIL] Missing generated CSS: {OUTPUT_CSS}\n"
                "Run `python tools/build_web_css.py` to generate it.",
                file=sys.stderr,
            )
            return 1
        if OUTPUT_CSS.read_bytes() != tmp_output.read_bytes():
            print(
                f"[FAIL] {OUTPUT_CSS} is out of date.\n"
                "Run `python tools/build_web_css.py` and commit the updated CSS.",
                file=sys.stderr,
            )
            return 1
    print(f"[OK] {OUTPUT_CSS} is up to date")
    return 0


def build_mode(binary: Path) -> int:
    OUTPUT_CSS.parent.mkdir(parents=True, exist_ok=True)
    run_tailwind(binary, OUTPUT_CSS)
    print(f"[OK] Wrote {OUTPUT_CSS}")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build static Tailwind CSS for gbdraw web UI.")
    parser.add_argument("--check", action="store_true", help="Verify generated CSS is up to date.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_inputs()
    binary = get_tailwind_binary()
    ensure_tailwind_binary(binary)
    return check_mode(binary) if args.check else build_mode(binary)


if __name__ == "__main__":
    raise SystemExit(main())
