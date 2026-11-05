#!/usr/bin/env python3
"""Build the Worker's Rust drawing engine into gbdraw/web/wasm/gbdraw-core/ (generated, gitignored)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    subprocess.run(
        ["cargo", "build", "--profile", "wasm", "--locked", "-p", "gbdraw-core-wasm", "--target", "wasm32-unknown-unknown"],
        cwd=ROOT / "rust",
        check=True,
    )
    target = ROOT / "gbdraw" / "web" / "wasm" / "gbdraw-core"
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "rust" / "target" / "wasm32-unknown-unknown" / "wasm" / "gbdraw_core.wasm", target / "gbdraw_core.wasm")


if __name__ == "__main__":
    main()
