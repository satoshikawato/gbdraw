from __future__ import annotations

import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PALETTE_ENTRYPOINT = Path("gbdraw/web/gallery/palettes/index.html")


def test_palette_explorer_entrypoint_is_tracked() -> None:
    entrypoint = REPO_ROOT / PALETTE_ENTRYPOINT
    assert entrypoint.is_file()

    result = subprocess.run(
        ["git", "ls-files", "--error-unmatch", PALETTE_ENTRYPOINT.as_posix()],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        "The palette explorer entrypoint is not tracked by Git, so clean deployments omit it. "
        "Check the repository's *.html ignore rules."
    )
