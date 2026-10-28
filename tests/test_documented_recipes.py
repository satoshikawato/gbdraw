from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pytest

from docs.recipes._scenario_support import destination_page


pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
MANIFEST = DOCS_ROOT / "scenarios" / "manifest.json"
LEGACY_ROUTERS = (DOCS_ROOT / "QUICKSTART.md",)
MARKDOWN_LINK_RE = re.compile(r"(?<!!)\[[^\]\n]+\]\(([^)\n]+)\)")


def _local_target(source: Path, raw_target: str) -> Path | None:
    target = raw_target.strip()
    if target.startswith("<"):
        target = target[1:].split(">", 1)[0]
    else:
        target = target.split(maxsplit=1)[0]
    if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.IGNORECASE) or target.startswith("//"):
        return None
    path_part = target.split("#", 1)[0].split("?", 1)[0]
    return source.resolve() if not path_part else (source.parent / path_part).resolve()


def _canonical_destinations_by_source() -> dict[Path, set[Path]]:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    destinations: dict[Path, set[Path]] = defaultdict(set)
    for chapter in manifest["scenarios"]:
        if "destination" not in chapter:
            continue
        destination = (REPO_ROOT / str(destination_page(chapter))).resolve()
        for source in chapter["sources"]:
            destinations[(REPO_ROOT / source).resolve()].add(destination)
    return destinations


def test_migrated_legacy_pages_are_short_routers() -> None:
    for path in LEGACY_ROUTERS:
        source = path.read_text(encoding="utf-8")
        assert len(source.splitlines()) <= 80, path.relative_to(REPO_ROOT)
        assert "```python" not in source
        assert "<!-- executable:" not in source
        assert "tests/test_inputs" not in source
        assert "wget " not in source


def _first_figures_block(path: Path) -> str:
    blocks = re.findall(r"```bash\n(.*?)```", path.read_text(encoding="utf-8"), re.S)
    blocks = [block for block in blocks if block.startswith("data=gbdraw/gbdraw/web/tutorial-data\n")]
    assert len(blocks) == 1, path.relative_to(REPO_ROOT)
    return blocks[0]


def test_readme_and_quickstart_first_figures_run_from_a_clean_directory(tmp_path: Path) -> None:
    """The README and the Quickstart give one block; it runs as written beside a checkout."""

    block = _first_figures_block(REPO_ROOT / "README.md")
    # The Quickstart's only code block is this one, so no unexecuted command slips in.
    quickstart = (DOCS_ROOT / "QUICKSTART.md").read_text(encoding="utf-8")
    assert re.findall(r"```([^\n]*)\n(.*?)```", quickstart, re.S) == [("bash", block)]
    (tmp_path / "gbdraw").symlink_to(REPO_ROOT, target_is_directory=True)
    # The documented `gbdraw` command, run from this checkout, not an installed gbdraw.
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    python = shlex.quote(sys.executable)
    shim = f"gbdraw() {{ {python} -c 'from gbdraw.cli import main; main()' \"$@\"; }}\n"
    imported = subprocess.run(
        ["bash", "-e", "-c", f"{python} -c 'import gbdraw; print(gbdraw.__file__)'"],
        cwd=tmp_path, env=env, capture_output=True, text=True, check=True,
    ).stdout.strip()
    assert Path(imported) == REPO_ROOT / "gbdraw" / "__init__.py"
    result = subprocess.run(
        ["bash", "-e", "-c", shim + block], cwd=tmp_path, env=env, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr[-2000:]
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "gbdraw", "human_mitochondrion.svg", "lambda_de3.svg"
    ]
    linear = (tmp_path / "lambda_de3.svg").read_text(encoding="utf-8")
    assert linear.count('data-gbdraw-match-id="comparison1_match') == 6
    assert "NC_012920.1" in (tmp_path / "human_mitochondrion.svg").read_text(encoding="utf-8")


def test_legacy_routers_link_every_manifest_owned_replacement() -> None:
    destinations_by_source = _canonical_destinations_by_source()
    for path in LEGACY_ROUTERS:
        expected = destinations_by_source.get(path.resolve(), set())
        if path.name == "QUICKSTART.md":
            expected = expected | {
                DOCS_ROOT / "TUTORIALS" / "first-circular-genome-diagram.md",
                DOCS_ROOT / "TUTORIALS" / "first-linear-genome-diagram.md",
            }
        assert expected, path.relative_to(REPO_ROOT)
        linked = {
            target
            for raw_target in MARKDOWN_LINK_RE.findall(path.read_text(encoding="utf-8"))
            if (target := _local_target(path, raw_target)) is not None
        }
        assert {target.resolve() for target in expected} <= linked, path.relative_to(REPO_ROOT)


def test_legacy_paths_are_never_canonical_chapter_destinations() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    destinations = {
        (REPO_ROOT / str(destination_page(chapter))).resolve()
        for chapter in manifest["scenarios"]
        if "destination" in chapter
    }
    assert destinations.isdisjoint(path.resolve() for path in LEGACY_ROUTERS)
