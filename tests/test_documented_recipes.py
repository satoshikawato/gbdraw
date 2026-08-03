from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
MANIFEST = DOCS_ROOT / "scenarios" / "manifest.json"
LEGACY_ROUTERS = (
    DOCS_ROOT / "QUICKSTART.md",
    *(sorted((DOCS_ROOT / "TUTORIALS").glob("[0-9]*_*.md"))),
    DOCS_ROOT / "GFF3_FASTA.md",
    DOCS_ROOT / "EXPORT.md",
    DOCS_ROOT / "PYTHON_API.md",
    DOCS_ROOT / "TYPED_API.md",
    DOCS_ROOT / "WORKFLOW_GUIDE.md",
)
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
    for chapter in manifest["chapters"]:
        destination = (REPO_ROOT / chapter["destination"]).resolve()
        for source in chapter["sources"]:
            destinations[(REPO_ROOT / source).resolve()].add(destination)
    return destinations


def test_migrated_legacy_pages_are_short_non_runnable_routers() -> None:
    for path in LEGACY_ROUTERS:
        source = path.read_text(encoding="utf-8")
        assert len(source.splitlines()) <= 80, path.relative_to(REPO_ROOT)
        assert "```bash" not in source
        assert "```python" not in source
        assert "<!-- executable:" not in source
        assert "tests/test_inputs" not in source
        assert "wget " not in source


def test_legacy_routers_link_every_manifest_owned_replacement() -> None:
    destinations_by_source = _canonical_destinations_by_source()
    for path in LEGACY_ROUTERS:
        expected = destinations_by_source.get(path.resolve(), set())
        if path.name == "QUICKSTART.md":
            expected = {
                DOCS_ROOT / "TUTORIALS" / "GUI" / "first-circular-genome-diagram.md",
                DOCS_ROOT / "TUTORIALS" / "CLI" / "first-circular-genome-diagram.md",
                DOCS_ROOT / "TUTORIALS" / "PYTHON" / "first-genome-diagram.md",
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
        (REPO_ROOT / chapter["destination"]).resolve()
        for chapter in manifest["chapters"]
    }
    assert destinations.isdisjoint(path.resolve() for path in LEGACY_ROUTERS)
