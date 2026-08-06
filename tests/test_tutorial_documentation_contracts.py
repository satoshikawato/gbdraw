from __future__ import annotations

import json
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
TUTORIAL_ROOT = DOCS_ROOT / "TUTORIALS"
TUTORIAL_INDEX = TUTORIAL_ROOT / "README.md"
MANIFEST = DOCS_ROOT / "scenarios" / "manifest.json"
MARKDOWN_LINK_RE = re.compile(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)")
FIRST_H1_RE = re.compile(r"^#\s+(.+?)\s*$", re.MULTILINE)


def _markdown_links(path: Path) -> list[tuple[str, str]]:
    return MARKDOWN_LINK_RE.findall(path.read_text(encoding="utf-8"))


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


def _tutorial_chapters() -> list[dict[str, object]]:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    return [chapter for chapter in manifest["chapters"] if chapter["role"] == "tutorial"]


def _public_docs() -> list[Path]:
    paths = [REPO_ROOT / "README.md"]
    for path in sorted(DOCS_ROOT.rglob("*.md")):
        relative = path.relative_to(DOCS_ROOT)
        if relative.parts[0] in {"internal", "skills"}:
            continue
        paths.append(path)
    return paths


def test_surface_indexes_list_each_canonical_tutorial_once_with_its_h1() -> None:
    chapters = _tutorial_chapters()
    indexed: list[Path] = []

    for surface in ("GUI", "CLI", "PYTHON"):
        index = TUTORIAL_ROOT / surface / "README.md"
        surface_chapters = [
            chapter
            for chapter in chapters
            if str(chapter["surface"]).casefold() == surface.casefold()
        ]
        links = [
            (label, target)
            for label, raw_target in _markdown_links(index)
            if (target := _local_target(index, raw_target))
            in {REPO_ROOT / str(chapter["destination"]) for chapter in surface_chapters}
        ]
        assert [target for _, target in links] == [
            REPO_ROOT / str(chapter["destination"]) for chapter in surface_chapters
        ]
        for (label, target), chapter in zip(links, surface_chapters, strict=True):
            heading = FIRST_H1_RE.search(target.read_text(encoding="utf-8"))
            assert heading is not None
            assert label == heading.group(1) == chapter["title"]
        indexed.extend(target for _, target in links)

    assert indexed == [REPO_ROOT / str(chapter["destination"]) for chapter in chapters]


def test_tutorial_root_routes_by_surface() -> None:
    surface_targets = [
        TUTORIAL_ROOT / "GUI" / "README.md",
        TUTORIAL_ROOT / "CLI" / "README.md",
        TUTORIAL_ROOT / "PYTHON" / "README.md",
    ]
    links = _markdown_links(TUTORIAL_INDEX)
    resolved = [_local_target(TUTORIAL_INDEX, target) for _, target in links]

    assert [target.resolve() for target in surface_targets] == [
        target for target in resolved if target in {path.resolve() for path in surface_targets}
    ]


def test_readme_and_docs_landing_route_through_the_tutorial_index() -> None:
    for source in (REPO_ROOT / "README.md", DOCS_ROOT / "DOCS.md"):
        targets = {
            target
            for _, raw_target in _markdown_links(source)
            if (target := _local_target(source, raw_target)) is not None
        }
        assert TUTORIAL_INDEX.resolve() in targets


def test_public_tutorial_labels_route_to_the_canonical_index() -> None:
    incorrect: list[str] = []
    for source in _public_docs():
        for label, raw_target in _markdown_links(source):
            if label.strip().strip("*_`") != "Tutorials":
                continue
            if _local_target(source, raw_target) != TUTORIAL_INDEX.resolve():
                incorrect.append(f"{source.relative_to(REPO_ROOT)} -> {raw_target}")
    assert incorrect == []


def test_tutorial_root_contains_only_the_index() -> None:
    assert {path.resolve() for path in TUTORIAL_ROOT.glob("*.md")} == {TUTORIAL_INDEX.resolve()}


def test_public_docs_do_not_expose_the_example_maintenance_script() -> None:
    public_docs = [*_public_docs(), REPO_ROOT / "examples" / "color_palette_examples.md"]
    offenders = [
        path.relative_to(REPO_ROOT).as_posix()
        for path in public_docs
        if "tools/reproduce_examples.py" in path.read_text(encoding="utf-8")
    ]
    assert offenders == []


def test_public_docs_do_not_expose_fixture_provenance_maintenance_details() -> None:
    forbidden = ("tutorial-fixture-provenance.md", "R-FIXTURE-01")
    offenders = [
        f"{path.relative_to(REPO_ROOT).as_posix()}: {token}"
        for path in _public_docs()
        for token in forbidden
        if token in path.read_text(encoding="utf-8")
    ]

    assert offenders == []
