from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
TUTORIAL_ROOT = DOCS_ROOT / "TUTORIALS"
CANONICAL_TUTORIAL_INDEX = TUTORIAL_ROOT / "README.md"
LEGACY_TUTORIAL_INDEX = TUTORIAL_ROOT / "TUTORIALS.md"

MARKDOWN_LINK_RE = re.compile(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)")
NUMBERED_TUTORIAL_RE = re.compile(r"(\d+)_.+\.md")
FIRST_H1_RE = re.compile(r"^#\s+(.+?)\s*$", re.MULTILINE)
DEVELOPER_DOC_MARKERS = (
    "IMPLEMENTATION_PLAN",
    "_AUDIT",
    "_OPERATION_",
    "_CAPTURE_PLAN",
)


def _markdown_links(path: Path) -> list[tuple[str, str]]:
    return MARKDOWN_LINK_RE.findall(path.read_text(encoding="utf-8"))


def _section(path: Path, heading: str) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf"^## {re.escape(heading)}\s*$\n(.*?)(?=^## |\Z)",
        text,
        re.MULTILINE | re.DOTALL,
    )
    assert match is not None, f"{path.name} has no {heading!r} section"
    return match.group(1)


def _local_target(source: Path, raw_target: str) -> Path | None:
    target = raw_target.strip()
    if target.startswith("<"):
        target = target[1:].split(">", 1)[0]
    else:
        target = target.split(maxsplit=1)[0]

    if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.IGNORECASE) or target.startswith(
        "//"
    ):
        return None

    path_part = target.split("#", 1)[0].split("?", 1)[0]
    if not path_part:
        return source.resolve()
    return (source.parent / path_part).resolve()


def _numbered_tutorials() -> list[Path]:
    tutorials: list[tuple[int, Path]] = []
    for path in TUTORIAL_ROOT.glob("[0-9]*_*.md"):
        match = NUMBERED_TUTORIAL_RE.fullmatch(path.name)
        if match is not None:
            tutorials.append((int(match.group(1)), path.resolve()))
    return [path for _, path in sorted(tutorials)]


def _public_docs() -> list[Path]:
    paths = [REPO_ROOT / "README.md"]
    for path in sorted(DOCS_ROOT.rglob("*.md")):
        relative = path.relative_to(DOCS_ROOT)
        if relative.parts[0] == "skills":
            continue
        if any(marker in path.name for marker in DEVELOPER_DOC_MARKERS):
            continue
        paths.append(path)
    return paths


def test_numbered_tutorials_are_contiguous_from_one_through_nine() -> None:
    candidates = sorted(TUTORIAL_ROOT.glob("[0-9]*_*.md"))
    malformed = [
        path.name
        for path in candidates
        if NUMBERED_TUTORIAL_RE.fullmatch(path.name) is None
    ]
    assert malformed == []

    numbers = sorted(int(path.name.split("_", 1)[0]) for path in candidates)
    assert numbers == list(range(1, 10))


def test_canonical_index_lists_each_numbered_tutorial_once_in_order() -> None:
    tutorials = _numbered_tutorials()
    tutorial_set = set(tutorials)
    indexed_links = [
        (label, target)
        for label, raw_target in _markdown_links(CANONICAL_TUTORIAL_INDEX)
        if (target := _local_target(CANONICAL_TUTORIAL_INDEX, raw_target))
        in tutorial_set
    ]

    assert [target for _, target in indexed_links] == tutorials

    for label, target in indexed_links:
        heading = FIRST_H1_RE.search(target.read_text(encoding="utf-8"))
        assert heading is not None, f"{target.name} has no H1"
        assert label == heading.group(1), f"{target.name} index label differs from its H1"


def test_readme_and_docs_landing_page_use_only_the_canonical_index() -> None:
    entry_points = (
        (REPO_ROOT / "README.md", "Documentation"),
        (DOCS_ROOT / "DOCS.md", "Tutorials"),
    )
    numbered_tutorials = set(_numbered_tutorials())

    for source, navigation_heading in entry_points:
        text = source.read_text(encoding="utf-8")
        targets = {
            target
            for _, raw_target in _markdown_links(source)
            if (target := _local_target(source, raw_target)) is not None
        }
        assert CANONICAL_TUTORIAL_INDEX.resolve() in targets
        assert LEGACY_TUTORIAL_INDEX.resolve() not in targets
        assert "TUTORIALS/TUTORIALS.md" not in text

        navigation_targets = {
            target
            for _, raw_target in MARKDOWN_LINK_RE.findall(
                _section(source, navigation_heading)
            )
            if (target := _local_target(source, raw_target)) is not None
        }
        assert CANONICAL_TUTORIAL_INDEX.resolve() in navigation_targets
        assert navigation_targets.isdisjoint(numbered_tutorials)


def test_tutorial_directory_has_one_index_and_no_compatibility_pages() -> None:
    expected = {CANONICAL_TUTORIAL_INDEX.resolve(), *_numbered_tutorials()}
    actual = {path.resolve() for path in TUTORIAL_ROOT.glob("*.md")}

    assert actual == expected


def test_public_tutorial_index_links_resolve_to_the_canonical_index() -> None:
    incorrect: list[str] = []
    canonical = CANONICAL_TUTORIAL_INDEX.resolve()
    legacy = LEGACY_TUTORIAL_INDEX.resolve()

    for source in _public_docs():
        for label, raw_target in _markdown_links(source):
            visible_label = label.strip().strip("*_`")
            target = _local_target(source, raw_target)
            if target == legacy or (visible_label == "Tutorials" and target != canonical):
                relative = source.relative_to(REPO_ROOT).as_posix()
                incorrect.append(f"{relative} -> {raw_target}")

    assert incorrect == []


def test_public_docs_do_not_expose_the_example_maintenance_script() -> None:
    public_docs = [
        *_public_docs(),
        REPO_ROOT / "examples" / "color_palette_examples.md",
    ]
    offenders = [
        path.relative_to(REPO_ROOT).as_posix()
        for path in public_docs
        if "tools/reproduce_examples.py" in path.read_text(encoding="utf-8")
    ]
    assert offenders == []
