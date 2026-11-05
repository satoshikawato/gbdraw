#!/usr/bin/env python3
"""Write the generated parts of the architecture description.

``LAYERS`` in ``tests/test_python_structure_ratchet.py`` and the first line of
each unit's docstring give the unit list in the "Project Structure" section of
``CLAUDE.md`` and the layer diagram in ``docs/internal/ARCHITECTURE.md``. Each
goes between its BEGIN GENERATED and END GENERATED markers.

``--check`` writes nothing. It fails when a generated part is out of date, and
when ``docs/internal/ARCHITECTURE.md`` names a ``gbdraw`` module, a name in
one, or a repository path that does not exist.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "gbdraw"
CLAUDE_MD = REPO_ROOT / "CLAUDE.md"
ARCHITECTURE_MD = REPO_ROOT / "docs" / "internal" / "ARCHITECTURE.md"
RATCHET = REPO_ROOT / "tests" / "test_python_structure_ratchet.py"

DOTTED_NAME = re.compile(r"(?<![\w./-])gbdraw(?:\.[A-Za-z_]\w*)+")
REPOSITORY_PATH = re.compile(r"(?<![\w./-])(?:gbdraw|tests|tools|docs)/[\w./-]*[\w/]")
UNITS_PER_LINE = 5


def read_layers() -> tuple[tuple[str, tuple[str, ...]], ...]:
    spec = importlib.util.spec_from_file_location("python_structure_ratchet", RATCHET)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.LAYERS


def _module_file(parts: list[str]) -> Path | None:
    base = REPO_ROOT.joinpath(*parts)
    for path in (base.with_suffix(".py"), base / "__init__.py"):
        if path.is_file():
            return path
    return None


def _unit_file(unit: str) -> Path | None:
    return _module_file(["gbdraw"] if unit == "__init__" else ["gbdraw", unit])


def unit_list(layers: tuple[tuple[str, tuple[str, ...]], ...]) -> str:
    """The units by layer, each with its docstring's first line; a unit not in the package yet is left out."""
    rows: list[tuple[str, str] | str] = []
    missing: list[str] = []
    for layer, units in layers:
        rows.append(layer)
        for unit in units:
            path = _unit_file(unit)
            if path is None:
                continue
            docstring = ast.get_docstring(ast.parse(path.read_text(encoding="utf-8")))
            if not docstring:
                missing.append(path.relative_to(REPO_ROOT).as_posix())
                continue
            name = f"{unit}/" if path.name == "__init__.py" and unit != "__init__" else f"{unit}.py"
            rows.append((name, docstring.strip().splitlines()[0]))
    if missing:
        raise SystemExit(f"units without a docstring: {', '.join(missing)}")
    width = max(len(row[0]) for row in rows if isinstance(row, tuple)) + 2
    lines = ["gbdraw/"]
    for row in rows:
        lines.append(f"  {row}" if isinstance(row, str) else f"    {row[0]:<{width}}{row[1]}")
    return "```text\n" + "\n".join(lines) + "\n```"


def layer_diagram(layers: tuple[tuple[str, tuple[str, ...]], ...]) -> str:
    """The layers from the top down; an arrow points from a layer to the one it may import from.

    Like ``unit_list``, it leaves out a unit that is not in the package yet.
    """
    lines = ["```mermaid", "flowchart TB"]
    for index, (layer, names) in reversed(list(enumerate(layers))):
        # Mermaid reads `__init__` as bold text, and its escapes differ between versions.
        units = ["package root" if unit == "__init__" else unit for unit in names if _unit_file(unit)]
        groups = [" · ".join(units[start : start + UNITS_PER_LINE]) for start in range(0, len(units), UNITS_PER_LINE)]
        lines.append(f'  L{index}["<b>{layer}</b><br/>{"<br/>".join(groups)}"]')
    lines.extend(f"  L{index} --> L{index - 1}" for index in range(len(layers) - 1, 0, -1))
    return "\n".join([*lines, "```"])


def replace_block(text: str, name: str, content: str, path: Path) -> str:
    begin = f"<!-- BEGIN GENERATED {name}: python tools/generate_architecture_docs.py -->"
    end = f"<!-- END GENERATED {name} -->"
    if text.count(begin) != 1 or text.count(end) != 1 or text.index(begin) > text.index(end):
        raise SystemExit(f"{path.relative_to(REPO_ROOT)}: expected one {begin!r} before one {end!r}")
    head, rest = text.split(begin)
    return f"{head}{begin}\n{content}\n{end}{rest.split(end, 1)[1]}"


def _defines(path: Path, names: list[str]) -> bool:
    """Whether the module defines ``names[0]`` at top level and, for a class, ``names[1]`` in its body."""
    body: list[ast.stmt] = ast.parse(path.read_text(encoding="utf-8")).body
    for name in names:
        found = None
        for node in body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                defined = [node.name]
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                defined = [alias.asname or alias.name for alias in node.names]
            elif isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                defined = [target.id for target in targets if isinstance(target, ast.Name)]
            else:
                continue
            if name in defined:
                found = node
        if found is None:
            return False
        body = found.body if isinstance(found, ast.ClassDef) else []
    return True


def unresolved_references(text: str) -> list[str]:
    """The ``gbdraw`` dotted names and repository paths in ``text`` that do not exist."""
    unresolved = []
    for match in sorted(set(DOTTED_NAME.findall(text))):
        parts = match.split(".")
        modules = ((end, _module_file(parts[:end])) for end in range(len(parts), 0, -1))
        end, module = next(((end, path) for end, path in modules if path is not None), (0, None))
        if module is None or not _defines(module, parts[end:]):
            unresolved.append(match)
    unresolved.extend(path for path in sorted(set(REPOSITORY_PATH.findall(text))) if not (REPO_ROOT / path).exists())
    return unresolved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="Fail instead of writing when a part is out of date.")
    args = parser.parse_args()
    layers = read_layers()
    blocks = ((CLAUDE_MD, "units", unit_list(layers)), (ARCHITECTURE_MD, "layers", layer_diagram(layers)))
    problems = []
    for path, name, content in blocks:
        text = path.read_text(encoding="utf-8")
        updated = replace_block(text, name, content, path)
        if updated == text:
            continue
        if args.check:
            problems.append(f"{path.relative_to(REPO_ROOT)}: the {name} part is out of date")
        else:
            path.write_text(updated, encoding="utf-8")
    problems.extend(
        f"{ARCHITECTURE_MD.relative_to(REPO_ROOT)} names {reference}, which does not exist "
        "(a name inside a class needs the module that defines the class)"
        for reference in unresolved_references(ARCHITECTURE_MD.read_text(encoding="utf-8"))
    )
    if problems:
        print("\n".join(problems), file=sys.stderr)
        print("Run python tools/generate_architecture_docs.py, and fix the names it reports.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
