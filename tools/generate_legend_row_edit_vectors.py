#!/usr/bin/env python3
"""Generate the shared Legend row edit vectors (Python's drawn rows, exact).

``tests/fixtures/legend_row_edit_vectors.json`` holds, for each case, the
Legend row edits as config leaves (``objects.legends.row_names``,
``hidden_rows``, ``row_order``), the Web editor operations that make the same
edits, the Legend Python draws without the edits (``unedited``) and with them
(``edited``), and the feature Legend groups of ``edited`` without geometry
(``expected``).

- ``tests/test_legend_row_edit_vectors.py`` redraws every case with Python
  (this generator's ``--check``): Python's drawn rows equal ``expected``.
- ``tests/web/legend-row-edit-vectors.test.mjs`` runs the Web's executor
  (``applyLegendOperations`` and ``nameLegendRows``, through Generate
  admission) over ``unedited`` with ``operations``: its rows equal
  ``expected``, and a second pass over ``edited`` changes nothing.

``expected`` leaves out ``transform``: Python lays the edited Legend out again,
while the executor keeps the rows' places until the Web lays it out. The Legend
layout vectors (``legend_layout_vectors.json``) pin the places.

Usage:
  python tools/generate_legend_row_edit_vectors.py            # redraw the vectors
  python tools/generate_legend_row_edit_vectors.py --check    # fail when stale
"""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import TypedDict

REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "tests" / "fixtures" / "legend_row_edit_vectors.json"
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import FeatureLocation, SeqFeature  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

import gbdraw  # noqa: E402
from gbdraw.config.models.objects import (  # noqa: E402
    LEGEND_HIDDEN_ROWS_PATH,
    LEGEND_ROW_NAMES_PATH,
    LEGEND_ROW_ORDER_PATH,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS  # noqa: E402

# Geometry the executor leaves for the Web's layout (see the module docstring).
GEOMETRY_ATTRIBUTES = frozenset({"transform"})
FEATURE_TYPES = ("CDS", "tRNA", "rRNA", "repeat_region", "misc_feature")
# misc_feature and repeat_region share one default swatch, so a rename of one onto the other merges the rows.
LAYOUT = ("CDS", "CDS", "tRNA", "rRNA", "repeat_region", "misc_feature", "CDS", "tRNA")


class LegendNode(TypedDict):
    """An element as the vectors compare it: tag, attributes without geometry, text of a leaf, children."""

    tag: str
    attributes: dict[str, str]
    text: str | None
    children: list[LegendNode]


class LegendGroupFacts(TypedDict):
    """One feature Legend group (the group the Web's executor edits) and its rows, in document order."""

    attributes: dict[str, str]
    rows: list[LegendNode]


@dataclass(frozen=True)
class Case:
    id: str
    mode: str
    row_names: tuple[tuple[str, str], ...] = ()
    hidden_rows: tuple[str, ...] = ()
    row_order: tuple[str, ...] = ()
    comparison: bool = False
    note: str = ""
    legend: str = "right"
    extra: dict[str, object] = field(default_factory=dict)


CIRCULAR_ROWS = ("CDS", "tRNA", "rRNA", "repeat_region", "misc_feature", "GC content", "GC skew (+)", "GC skew (-)")
LINEAR_ROWS = ("CDS", "tRNA", "rRNA", "repeat_region", "misc_feature")

CASES: tuple[Case, ...] = (
    Case("circular-rename-plain", "circular", row_names=(("CDS", "Protein-coding genes"),)),
    Case("circular-rename-merges", "circular", row_names=(("misc_feature", "repeat_region"),),
         note="same name and swatch: the later row is drawn hidden, merged into the earlier"),
    Case("circular-rename-numbers", "circular", row_names=(("tRNA", "RNA genes"), ("rRNA", "RNA genes")),
         note="same name, other swatches: the names are numbered"),
    Case("circular-delete-one", "circular", hidden_rows=("tRNA",)),
    Case("circular-delete-all-but-one", "circular", hidden_rows=CIRCULAR_ROWS[1:]),
    Case("circular-delete-merged", "circular", row_names=(("misc_feature", "repeat_region"),),
         hidden_rows=("repeat_region",), note="the row a rename merged into is deleted: the renamed row shows"),
    Case("circular-order-move", "circular",
         row_order=("rRNA", "CDS", "tRNA", "repeat_region", "misc_feature", "GC content", "GC skew (+)", "GC skew (-)")),
    Case("circular-order-reverse", "circular", row_order=tuple(reversed(CIRCULAR_ROWS))),
    Case("circular-order-sort-a-z", "circular", row_order=tuple(sorted(CIRCULAR_ROWS, key=str.casefold))),
    Case("circular-order-partial", "circular", row_order=("misc_feature", "CDS"),
         note="unlisted rows follow the listed ones in their order"),
    Case("circular-combined", "circular",
         row_names=(("CDS", "Genes"), ("rRNA", "Ribosomal RNA")), hidden_rows=("tRNA", "GC skew (-)"),
         row_order=("rRNA", "GC content", "CDS", "tRNA", "repeat_region", "misc_feature", "GC skew (+)", "GC skew (-)")),
    Case("circular-missing-key", "circular", row_names=(("No such row", "Renamed"),), hidden_rows=("No such row",),
         row_order=("No such row", "tRNA"), note="a key the table lacks is ignored"),
    Case("circular-comparison-ring", "circular", comparison=True,
         row_names=(("CDS", "Genes"),), hidden_rows=("rRNA",), row_order=("tRNA", "CDS"),
         note="the ring block keeps its rows; the edits reach the solid rows"),
    Case("linear-rename-plain", "linear", row_names=(("CDS", "Protein-coding genes"),)),
    Case("linear-rename-merges", "linear", row_names=(("misc_feature", "repeat_region"),)),
    Case("linear-rename-numbers", "linear", row_names=(("tRNA", "RNA genes"), ("rRNA", "RNA genes"))),
    Case("linear-delete-one", "linear", hidden_rows=("rRNA",)),
    Case("linear-delete-all-but-one", "linear", hidden_rows=LINEAR_ROWS[1:]),
    Case("linear-order-reverse", "linear", row_order=tuple(reversed(LINEAR_ROWS))),
    Case("linear-order-sort-a-z", "linear", row_order=tuple(sorted(LINEAR_ROWS, key=str.casefold))),
    Case("linear-combined", "linear", row_names=(("CDS", "Genes"),), hidden_rows=("misc_feature",),
         row_order=("repeat_region", "rRNA", "tRNA", "CDS", "misc_feature")),
    Case("linear-combined-top", "linear", legend="top", row_names=(("tRNA", "Transfer RNA"),), hidden_rows=("CDS",),
         row_order=tuple(reversed(LINEAR_ROWS)), note="a horizontal Legend"),
    Case("linear-pairwise-comparison", "linear", comparison=True, row_names=(("CDS", "Genes"),),
         hidden_rows=("tRNA",), row_order=("rRNA", "CDS")),
)


def _record(name: str) -> SeqRecord:
    rng = random.Random(name)
    record = SeqRecord(Seq("".join(rng.choice("ACGT") for _ in range(6000))), id=name, name=name)
    record.annotations["molecule_type"] = "DNA"
    for index, feature_type in enumerate(LAYOUT):
        location = FeatureLocation(100 + index * 700, 600 + index * 700, strand=1 if index % 2 else -1)
        record.features.append(SeqFeature(location, type=feature_type, qualifiers={"product": [f"p{index}"]}))
    return record


def _hits(query: str, subject: str) -> pd.DataFrame:
    rows = [
        [query, subject, identity, 400, 0, 0, 1 + 1500 * index, 400 + 1500 * index,
         1 + 1500 * index, 400 + 1500 * index, 1e-50, 500.0]
        for index, identity in enumerate((60.0, 80.0, 95.0))
    ]
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def _overrides(case: Case, edited: bool) -> dict[str, object]:
    if not edited:
        return {}
    leaves: dict[str, object] = {}
    if case.row_names:
        leaves[LEGEND_ROW_NAMES_PATH] = [list(pair) for pair in case.row_names]
    if case.hidden_rows:
        leaves[LEGEND_HIDDEN_ROWS_PATH] = list(case.hidden_rows)
    if case.row_order:
        leaves[LEGEND_ROW_ORDER_PATH] = list(case.row_order)
    return leaves


def _draw(case: Case, edited: bool, workdir: Path) -> str:
    features = gbdraw.FeatureOptions(types=list(FEATURE_TYPES))
    overrides = _overrides(case, edited)
    if case.mode == "circular":
        record = _record("recA")
        rings = None
        if case.comparison:
            rings = gbdraw.ComparisonRingOptions(
                tracks=[gbdraw.ComparisonRingTrackOptions(source=_hits("hitA", "recA"))], reference="subject"
            )
        options = gbdraw.CircularOptions(
            features=features, legend=case.legend, config_overrides=overrides, comparison_rings=rings,
            thresholds=gbdraw.Thresholds(identity=50, bitscore=0),
        )
        return gbdraw.draw_circular([record], options=options).to_svg()
    records = [_record("recA")]
    comparisons = gbdraw.LinearComparisonOptions()
    if case.comparison:
        records.append(_record("recB"))
        table = workdir / "pair.tsv"
        _hits("recA", "recB").to_csv(table, sep="\t", header=False, index=False)
        comparisons = gbdraw.LinearComparisonOptions(blast_files=[str(table)])
    options = gbdraw.LinearOptions(
        features=features, legend=case.legend, config_overrides=overrides, comparisons=comparisons,
        thresholds=gbdraw.Thresholds(identity=50, bitscore=0),
    )
    return gbdraw.draw_linear(records, options=options).to_svg()


_G_TAG = re.compile(r"<g\b[^>]*?(/?)>|</g>")


def legend_fragment(svg: str) -> str:
    """The ``#legend`` group of a drawn SVG as Python wrote it, in an ``<svg>`` of its own."""

    start = re.search(r'<g\b[^>]*\bid="legend"[^>]*>', svg)
    if start is None:
        raise ValueError("The drawing has no Legend.")
    depth = 0
    for match in _G_TAG.finditer(svg, start.start()):
        if match.group(0) == "</g>":
            depth -= 1
        elif not match.group(1):
            depth += 1
        if depth == 0:
            return f'<svg xmlns="http://www.w3.org/2000/svg">{svg[start.start():match.end()]}</svg>'
    raise ValueError("The Legend group is not closed.")


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _node(element: ET.Element) -> LegendNode:
    children = list(element)
    return {
        "tag": _local(element.tag),
        "attributes": {name: value for name, value in element.attrib.items() if name not in GEOMETRY_ATTRIBUTES},
        "text": None if children else (element.text or ""),
        "children": [_node(child) for child in children],
    }


def _by_id(root: ET.Element, element_id: str) -> ET.Element | None:
    return next((element for element in root.iter() if element.get("id") == element_id), None)


def feature_legend_groups(fragment: str) -> list[LegendGroupFacts]:
    """The feature Legend groups of ``fragment`` as the Web finds them (``getAllFeatureLegendGroups``), each
    with its rows: the direct child rows, else every row inside it."""

    root = ET.fromstring(fragment)
    legend = _by_id(root, "legend")
    if legend is None:
        return []
    def inner(outer: ET.Element, element_id: str) -> ET.Element:
        found = _by_id(outer, element_id)
        return outer if found is None else found

    horizontal, vertical = _by_id(legend, "legend_horizontal"), _by_id(legend, "legend_vertical")
    if horizontal is not None and vertical is not None:
        groups = [inner(horizontal, "feature_legend_h"), inner(vertical, "feature_legend_v")]
    else:
        groups = [inner(legend, "feature_legend")]
    facts: list[LegendGroupFacts] = []
    for group in groups:
        direct = [child for child in group if _local(child.tag) == "g" and "data-legend-key" in child.attrib]
        rows = direct or [element for element in group.iter() if _local(element.tag) == "g"
                          and "data-legend-key" in element.attrib and element is not group]
        facts.append({
            "attributes": {name: value for name, value in group.attrib.items() if name not in GEOMETRY_ATTRIBUTES},
            "rows": [_node(row) for row in rows],
        })
    return facts


def web_operations(case: Case) -> dict[str, list[dict[str, object]]]:
    """The Web editor operations that make the case's edits (the projection Generate's leaves come from)."""

    return {
        "legendRenames": [{"from": key, "to": name, "allowMissing": True} for key, name in case.row_names],
        "legendDeletes": [{"caption": key, "allowMissing": True} for key in case.hidden_rows],
        "legendOrder": [{"captions": list(case.row_order)}] if case.row_order else [],
    }


def vector(case: Case, workdir: Path) -> dict[str, object]:
    unedited = legend_fragment(_draw(case, False, workdir))
    edited = legend_fragment(_draw(case, True, workdir))
    return {
        "id": case.id,
        "mode": case.mode,
        "note": case.note,
        "edits": {
            "row_names": [list(pair) for pair in case.row_names],
            "hidden_rows": list(case.hidden_rows),
            "row_order": list(case.row_order),
        },
        "operations": web_operations(case),
        "unedited": unedited,
        "edited": edited,
        "expected": feature_legend_groups(edited),
    }


def render() -> str:
    # One case per line keeps a regeneration's diff readable case by case.
    with tempfile.TemporaryDirectory() as directory:
        cases = [vector(case, Path(directory)) for case in CASES]
    compact = lambda value: json.dumps(value, ensure_ascii=False, separators=(",", ":"))  # noqa: E731
    lines = ["{", f' "description": {compact(__doc__.strip().splitlines()[0])},', ' "cases": [']
    lines.extend(f"  {compact(case)}{',' if index < len(cases) - 1 else ''}" for index, case in enumerate(cases))
    lines.extend([" ]", "}"])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Fail if the vectors are stale.")
    args = parser.parse_args()
    expected = render()
    if args.check:
        if not TARGET.is_file() or TARGET.read_text(encoding="utf-8") != expected:
            print("Legend row edit vectors are stale. Run: python tools/generate_legend_row_edit_vectors.py")
            return 1
        return 0
    TARGET.write_text(expected, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
