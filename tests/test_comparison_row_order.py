"""G-F(5) (Web GUI audit 2026-09-30, R5): the drawn comparison does not depend
on the row order of its BLAST table.

A table and the same rows in another order draw the same matches with the same
IDs and stacking. A swapped query/subject table is rejected
(tests/test_comparison_tables.py). The last test shows that the check fails when
matches are drawn in input order.
"""

from __future__ import annotations

import random

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import LinearComparison
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.render.groups.linear import pairwise_match

# Distinct spans, both orientations, overlapping matches with different identities.
_ROWS = [
    ("Q1", "S1", 99.0, 900, 0, 0, 101, 1000, 201, 1100, 1e-100, 1600.0),
    ("Q1", "S1", 82.5, 400, 0, 0, 1201, 1600, 1900, 1501, 1e-40, 500.0),
    ("Q1", "S1", 91.0, 300, 0, 0, 1701, 2000, 2101, 2400, 1e-60, 520.0),
    ("Q1", "S1", 75.0, 200, 0, 0, 2201, 2400, 2800, 2601, 1e-20, 210.0),
    ("Q1", "S1", 88.0, 150, 0, 0, 2451, 2600, 151, 300, 1e-30, 260.0),
    ("Q1", "S1", 95.0, 500, 0, 0, 601, 1100, 701, 1200, 1e-80, 880.0),
]


def _records() -> list[SeqRecord]:
    records = []
    for seed, record_id in enumerate(("Q1", "S1")):
        generator = random.Random(seed)
        record = SeqRecord(Seq("".join(generator.choice("ACGT") for _ in range(3000))), id=record_id)
        record.annotations["molecule_type"] = "DNA"
        records.append(record)
    return records


def _render(rows: list[tuple]) -> str:
    table = pd.DataFrame(rows, columns=COMPARISON_COLUMNS)
    return assemble_linear_diagram_from_records(
        _records(),
        cfg=apply_config_overrides(
            None,
            {"labels.linear.scope": "none", "canvas.show_gc": False, "canvas.show_skew": False},
        ),
        linear_comparisons=[LinearComparison(0, 1, table)],
        legend="none",
    ).tostring()


def _reorderings() -> list[list[tuple]]:
    orders = [list(reversed(_ROWS))]
    for seed in (1, 2):
        shuffled = list(_ROWS)
        random.Random(seed).shuffle(shuffled)
        orders.append(shuffled)
    return orders


@pytest.mark.linear
def test_reordered_comparison_rows_draw_the_same_matches() -> None:
    expected = _render(_ROWS)
    assert expected.count("data-qstart=") == len(_ROWS)
    for rows in _reorderings():
        assert _render(rows) == expected


@pytest.mark.linear
def test_row_order_check_fails_when_matches_are_drawn_in_input_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(pairwise_match, "_match_draw_order_key", lambda row: 0)
    expected = _render(_ROWS)
    assert any(_render(rows) != expected for rows in _reorderings())
