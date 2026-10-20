"""Read the canonical comparison table of a converted protein pair."""

from __future__ import annotations

import io
from collections.abc import Mapping
from typing import Any

import pandas as pd

from gbdraw.layout.linear_multi_record import LinearRecordPlacement, RecordKey


def converted_pair_rows(pair: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The rows of a converted pair's ``canonicalTsv``, read as the request decoder
    reads a comparison table (pandas dtypes; an empty cell is NaN). The pair
    carries no ``rows`` copy (D-66)."""

    assert "rows" not in pair
    return pd.read_csv(io.StringIO(pair["canonicalTsv"]), sep="\t").to_dict(orient="records")


def ribbon_placement(record_index: int, *, px_per_bp: float, x: float = 0.0) -> LinearRecordPlacement:
    """A record's placement for comparison ribbons drawn outside the assembler: record ``n`` on
    row ``n`` at ``x``, scaled by ``px_per_bp``. Ribbons place their ends only through it."""

    return LinearRecordPlacement(
        record_index=record_index,
        record_key=RecordKey(f"record-{record_index + 1}"),
        row=record_index,
        column=0,
        x=x,
        axis_y=0.0,
        sequence_width=0.0,
        left_inset=0.0,
        right_inset=0.0,
        top_extent=0.0,
        bottom_extent=0.0,
        comparison_top_y=0.0,
        comparison_bottom_y=0.0,
        px_per_bp=px_per_bp,
    )
