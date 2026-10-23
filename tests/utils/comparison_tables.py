"""Read the canonical comparison table of a converted protein pair."""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from types import SimpleNamespace
from typing import Any

import pandas as pd

from gbdraw.layout.linear_multi_record import LinearRecordPlacement, RecordKey
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup


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


def ribbon_config(**overrides: Any) -> SimpleNamespace:
    """The ``BlastMatchConfigurator`` attributes ``PairWiseMatchGroup`` reads, with ``overrides``."""

    values: dict[str, Any] = {
        "fill_color": "#d3d3d3",
        "identity": 0,
        "min_color": "#ffffff",
        "max_color": "#000000",
        "fill_opacity": 0.75,
        "stroke_color": "none",
        "stroke_width": 0,
        "collinearity_orientation_colors": {},
        "collinearity_orientation_min_colors": {},
        "match_style": "ribbon",
        "curve_tension": 0.5,
        "sequence_length_dict": {},
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def ribbon_group(
    records: Sequence[Any],
    comparison: pd.DataFrame | None = None,
    *,
    px_per_bp: float = 1.0,
    height: float = 40.0,
    query_index: int = 0,
    subject_index: int = 1,
    subject_x: float = 0.0,
    feature_dom_index: Any = None,
    **config: Any,
) -> PairWiseMatchGroup:
    """A ``PairWiseMatchGroup`` built as the builder builds it, between record ``query_index`` and
    ``subject_index`` (``ribbon_placement``), ribbon ends at 0 and ``height``, drawing ``comparison``
    (no rows by default); ``config`` overrides ``ribbon_config``."""

    return PairWiseMatchGroup(
        {},
        comparison if comparison is not None else pd.DataFrame(columns=["identity", "qstart", "qend", "sstart", "send"]),
        height,
        1,
        ribbon_config(**config),
        records,
        query_y=0.0,
        subject_y=height,
        feature_dom_index=feature_dom_index,
        record_transforms=None,
        query_placement=ribbon_placement(query_index, px_per_bp=px_per_bp),
        subject_placement=ribbon_placement(subject_index, px_per_bp=px_per_bp, x=subject_x),
    )
