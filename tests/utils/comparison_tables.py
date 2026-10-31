"""Read the canonical comparison table of a converted protein pair."""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from types import SimpleNamespace
from typing import Any

import pandas as pd

from gbdraw.comparisons.identity_scale import identity_color_domain
from gbdraw.layout.linear_multi_record import LinearComparisonEndpoint
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup


def converted_pair_rows(pair: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The rows of a converted pair's ``canonicalTsv``, read as the request decoder
    reads a comparison table (pandas dtypes; an empty cell is NaN). The pair
    carries no ``rows`` copy (D-66)."""

    assert "rows" not in pair
    return pd.read_csv(io.StringIO(pair["canonicalTsv"]), sep="\t").to_dict(orient="records")


def ribbon_placement(
    record_index: int, *, px_per_bp: float, x: float = 0.0, anchor_y: float = 0.0
) -> LinearComparisonEndpoint:
    """A record's end of comparison ribbons drawn outside the assembler: record ``n`` on row ``n`` at
    ``x``, scaled by ``px_per_bp``, its ribbon ends at ``anchor_y``. Ribbons place their ends only
    through it."""

    return LinearComparisonEndpoint(
        record_index=record_index,
        row=record_index,
        x=x,
        px_per_bp=px_per_bp,
        anchor_y=anchor_y,
        attachment_offsets_px={},
    )


def ribbon_config(**overrides: Any) -> SimpleNamespace:
    """The ``BlastMatchConfigurator`` attributes ``PairWiseMatchGroup`` reads, with ``overrides``."""

    values: dict[str, Any] = {
        "fill_color": "#d3d3d3",
        "identity": 0,
        "identity_domain": identity_color_domain(0, None),
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
        1,
        ribbon_config(**config),
        records,
        feature_dom_index=feature_dom_index,
        record_transforms=None,
        query_placement=ribbon_placement(query_index, px_per_bp=px_per_bp),
        subject_placement=ribbon_placement(subject_index, px_per_bp=px_per_bp, x=subject_x, anchor_y=height),
    )
