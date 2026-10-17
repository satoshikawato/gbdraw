"""Read the canonical comparison table of a converted protein pair."""

from __future__ import annotations

import io
from collections.abc import Mapping
from typing import Any

import pandas as pd


def converted_pair_rows(pair: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The rows of a converted pair's ``canonicalTsv``, read as the request decoder
    reads a comparison table (pandas dtypes; an empty cell is NaN). The pair
    carries no ``rows`` copy (D-66)."""

    assert "rows" not in pair
    return pd.read_csv(io.StringIO(pair["canonicalTsv"]), sep="\t").to_dict(orient="records")
