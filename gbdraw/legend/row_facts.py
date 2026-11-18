"""Per-Result Legend row facts: which rows the Legend drew, which the draft removed.

A Legend style or rename is stored against a row key (a feature type, an
``other ...`` row, a color rule's key, or a track's key), never against the name
the Legend draws. The draft can remove that row without retiring the preference:
hiding every feature of the row, recaptioning them with a color rule, or giving a
track the name and swatch of another row leaves Python with nothing to draw under
that key. So does a record selection that draws a record OFF: its rows wait for
it. The Web admits such a Generate only when Python reports the key as
*suppressed*; a key that is neither drawn nor suppressed is stale. Python owns the
answer because it owns the row derivation (``prepare_legend_table``).
"""

from __future__ import annotations

import copy
from collections.abc import Sequence
from functools import cached_property
from typing import NamedTuple, TypedDict, TypeVar

from Bio.SeqRecord import SeqRecord

from gbdraw.features.colors import precompute_used_color_rules
from gbdraw.io.genome import filter_features_by_type
from gbdraw.legend.sources import LegendSources
from gbdraw.legend.table import (
    LegendRowKey,
    LegendTable,
    _key_of_name,
    _other_feature_legend_key,
    prepare_legend_table,
)
from gbdraw.render_stages import render_stage

LEGEND_ROW_FACTS_ATTR = "_gbdraw_legend_row_facts"
_Canvas = TypeVar("_Canvas")


class LeftOutRecord(NamedTuple):
    """A record of a drawing's files that its record selection leaves out, as the record load parsed it (shared:
    read it, never change it)."""

    record: SeqRecord
    # A GFF3 record is drawn with its features flattened to the types the request draws, as the load does.
    gff3: bool
    # Those types; None keeps every type (still flattened).
    gff_types: frozenset[str] | None

    def as_loaded(self) -> SeqRecord:
        """The record with the features the request would load (no feature identity: rows need none)."""

        return filter_features_by_type(self.record, self.gff_types) if self.gff3 else self.record


class LegendRowFacts:
    """The Legend rows one drawing drew, and the feature rows its draft removed.

    ``suppressed`` is derived on first read, so a render that never asks pays nothing.
    """

    def __init__(
        self,
        *,
        records: Sequence[SeqRecord],
        sources: LegendSources,
        legend_table: LegendTable,
    ) -> None:
        self._records = tuple(records)
        # The records of the drawing's files its record selection leaves out (``cover_records_left_out``).
        self._left_out: tuple[LeftOutRecord, ...] = ()
        self._sources = sources
        # A row a Legend row edit hides is drawn too: the Result holds it, hidden.
        self.drawn: tuple[LegendRowKey, ...] = legend_table.drawn_keys()
        self._shared = legend_table.shared_keys

    @cached_property
    def suppressed(self) -> tuple[LegendRowKey, ...]:
        """Rows the records can name that the draft did not draw.

        The rows are those of :func:`producible_legend_rows` for the drawn
        records and the records left out, plus the keys of the sources whose
        row another source of the same name and swatch draws.
        """

        producible = render_stage(
            "legend_table",
            producible_drawing_rows,
            self._records,
            self._left_out,
            self._sources,
        )
        rows = dict.fromkeys([*producible, *self._shared])
        drawn = set(self.drawn)
        return tuple(row for row in rows if row not in drawn)


def producible_drawing_rows(
    records: tuple[SeqRecord, ...],
    left_out: tuple[LeftOutRecord, ...],
    sources: LegendSources,
) -> tuple[LegendRowKey, ...]:
    """The rows of :func:`producible_legend_rows` for the drawn records and the records left out (render stage
    ``legend_table``): a left-out GFF3 record is flattened here, so a stage hit does no feature work."""

    return producible_legend_rows((*records, *(entry.as_loaded() for entry in left_out)), sources)


def producible_legend_rows(
    records: tuple[SeqRecord, ...],
    sources: LegendSources,
) -> tuple[LegendRowKey, ...]:
    """Every feature row the records can name (render stage ``legend_table``).

    The rows of every feature type of the records (hidden or not selected
    included): the type row, its ``other`` row, and every row
    ``prepare_legend_table`` yields when no feature is hidden, plus the GC
    content and GC skew rows (a switched-off track draws none).
    """

    feature_types = list(
        dict.fromkeys(
            str(feature.type) for record in records for feature in record.features
        )
    )
    features = sources.features
    # Hidden features count, so no feature override is read: a record left out has none.
    used_rules, default_used = precompute_used_color_rules(
        records,
        features.specific_color_rules,
        features.default_color_map,
        set(features.selected_features),
        feature_visibility_rules=features.feature_visibility_rules,
        include_hidden=True,
    )
    producible = prepare_legend_table(
        sources.gc,
        sources.skew,
        features,
        feature_types,
        used_color_rules=used_rules,
        default_used_features=default_used,
        show_gc=True,
        show_skew=True,
        show_depth=False,
    )
    return tuple(
        dict.fromkeys(
            [
                *producible,
                *map(_key_of_name, feature_types),
                *(_other_feature_legend_key(feature_type) for feature_type in feature_types),
            ]
        )
    )


def attach_legend_row_facts(canvas: object, facts: LegendRowFacts) -> None:
    """Record the facts on the drawing, as the track-slot geometry is recorded."""

    setattr(canvas, LEGEND_ROW_FACTS_ATTR, facts)


def cover_records_left_out(canvas: _Canvas, records: Sequence[LeftOutRecord]) -> _Canvas:
    """Let the facts on the drawing name the rows of ``records``, the records of its files that the record
    selection leaves out (the records a Web drawing draws OFF), as rows the draft did not draw; the drawn rows
    and the drawing stay as they are. Returns ``canvas``."""

    facts = getattr(canvas, LEGEND_ROW_FACTS_ATTR, None)
    if isinstance(facts, LegendRowFacts) and records:
        covered = copy.copy(facts)
        covered.__dict__.pop("suppressed", None)
        covered._left_out = tuple(records)
        attach_legend_row_facts(canvas, covered)
    return canvas


class LegendRowFactsPayload(TypedDict):
    """One Result's Legend row facts in a Web render's run metadata (``metadata.legendRows``); not saved in a
    Session."""

    resultIndex: int
    resultName: str
    drawn: list[LegendRowKey]
    suppressed: list[LegendRowKey]


def collect_legend_row_facts(
    canvas: object,
    *,
    result_index: int,
    result_name: str,
) -> LegendRowFactsPayload:
    """Serialize one Result's Legend row facts for the browser run metadata."""

    facts = getattr(canvas, LEGEND_ROW_FACTS_ATTR, None)
    return {
        "resultIndex": int(result_index),
        "resultName": str(result_name),
        "drawn": list(facts.drawn) if isinstance(facts, LegendRowFacts) else [],
        "suppressed": list(facts.suppressed) if isinstance(facts, LegendRowFacts) else [],
    }


__all__ = [
    "LeftOutRecord",
    "LegendRowFacts",
    "LegendRowFactsPayload",
    "attach_legend_row_facts",
    "collect_legend_row_facts",
    "cover_records_left_out",
    "producible_drawing_rows",
    "producible_legend_rows",
]
