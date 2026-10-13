"""Per-Result Legend row facts: which rows the Legend drew, which the draft removed.

A Legend style or rename is stored against a row key (a feature type, an
``other ...`` row, a color rule's key, or a track's key), never against the name
the Legend draws. The draft can remove that row without retiring the preference:
hiding every feature of the row, recaptioning them with a color rule, or giving a
track the name and swatch of another row leaves Python with nothing to draw under
that key. The Web admits such a Generate only when Python reports the key as
*suppressed*; a key that is neither drawn nor suppressed is stale. Python owns the
answer because it owns the row derivation (``prepare_legend_table``).
"""

from __future__ import annotations

from collections.abc import Sequence
from functools import cached_property
from typing import TYPE_CHECKING, Any

from Bio.SeqRecord import SeqRecord

from gbdraw.features.colors import precompute_used_color_rules
from gbdraw.legend.table import (
    LegendRowKey,
    LegendTable,
    _key_of_name,
    _other_feature_legend_key,
    prepare_legend_table,
)

if TYPE_CHECKING:
    from gbdraw.configurators import (
        FeatureDrawingConfigurator,
        GcContentConfigurator,
        GcSkewConfigurator,
    )

LEGEND_ROW_FACTS_ATTR = "_gbdraw_legend_row_facts"


class LegendRowFacts:
    """The Legend rows one drawing drew, and the feature rows its draft removed.

    ``suppressed`` is derived on first read, so a render that never asks pays nothing.
    """

    def __init__(
        self,
        *,
        records: Sequence[SeqRecord],
        feature_config: FeatureDrawingConfigurator,
        gc_config: GcContentConfigurator,
        skew_config: GcSkewConfigurator,
        legend_table: LegendTable,
    ) -> None:
        self._records = tuple(records)
        self._feature_config = feature_config
        self._gc_config = gc_config
        self._skew_config = skew_config
        self.drawn: tuple[LegendRowKey, ...] = tuple(legend_table)
        self._shared = legend_table.shared_keys

    @cached_property
    def suppressed(self) -> tuple[LegendRowKey, ...]:
        """Rows the records can name that the draft did not draw.

        The rows are those of every feature type of the records (hidden or not
        selected included): the type row, its ``other`` row, and every row
        ``prepare_legend_table`` yields when no feature is hidden, plus the GC
        content and GC skew rows (a switched-off track draws none), plus the
        keys of the sources whose row another source of the same name and
        swatch draws.
        """

        feature_config = self._feature_config
        feature_types = list(
            dict.fromkeys(
                str(feature.type) for record in self._records for feature in record.features
            )
        )
        used_rules, default_used = precompute_used_color_rules(
            self._records,
            feature_config.specific_color_rules,
            feature_config.default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
            record_features=feature_config.record_features,
            include_hidden=True,
        )
        producible = prepare_legend_table(
            self._gc_config,
            self._skew_config,
            feature_config,
            feature_types,
            used_color_rules=used_rules,
            default_used_features=default_used,
            show_gc=True,
            show_skew=True,
            show_depth=False,
        )
        rows: dict[LegendRowKey, None] = dict.fromkeys(
            [
                *producible,
                *map(_key_of_name, feature_types),
                *(_other_feature_legend_key(feature_type) for feature_type in feature_types),
                *self._shared,
            ]
        )
        drawn = set(self.drawn)
        return tuple(row for row in rows if row not in drawn)


def attach_legend_row_facts(canvas: object, facts: LegendRowFacts) -> None:
    """Record the facts on the drawing, as the track-slot geometry is recorded."""

    setattr(canvas, LEGEND_ROW_FACTS_ATTR, facts)


def collect_legend_row_facts(
    canvas: object,
    *,
    result_index: int,
    result_name: str,
) -> dict[str, Any]:
    """Serialize one Result's Legend row facts for the browser run metadata."""

    facts = getattr(canvas, LEGEND_ROW_FACTS_ATTR, None)
    return {
        "resultIndex": int(result_index),
        "resultName": str(result_name),
        "drawn": list(facts.drawn) if isinstance(facts, LegendRowFacts) else [],
        "suppressed": list(facts.suppressed) if isinstance(facts, LegendRowFacts) else [],
    }


__all__ = [
    "LegendRowFacts",
    "attach_legend_row_facts",
    "collect_legend_row_facts",
]
