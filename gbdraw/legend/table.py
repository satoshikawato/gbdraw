#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Callable, Hashable, Iterator, Mapping, MutableMapping
from typing import Iterable, List, NamedTuple, Optional, Sequence, Set, Tuple

from pandas import DataFrame

from gbdraw.core.color import (
    DEFAULT_COLLINEAR_ORIENTATION_COLORS,
    DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
    normalize_hex_color,
)
from gbdraw.exceptions import ValidationError


LegendEntry = Mapping[str, object]


def _unique_legend_key(reserved: Iterable[str], preferred: str) -> str:
    """A row key no reserved key equals: ``preferred``, else ``preferred (2)``, ``(3)``, ..."""

    taken = set(reserved)
    if preferred not in taken:
        return preferred
    suffix = 2
    while f"{preferred} ({suffix})" in taken:
        suffix += 1
    return f"{preferred} ({suffix})"


def _legend_fill_identity(color: object) -> str:
    """Return the color identity used to compare two legend colors (a name equals its hex)."""

    from gbdraw.io.colors import resolve_color_to_hex

    text = str(color or "").strip()
    if text.lower() == "none":
        return "none"
    try:
        return normalize_hex_color(resolve_color_to_hex(text))
    except (ValueError, ValidationError):
        return text.lower()


def legend_swatch(entry: LegendEntry) -> tuple[object, ...]:
    """What two legend rows compare to be "the same swatch".

    A solid row is its fill and its stroke (color and width; a stroke of no
    width or no color is no stroke); a hatch is part of the fill (its pattern
    paint). A gradient row is its two end colors.
    """

    if entry.get("type") == "gradient":
        return (
            "gradient",
            _legend_fill_identity(entry.get("min_color")),
            _legend_fill_identity(entry.get("max_color")),
        )
    try:
        width = float(str(entry.get("width") or 0))
    except ValueError:
        width = 0.0
    stroke = _legend_fill_identity(entry.get("stroke") or "none")
    if width <= 0 or stroke == "none":
        return ("solid", _legend_fill_identity(entry.get("fill")), "none", 0.0)
    return ("solid", _legend_fill_identity(entry.get("fill")), stroke, width)


def drawn_legend_names(rows: Sequence[tuple[str, Hashable]]) -> list[str]:
    """THE naming rule: the name each legend row is drawn with.

    ``rows`` are the rows being drawn, in legend order, each as its name and
    its swatch. Rows of one name and one swatch are one row and share a drawn
    name. A name left with one row is drawn as written. A name left with two
    or more rows is drawn ``name (1)``, ``name (2)``, ... in legend order; a
    number skips a drawn name that another row already has.
    """

    swatches: dict[str, list[Hashable]] = {}
    for name, swatch in rows:
        group = swatches.setdefault(name, [])
        if swatch not in group:
            group.append(swatch)
    taken = {name for name, group in swatches.items() if len(group) == 1}
    numbered: dict[tuple[str, Hashable], str] = {}
    for name, group in swatches.items():
        if len(group) == 1:
            continue
        number = 0
        for swatch in group:
            number += 1
            while f"{name} ({number})" in taken:
                number += 1
            numbered[name, swatch] = f"{name} ({number})"
            taken.add(numbered[name, swatch])
    return [numbered.get((name, swatch), name) for name, swatch in rows]


class LegendRowNames(NamedTuple):
    """A legend row's name as written, and the name the legend draws."""

    name: str
    drawn: str


class _LegendRow(NamedTuple):
    name: str
    swatch: tuple[object, ...]
    entry: LegendEntry


class LegendTable(Mapping[str, LegendEntry]):
    """The legend rows of one drawing, in legend order, by row key.

    A row has three separate facts. Its *key* is its identity (the Web stores
    row edits and Sessions under it, ``data-legend-key``); it is never drawn.
    Its *name* is what the user or the renderer wrote and is never rewritten.
    Its *drawn name* is derived from the names and swatches of the rows being
    drawn (``drawn_legend_names``). Read as a mapping, the table gives each
    key's entry (the swatch properties the legend draws).
    """

    def __init__(self, rows: Mapping[object, LegendEntry] | None = None) -> None:
        self._rows: dict[str, _LegendRow] = {}
        self._shared_keys: list[str] = []
        for key, entry in (rows or {}).items():
            self.add(str(key), entry)

    def add(self, name: str, entry: LegendEntry, *, key: str | None = None) -> None:
        """Add a row under its key (its name unless given); one of the same name and swatch is shared."""

        self._add(name, entry, name if key is None else key)

    def add_source(self, name: str, entry: LegendEntry) -> None:
        """Add the row of a track or annotation source, which has no key of its own.

        Its key is its name, numbered ``(2)``, ``(3)``, ... after the keys taken
        before it. A source that shares a row still takes its key, so the keys
        of the other sources do not depend on which rows are shared.
        """

        self._add(name, entry, _unique_legend_key([*self._rows, *self._shared_keys], name))

    def _add(self, name: str, entry: LegendEntry, key: str) -> None:
        swatch = legend_swatch(entry)
        if any(row.name == name and row.swatch == swatch for row in self._rows.values()):
            if key not in self._rows:
                self._shared_keys.append(key)
            return
        self._rows[key] = _LegendRow(name, swatch, entry)

    def discard(self, key: str) -> None:
        self._rows.pop(key, None)

    @property
    def shared_keys(self) -> tuple[str, ...]:
        """Keys of the sources whose row another source of the same name and swatch draws."""

        return tuple(key for key in dict.fromkeys(self._shared_keys) if key not in self._rows)

    def names(self) -> dict[str, LegendRowNames]:
        """Each row's name and drawn name, by key, in legend order."""

        rows = list(self._rows.values())
        drawn = drawn_legend_names([(row.name, row.swatch) for row in rows])
        return {
            key: LegendRowNames(row.name, text)
            for (key, row), text in zip(self._rows.items(), drawn)
        }

    def __getitem__(self, key: str) -> LegendEntry:
        return self._rows[key].entry

    def __iter__(self) -> Iterator[str]:
        return iter(self._rows)

    def __len__(self) -> int:
        return len(self._rows)


def legend_row_names(legend_table: Mapping[object, LegendEntry]) -> dict[str, LegendRowNames]:
    """The names the legend layouts measure and the legend drawers write, by row key.

    A plain mapping (a test's table) is read as rows named by their keys.
    """

    table = legend_table if isinstance(legend_table, LegendTable) else LegendTable(legend_table)
    return table.names()


def write_legend_row_names(
    legend_table: Mapping[object, LegendEntry] | None,
) -> Callable[[MutableMapping[str, object], str], str]:
    """The legend drawers' one writer of a row's key, name, and drawn name.

    The returned function writes ``data-legend-key`` and ``data-legend-name``
    on a row group's attributes and returns the name the legend draws.
    """

    names = legend_row_names(legend_table or {})

    def write(attribs: MutableMapping[str, object], key: str) -> str:
        attribs["data-legend-key"] = key
        attribs["data-legend-name"] = names[key].name
        return names[key].drawn

    return write


def rule_legend_keys(
    rules: Iterable[tuple[object, object]],
    renderer_fills: Mapping[str, str],
) -> dict[tuple[str, str], str]:
    """Row keys for the captions of the color table, by caption and fill identity.

    A rule row's key is its caption. A caption written with several colors in
    the table, or naming a renderer row of another color (``renderer_fills``),
    takes ``caption [fill]`` instead, made unique against every caption and
    renderer row. The allocation reads the whole table, never its row order.
    The key is an identity only: the legend draws the caption
    (``drawn_legend_names``).
    """

    fills: dict[str, set[str]] = {}
    for caption, color in rules:
        if str(caption or "").strip():
            fills.setdefault(str(caption), set()).add(_legend_fill_identity(color))
    reserved = {*renderer_fills, *fills}
    keys: dict[tuple[str, str], str] = {}
    for caption in sorted(fills):
        renderer_fill = renderer_fills.get(caption)
        for fill in sorted(fills[caption]):
            if len(fills[caption]) > 1 or renderer_fill not in (None, fill):
                key = _unique_legend_key(reserved, f"{caption} [{fill}]")
                reserved.add(key)
            else:
                key = caption
            keys[caption, fill] = key
    return keys


def _normalize_collinearity_legend_color_mode(value: object) -> str | None:
    # Imported here: the legend layouts read this module's row names, and
    # gbdraw.analysis imports the configurators that import those layouts.
    from gbdraw.analysis.collinearity import normalize_collinearity_color_mode

    if not str(value or "").strip():
        return None
    try:
        return normalize_collinearity_color_mode(str(value))
    except ValidationError:
        return None


def collect_comparison_collinearity_color_modes(
    comparison_dataframes: Sequence[DataFrame] | None,
) -> set[str]:
    """Collect recognized collinearity color modes from comparison metadata."""

    color_modes: set[str] = set()
    if not comparison_dataframes:
        return color_modes
    for frame in comparison_dataframes:
        if "collinearity_color_mode" not in frame.columns:
            continue
        for value in frame["collinearity_color_mode"].dropna():
            normalized = _normalize_collinearity_legend_color_mode(value)
            if normalized is not None:
                color_modes.add(normalized)
    return color_modes


def configure_pairwise_identity_legend_from_comparisons(
    blast_config,
    comparison_dataframes: Sequence[DataFrame] | None,
    *,
    additional_color_modes: Iterable[object] | None = None,
) -> set[str]:
    """Configure pairwise identity legend labels from collinearity metadata."""

    color_modes = collect_comparison_collinearity_color_modes(comparison_dataframes)
    for value in additional_color_modes or ():
        normalized = _normalize_collinearity_legend_color_mode(value)
        if normalized is not None:
            color_modes.add(normalized)

    if not color_modes:
        return color_modes

    if hasattr(blast_config, "pairwise_identity_legend_entries"):
        delattr(blast_config, "pairwise_identity_legend_entries")
    if hasattr(blast_config, "pairwise_identity_legend_label"):
        delattr(blast_config, "pairwise_identity_legend_label")

    blast_config.hide_pairwise_identity_legend = color_modes == {"orientation"}
    if "orientation_identity" in color_modes:
        orientation_min_colors = getattr(
            blast_config,
            "collinearity_orientation_min_colors",
            DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
        )
        orientation_colors = getattr(
            blast_config,
            "collinearity_orientation_colors",
            DEFAULT_COLLINEAR_ORIENTATION_COLORS,
        )
        blast_config.pairwise_identity_legend_entries = [
            {
                "label": "Collinear",
                "min_color": orientation_min_colors.get(
                    "plus", DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS["plus"]
                ),
                "max_color": orientation_colors.get(
                    "plus", DEFAULT_COLLINEAR_ORIENTATION_COLORS["plus"]
                ),
            },
            {
                "label": "Inverted",
                "min_color": orientation_min_colors.get(
                    "minus", DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS["minus"]
                ),
                "max_color": orientation_colors.get(
                    "minus", DEFAULT_COLLINEAR_ORIENTATION_COLORS["minus"]
                ),
            },
        ]
    elif "average_identity" in color_modes:
        blast_config.pairwise_identity_legend_label = "Average identity"

    return color_modes


def _other_feature_legend_key(feature_type: str) -> str:
    return "other proteins" if feature_type == "CDS" else f"other {feature_type}s"


def _resolve_legend_feature_color(default_colors: DataFrame, feature_type: str) -> str:
    matching_rows = default_colors[default_colors["feature_type"] == feature_type]
    if not matching_rows.empty:
        return str(matching_rows["color"].values[0])

    default_rows = default_colors[default_colors["feature_type"] == "default"]
    if not default_rows.empty:
        return str(default_rows["color"].values[0])

    return "#d3d3d3"


def _generated_legend_fills(
    feature_specific_colors: dict,
    features_present: Sequence[str],
    default_colors: DataFrame,
    *,
    used_color_rules: Optional[Set[Tuple[str, str]]],
    default_used_features: Set[str],
    gc_config,
    skew_config,
    depth_config,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
) -> dict[str, str]:
    """Return the renderer-generated legend rows and their fill identities."""

    fills: dict[str, str] = {}
    for feature_type in features_present:
        feature_color = _legend_fill_identity(
            _resolve_legend_feature_color(default_colors, feature_type)
        )
        if feature_type not in feature_specific_colors:
            fills[feature_type] = feature_color
        elif used_color_rules is None:
            fills[_other_feature_legend_key(feature_type)] = feature_color
        elif any(
            (caption, color) in used_color_rules
            for caption, color in feature_specific_colors[feature_type]
            if str(caption or "").strip()
        ):
            if feature_type in default_used_features:
                fills[_other_feature_legend_key(feature_type)] = feature_color
        else:
            fills[feature_type] = feature_color
    dinucleotide = gc_config.dinucleotide
    if show_depth and depth_config is not None:
        fills["Depth"] = _legend_fill_identity(depth_config.fill_color)
    for show, config, name in ((show_gc, gc_config, "content"), (show_skew, skew_config, "skew")):
        if not show:
            continue
        if config.high_fill_color == config.low_fill_color:
            fills[f"{dinucleotide} {name}"] = _legend_fill_identity(config.high_fill_color)
        else:
            fills[f"{dinucleotide} {name} (+)"] = _legend_fill_identity(config.high_fill_color)
            fills[f"{dinucleotide} {name} (-)"] = _legend_fill_identity(config.low_fill_color)
    return fills


def prepare_legend_table(
    gc_config,
    skew_config,
    feature_config,
    features_present: List[str],
    blast_config=None,
    has_blast: bool = False,
    used_color_rules: Optional[Set[Tuple[str, str]]] = None,
    default_used_features: Optional[Set[str]] = None,
    depth_config=None,
    *,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
) -> LegendTable:
    """
    Prepare the legend table for the diagram.

    Args:
        used_color_rules: Optional set of (caption, color) tuples that were actually
            matched during feature creation. If provided, only these rules will be
            included in the legend. If None, all rules for present feature types
            will be included (legacy behavior).
        default_used_features: Optional set of feature types that fell back to
            default colors. Used to decide whether to include "other X" entries
            when specific rules exist.
    """
    legend_table = LegendTable()
    color_table: Optional[DataFrame] = feature_config.color_table
    default_colors: DataFrame = feature_config.default_colors
    block_stroke_color: str = feature_config.block_stroke_color
    block_stroke_width: float = feature_config.block_stroke_width
    gc_stroke_color: str = gc_config.stroke_color
    gc_stroke_width: float = gc_config.stroke_width
    gc_high_fill_color: str = gc_config.high_fill_color
    gc_low_fill_color: str = gc_config.low_fill_color
    skew_high_fill_color: str = skew_config.high_fill_color
    skew_low_fill_color: str = skew_config.low_fill_color
    skew_stroke_color: str = skew_config.stroke_color
    skew_stroke_width: float = skew_config.stroke_width
    dinucleotide = gc_config.dinucleotide
    feature_specific_colors: dict[str, list[tuple[str, str]]] = {}
    default_used_features = default_used_features or set()
    if color_table is not None and not color_table.empty:
        for _, row in color_table.iterrows():
            feature_type = row["feature_type"]
            if feature_type not in feature_specific_colors:
                feature_specific_colors[feature_type] = []
            feature_specific_colors[feature_type].append((row["caption"], row["color"]))
    # The renderer's own rows (feature types, "other ..." and numeric tracks)
    # keep their keys; a rule row takes its key from the whole color table and
    # these rows (``rule_legend_keys``). The legend draws names, never keys.
    generated_fills = _generated_legend_fills(
        feature_specific_colors,
        features_present,
        default_colors,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=show_depth,
    )
    rule_keys = rule_legend_keys(
        (entry for entries in feature_specific_colors.values() for entry in entries),
        generated_fills,
    )
    for selected_feature in features_present:
        if selected_feature in feature_specific_colors.keys():
            has_matching_rules = False
            for entry in feature_specific_colors[selected_feature]:
                specific_caption = entry[0]
                if not str(specific_caption or "").strip():
                    continue
                specific_fill_color = entry[1]
                # Only add to legend if this rule was actually used (or if used_color_rules not provided)
                if used_color_rules is None or (specific_caption, specific_fill_color) in used_color_rules:
                    has_matching_rules = True
                    legend_table.add(
                        str(specific_caption),
                        {
                            "type": "solid",
                            "fill": specific_fill_color,
                            "stroke": block_stroke_color,
                            "width": block_stroke_width,
                        },
                        key=rule_keys[str(specific_caption), _legend_fill_identity(specific_fill_color)],
                    )
            # Only add "other X" entry if at least one specific rule was used
            if has_matching_rules:
                allow_other_entry = False
                if used_color_rules is None:
                    allow_other_entry = True
                elif selected_feature in default_used_features:
                    allow_other_entry = True
                new_selected_key_name = _other_feature_legend_key(selected_feature)
                if allow_other_entry:
                    feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
                    legend_table.add(new_selected_key_name, {
                        "type": "solid",
                        "fill": feature_fill_color,
                        "stroke": block_stroke_color,
                        "width": block_stroke_width,
                    })
            elif used_color_rules is None:
                # Legacy behavior: add "other X" even if no rules matched
                new_selected_key_name = _other_feature_legend_key(selected_feature)
                feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
                legend_table.add(new_selected_key_name, {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                })
            else:
                # used_color_rules is provided but no specific rules matched
                # Just add the feature type with default color
                feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
                legend_table.add(selected_feature, {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                })
        else:
            feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
            legend_table.add(selected_feature, {
                "type": "solid",
                "fill": feature_fill_color,
                "stroke": block_stroke_color,
                "width": block_stroke_width,
            })
    if show_depth:
        legend_table.add("Depth", {
            "type": "solid",
            "fill": depth_config.fill_color,
            "stroke": depth_config.stroke_color,
            "width": depth_config.stroke_width,
        })
    if show_gc:
        if gc_high_fill_color == gc_low_fill_color:
            legend_table.add(f"{dinucleotide} content", {
                "type": "solid",
                "fill": gc_high_fill_color,
                "stroke": gc_stroke_color,
                "width": gc_stroke_width,
            })
        else:
            legend_table.add(f"{dinucleotide} content (+)", {
                "type": "solid",
                "fill": gc_high_fill_color,
                "stroke": gc_stroke_color,
                "width": gc_stroke_width,
            })
            legend_table.add(f"{dinucleotide} content (-)", {
                "type": "solid",
                "fill": gc_low_fill_color,
                "stroke": gc_stroke_color,
                "width": gc_stroke_width,
            })
    if show_skew:
        if skew_high_fill_color == skew_low_fill_color:
            legend_table.add(f"{dinucleotide} skew", {
                "type": "solid",
                "fill": skew_high_fill_color,
                "stroke": skew_stroke_color,
                "width": skew_stroke_width,
            })
        else:
            legend_table.add(f"{dinucleotide} skew (+)", {
                "type": "solid",
                "fill": skew_high_fill_color,
                "stroke": skew_stroke_color,
                "width": skew_stroke_width,
            })
            legend_table.add(f"{dinucleotide} skew (-)", {
                "type": "solid",
                "fill": skew_low_fill_color,
                "stroke": skew_stroke_color,
                "width": skew_stroke_width,
            })
    if has_blast and blast_config and not bool(getattr(blast_config, "hide_pairwise_identity_legend", False)):
        identity_legend_entries = getattr(blast_config, "pairwise_identity_legend_entries", None)
        if identity_legend_entries:
            for entry in identity_legend_entries:
                label = str(entry["label"])
                min_color = entry["min_color"]
                if float(blast_config.identity) >= 100:
                    min_color = entry["max_color"]
                legend_table.add(label, {
                    "type": "gradient",
                    "min_color": min_color,
                    "max_color": entry["max_color"],
                    "stroke": "none",
                    "width": 0,
                    "min_value": blast_config.identity,
                })
        else:
            identity_legend_label = str(
                getattr(blast_config, "pairwise_identity_legend_label", "Pairwise match identity")
            )
            min_color = blast_config.min_color
            if float(blast_config.identity) >= 100:
                min_color = blast_config.max_color
            legend_table.add(identity_legend_label, {
                "type": "gradient",
                "min_color": min_color,
                "max_color": blast_config.max_color,
                "stroke": "none",
                "width": 0,
                "min_value": blast_config.identity,
            })
    return legend_table


__all__ = [
    "LegendRowNames",
    "LegendTable",
    "collect_comparison_collinearity_color_modes",
    "configure_pairwise_identity_legend_from_comparisons",
    "drawn_legend_names",
    "legend_row_names",
    "legend_swatch",
    "prepare_legend_table",
    "rule_legend_keys",
    "write_legend_row_names",
]

