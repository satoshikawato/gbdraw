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
        self._sources: dict[str, int] = {}
        self._shared_keys: list[str] = []
        self._default_tracks: dict[str, str] = {}
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

    def add_default_track(self, key: str, entry: LegendEntry) -> None:
        """Add the renderer's own row of a numeric track, which explicit slots replace (``discard_default_track``)."""

        self._default_tracks[key] = self._add(key, entry, key)

    def discard_default_track(self, key: str) -> None:
        """Take a default track's row back; a row another source shares stays for that source (OV-588)."""

        row_key = self._default_tracks.pop(key, None)
        if row_key is None:
            return
        self._sources[row_key] -= 1
        if not self._sources[row_key]:
            del self._rows[row_key], self._sources[row_key]
        if key in self._shared_keys:
            self._shared_keys.remove(key)

    def _add(self, name: str, entry: LegendEntry, key: str) -> str:
        """Add a source's row and return the key of the row that draws it."""

        swatch = legend_swatch(entry)
        for row_key, row in self._rows.items():
            if row.name == name and row.swatch == swatch:
                self._sources[row_key] += 1
                if key not in self._rows:
                    self._shared_keys.append(key)
                return row_key
        self._rows[key] = _LegendRow(name, swatch, entry)
        self._sources[key] = 1
        return key

    @property
    def shared_keys(self) -> tuple[str, ...]:
        """Keys of the sources whose row another source of the same name and swatch draws."""

        return tuple(key for key in dict.fromkeys(self._shared_keys) if key not in self._rows)

    def names(self) -> dict[str, LegendRowNames]:
        """Each row's name and drawn name, by key.

        The names are numbered in the order the legends draw the rows: the
        solid rows, then the gradient rows (the block beside or below them),
        each in table order. The Web app numbers the rows of a Result in the
        order it shows them, so both give one answer.
        """

        solid_first = sorted(self._rows.items(), key=lambda item: item[1].swatch[0] == "gradient")
        drawn = drawn_legend_names([(row.name, row.swatch) for _, row in solid_first])
        texts = {key: text for (key, _), text in zip(solid_first, drawn)}
        return {key: LegendRowNames(row.name, texts[key]) for key, row in self._rows.items()}

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

    The returned function writes ``data-legend-key`` on a row group's
    attributes, and ``data-legend-name`` where the row's name is not its key (a
    reader takes the key as the name of a row without it, so a drawing whose
    names are its keys keeps its bytes), and returns the name the legend draws.
    """

    names = legend_row_names(legend_table or {})

    def write(attribs: MutableMapping[str, object], key: str) -> str:
        attribs["data-legend-key"] = key
        if names[key].name != key:
            attribs["data-legend-name"] = names[key].name
        return names[key].drawn

    return write


def rule_legend_keys(
    rules: Iterable[tuple[object, object]],
    renderer_swatches: Mapping[str, tuple[object, ...]],
    rule_swatch: Callable[[str], tuple[object, ...]],
) -> dict[tuple[str, str], str]:
    """Row keys for the captions of the color table, by caption and fill identity.

    A rule row's key is its caption. A caption written with several colors in
    the table, or naming a renderer row of another swatch (``renderer_swatches``
    against ``rule_swatch``, the swatch a rule row of a fill is drawn with),
    takes ``caption [fill]`` instead, made unique against every caption and
    renderer row. The allocation reads the whole table, never its row order.
    The key is an identity only: the legend draws the caption
    (``drawn_legend_names``).
    """

    fills: dict[str, set[str]] = {}
    for caption, color in rules:
        if str(caption or "").strip():
            fills.setdefault(str(caption), set()).add(_legend_fill_identity(color))
    reserved = {*renderer_swatches, *fills}
    keys: dict[tuple[str, str], str] = {}
    for caption in sorted(fills):
        renderer_swatch = renderer_swatches.get(caption)
        for fill in sorted(fills[caption]):
            if len(fills[caption]) > 1 or renderer_swatch not in (None, rule_swatch(fill)):
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


def _solid_legend_entry(fill: object, stroke: object, width: object) -> dict[str, object]:
    return {"type": "solid", "fill": fill, "stroke": stroke, "width": width}


def _generated_legend_rows(
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
    block_stroke: tuple[object, object] = ("none", 0),
) -> dict[str, dict[str, object]]:
    """The renderer's own legend rows by key, as drawn: feature types, "other ...", and numeric tracks."""

    rows: dict[str, dict[str, object]] = {}
    for feature_type in features_present:
        captioned = [
            (caption, color)
            for caption, color in feature_specific_colors.get(feature_type, ())
            if str(caption or "").strip()
        ]
        if feature_type not in feature_specific_colors:
            key: str | None = feature_type
        elif used_color_rules is None:
            key = _other_feature_legend_key(feature_type)
        elif any(rule in used_color_rules for rule in captioned):
            key = _other_feature_legend_key(feature_type) if feature_type in default_used_features else None
        else:
            key = feature_type
        if key is not None:
            rows[key] = _solid_legend_entry(
                _resolve_legend_feature_color(default_colors, feature_type), *block_stroke
            )
    dinucleotide = gc_config.dinucleotide
    if show_depth and depth_config is not None:
        rows["Depth"] = _solid_legend_entry(
            depth_config.fill_color,
            getattr(depth_config, "stroke_color", "none"),
            getattr(depth_config, "stroke_width", 0),
        )
    for show, config, name in ((show_gc, gc_config, "content"), (show_skew, skew_config, "skew")):
        if not show:
            continue
        stroke = (getattr(config, "stroke_color", "none"), getattr(config, "stroke_width", 0))
        if config.high_fill_color == config.low_fill_color:
            rows[f"{dinucleotide} {name}"] = _solid_legend_entry(config.high_fill_color, *stroke)
        else:
            rows[f"{dinucleotide} {name} (+)"] = _solid_legend_entry(config.high_fill_color, *stroke)
            rows[f"{dinucleotide} {name} (-)"] = _solid_legend_entry(config.low_fill_color, *stroke)
    return rows


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
    block_stroke = (feature_config.block_stroke_color, feature_config.block_stroke_width)
    feature_specific_colors: dict[str, list[tuple[str, str]]] = {}
    if color_table is not None and not color_table.empty:
        for _, row in color_table.iterrows():
            feature_specific_colors.setdefault(row["feature_type"], []).append((row["caption"], row["color"]))
    # The renderer's own rows (feature types, "other ..." and numeric tracks)
    # keep their keys; a rule row takes its key from the whole color table and
    # these rows (``rule_legend_keys``). The legend draws names, never keys.
    generated = _generated_legend_rows(
        feature_specific_colors,
        features_present,
        feature_config.default_colors,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features or set(),
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=show_depth,
        block_stroke=block_stroke,
    )
    rule_keys = rule_legend_keys(
        (entry for entries in feature_specific_colors.values() for entry in entries),
        {key: legend_swatch(entry) for key, entry in generated.items()},
        lambda fill: legend_swatch(_solid_legend_entry(fill, *block_stroke)),
    )
    for selected_feature in features_present:
        for specific_caption, specific_fill_color in feature_specific_colors.get(selected_feature, ()):
            if not str(specific_caption or "").strip():
                continue
            # Only a rule a feature took has a row (every rule when the used rules are not given).
            if used_color_rules is None or (specific_caption, specific_fill_color) in used_color_rules:
                legend_table.add(
                    str(specific_caption),
                    _solid_legend_entry(specific_fill_color, *block_stroke),
                    key=rule_keys[str(specific_caption), _legend_fill_identity(specific_fill_color)],
                )
        for key in (selected_feature, _other_feature_legend_key(selected_feature)):
            if key in generated:
                legend_table.add(key, generated.pop(key))
    for key, entry in generated.items():
        legend_table.add_default_track(key, entry)
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

