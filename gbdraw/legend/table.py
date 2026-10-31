#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Callable, Container, Hashable, Iterator, Mapping, MutableMapping
from typing import Iterable, List, NamedTuple, NewType, Optional, Sequence, Set, Tuple

from pandas import DataFrame

from gbdraw.core.color import (
    DEFAULT_COLLINEAR_ORIENTATION_COLORS,
    DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
    normalize_hex_color,
)
from gbdraw.exceptions import ValidationError
from gbdraw.legend.sources import (
    ComparisonLegendSource,
    DepthLegendSource,
    FeatureLegendSource,
    NumericLegendSource,
)


LegendEntry = Mapping[str, object]

# A legend row's identity (``data-legend-key``), such as ``name [fill]`` or ``label (2)``; never drawn.
LegendRowKey = NewType("LegendRowKey", str)
# A legend row's name as the user or the renderer wrote it; never a key, never rewritten.
LegendRowName = NewType("LegendRowName", str)
# The name the legend draws for a row: its name, numbered when rows of other swatches share it.
DrawnLegendName = NewType("DrawnLegendName", str)


def _key_of_name(name: str) -> LegendRowKey:
    """The key of a row whose key is its name: a renderer row, and a row no other row contests."""

    return LegendRowKey(name)


def _unique_legend_key(reserved: Iterable[str], preferred: str) -> LegendRowKey:
    """A row key no reserved key equals: ``preferred``, else ``preferred (2)``, ``(3)``, ..."""

    taken = set(reserved)
    if preferred not in taken:
        return LegendRowKey(preferred)
    suffix = 2
    while f"{preferred} ({suffix})" in taken:
        suffix += 1
    return LegendRowKey(f"{preferred} ({suffix})")


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
    paint). A gradient row is its two end colors. A similarity-ring row is its
    ring (``track_index``), so every drawn ring keeps a row of its own.
    """

    if entry.get("type") == "ring":
        return ("ring", int(str(entry.get("track_index"))))
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


def drawn_legend_names(rows: Sequence[tuple[LegendRowName, Hashable]]) -> list[DrawnLegendName]:
    """THE naming rule: the name each legend row is drawn with.

    ``rows`` are the rows being drawn, in legend order, each as its name and
    its swatch. Rows of one name and one swatch are one row and share a drawn
    name. A name left with one row is drawn as written. A name left with two
    or more rows is drawn ``name (1)``, ``name (2)``, ... in legend order; a
    number skips a drawn name that another row already has.
    """

    swatches: dict[LegendRowName, list[Hashable]] = {}
    for name, swatch in rows:
        group = swatches.setdefault(name, [])
        if swatch not in group:
            group.append(swatch)
    taken: set[str] = {name for name, group in swatches.items() if len(group) == 1}
    numbered: dict[tuple[LegendRowName, Hashable], DrawnLegendName] = {}
    for name, group in swatches.items():
        if len(group) == 1:
            continue
        number = 0
        for swatch in group:
            number += 1
            while f"{name} ({number})" in taken:
                number += 1
            numbered[name, swatch] = DrawnLegendName(f"{name} ({number})")
            taken.add(numbered[name, swatch])
    return [numbered.get((name, swatch), DrawnLegendName(name)) for name, swatch in rows]


class LegendRowNames(NamedTuple):
    """A legend row's name as written, and the name the legend draws."""

    name: LegendRowName
    drawn: DrawnLegendName


class _LegendRow(NamedTuple):
    name: LegendRowName
    swatch: tuple[object, ...]
    entry: LegendEntry


class LegendTable(Mapping[LegendRowKey, LegendEntry]):
    """The legend rows of one drawing, in legend order, by row key.

    A row has three separate facts. Its *key* is its identity (the Web stores
    row edits and Sessions under it, ``data-legend-key``); it is never drawn.
    Its *name* is what the user or the renderer wrote and is never rewritten.
    Its *drawn name* is derived from the names and swatches of the rows being
    drawn (``drawn_legend_names``). Read as a mapping, the table gives each
    key's entry (the swatch properties the legend draws).
    """

    def __init__(self, rows: Mapping[object, LegendEntry] | None = None) -> None:
        self._rows: dict[LegendRowKey, _LegendRow] = {}
        self._sources: dict[LegendRowKey, int] = {}
        self._shared_keys: list[LegendRowKey] = []
        self._default_tracks: dict[LegendRowKey, LegendRowKey] = {}
        for key, entry in (rows or {}).items():
            self.add(str(key), entry)

    def add(self, name: str, entry: LegendEntry, *, key: LegendRowKey | None = None) -> None:
        """Add a row under its key (its name unless given); one of the same name and swatch is shared."""

        self._add(name, entry, _key_of_name(name) if key is None else key)

    def add_source(self, name: str, entry: LegendEntry) -> None:
        """Add the row of a track or annotation source, which has no key of its own.

        Its key is its name, numbered ``(2)``, ``(3)``, ... after the keys taken
        before it. A source that shares a row still takes its key, so the keys
        of the other sources do not depend on which rows are shared.
        """

        self._add(name, entry, _unique_legend_key([*self._rows, *self._shared_keys], name))

    def add_numeric_source(self, name: str, swatch: NumericLegendSource) -> None:
        """Add a content or skew slot's rows: one when its fills are equal, else ``(+)`` and ``(-)``."""

        stroke = {"stroke": swatch.stroke_color, "width": swatch.stroke_width}
        if swatch.high_fill_color == swatch.low_fill_color:
            self.add_source(name, {"type": "solid", "fill": swatch.high_fill_color, **stroke})
        else:
            self.add_source(f"{name} (+)", {"type": "solid", "fill": swatch.high_fill_color, **stroke})
            self.add_source(f"{name} (-)", {"type": "solid", "fill": swatch.low_fill_color, **stroke})

    def add_default_track(self, name: str, entry: LegendEntry) -> None:
        """Add the renderer's own row of a numeric track, which explicit slots replace (``discard_default_track``)."""

        key = _key_of_name(name)
        self._default_tracks[key] = self._add(name, entry, key)

    def discard_default_track(self, name: str) -> None:
        """Take a default track's row back; a row another source shares stays for that source (OV-588)."""

        key = _key_of_name(name)
        row_key = self._default_tracks.pop(key, None)
        if row_key is None:
            return
        self._sources[row_key] -= 1
        if not self._sources[row_key]:
            del self._rows[row_key], self._sources[row_key]
        if key in self._shared_keys:
            self._shared_keys.remove(key)

    def _add(self, name: str, entry: LegendEntry, key: LegendRowKey) -> LegendRowKey:
        """Add a source's row and return the key of the row that draws it."""

        swatch = legend_swatch(entry)
        for row_key, row in self._rows.items():
            if row.name == name and row.swatch == swatch:
                self._sources[row_key] += 1
                if key not in self._rows:
                    self._shared_keys.append(key)
                return row_key
        self._rows[key] = _LegendRow(LegendRowName(name), swatch, entry)
        self._sources[key] = 1
        return key

    @property
    def shared_keys(self) -> tuple[LegendRowKey, ...]:
        """Keys of the sources whose row another source of the same name and swatch draws."""

        return tuple(key for key in dict.fromkeys(self._shared_keys) if key not in self._rows)

    def names(self) -> dict[LegendRowKey, LegendRowNames]:
        """Each row's name and drawn name, by key.

        The names are numbered in the order the legends draw the rows: the
        solid rows, then the gradient and ring rows (the block beside or below
        them), each in table order. The Web app numbers the rows of a Result in the
        order it shows them, so both give one answer.
        """

        solid_first = sorted(self._rows.items(), key=lambda item: item[1].swatch[0] != "solid")
        drawn = drawn_legend_names([(row.name, row.swatch) for _, row in solid_first])
        texts = {key: text for (key, _), text in zip(solid_first, drawn)}
        return {key: LegendRowNames(row.name, texts[key]) for key, row in self._rows.items()}

    def __getitem__(self, key: LegendRowKey) -> LegendEntry:
        return self._rows[key].entry

    def __iter__(self) -> Iterator[LegendRowKey]:
        return iter(self._rows)

    def __len__(self) -> int:
        return len(self._rows)


def as_legend_table(legend_table: LegendTable | Mapping[object, LegendEntry]) -> LegendTable:
    """The table itself; a plain mapping (a test's table) is read as rows named by their keys."""

    return legend_table if isinstance(legend_table, LegendTable) else LegendTable(legend_table)


def legend_row_names(legend_table: LegendTable | Mapping[object, LegendEntry]) -> dict[LegendRowKey, LegendRowNames]:
    """The names the legend layouts measure and the legend drawers write, by row key."""

    return as_legend_table(legend_table).names()


def write_legend_row_names(
    legend_table: LegendTable | None,
) -> Callable[[MutableMapping[str, object], LegendRowKey], DrawnLegendName]:
    """The legend drawers' one writer of a row's key, name, and drawn name.

    The returned function writes ``data-legend-key`` on a row group's
    attributes, and ``data-legend-name`` where the row's name is not its key (a
    reader takes the key as the name of a row without it, so a drawing whose
    names are its keys keeps its bytes), and returns the name the legend draws.
    """

    names = legend_row_names(legend_table if legend_table is not None else LegendTable())

    def write(attribs: MutableMapping[str, object], key: LegendRowKey) -> DrawnLegendName:
        attribs["data-legend-key"] = key
        if names[key].name != key:
            attribs["data-legend-name"] = names[key].name
        return names[key].drawn

    return write


def _rule_fill_key(caption: str, fill: str) -> str:
    """The key of a rule row that its fill tells apart: ``caption [fill]``."""

    return f"{caption} [{fill}]"


def saved_rule_caption(caption: str, color: object) -> str:
    """The caption of a color rule as a saved color table or Session is read.

    Builds from first-parent ``main`` fe6861f0 until the naming rule
    (``drawn_legend_names``) saved the key ``caption [fill]`` as the caption of
    a rule whose caption was written with several colors, in the Session and
    in the downloaded color table. A caption that ends in the fill identity of
    its own rule is read without it, once; one that would be left blank is read
    as written (PD-OI-042 scenario revision 2, accepted residual risk).
    """

    name = caption.removesuffix(_rule_fill_key("", _legend_fill_identity(color)))
    return name if name.strip() else caption


def saved_rule_captions(color_table: DataFrame) -> DataFrame:
    """A saved color table with each caption read by ``saved_rule_caption``.

    Returns ``color_table`` itself when no caption changes.
    """

    if "caption" not in color_table:
        return color_table
    captions = [
        saved_rule_caption(caption, color)
        for caption, color in zip(color_table["caption"], color_table["color"])
    ]
    if captions == color_table["caption"].tolist():
        return color_table
    return color_table.assign(caption=captions)


def rule_legend_keys(
    rules: Iterable[tuple[object, object]],
    renderer_swatches: Mapping[str, tuple[object, ...]],
    rule_swatch: Callable[[str], tuple[object, ...]],
) -> dict[tuple[str, str], LegendRowKey]:
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
    keys: dict[tuple[str, str], LegendRowKey] = {}
    for caption in sorted(fills):
        renderer_swatch = renderer_swatches.get(caption)
        for fill in sorted(fills[caption]):
            if len(fills[caption]) > 1 or renderer_swatch not in (None, rule_swatch(fill)):
                key = _unique_legend_key(reserved, _rule_fill_key(caption, fill))
                reserved.add(key)
            else:
                key = _key_of_name(caption)
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


class PairwiseIdentityGradient(NamedTuple):
    """One identity gradient row of the pairwise match legend."""

    label: str
    min_color: str
    max_color: str


class PairwiseIdentityLegend(NamedTuple):
    """How the pairwise match identity legend is drawn.

    ``gradients`` replaces the single ``label`` row when collinear and inverted
    matches have their own colors; ``hidden`` drops the rows.
    """

    hidden: bool = False
    label: str = "Pairwise match identity"
    gradients: tuple[PairwiseIdentityGradient, ...] = ()


def pairwise_identity_legend_from_comparisons(
    blast_config: object,
    comparison_dataframes: Sequence[DataFrame] | None,
    *,
    additional_color_modes: Iterable[object] | None = None,
) -> PairwiseIdentityLegend | None:
    """The pairwise identity legend that collinearity color modes ask for.

    ``None`` when the comparisons carry no recognized color mode, so the
    caller keeps the legend it had.
    """

    color_modes = collect_comparison_collinearity_color_modes(comparison_dataframes)
    for value in additional_color_modes or ():
        normalized = _normalize_collinearity_legend_color_mode(value)
        if normalized is not None:
            color_modes.add(normalized)
    if not color_modes:
        return None
    hidden = color_modes == {"orientation"}
    if "orientation_identity" in color_modes:
        min_colors = getattr(
            blast_config,
            "collinearity_orientation_min_colors",
            DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
        )
        max_colors = getattr(
            blast_config,
            "collinearity_orientation_colors",
            DEFAULT_COLLINEAR_ORIENTATION_COLORS,
        )
        return PairwiseIdentityLegend(
            hidden=hidden,
            gradients=tuple(
                PairwiseIdentityGradient(
                    label,
                    min_colors.get(strand, DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS[strand]),
                    max_colors.get(strand, DEFAULT_COLLINEAR_ORIENTATION_COLORS[strand]),
                )
                for label, strand in (("Collinear", "plus"), ("Inverted", "minus"))
            ),
        )
    if "average_identity" in color_modes:
        return PairwiseIdentityLegend(hidden=hidden, label="Average identity")
    return PairwiseIdentityLegend(hidden=hidden)


def _other_feature_legend_key(feature_type: str) -> LegendRowKey:
    return _key_of_name("other proteins" if feature_type == "CDS" else f"other {feature_type}s")


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


def feature_type_rows(
    feature_specific_colors: Mapping[str, Sequence[tuple[str, str]]],
    features_present: Sequence[str],
    used_color_rules: Optional[Set[Tuple[str, str]]],
    default_used_features: Set[str],
) -> list[tuple[str, LegendRowKey]]:
    """The renderer's own row of each drawn feature type that has one, by type:
    the type itself without rules of its type, else "other ..." while some of
    its features keep the default color (always without the used rules), and
    the type again when none of its captioned rules is used."""

    keys: list[tuple[str, LegendRowKey]] = []
    for feature_type in features_present:
        captioned = [
            (caption, color)
            for caption, color in feature_specific_colors.get(feature_type, ())
            if str(caption or "").strip()
        ]
        if feature_type not in feature_specific_colors:
            key: LegendRowKey | None = _key_of_name(feature_type)
        elif used_color_rules is None:
            key = _other_feature_legend_key(feature_type)
        elif any(rule in used_color_rules for rule in captioned):
            key = _other_feature_legend_key(feature_type) if feature_type in default_used_features else None
        else:
            key = _key_of_name(feature_type)
        if key is not None:
            keys.append((feature_type, key))
    return keys


def feature_legend_row_order(
    feature_specific_colors: Mapping[str, Sequence[tuple[str, str]]],
    features_present: Sequence[str],
    used_color_rules: Optional[Set[Tuple[str, str]]],
    renderer_rows: Container[str],
) -> Iterator[tuple[str, str, str] | LegendRowKey]:
    """THE order of the feature rows of a legend: for each drawn feature type,
    in the order the types first appear in the drawing (``features_present``),
    the captioned rules of that type in table order, each as its type, caption
    and color (a rule the features did not take has no row, unless the used
    rules are not given), then the renderer's own row of the type (its type,
    or "other ...") when ``renderer_rows`` holds it. Each renderer row comes
    once. The legend adds the rows in this order; rows of one name and one
    swatch are one row, at the first (``LegendTable``)."""

    placed: set[str] = set()
    for feature_type in features_present:
        for caption, color in feature_specific_colors.get(feature_type, ()):
            if str(caption or "").strip() and (used_color_rules is None or (caption, color) in used_color_rules):
                yield feature_type, str(caption), color
        for key in (_key_of_name(feature_type), _other_feature_legend_key(feature_type)):
            if key in renderer_rows and key not in placed:
                placed.add(key)
                yield key


def _generated_legend_rows(
    feature_specific_colors: dict,
    features_present: Sequence[str],
    default_colors: DataFrame,
    *,
    used_color_rules: Optional[Set[Tuple[str, str]]],
    default_used_features: Set[str],
    gc: NumericLegendSource,
    skew: NumericLegendSource,
    depth: DepthLegendSource | None,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
    block_stroke: tuple[object, object] = ("none", 0),
) -> dict[str, dict[str, object]]:
    """The renderer's own legend rows by key, as drawn: feature types, "other ...", and numeric tracks."""

    rows: dict[str, dict[str, object]] = {
        key: _solid_legend_entry(_resolve_legend_feature_color(default_colors, feature_type), *block_stroke)
        for feature_type, key in feature_type_rows(
            feature_specific_colors, features_present, used_color_rules, default_used_features
        )
    }
    dinucleotide = gc.dinucleotide
    if show_depth and depth is not None:
        rows["Depth"] = _solid_legend_entry(depth.fill_color, depth.stroke_color, depth.stroke_width)
    for show, source, name in ((show_gc, gc, "content"), (show_skew, skew, "skew")):
        if not show:
            continue
        stroke = (source.stroke_color, source.stroke_width)
        if source.high_fill_color == source.low_fill_color:
            rows[f"{dinucleotide} {name}"] = _solid_legend_entry(source.high_fill_color, *stroke)
        else:
            rows[f"{dinucleotide} {name} (+)"] = _solid_legend_entry(source.high_fill_color, *stroke)
            rows[f"{dinucleotide} {name} (-)"] = _solid_legend_entry(source.low_fill_color, *stroke)
    return rows


def prepare_legend_table(
    gc: NumericLegendSource,
    skew: NumericLegendSource,
    features: FeatureLegendSource,
    features_present: List[str],
    *,
    comparison: ComparisonLegendSource | None = None,
    used_color_rules: Optional[Set[Tuple[str, str]]] = None,
    default_used_features: Optional[Set[str]] = None,
    depth: DepthLegendSource | None = None,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
    pairwise_identity_legend: PairwiseIdentityLegend = PairwiseIdentityLegend(),
) -> LegendTable:
    """
    Prepare the legend table for the diagram.

    Args:
        comparison: The pairwise identity gradient; ``None`` without comparisons.
        used_color_rules: Optional set of (caption, color) tuples that were actually
            matched during feature creation. If provided, only these rules will be
            included in the legend. If None, all rules for present feature types
            will be included (legacy behavior).
        default_used_features: Optional set of feature types that fell back to
            default colors. Used to decide whether to include "other X" entries
            when specific rules exist.
    """
    legend_table = LegendTable()
    color_table: Optional[DataFrame] = features.color_table
    block_stroke = (features.block_stroke_color, features.block_stroke_width)
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
        features.default_colors,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features or set(),
        gc=gc,
        skew=skew,
        depth=depth,
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
    for row in feature_legend_row_order(feature_specific_colors, features_present, used_color_rules, generated):
        if isinstance(row, tuple):
            _, caption, fill = row
            legend_table.add(
                caption,
                _solid_legend_entry(fill, *block_stroke),
                key=rule_keys[caption, _legend_fill_identity(fill)],
            )
        else:
            legend_table.add(row, generated.pop(row))
    for key, entry in generated.items():
        legend_table.add_default_track(key, entry)
    if comparison is not None and not pairwise_identity_legend.hidden:
        gradients = pairwise_identity_legend.gradients or (
            PairwiseIdentityGradient(
                pairwise_identity_legend.label,
                comparison.min_color,
                comparison.max_color,
            ),
        )
        for gradient in gradients:
            # A source row: a rule captioned like it keeps its own row (OV-646).
            legend_table.add_source(gradient.label, {
                "type": "gradient",
                "min_color": (
                    gradient.max_color if comparison.domain.minimum >= comparison.domain.maximum else gradient.min_color
                ),
                "max_color": gradient.max_color,
                "stroke": "none",
                "width": 0,
                "min_value": comparison.domain.minimum,
                "max_value": comparison.domain.maximum,
            })
    return legend_table


__all__ = [
    "DrawnLegendName",
    "LegendRowKey",
    "LegendRowName",
    "LegendRowNames",
    "LegendTable",
    "PairwiseIdentityGradient",
    "PairwiseIdentityLegend",
    "as_legend_table",
    "collect_comparison_collinearity_color_modes",
    "drawn_legend_names",
    "feature_legend_row_order",
    "feature_type_rows",
    "legend_row_names",
    "legend_swatch",
    "pairwise_identity_legend_from_comparisons",
    "prepare_legend_table",
    "rule_legend_keys",
    "write_legend_row_names",
]

