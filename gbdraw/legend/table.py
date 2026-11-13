#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import json
from collections.abc import Callable, Hashable, Iterator, Mapping, MutableMapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, List, NamedTuple, NewType, Optional, Sequence, Set, Tuple, TypeVar

from pandas import DataFrame

from gbdraw.core.color import (
    DEFAULT_COLLINEAR_ORIENTATION_COLORS,
    DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
    normalize_hex_color,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import RuleUse, build_default_color_map, feature_type_color, legend_rule_use
from gbdraw.svg.result_records import LEGEND_MERGED_ROW, result_base_attribute
from gbdraw.legend.sources import (
    ComparisonLegendSource,
    DepthLegendSource,
    FeatureDensityLegendSource,
    FeatureLegendSource,
    NumericLegendSource,
)
from gbdraw.tracks.parsing import SlotView, feature_density_slot_track_id, slot_legend_label

if TYPE_CHECKING:
    from gbdraw.config.models.objects import ObjectsLegendsConfig


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


_Placed = TypeVar("_Placed")


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
        # The request's Legend row edits (``apply_legend_row_edits``): the rows drawn hidden, each with the key of
        # the row it merges into (None: a hidden row), the names the edits give, the text each row had before the
        # edits, and Python's order of the feature rows when the edits reorder them.
        self._hidden: dict[LegendRowKey, LegendRowKey | None] = {}
        self._renamed: dict[LegendRowKey, LegendRowName] = {}
        self._unedited_texts: dict[LegendRowKey, DrawnLegendName] = {}
        self._python_order: tuple[LegendRowKey, ...] | None = None
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

    def _name(self, key: LegendRowKey) -> LegendRowName:
        """A row's name: the one an edit gives it, else the one written."""

        return self._renamed.get(key, self._rows[key].name)

    def names(self) -> dict[LegendRowKey, LegendRowNames]:
        """Each shown row's name and drawn name, by key.

        The names are numbered in the order the legends draw the rows: the
        solid rows, then the gradient and ring rows (the block beside or below
        them), each in table order. The Web app numbers the rows of a Result in the
        order it shows them, so both give one answer.
        """

        shown = [(key, row) for key, row in self._rows.items() if key not in self._hidden]
        solid_first = sorted(shown, key=lambda item: item[1].swatch[0] != "solid")
        drawn = drawn_legend_names([(self._name(key), row.swatch) for key, row in solid_first])
        texts = {key: text for (key, _), text in zip(solid_first, drawn)}
        return {key: LegendRowNames(self._name(key), texts[key]) for key, _ in shown}

    def ring_names(self) -> dict[int, DrawnLegendName]:
        """Each similarity ring's drawn name (``names``) by its ring (``track_index``), for its ring group."""

        names = self.names()
        return {
            int(str(row.entry["track_index"])): names[key].drawn
            for key, row in self._rows.items()
            if row.swatch[0] == "ring" and key in names
        }

    def drawn_keys(self) -> tuple[LegendRowKey, ...]:
        """The keys of every row the legend draws, the hidden ones included, in the order it draws them."""

        return tuple(self._rows)

    def hidden_entry(self, key: LegendRowKey) -> LegendEntry | None:
        """The entry of a row the legend draws hidden (a deleted or merged row), else None."""

        return self._rows[key].entry if key in self._hidden else None

    def feature_rows(self, placed: Mapping[LegendRowKey, _Placed]) -> list[tuple[LegendRowKey, LegendEntry, _Placed | None]]:
        """The feature rows a legend drawer writes, in the order it writes them: each row of ``placed`` (the shown
        rows its layout placed) with its placement, and each hidden row, which no layout places, at the place of the
        next shown row (else of the last one; None when no row is shown). The Web's executor orders a group's rows
        by where they are, so a hidden row there keeps a second pass of the order a no-op."""

        keys = [key for key in self._rows if key in placed or key in self._hidden]
        rows: list[tuple[LegendRowKey, LegendEntry, _Placed | None]] = []
        for index, key in enumerate(keys):
            place = placed.get(key)
            if place is None:
                later = (placed[other] for other in keys[index + 1:] if other in placed)
                earlier = (placed[other] for other in reversed(keys[:index]) if other in placed)
                place = next(later, None)
                if place is None:
                    place = next(earlier, None)
            rows.append((key, self._rows[key].entry, place))
        return rows

    def order_record(self) -> str | None:
        """``LEGEND_ORDER_RECORD`` of a feature Legend group: Python's order of its rows when the edits reorder
        them, as the Web's executor writes it (JSON of the keys), else None."""

        if self._python_order is None:
            return None
        return json.dumps(list(self._python_order), ensure_ascii=False, separators=(",", ":"))

    def _row_attributes(self, key: LegendRowKey) -> dict[str, str]:
        """The attributes a row's group carries: its key, its name where it is not the key, and the records of an
        edit as the Web's executor writes them (a rename keeps Python's key; a hidden row keeps Python's shown
        ``display``; a merged row names the row it merges into)."""

        row = self._rows[key]
        attributes = {"data-legend-key": str(self._renamed.get(key, key))}
        if key in self._renamed:
            attributes[result_base_attribute("data-legend-key")] = key
        if row.name != key:
            attributes["data-legend-name"] = row.name
        if key in self._hidden:
            merged_into = self._hidden[key]
            attributes["display"] = "none"
            if merged_into is None:
                attributes[result_base_attribute("display")] = ""
            else:
                attributes[LEGEND_MERGED_ROW] = merged_into
        return attributes

    def __getitem__(self, key: LegendRowKey) -> LegendEntry:
        if key in self._hidden:
            raise KeyError(key)
        return self._rows[key].entry

    def __iter__(self) -> Iterator[LegendRowKey]:
        return (key for key in self._rows if key not in self._hidden)

    def __len__(self) -> int:
        return len(self._rows) - len(self._hidden)


def as_legend_table(legend_table: LegendTable | Mapping[object, LegendEntry]) -> LegendTable:
    """The table itself; a plain mapping (a test's table) is read as rows named by their keys."""

    return legend_table if isinstance(legend_table, LegendTable) else LegendTable(legend_table)


def draws_legend(legend_table: LegendTable | Mapping[object, LegendEntry]) -> bool:
    """Whether a Legend is drawn for the table: it has a row, shown or hidden.

    A Legend whose rows a request's Legend row edits all hide is drawn with those rows and their records, as an
    item that takes no room in the composition (D-36): the Result keeps the rows, so the Web can restore them.
    """

    return bool(as_legend_table(legend_table).drawn_keys())


def legend_row_names(legend_table: LegendTable | Mapping[object, LegendEntry]) -> dict[LegendRowKey, LegendRowNames]:
    """The names the legend layouts measure and the legend drawers write, by row key."""

    return as_legend_table(legend_table).names()


def write_legend_row_names(
    legend_table: LegendTable | None,
) -> Callable[[MutableMapping[str, object], LegendRowKey], DrawnLegendName]:
    """The legend drawers' one writer of a row's key, name, records, and drawn name.

    The returned function writes ``data-legend-key`` on a row group's
    attributes, and ``data-legend-name`` where the row's name is not its key (a
    reader takes the key as the name of a row without it, so a drawing whose
    names are its keys keeps its bytes), and the records of the request's row
    edits (``LegendTable._row_attributes``), and returns the text the legend
    draws: the drawn name of a shown row, and the text a hidden row had before
    the edits (the Web's executor leaves a hidden row's text as it was).
    """

    table = as_legend_table(legend_table if legend_table is not None else LegendTable())
    names = table.names()

    def write(attribs: MutableMapping[str, object], key: LegendRowKey) -> DrawnLegendName:
        attribs.update(table._row_attributes(key))
        return names[key].drawn if key in names else table._unedited_texts[key]

    return write


@dataclass(frozen=True)
class LegendRowEdits:
    """A request's Legend row edits, each by row key (``data-legend-key``, ``LegendRowKey``): the name a row is
    given (``objects.legends.row_names``), the rows hidden (``hidden_rows``), and the rows' order (``row_order``;
    empty: Python's). A render stage keys on it."""

    names: tuple[tuple[LegendRowKey, LegendRowName], ...] = ()
    hidden: tuple[LegendRowKey, ...] = ()
    order: tuple[LegendRowKey, ...] = ()

    @classmethod
    def of_config(cls, legends: ObjectsLegendsConfig) -> LegendRowEdits:
        return cls(
            names=tuple((LegendRowKey(key), LegendRowName(name)) for key, name in legends.row_names),
            hidden=tuple(LegendRowKey(key) for key in legends.hidden_rows),
            order=tuple(LegendRowKey(key) for key in legends.row_order),
        )

    def __bool__(self) -> bool:
        return bool(self.names or self.hidden or self.order)


def apply_legend_row_edits(legend_table: LegendTable, edits: LegendRowEdits) -> LegendTable:
    """THE one place Python applies a request's Legend row edits; returns the table the legend measures and draws.

    The edits reach the feature rows (the solid rows; the gradient and ring rows are drawn in their own block, which
    no edit reaches) and are applied as the Web's Result executor applies them to Python's drawing of the request
    without them (``applyLegendOperations`` and ``nameLegendRows`` in ``services/svg-result-ingestion.js``; vector
    ``tests/fixtures/legend_row_edit_vectors.json``):

    1. A rename gives the row of a key its name; the row shows that name as its key and keeps Python's key as a
       record. The drawn names then follow the naming rule (``drawn_legend_names``).
    2. A hidden row stays in the drawing, hidden, with a record that Python drew it shown; it is not measured.
    3. The order puts the rows of the listed keys first, in the listed order, then the other rows in their order;
       the feature Legend group records Python's order when the rows move.
    4. Then the shown rows of one name and one swatch are one row (D-28): a later one is hidden and names the
       earlier one's key. A hidden row keeps the text it had before the edits.

    A key the table lacks is ignored.
    """

    if not edits:
        return legend_table
    solid = [key for key, row in legend_table._rows.items() if row.swatch[0] == "solid"]
    edited = LegendTable()
    edited._sources = dict(legend_table._sources)
    edited._shared_keys = list(legend_table._shared_keys)
    edited._default_tracks = dict(legend_table._default_tracks)
    # The executor ranks a row by its key among the listed keys, trimmed; an unlisted row comes after the listed
    # ones, in its place.
    rank: dict[str, int] = {}
    for key in edits.order:
        rank.setdefault(key.strip(), len(rank))
    ordered = sorted(solid, key=lambda key: rank.get(key, len(rank)))
    rows = dict(legend_table._rows)
    if ordered != solid:
        rows = {key: legend_table._rows[key] for key in ordered}
        rows.update({key: row for key, row in legend_table._rows.items() if key not in rows})
        edited._python_order = tuple(solid)
    edited._rows = rows
    edited._unedited_texts = {key: names.drawn for key, names in legend_table.names().items()}
    edited._renamed = {
        key: name for key, name in dict(edits.names).items() if key in solid and name != key
    }
    edited._hidden = {key: None for key in edits.hidden if key in solid}
    seen: dict[tuple[LegendRowName, tuple[object, ...]], LegendRowKey] = {}
    for key in ordered:
        if key in edited._hidden:
            continue
        identity = (edited._name(key), rows[key].swatch)
        if identity in seen:
            edited._hidden[key] = seen[identity]
        else:
            seen[identity] = key
    return edited


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


def _solid_legend_entry(fill: object, stroke: object, width: object) -> dict[str, object]:
    return {"type": "solid", "fill": fill, "stroke": stroke, "width": width}


# A color table row as the Legend reads it: the feature type it is for (``*``:
# every type), its caption, and its color.
LegendRule = tuple[str, str, str]


@dataclass(frozen=True)
class FeatureLegendRows:
    """THE owner of the feature rows of a Legend and of the row each drawn feature is in.

    ``rules`` are the color table rows, in table order. A rule for every type
    (``*``) is read as the same rule written for each type (OV-3460). ``present``
    are the drawn types in the order they first appear in the drawing; ``used``
    and ``default_types`` are what ``legend_rule_use`` reads of the drawn
    features (``used`` None: every captioned rule counts as used, for callers
    that do not read the features). The Web's twin is ``featureLegendRows``
    (gbdraw/web/js/services/specific-color-rules.js); both answer
    tests/fixtures/legend_generated_rows.json (``feature_rows``).
    """

    rules: tuple[LegendRule, ...]
    present: tuple[str, ...]
    used: frozenset[RuleUse] | None
    default_types: frozenset[str]

    @classmethod
    def of(cls, rules: Iterable[LegendRule], drawn: Sequence[tuple[str, RuleUse | None]]) -> FeatureLegendRows:
        """The rows of the drawn features, each given by its type and the rule it takes (None: no rule)."""

        used, default_types = legend_rule_use(drawn)
        return cls(
            tuple(rules),
            tuple(dict.fromkeys(feature_type for feature_type, _ in drawn)),
            frozenset(used),
            frozenset(default_types),
        )

    def type_rules(self, feature_type: str) -> list[LegendRule]:
        """The captioned rules of a type, its own and those for every type, in table order."""

        return [rule for rule in self.rules if rule[0] in (feature_type, "*") and str(rule[1] or "").strip()]

    def _drawn(self, rule: LegendRule) -> bool:
        return self.used is None or (rule[1], rule[2]) in self.used

    def rules_draw(self, feature_type: str) -> bool:
        """Whether a captioned rule of the type is used: its features of no rule then form an "other" row."""

        return any(self._drawn(rule) for rule in self.type_rules(feature_type))

    def type_row(self, feature_type: str) -> LegendRowKey | None:
        """The renderer's own row of a drawn type, if it has one: the type itself
        without rules of its type; else "other ..." while some of its features
        keep the type's color (always without the used rules); and the type
        again when none of its captioned rules is used."""

        if not any(rule[0] in (feature_type, "*") for rule in self.rules):
            return _key_of_name(feature_type)
        if self.used is None:
            return _other_feature_legend_key(feature_type)
        if self.rules_draw(feature_type):
            return _other_feature_legend_key(feature_type) if feature_type in self.default_types else None
        return _key_of_name(feature_type)

    def row_of(self, feature_type: str, winner: RuleUse | None) -> LegendRowKey | RuleUse | None:
        """The row a feature of the type is in, by the rule it takes (None: no
        rule): a captioned rule's row, as (caption, color), whose key the whole
        table allocates (``rule_legend_keys``); else its type's row or its
        "other" row; no row for an uncaptioned rule once captioned rules of the
        type are used. A feature that is not drawn gets the row it would join."""

        if winner is not None and str(winner[0] or "").strip():
            return winner
        if not self.rules_draw(feature_type):
            return _key_of_name(feature_type)
        return _other_feature_legend_key(feature_type) if winner is None else None

    def rows(self) -> Iterator[LegendRule | LegendRowKey]:
        """THE order of the feature rows of a legend: for each drawn type, in
        ``present`` order, its captioned rules in table order (a rule the
        features did not take has no row), then its renderer row. A rule row
        comes as its rule; each renderer row comes once. The legend adds the
        rows in this order; rows of one name and one swatch are one row, at the
        first (``LegendTable``)."""

        renderer_rows = {self.type_row(feature_type) for feature_type in self.present}
        placed: set[str] = set()
        for feature_type in self.present:
            yield from (rule for rule in self.type_rules(feature_type) if self._drawn(rule))
            for key in (_key_of_name(feature_type), _other_feature_legend_key(feature_type)):
                if key in renderer_rows and key not in placed:
                    placed.add(key)
                    yield key


def _generated_legend_rows(
    feature_rows: FeatureLegendRows,
    default_colors: DataFrame,
    *,
    gc: NumericLegendSource,
    skew: NumericLegendSource,
    depth: DepthLegendSource | None,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
    block_stroke: tuple[object, object] = ("none", 0),
) -> dict[str, dict[str, object]]:
    """The renderer's own legend rows by key, as drawn: feature types, "other ...", and numeric tracks."""

    # A type's row shows the color its features are drawn in (one owner, OV-3114).
    type_colors = build_default_color_map(default_colors)
    rows: dict[str, dict[str, object]] = {}
    for feature_type in feature_rows.present:
        key = feature_rows.type_row(feature_type)
        if key is not None:
            rows[key] = _solid_legend_entry(feature_type_color(type_colors, feature_type), *block_stroke)
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
    rules: list[LegendRule] = []
    if color_table is not None and not color_table.empty:
        rules = [(row["feature_type"], row["caption"], row["color"]) for _, row in color_table.iterrows()]
    feature_rows = FeatureLegendRows(
        tuple(rules),
        tuple(features_present),
        None if used_color_rules is None else frozenset(used_color_rules),
        frozenset(default_used_features or ()),
    )
    # The renderer's own rows (feature types, "other ..." and numeric tracks)
    # keep their keys; a rule row takes its key from the whole color table and
    # these rows (``rule_legend_keys``). The legend draws names, never keys.
    generated = _generated_legend_rows(
        feature_rows,
        features.default_colors,
        gc=gc,
        skew=skew,
        depth=depth,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=show_depth,
        block_stroke=block_stroke,
    )
    rule_keys = rule_legend_keys(
        ((caption, color) for _, caption, color in rules),
        {key: legend_swatch(entry) for key, entry in generated.items()},
        lambda fill: legend_swatch(_solid_legend_entry(fill, *block_stroke)),
    )
    for row in feature_rows.rows():
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


def sync_feature_density_legend_entries(
    out: LegendTable,
    sources: Sequence[FeatureDensityLegendSource],
    *,
    slots: Iterable[SlotView] | None,
) -> LegendTable:
    """One row per drawn ``feature_density`` slot, in slot order, named by its slot ``legend_label`` or its track's
    caption (both modes): a heatmap's gradient from 0 to its scale's top, a histogram's swatch."""

    by_id = {source.track_id: source for source in sources}
    for slot in slots or ():
        if not slot.enabled or slot.renderer != "feature_density":
            continue
        source = by_id.get(feature_density_slot_track_id(slot.id, slot.params))
        if source is None:
            continue
        label = slot_legend_label(slot, source.caption)
        if source.style == "heatmap":
            out.add_source(label, {
                "type": "gradient",
                "min_color": source.low_color,
                "max_color": source.color,
                "stroke": "none",
                "width": 0,
                "min_value": 0,
                "max_value": source.max_value,
                "min_label": "0",
                "max_label": source.max_text,
            })
        else:
            out.add_source(label, {
                "type": "solid",
                "fill": source.color,
                "stroke": source.stroke_color,
                "width": source.stroke_width,
            })
    return out


__all__ = [
    "DrawnLegendName",
    "LegendRowEdits",
    "LegendRowKey",
    "LegendRowName",
    "LegendRowNames",
    "LegendTable",
    "apply_legend_row_edits",
    "draws_legend",
    "PairwiseIdentityGradient",
    "sync_feature_density_legend_entries",
    "PairwiseIdentityLegend",
    "as_legend_table",
    "collect_comparison_collinearity_color_modes",
    "drawn_legend_names",
    "FeatureLegendRows",
    "LegendRule",
    "legend_row_names",
    "legend_swatch",
    "pairwise_identity_legend_from_comparisons",
    "prepare_legend_table",
    "rule_legend_keys",
    "write_legend_row_names",
]

