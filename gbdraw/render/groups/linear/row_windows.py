"""The drawing of a Linear row of record windows (PD-OI-122): its heading (one name and the numbered list of its
windows) and its marks (window numbers, edge bars and dashes, "//" with the length left out).

Every value comes from ``gbdraw.layout.linear_row_windows``; this module only draws it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from Bio.SeqRecord import SeqRecord

from gbdraw.svg.elements import Circle, Group, Line, Rect, Text, TSpan

from ....core.text import create_text_element, drawn_text_box, measure_mixed_text_line, parse_mixed_content_text
from ....layout.linear_row_windows import (
    LIST_LABEL_SEPARATOR,
    EdgeMark,
    LinearRowWindows,
    RowWindowLine,
    RowWindowStyle,
)
from ....layout.spatial import Aabb, union_aabbs
from .definition import DefinitionGroup, LinearDefinitionStyle, DefinitionLine

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation

ROW_MARK_ATTRIBUTE = "data-gbdraw-row-mark"
WINDOW_NUMBER_ATTRIBUTE = "data-gbdraw-window-number"


def _px(value: float) -> float:
    return round(float(value), 3) + 0.0


def circled_number(number: int, cx: float, cy: float, style: RowWindowStyle) -> Group:
    """A window's number in a circle, centred at ``(cx, cy)``."""

    group = Group()
    group.attribs[ROW_MARK_ATTRIBUTE] = "number"
    group.attribs[WINDOW_NUMBER_ATTRIBUTE] = str(number)
    group.add(
        Circle(
            center=(_px(cx), _px(cy)),
            r=_px(style.circle_radius),
            fill="white",
            stroke=style.text_fill,
            stroke_width=_px(style.circle_stroke),
        )
    )
    group.add(
        create_text_element(
            str(number), _px(cx), _px(cy), _px(style.number_font_size), "bold", style.font_family, fill=style.text_fill
        )
    )
    return group


class RowWindowHeadingGroup(DefinitionGroup):
    """A windowed row's heading at its left: the row name, then one line per window (its number, the record's own
    label when it differs from the name, its record ID and span). The name line is centred on the group's origin
    and the list hangs below it, a left-aligned block placed by the text anchor."""

    def __init__(
        self,
        record: SeqRecord,
        length_param: str,
        *,
        style: LinearDefinitionStyle,
        windows: LinearRowWindows,
        window_style: RowWindowStyle,
        text_anchor: str,
        group_id: str,
        record_index: int,
        record_presentation: RecordPresentation | None,
    ) -> None:
        self._windows = windows
        self._window_style = window_style
        super().__init__(
            record,
            length_param,
            style=style,
            text_anchor=text_anchor,
            group_id=group_id,
            line_kinds=None,
            record_index=record_index,
            definition_part="row",
            record_presentation=record_presentation,
        )

    def _list_parts(self, line: RowWindowLine) -> list[dict[str, str | bool | None]]:
        label = parse_mixed_content_text(line.label) if line.label else []
        return [*label, {"text": (LIST_LABEL_SEPARATOR if line.label else "") + line.text, "italic": False}]

    def build_definition_lines(self) -> list[DefinitionLine]:
        lines = [line for line in super().build_definition_lines() if line.kind == "name"]
        style = self._window_style
        for window_line in self._windows.lines:
            parts = self._list_parts(window_line)
            width, height = measure_mixed_text_line(
                [(str(part.get("text") or ""), "italic" if part.get("italic") else "normal") for part in parts],
                style.font_family,
                style.list_font_size,
                font_weight="bold" if style.list_font_weight == "bold" else "normal",
            )
            height = max(height, drawn_text_box("0", style.font_family, style.list_font_size)[1])
            lines.append(
                DefinitionLine(
                    kind="window",
                    text="".join(str(part.get("text") or "") for part in parts),
                    width=style.number_room + width,
                    height=max(height, 2.0 * (style.circle_radius + style.circle_stroke)),
                    parts=parts,
                )
            )
        return lines

    def calculate_start_coordinates(self) -> None:
        super().calculate_start_coordinates()
        if self.definition_lines:
            # The first line, not the block, is centred on the origin.
            shift = -self.definition_lines[0].y
            for line in self.definition_lines:
                line.y += shift

    def _list_left(self) -> float:
        block = max((line.width for line in self.definition_lines if line.kind == "window"), default=0.0)
        anchor = str(self.linear_text_anchor).strip().lower()
        return self.text_x if anchor == "start" else self.text_x - (block if anchor == "end" else 0.5 * block)

    def calculate_local_bounds(self) -> Aabb:
        name_bounds = super().calculate_local_bounds()
        left = self._list_left()
        boxes = [
            Aabb(left, line.y - 0.5 * line.height, left + line.width, line.y + 0.5 * line.height)
            for line in self.definition_lines
            if line.kind == "window"
        ]
        names = [name_bounds] if any(line.kind == "name" for line in self.definition_lines) else []
        merged = union_aabbs([*names, *boxes])
        return merged if merged is not None else Aabb(0.0, 0.0, 0.0, 0.0)

    def add_elements_to_group(self) -> None:
        style = self._window_style
        left = self._list_left()
        number = 0
        for line in self.definition_lines:
            if line.kind != "window":
                self.definition_group.add(self.create_name_text(line, line.y))
                continue
            window_line = self._windows.lines[number]
            number += 1
            self.definition_group.add(
                circled_number(window_line.number, left + style.circle_radius, line.y, style)
            )
            text = Text(
                "",
                insert=(_px(left + style.number_room), _px(line.y)),
                fill=style.list_fill,
                font_size=_px(style.list_font_size),
                font_weight=style.list_font_weight,
                font_family=style.font_family,
                text_anchor="start",
                dominant_baseline="middle",
            )
            text.attribs["data-definition-line-kind"] = "window"
            text.attribs[WINDOW_NUMBER_ATTRIBUTE] = str(window_line.number)
            for part in line.parts or []:
                part_text = part.get("text")
                if part_text:
                    text.add(TSpan(str(part_text), font_style="italic") if part.get("italic") else TSpan(str(part_text)))
            self.definition_group.add(text)


@dataclass(frozen=True)
class RowWindowPlace:
    """Where one unit of a windowed row is drawn: its left x and width, its axis y, and the top of its feature
    lanes from the axis (canvas coordinates without the horizontal offset)."""

    record_index: int
    x: float
    width: float
    axis_y: float
    feature_top_y: float


class RowWindowMarksGroup:
    """The marks of one windowed row: each window's number at its upper left, a bar or a dashed axis at each
    edge, and "//" (with the length left out above it) between adjacent windows of one record."""

    def __init__(
        self, windows: LinearRowWindows, places: tuple[RowWindowPlace, ...], style: RowWindowStyle, group_id: str
    ) -> None:
        self.group = Group(id=group_id)
        self.group.attribs["data-gbdraw-role"] = "row-window-marks"
        by_index = {place.record_index: place for place in places}
        for record_index in windows.record_indices:
            place = by_index[record_index]
            left, right = windows.edges_of(record_index)
            self._edge(place.x, -1.0, place.axis_y, left, style)
            self._edge(place.x + place.width, 1.0, place.axis_y, right, style)
            self.group.add(
                circled_number(
                    windows.number_of(record_index),
                    place.x - style.circle_gap - style.circle_radius,
                    place.axis_y + place.feature_top_y,
                    style,
                )
            )
        for gap in windows.gaps:
            left_place, right_place = by_index[gap.left_record_index], by_index[gap.right_record_index]
            self._gap(
                0.5 * (left_place.x + left_place.width + right_place.x),
                left_place.axis_y,
                left_place.axis_y + min(left_place.feature_top_y, right_place.feature_top_y),
                gap.text,
                style,
            )

    def _edge(self, x: float, side: float, axis_y: float, mark: EdgeMark, style: RowWindowStyle) -> None:
        if mark == "bar":
            bar_x = x + side * style.bar_offset
            line = Line(
                start=(_px(bar_x), _px(axis_y - style.bar_half)),
                end=(_px(bar_x), _px(axis_y + style.bar_half)),
                stroke=style.axis_stroke_color,
                stroke_width=_px(style.bar_stroke),
            )
        else:
            on, off = style.dash_pattern
            line = Line(
                start=(_px(x), _px(axis_y)),
                end=(_px(x + side * style.dash_length), _px(axis_y)),
                stroke=style.axis_stroke_color,
                stroke_width=_px(style.axis_stroke_width),
                stroke_dasharray=f"{_px(on)} {_px(off)}",
            )
        line.attribs[ROW_MARK_ATTRIBUTE] = mark
        self.group.add(line)

    def _gap(self, cx: float, axis_y: float, text_y: float, text: str, style: RowWindowStyle) -> None:
        gap = Group()
        gap.attribs[ROW_MARK_ATTRIBUTE] = "gap"
        for offset in (-style.slash_spacing, style.slash_spacing):
            gap.add(
                Line(
                    start=(_px(cx + offset - style.slash_slant), _px(axis_y + style.slash_half)),
                    end=(_px(cx + offset + style.slash_slant), _px(axis_y - style.slash_half)),
                    stroke=style.slash_color,
                    stroke_width=_px(style.slash_stroke),
                    stroke_linecap="round",
                )
            )
        if text:
            width, height = style.omitted_text_box(text)
            gap.add(
                Rect(
                    insert=(_px(cx - 0.5 * width - style.pad), _px(text_y - 0.5 * height)),
                    size=(_px(width + 2.0 * style.pad), _px(height)),
                    fill="white",
                    fill_opacity=0.9,
                )
            )
            gap.add(
                create_text_element(
                    text, _px(cx), _px(text_y), _px(style.list_font_size), "bold", style.font_family,
                    fill=style.text_fill,
                )
            )
        self.group.add(gap)

    def get_group(self) -> Group:
        return self.group


__all__ = [
    "RowWindowHeadingGroup",
    "RowWindowMarksGroup",
    "RowWindowPlace",
    "circled_number",
]
