"""TSV input for annotation sets."""

from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Literal

import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector

from .models import (
    AnnotationSet,
    CoordinateSpan,
    FeatureSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
    parse_feature_selector,
)


ANNOTATION_TABLE_REQUIRED_COLUMNS = frozenset({"set_id", "id", "mark"})
ANNOTATION_TABLE_COORDINATE_COLUMNS = frozenset(
    {"record", "start", "end", "coordinate_space", "wraps_origin", "out_of_bounds"}
)
ANNOTATION_TABLE_FEATURE_COLUMNS = frozenset(
    {"record", "feature_selector", "envelope", "circular_path"}
)
ANNOTATION_TABLE_OPTIONAL_COLUMNS = frozenset(
    {
        "label",
        "lane",
        "legend_label",
        "stroke",
        "stroke_width",
        "stroke_dasharray",
        "line_cap",
        "fill",
        "fill_opacity",
        "hatch_angle",
        "hatch_spacing",
        "hatch_color",
        "hatch_width",
        "hatch_cross",
        "label_color",
        "label_font_size",
        "label_orientation",
        "label_position",
        "label_offset",
    }
)
ANNOTATION_TABLE_COLUMNS = (
    ANNOTATION_TABLE_REQUIRED_COLUMNS
    | ANNOTATION_TABLE_COORDINATE_COLUMNS
    | ANNOTATION_TABLE_FEATURE_COLUMNS
    | ANNOTATION_TABLE_OPTIONAL_COLUMNS
)


def _text(row: pd.Series, column: str) -> str:
    value = row.get(column, "")
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _error(row_number: int, column: str, message: str) -> ValidationError:
    return ValidationError(f"Annotation table row {row_number}, column {column!r}: {message}")


def _bool(text: str, *, row_number: int, column: str, default: bool = False) -> bool:
    if not text:
        return default
    normalized = text.lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise _error(row_number, column, "expected true/false")


def _number(text: str, *, row_number: int, column: str, integer: bool = False):
    try:
        return int(text) if integer else float(text)
    except (TypeError, ValueError) as exc:
        kind = "integer" if integer else "number"
        raise _error(row_number, column, f"expected a {kind}, got {text!r}") from exc


def _style(row: pd.Series, row_number: int) -> RegionAnnotationStyle | None:
    style_columns = {
        "stroke",
        "stroke_width",
        "stroke_dasharray",
        "line_cap",
        "fill",
        "fill_opacity",
        "label_color",
        "label_font_size",
        "label_orientation",
        "label_position",
        "label_offset",
    }
    hatch_columns = {"hatch_angle", "hatch_spacing", "hatch_color", "hatch_width", "hatch_cross"}
    if not any(_text(row, column) for column in style_columns | hatch_columns):
        return None

    defaults = RegionAnnotationStyle()
    kwargs = {field.name: getattr(defaults, field.name) for field in fields(RegionAnnotationStyle)}
    if _text(row, "stroke"):
        kwargs["stroke"] = _text(row, "stroke")
    if _text(row, "stroke_width"):
        kwargs["stroke_width"] = _number(_text(row, "stroke_width"), row_number=row_number, column="stroke_width")
    if _text(row, "stroke_dasharray"):
        raw = _text(row, "stroke_dasharray").replace(",", " ")
        kwargs["stroke_dasharray"] = tuple(
            _number(item, row_number=row_number, column="stroke_dasharray")
            for item in raw.split()
        )
    if _text(row, "line_cap"):
        kwargs["line_cap"] = _text(row, "line_cap").lower()
    if "fill" in row.index and _text(row, "fill"):
        kwargs["fill"] = _text(row, "fill")
    if _text(row, "fill_opacity"):
        kwargs["fill_opacity"] = _number(_text(row, "fill_opacity"), row_number=row_number, column="fill_opacity")
    if _text(row, "label_color"):
        kwargs["label_color"] = _text(row, "label_color")
    if _text(row, "label_font_size"):
        kwargs["label_font_size"] = _number(_text(row, "label_font_size"), row_number=row_number, column="label_font_size")
    if _text(row, "label_orientation"):
        kwargs["label_orientation"] = _text(row, "label_orientation").lower()
    if _text(row, "label_position"):
        kwargs["label_position"] = _text(row, "label_position").lower()
    if _text(row, "label_offset"):
        kwargs["label_offset"] = _number(_text(row, "label_offset"), row_number=row_number, column="label_offset")

    if any(_text(row, column) for column in hatch_columns):
        hatch_defaults = HatchStyle()
        kwargs["hatch"] = HatchStyle(
            angle=(
                _number(_text(row, "hatch_angle"), row_number=row_number, column="hatch_angle")
                if _text(row, "hatch_angle")
                else hatch_defaults.angle
            ),
            spacing=(
                _number(_text(row, "hatch_spacing"), row_number=row_number, column="hatch_spacing")
                if _text(row, "hatch_spacing")
                else hatch_defaults.spacing
            ),
            color=_text(row, "hatch_color") or hatch_defaults.color,
            width=(
                _number(_text(row, "hatch_width"), row_number=row_number, column="hatch_width")
                if _text(row, "hatch_width")
                else hatch_defaults.width
            ),
            cross=_bool(
                _text(row, "hatch_cross"),
                row_number=row_number,
                column="hatch_cross",
                default=hatch_defaults.cross,
            ),
        )
    try:
        return RegionAnnotationStyle(**kwargs)
    except ValidationError as exc:
        raise ValidationError(f"Annotation table row {row_number}: {exc}") from exc


def _target(
    row: pd.Series,
    row_number: int,
    *,
    mode: Literal["circular", "linear"] | None,
):
    has_coordinates = bool(_text(row, "start") or _text(row, "end"))
    has_features = bool(_text(row, "feature_selector"))
    if has_coordinates == has_features:
        raise _error(
            row_number,
            "start/feature_selector",
            "provide exactly one coordinate target or feature_selector target",
        )
    record_text = _text(row, "record")
    try:
        record = parse_record_selector(record_text or None)
    except ValueError as exc:
        raise _error(row_number, "record", str(exc)) from exc

    if has_coordinates:
        if not _text(row, "start") or not _text(row, "end"):
            raise _error(row_number, "start/end", "both start and end are required")
        wraps = _bool(
            _text(row, "wraps_origin"),
            row_number=row_number,
            column="wraps_origin",
        )
        if mode == "linear" and wraps:
            raise _error(row_number, "wraps_origin", "origin-spanning targets are circular-only")
        return CoordinateSpan(
            record=record,
            start=_number(_text(row, "start"), row_number=row_number, column="start", integer=True),
            end=_number(_text(row, "end"), row_number=row_number, column="end", integer=True),
            coordinate_space=(_text(row, "coordinate_space") or "source").lower(),
            wraps_origin=wraps,
            out_of_bounds=(_text(row, "out_of_bounds") or "clip").lower(),
        )

    raw_selectors = _text(row, "feature_selector")
    selector_parts = [item.strip() for item in raw_selectors.split(";") if item.strip()]
    return FeatureSpan(
        record=record,
        selectors=tuple(parse_feature_selector(item) for item in selector_parts),
        envelope=(_text(row, "envelope") or "outer_bounds").lower(),
        circular_path=(_text(row, "circular_path") or "shortest").lower(),
    )


def annotation_sets_from_dataframe(
    dataframe: DataFrame,
    *,
    mode: Literal["circular", "linear"] | None = None,
) -> tuple[AnnotationSet, ...]:
    """Validate a materialized table and preserve first set appearance order."""

    if not isinstance(dataframe, DataFrame):
        raise ValidationError("Annotation table must be a pandas DataFrame.")
    columns = {str(column).strip() for column in dataframe.columns}
    missing = ANNOTATION_TABLE_REQUIRED_COLUMNS - columns
    if missing:
        raise ValidationError(f"Annotation table is missing required columns: {', '.join(sorted(missing))}.")
    unknown = columns - ANNOTATION_TABLE_COLUMNS
    if unknown:
        raise ValidationError(f"Annotation table has unknown columns: {', '.join(sorted(unknown))}.")

    normalized = dataframe.copy()
    normalized.columns = [str(column).strip() for column in dataframe.columns]
    grouped: dict[str, list[RegionAnnotation]] = {}
    for row_index, row in normalized.iterrows():
        row_number = int(row_index) + 2 if isinstance(row_index, int) else len(sum(grouped.values(), [])) + 2
        set_id = _text(row, "set_id")
        annotation_id = _text(row, "id")
        if not set_id:
            raise _error(row_number, "set_id", "value cannot be empty")
        if not annotation_id:
            raise _error(row_number, "id", "value cannot be empty")
        mark = _text(row, "mark").lower()
        if not mark:
            raise _error(row_number, "mark", "value cannot be empty")
        lane_text = _text(row, "lane")
        try:
            annotation = RegionAnnotation(
                id=annotation_id,
                target=_target(row, row_number, mode=mode),
                label=_text(row, "label"),
                mark=mark,
                lane=(
                    _number(lane_text, row_number=row_number, column="lane", integer=True)
                    if lane_text
                    else None
                ),
                style=_style(row, row_number),
                legend_label=_text(row, "legend_label") or None,
            )
        except ValidationError as exc:
            if str(exc).startswith("Annotation table row"):
                raise
            raise ValidationError(f"Annotation table row {row_number}: {exc}") from exc
        grouped.setdefault(set_id, []).append(annotation)

    return tuple(AnnotationSet(id=set_id, annotations=tuple(items)) for set_id, items in grouped.items())


def read_annotation_table(
    path: str | Path,
    *,
    mode: Literal["circular", "linear"] | None = None,
) -> tuple[AnnotationSet, ...]:
    """Read a UTF-8 TSV annotation table into typed annotation sets."""

    try:
        dataframe = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        raise ValidationError(f"Unable to read annotation table {str(path)!r}: {exc}") from exc
    return annotation_sets_from_dataframe(dataframe, mode=mode)


__all__ = [
    "ANNOTATION_TABLE_COLUMNS",
    "ANNOTATION_TABLE_REQUIRED_COLUMNS",
    "annotation_sets_from_dataframe",
    "read_annotation_table",
]
