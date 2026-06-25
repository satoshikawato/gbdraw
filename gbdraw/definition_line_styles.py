#!/usr/bin/env python
# coding: utf-8

"""Parsing and normalization helpers for linear definition line styles."""

from __future__ import annotations

import math
from typing import Iterable, Mapping

DEFINITION_LINE_KINDS: tuple[str, ...] = ("name", "replicon", "accession", "length")
DEFINITION_LINE_STYLE_PROPERTIES: tuple[str, ...] = ("font_size", "font_weight", "fill")

_LINE_KEY_ALIASES = {
    "name": "name",
    "species": "name",
    "record_label": "name",
    "record-label": "name",
    "replicon": "replicon",
    "accession": "accession",
    "length": "length",
    "coordinates": "length",
}

_PROPERTY_ALIASES = {
    "size": "font_size",
    "font_size": "font_size",
    "font-size": "font_size",
    "weight": "font_weight",
    "font_weight": "font_weight",
    "font-weight": "font_weight",
    "color": "fill",
    "fill": "fill",
}

_NAMED_FONT_WEIGHTS = {"normal", "bold", "lighter", "bolder"}
_NUMERIC_FONT_WEIGHTS = {str(weight) for weight in range(100, 1000, 100)}


def normalize_definition_line_key(value: object) -> str:
    key = str(value or "").strip().lower()
    normalized = _LINE_KEY_ALIASES.get(key)
    if normalized is None:
        allowed = ", ".join(sorted(_LINE_KEY_ALIASES))
        raise ValueError(f"unknown definition line key '{value}': expected one of {allowed}")
    return normalized


def normalize_definition_line_style_property(value: object) -> str:
    key = str(value or "").strip().lower()
    normalized = _PROPERTY_ALIASES.get(key)
    if normalized is None:
        allowed = ", ".join(sorted(_PROPERTY_ALIASES))
        raise ValueError(f"unknown definition line style property '{value}': expected one of {allowed}")
    return normalized


def normalize_definition_line_font_size(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "auto", "none", "null", "default"}:
            return None
        value = normalized
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValueError("definition line font_size must be a positive finite number")
    return parsed


def normalize_definition_line_font_weight(value: object) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip().lower()
    if normalized in {"", "auto", "none", "null", "default"}:
        return None
    if normalized in _NAMED_FONT_WEIGHTS:
        return normalized
    try:
        numeric = float(normalized)
    except ValueError as exc:
        raise ValueError(
            "definition line font_weight must be normal, bold, lighter, bolder, or 100 through 900"
        ) from exc
    if numeric.is_integer() and str(int(numeric)) in _NUMERIC_FONT_WEIGHTS:
        return str(int(numeric))
    raise ValueError("definition line font_weight must be normal, bold, lighter, bolder, or 100 through 900")


def normalize_definition_line_fill(value: object, *, allow_empty: bool = True) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if normalized:
        return normalized
    if allow_empty:
        return None
    raise ValueError("definition line fill/color must not be empty")


def split_definition_line_style_properties(text: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    depth = 0
    quote: str | None = None
    escaped = False

    for char in str(text):
        if escaped:
            current.append(char)
            escaped = False
            continue
        if quote is not None:
            current.append(char)
            if char == "\\":
                escaped = True
            elif char == quote:
                quote = None
            continue
        if char in {"'", '"'}:
            quote = char
            current.append(char)
            continue
        if char == "(":
            depth += 1
            current.append(char)
            continue
        if char == ")":
            depth = max(0, depth - 1)
            current.append(char)
            continue
        if char == "," and depth == 0:
            part = "".join(current).strip()
            if not part:
                raise ValueError("definition line style property must not be empty")
            parts.append(part)
            current = []
            continue
        current.append(char)

    part = "".join(current).strip()
    if not part:
        raise ValueError("definition line style property must not be empty")
    parts.append(part)
    return parts


def parse_definition_line_style_assignment(raw: str) -> tuple[str, dict[str, object]]:
    text = str(raw).strip()
    line_key_raw, separator, properties_raw = text.partition(":")
    if not separator:
        raise ValueError("definition line style must be in LINE:PROPERTY=VALUE format")
    line_key = normalize_definition_line_key(line_key_raw)
    if not properties_raw.strip():
        raise ValueError("definition line style properties must not be empty")

    style: dict[str, object] = {}
    for property_assignment in split_definition_line_style_properties(properties_raw):
        property_raw, property_separator, value_raw = property_assignment.partition("=")
        if not property_separator:
            raise ValueError("definition line style properties must be KEY=VALUE pairs")
        property_key = normalize_definition_line_style_property(property_raw)
        value_text = value_raw.strip()
        if not value_text:
            raise ValueError(f"definition line style {property_key} value must not be empty")
        if property_key == "font_size":
            value: object = normalize_definition_line_font_size(value_text)
        elif property_key == "font_weight":
            value = normalize_definition_line_font_weight(value_text)
        else:
            value = normalize_definition_line_fill(value_text, allow_empty=False)
        style[property_key] = value

    return line_key, style


def parse_definition_line_style_overrides(
    assignments: Iterable[str] | None,
) -> dict[str, dict[str, object]]:
    overrides: dict[str, dict[str, object]] = {}
    if assignments is None:
        return overrides
    for raw in assignments:
        line_key, style = parse_definition_line_style_assignment(raw)
        overrides.setdefault(line_key, {}).update(style)
    return overrides


def normalize_definition_line_style_mapping(
    line_styles: Mapping[str, Mapping[str, object]] | None,
) -> dict[str, dict[str, object]]:
    normalized: dict[str, dict[str, object]] = {}
    if line_styles is None:
        return normalized
    for line_key_raw, raw_style in line_styles.items():
        line_key = normalize_definition_line_key(line_key_raw)
        if not isinstance(raw_style, Mapping):
            raise ValueError("definition line style entries must be mappings")
        style: dict[str, object] = {}
        for property_raw, raw_value in raw_style.items():
            property_key = normalize_definition_line_style_property(property_raw)
            if property_key == "font_size":
                value = normalize_definition_line_font_size(raw_value)
            elif property_key == "font_weight":
                value = normalize_definition_line_font_weight(raw_value)
            else:
                value = normalize_definition_line_fill(raw_value)
            if value is not None:
                style[property_key] = value
        normalized.setdefault(line_key, {}).update(style)
    return normalized


__all__ = [
    "DEFINITION_LINE_KINDS",
    "DEFINITION_LINE_STYLE_PROPERTIES",
    "normalize_definition_line_fill",
    "normalize_definition_line_font_size",
    "normalize_definition_line_font_weight",
    "normalize_definition_line_key",
    "normalize_definition_line_style_mapping",
    "normalize_definition_line_style_property",
    "parse_definition_line_style_assignment",
    "parse_definition_line_style_overrides",
    "split_definition_line_style_properties",
]
