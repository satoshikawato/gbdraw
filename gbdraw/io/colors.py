#!/usr/bin/env python
# coding: utf-8

import logging
import re
import sys

# tomllib is available in Python 3.11+; use tomli as fallback for 3.10
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from importlib import resources
from typing import Mapping, Optional

import pandas as pd
from pandas import DataFrame

from ..core.color import normalize_hex_color
from ..exceptions import InputFileError, ParseError, ValidationError
from .table_text import read_literal_table

logger = logging.getLogger(__name__)

_COLOR_NAME_MAP = {
    "aliceblue": "#F0F8FF",
    "antiquewhite": "#FAEBD7",
    "aqua": "#00FFFF",
    "aquamarine": "#7FFFD4",
    "azure": "#F0FFFF",
    "beige": "#F5F5DC",
    "bisque": "#FFE4C4",
    "black": "#000000",
    "blanchedalmond": "#FFEBCD",
    "blue": "#0000FF",
    "blueviolet": "#8A2BE2",
    "brown": "#A52A2A",
    "burlywood": "#DEB887",
    "cadetblue": "#5F9EA0",
    "chartreuse": "#7FFF00",
    "chocolate": "#D2691E",
    "coral": "#FF7F50",
    "cornflowerblue": "#6495ED",
    "cornsilk": "#FFF8DC",
    "crimson": "#DC143C",
    "cyan": "#00FFFF",
    "darkblue": "#00008B",
    "darkcyan": "#008B8B",
    "darkgoldenrod": "#B8860B",
    "darkgray": "#A9A9A9",
    "darkgreen": "#006400",
    "darkgrey": "#A9A9A9",
    "darkkhaki": "#BDB76B",
    "darkmagenta": "#8B008B",
    "darkolivegreen": "#556B2F",
    "darkorange": "#FF8C00",
    "darkorchid": "#9932CC",
    "darkred": "#8B0000",
    "darksalmon": "#E9967A",
    "darkseagreen": "#8FBC8F",
    "darkslateblue": "#483D8B",
    "darkslategray": "#2F4F4F",
    "darkslategrey": "#2F4F4F",
    "darkturquoise": "#00CED1",
    "darkviolet": "#9400D3",
    "deeppink": "#FF1493",
    "deepskyblue": "#00BFFF",
    "dimgray": "#696969",
    "dimgrey": "#696969",
    "dodgerblue": "#1E90FF",
    "firebrick": "#B22222",
    "floralwhite": "#FFFAF0",
    "forestgreen": "#228B22",
    "fuchsia": "#FF00FF",
    "gainsboro": "#DCDCDC",
    "ghostwhite": "#F8F8FF",
    "gold": "#FFD700",
    "goldenrod": "#DAA520",
    "gray": "#808080",
    "grey": "#808080",
    "green": "#008000",
    "greenyellow": "#ADFF2F",
    "honeydew": "#F0FFF0",
    "hotpink": "#FF69B4",
    "indianred": "#CD5C5C",
    "indigo": "#4B0082",
    "ivory": "#FFFFF0",
    "khaki": "#F0E68C",
    "lavender": "#E6E6FA",
    "lavenderblush": "#FFF0F5",
    "lawngreen": "#7CFC00",
    "lemonchiffon": "#FFFACD",
    "lightblue": "#ADD8E6",
    "lightcoral": "#F08080",
    "lightcyan": "#E0FFFF",
    "lightgoldenrodyellow": "#FAFAD2",
    "lightgray": "#D3D3D3",
    "lightgreen": "#90EE90",
    "lightgrey": "#D3D3D3",
    "lightpink": "#FFB6C1",
    "lightsalmon": "#FFA07A",
    "lightseagreen": "#20B2AA",
    "lightskyblue": "#87CEFA",
    "lightslategray": "#778899",
    "lightslategrey": "#778899",
    "lightsteelblue": "#B0C4DE",
    "lightyellow": "#FFFFE0",
    "lime": "#00FF00",
    "limegreen": "#32CD32",
    "linen": "#FAF0E6",
    "magenta": "#FF00FF",
    "maroon": "#800000",
    "mediumaquamarine": "#66CDAA",
    "mediumblue": "#0000CD",
    "mediumorchid": "#BA55D3",
    "mediumpurple": "#9370DB",
    "mediumseagreen": "#3CB371",
    "mediumslateblue": "#7B68EE",
    "mediumspringgreen": "#00FA9A",
    "mediumturquoise": "#48D1CC",
    "mediumvioletred": "#C71585",
    "midnightblue": "#191970",
    "mintcream": "#F5FFFA",
    "mistyrose": "#FFE4E1",
    "moccasin": "#FFE4B5",
    "navajowhite": "#FFDEAD",
    "navy": "#000080",
    "oldlace": "#FDF5E6",
    "olive": "#808000",
    "olivedrab": "#6B8E23",
    "orange": "#FFA500",
    "orangered": "#FF4500",
    "orchid": "#DA70D6",
    "palegoldenrod": "#EEE8AA",
    "palegreen": "#98FB98",
    "paleturquoise": "#AFEEEE",
    "palevioletred": "#DB7093",
    "papayawhip": "#FFEFD5",
    "peachpuff": "#FFDAB9",
    "peru": "#CD853F",
    "pink": "#FFC0CB",
    "plum": "#DDA0DD",
    "powderblue": "#B0E0E6",
    "purple": "#800080",
    "rebeccapurple": "#663399",
    "red": "#FF0000",
    "rosybrown": "#BC8F8F",
    "royalblue": "#4169E1",
    "saddlebrown": "#8B4513",
    "salmon": "#FA8072",
    "sandybrown": "#F4A460",
    "seagreen": "#2E8B57",
    "seashell": "#FFF5EE",
    "sienna": "#A0522D",
    "silver": "#C0C0C0",
    "skyblue": "#87CEEB",
    "slateblue": "#6A5ACD",
    "slategray": "#708090",
    "slategrey": "#708090",
    "snow": "#FFFAFA",
    "springgreen": "#00FF7F",
    "steelblue": "#4682B4",
    "tan": "#D2B48C",
    "teal": "#008080",
    "thistle": "#D8BFD8",
    "tomato": "#FF6347",
    "turquoise": "#40E0D0",
    "violet": "#EE82EE",
    "wheat": "#F5DEB3",
    "white": "#FFFFFF",
    "whitesmoke": "#F5F5F5",
    "yellow": "#FFFF00",
    "yellowgreen": "#9ACD32",
}


def named_color_hex(color_name: str) -> str | None:
    """The hex code of an SVG/CSS color name (any case), or ``None``."""

    return _COLOR_NAME_MAP.get(color_name.lower())


_HEX_COLOR = re.compile(r"#(?:[0-9A-Fa-f]{3,4}|[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})")
_COLOR_FUNCTION = re.compile(r"(rgba?|hsla?)\((.*)\)", re.IGNORECASE | re.DOTALL)
_CSS_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_CSS_NUMBER_OR_PERCENT = re.compile(rf"{_CSS_NUMBER}%?")
_CSS_HUE = re.compile(rf"{_CSS_NUMBER}(?:deg|grad|rad|turn)?", re.IGNORECASE)
_PAINT_KEYWORDS = frozenset({"none", "transparent"})
# svgwrite's paint type takes these too, but the embedding document decides
# their color, so a user color is never one of them (OV-272).
_DOCUMENT_PAINT_KEYWORDS = frozenset({"currentcolor", "inherit"})
USER_COLOR_FORMS = "none, an SVG color name, #RGB, #RRGGBB, rgb(), or hsl()"


def _is_css_color_function(text: str) -> bool:
    match = _COLOR_FUNCTION.fullmatch(text)
    if match is None:
        return False
    arguments = match.group(2).strip()
    if "," in arguments:
        parts = [part.strip() for part in arguments.split(",")]
        if "/" in arguments or len(parts) not in (3, 4):
            return False
        channels, alpha = parts[:3], parts[3:]
    else:
        channel_text, _slash, alpha_text = arguments.partition("/")
        channels = channel_text.split()
        alpha = alpha_text.split() if _slash else []
        if len(channels) != 3 or (_slash and len(alpha) != 1):
            return False
    first = _CSS_HUE if match.group(1).lower().startswith("hsl") else _CSS_NUMBER_OR_PERCENT
    return bool(first.fullmatch(channels[0])) and all(
        _CSS_NUMBER_OR_PERCENT.fullmatch(part) for part in (*channels[1:], *alpha)
    )


def is_user_color(value: object) -> bool:
    """Whether ``value`` is a documented user color.

    ``none``, ``transparent`` and color names in any letter case,
    #RGB/#RGBA/#RRGGBB/#RRGGBBAA, and rgb()/rgba()/hsl()/hsla(): the forms the
    web app checks with the same vectors (tests/fixtures/default_color_domain.json).
    Not ``currentColor`` or ``inherit`` (OV-272), and not the other values of
    svgwrite's ``paint`` type: ``url()`` references, ``icc-color()``, and an
    empty value (Owner decision 2026-10-10, D-38).
    """

    if not isinstance(value, str):
        return False
    text = value.strip()
    lowered = text.lower()
    return lowered not in _DOCUMENT_PAINT_KEYWORDS and (
        lowered in _PAINT_KEYWORDS
        or lowered in _COLOR_NAME_MAP
        or _HEX_COLOR.fullmatch(text) is not None
        or _is_css_color_function(text)
    )


def check_user_color(
    value: object,
    *,
    where: str,
    diagnostic: Mapping[str, object],
) -> None:
    """Raise ``ValidationError`` unless ``value`` is a documented user color."""

    if not is_user_color(value):
        raise ValidationError(
            f"Invalid color {value!r} {where}. Use {USER_COLOR_FORMS}.",
            diagnostic=diagnostic,
        )


def resolve_color_to_hex(color_str: str) -> str:
    if not isinstance(color_str, str):
        raise ValidationError(f"Invalid color value (not a string): {color_str}.")

    if color_str.startswith("#"):
        check_user_color(
            color_str,
            where="(hex colors have 3, 4, 6, or 8 digits)",
            diagnostic={"code": "INPUT_INVALID", "reason": "COLOR"},
        )
        return color_str

    hex_code = named_color_hex(color_str)

    if hex_code:
        return hex_code

    raise ValidationError(
        f"Unknown color name: {color_str}. Please use a valid SVG color name or hex code."
    )


def load_default_colors(
    user_defined_default_colors: str,
    palette: str = "default",
    load_comparison: bool = False,
) -> DataFrame:
    column_names = ["feature_type", "color"]

    # ── 1) Load TOML
    try:
        toml_path = resources.files("gbdraw.data").joinpath("color_palettes.toml")
        with toml_path.open("rb") as fh:
            palettes_dict = tomllib.load(fh)
    except Exception as exc:
        logger.error(f"ERROR: failed to read colour_palettes.toml – {exc}")
        raise ParseError(f"Failed to read color_palettes.toml: {exc}") from exc

    if palette not in palettes_dict:
        logger.warning(f"Palette '{palette}' not found; using [default]")
        palette_dict = palettes_dict.get("default", {})
    else:
        palette_dict = palettes_dict[palette]

    default_colors = pd.DataFrame(palette_dict.items(), columns=column_names).set_index(
        "feature_type"
    )

    # ── 2) Apply user TSV overrides
    if user_defined_default_colors:
        try:
            user_df = (
                read_literal_table(
                    user_defined_default_colors,
                    names=column_names,
                    label="default colors file",
                    engine="c",
                ).set_index("feature_type")
            )
            # A `feature_type<TAB>color` header row (the Web writes one) names the columns.
            header = (user_df.index.astype(str).str.strip().str.lower() == "feature_type") & (
                user_df["color"].astype(str).str.strip().str.lower() == "color"
            )
            user_df = user_df[~header]
            # Drop rows with a missing or blank colour cell, as the web import does
            missing = user_df["color"].isna() | user_df["color"].astype(str).str.strip().eq("")
            if missing.any():
                for ft in user_df[missing].index.tolist():
                    logger.warning(
                        f"WARNING: colour missing for feature '{ft}' "
                        f"in '{user_defined_default_colors}' – "
                        "keeping built-in value."
                    )
                user_df = user_df[~missing]
            if load_comparison:  # if load_comparison is true, replace color names with hex codes
                for idx, row in user_df.iterrows():
                    resolved_color = resolve_color_to_hex(row["color"])
                    user_df.at[idx, "color"] = resolved_color

            default_colors = user_df.combine_first(default_colors)
            logger.info(f"User overrides applied: {user_defined_default_colors}")

        except FileNotFoundError:
            logger.error(
                f"ERROR: override file '{user_defined_default_colors}' not found"
            )
            raise InputFileError(
                f"Override file '{user_defined_default_colors}' not found"
            )
        except ParseError:
            raise
        except Exception as exc:
            logger.error(
                f"ERROR: failed to read '{user_defined_default_colors}' – {exc}"
            )
            raise ParseError(
                f"Failed to read '{user_defined_default_colors}': {exc}"
            ) from exc

    # ── 3) Return tidy DataFrame (index reset for downstream code)
    return default_colors.reset_index()


def _is_specific_table_color(color: str) -> bool:
    """Specific-color table domain: none, an SVG color name, #RGB, or #RRGGBB."""
    text = color.strip()
    if text.lower() == "none" or text.lower() in _COLOR_NAME_MAP:
        return True
    try:
        normalize_hex_color(text)
    except ValueError:
        return False
    return True


def read_color_table(color_table_file: str) -> Optional[DataFrame]:
    required_cols = ["feature_type", "qualifier_key", "value", "color", "caption"]
    mandatory_cols = ["feature_type", "qualifier_key", "value", "color"]

    # If user did not supply -t, just skip and return None
    if not color_table_file:
        return None

    try:
        df = read_literal_table(
            color_table_file,
            names=required_cols,
            label="color table",
            keep_default_na=False,  # "None", "NA" and "null" are values, not blanks
            na_values=[""],
        )
    except ParseError:
        raise
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in '{color_table_file}': {e}")
        raise ParseError(f"Malformed line in '{color_table_file}': {e}") from e
    except FileNotFoundError as e:
        logger.error(f"ERROR: Color table file not found: {e}")
        raise InputFileError(f"Color table file not found: {e}") from e
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{color_table_file}': {e}")
        raise ParseError(f"Failed to read '{color_table_file}': {e}") from e

    df["caption"] = df["caption"].fillna("")

    # Check for any rows with missing required values and error out if found
    null_rows = df[df[mandatory_cols].isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in mandatory_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{color_table_file}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        raise ValidationError(
            f"Missing values in '{color_table_file}'. See log for details."
        )

    for idx, color in df["color"].items():
        if not _is_specific_table_color(color):
            raise ValidationError(
                f"Invalid color {color!r} in '{color_table_file}' at line {idx + 1}. "
                "Use none, an SVG color name, #RGB, or #RRGGBB.",
                diagnostic={"code": "TABLE_INVALID", "field": "color", "reason": "COLOR", "row": idx + 1},
            )

    return df


__all__ = [
    "USER_COLOR_FORMS",
    "check_user_color",
    "is_user_color",
    "load_default_colors",
    "named_color_hex",
    "read_color_table",
    "resolve_color_to_hex",
]


