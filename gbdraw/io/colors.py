#!/usr/bin/env python
# coding: utf-8

import logging
import sys

# tomllib is available in Python 3.11+; use tomli as fallback for 3.10
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from importlib import resources
from typing import Optional

import pandas as pd
from pandas import DataFrame

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
    "red": "#FF0000",
    "rosybrown": "#BC8F8F",
    "royalblue": "#4169E1",
    "saddlebrown": "#8B4513",
    "salmon": "#FA8072",
    "sandybrown": "#F4A460",
    "seagreen": "#2E8B57",
    "seashell": "#2E8B57",
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


def resolve_color_to_hex(color_str: str) -> str:
    if not isinstance(color_str, str):
        logger.error(f"Invalid color value (not a string): {color_str}.")
        sys.exit(1)

    if color_str.startswith("#"):
        return color_str

    hex_code = _COLOR_NAME_MAP.get(color_str.lower())

    if hex_code:
        return hex_code

    logger.error(
        f"Unknown color name: {color_str}. Please use a valid SVG color name or hex code."
    )
    sys.exit(1)


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
        raise

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
                pd.read_csv(
                    user_defined_default_colors,
                    sep="\t",
                    names=column_names,
                    header=None,
                    dtype=str,
                ).set_index("feature_type")
            )
            # Drop rows with missing colour cells
            missing = user_df["color"].isna()
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
            sys.exit(1)
        except Exception as exc:
            logger.error(
                f"ERROR: failed to read '{user_defined_default_colors}' – {exc}"
            )
            sys.exit(1)

    # ── 3) Return tidy DataFrame (index reset for downstream code)
    return default_colors.reset_index()


def read_color_table(color_table_file: str) -> Optional[DataFrame]:
    required_cols = ["feature_type", "qualifier_key", "value", "color", "caption"]

    # If user did not supply -t, just skip and return None
    if not color_table_file:
        return None

    try:
        df = pd.read_csv(
            color_table_file,
            sep="\t",
            header=None,
            names=required_cols,
            dtype=str,
            on_bad_lines="error",  # raise on any row with wrong number of fields
            engine="python",  # required for on_bad_lines
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in '{color_table_file}': {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Color table file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{color_table_file}': {e}")
        sys.exit(1)

    # Check for any rows with missing values and error out if found
    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{color_table_file}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    return df


__all__ = [
    "load_default_colors",
    "read_color_table",
    "resolve_color_to_hex",
]


