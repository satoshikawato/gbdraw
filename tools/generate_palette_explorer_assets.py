#!/usr/bin/env python
"""Build the Circular palette explorer assets from one default-palette SVG."""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PALETTES_PATH = PROJECT_ROOT / "gbdraw" / "data" / "color_palettes.toml"
OUTPUT_DIR = PROJECT_ROOT / "gbdraw" / "web" / "gallery" / "palettes"
SVG_OUTPUT_PATH = OUTPUT_DIR / "circular.svg"
JSON_OUTPUT_PATH = OUTPUT_DIR / "palettes.json"

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
XLINK_NAMESPACE = "http://www.w3.org/1999/xlink"
PALETTE_ROLES = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "repeat_region",
    "skew_high",
    "skew_low",
    "gc_content",
)
REQUIRED_SVG_ROLES = frozenset(
    {"CDS", "rRNA", "tRNA", "tmRNA", "repeat_region", "skew_high", "skew_low", "gc_content"}
)
COLOR_ATTRIBUTES = ("fill", "stroke", "stop-color")
HEX_COLOR_RE = re.compile(r"^#[0-9a-f]{6}$", re.IGNORECASE)


def load_palettes() -> tuple[str, dict[str, dict[str, str]]]:
    with PALETTES_PATH.open("rb") as handle:
        source = tomllib.load(handle)
    title = str(source.pop("title", "gbdraw color palettes"))
    palettes = {
        str(name): {str(key): str(value) for key, value in values.items()}
        for name, values in source.items()
        if isinstance(values, dict)
    }
    return title, palettes


def palette_payload() -> dict[str, object]:
    title, palettes = load_palettes()
    return {
        "title": title,
        "defaultPalette": "default",
        "roles": list(PALETTE_ROLES),
        "palettes": palettes,
    }


def render_palette_json() -> str:
    return json.dumps(palette_payload(), indent=2, ensure_ascii=False) + "\n"


def normalize_hex_color(value: object) -> str | None:
    text = str(value or "").strip().lower()
    return text if HEX_COLOR_RE.fullmatch(text) else None


def default_role_colors(palettes: dict[str, dict[str, str]]) -> dict[str, str]:
    default = palettes.get("default", {})
    colors: dict[str, str] = {}
    for role in PALETTE_ROLES:
        color = normalize_hex_color(default.get(role))
        if color:
            colors[role] = color
    return colors


def annotate_svg(source_path: Path, output_path: Path) -> dict[str, int]:
    _, palettes = load_palettes()
    colors_by_role = default_role_colors(palettes)
    roles_by_color = {color: role for role, color in colors_by_role.items()}

    ET.register_namespace("", SVG_NAMESPACE)
    ET.register_namespace("xlink", XLINK_NAMESPACE)
    tree = ET.parse(source_path)
    root = tree.getroot()
    root.set("data-palette-template", "circular")
    root.set("role", "img")
    root.set("aria-label", "Circular genome diagram using the selected gbdraw color palette")

    counts = {role: 0 for role in PALETTE_ROLES}
    for element in root.iter():
        for attribute in COLOR_ATTRIBUTES:
            color = normalize_hex_color(element.get(attribute))
            role = roles_by_color.get(color or "")
            if not role:
                continue
            element.set("data-palette-key", role)
            element.set("data-palette-attribute", attribute)
            counts[role] += 1
            break

    missing = sorted(REQUIRED_SVG_ROLES.difference(role for role, count in counts.items() if count))
    if missing:
        raise ValueError(f"Source SVG does not expose required default-palette colors: {', '.join(missing)}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    return counts


def validate_svg(path: Path) -> dict[str, int]:
    tree = ET.parse(path)
    root = tree.getroot()
    if root.get("data-palette-template") != "circular":
        raise ValueError(f"{path} is not marked as a Circular palette template")
    counts = {role: 0 for role in PALETTE_ROLES}
    for element in root.iter():
        role = element.get("data-palette-key")
        attribute = element.get("data-palette-attribute")
        if role in counts and attribute in COLOR_ATTRIBUTES:
            counts[role] += 1
    missing = sorted(REQUIRED_SVG_ROLES.difference(role for role, count in counts.items() if count))
    if missing:
        raise ValueError(f"Palette template is missing semantic colors: {', '.join(missing)}")
    return counts


def check_assets() -> int:
    expected_json = render_palette_json()
    if not JSON_OUTPUT_PATH.exists() or JSON_OUTPUT_PATH.read_text(encoding="utf-8") != expected_json:
        print(f"Out of date: {JSON_OUTPUT_PATH.relative_to(PROJECT_ROOT)}")
        return 1
    if not SVG_OUTPUT_PATH.exists():
        print(f"Missing: {SVG_OUTPUT_PATH.relative_to(PROJECT_ROOT)}")
        return 1
    try:
        validate_svg(SVG_OUTPUT_PATH)
    except (ET.ParseError, ValueError) as exc:
        print(f"Invalid palette SVG: {exc}")
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-svg",
        type=Path,
        help="Default-palette Circular SVG to convert into the semantic palette template.",
    )
    parser.add_argument("--check", action="store_true", help="Validate the committed explorer assets.")
    args = parser.parse_args()

    if args.check:
        return check_assets()
    if args.source_svg is None:
        parser.error("--source-svg is required unless --check is used")
    if not args.source_svg.is_file():
        parser.error(f"source SVG does not exist: {args.source_svg}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    JSON_OUTPUT_PATH.write_text(render_palette_json(), encoding="utf-8")
    counts = annotate_svg(args.source_svg, SVG_OUTPUT_PATH)
    annotated = sum(counts.values())
    print(f"Wrote {JSON_OUTPUT_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {SVG_OUTPUT_PATH.relative_to(PROJECT_ROOT)} ({annotated} palette-aware elements)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
