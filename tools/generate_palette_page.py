#!/usr/bin/env python
"""Generate the compact palette reference from color_palettes.toml."""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PALETTES_PATH = PROJECT_ROOT / "gbdraw" / "data" / "color_palettes.toml"
OUTPUT_PATH = PROJECT_ROOT / "examples" / "color_palette_examples.md"
REPRESENTATIVE_PALETTES = ("default", "ajisai", "soft_pastels")


def _cell(palette: dict[str, str], key: str) -> str:
    color = palette.get(key, "—")
    if color == "—":
        return color
    return f'<span style="color:{color}">■</span> `{color}`'


def render_palette_page() -> str:
    with PALETTES_PATH.open("rb") as handle:
        source = tomllib.load(handle)
    palettes = {name: values for name, values in source.items() if name != "title"}

    lines = [
        "# Color palette reference",
        "",
        "This page is generated from [`gbdraw/data/color_palettes.toml`](../gbdraw/data/color_palettes.toml) by `python tools/generate_palette_page.py`. Choose a palette with `-p/--palette <name>`. The TOML file is the complete source of truth for every semantic color.",
        "",
        f"gbdraw currently provides **{len(palettes)} palettes**.",
        "",
        "## Interactive Circular preview",
        "",
        "Open the [Circular Palette Explorer](https://gbdraw.app/gallery/palettes/) to choose a palette and recolor one full-size SVG immediately. The explorer keeps the diagram geometry fixed and loads the palette values generated from the TOML source.",
        "",
        "## Representative full-size examples",
        "",
    ]
    for name in REPRESENTATIVE_PALETTES:
        lines.extend(
            [
                f"### `{name}`",
                "",
                f"![Circular diagram using the {name} palette](./AP027078_tuckin_separate_strands_{name}.svg)",
                "",
                f"![Linear comparison using the {name} palette](./hepatoplasmataceae_{name}.svg)",
                "",
            ]
        )

    lines.extend(
        [
            "## Primary color values",
            "",
            "The table keeps the most frequently compared roles visible. Open the TOML source for `tmRNA`, `ncRNA`, repeat, default, and collinear-block values.",
            "",
            "| Palette | CDS | rRNA | tRNA | Skew high | Skew low | GC content | Match maximum |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for name, palette in palettes.items():
        lines.append(
            "| "
            + " | ".join(
                (
                    f"`{name}`",
                    _cell(palette, "CDS"),
                    _cell(palette, "rRNA"),
                    _cell(palette, "tRNA"),
                    _cell(palette, "skew_high"),
                    _cell(palette, "skew_low"),
                    _cell(palette, "gc_content"),
                    _cell(palette, "pairwise_match_max"),
                )
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Reproduce the images",
            "",
            "```bash",
            "python tools/reproduce_examples.py --figure palette_circular_default",
            "python tools/generate_palette_explorer_assets.py --source-svg _reproduced/examples/AP027078_tuckin_separate_strands_default.svg",
            "```",
            "",
            "The explorer stores one Circular SVG plus the TOML-derived palette definitions. Only the three representative Circular/Linear pairs above are retained as standalone examples.",
            "",
            "## Accessibility",
            "",
            "Palette names describe appearance, not a guaranteed accessibility grade. Check contrast against the final background and at the intended print size. Do not encode a biological distinction by color alone: combine color with labels, track position, feature shape, stroke, or another redundant cue. For pairwise matches and quantitative tracks, verify that both endpoints and intermediate values remain distinguishable in the final export medium.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit non-zero when the committed page differs from generated content.",
    )
    args = parser.parse_args()
    rendered = render_palette_page()
    if args.check:
        return 0 if OUTPUT_PATH.read_text(encoding="utf-8") == rendered else 1
    OUTPUT_PATH.write_text(rendered, encoding="utf-8")
    print(f"Wrote {OUTPUT_PATH.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
