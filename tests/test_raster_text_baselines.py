"""CairoSVG exports keep text where browsers draw the SVG (B8)."""

from __future__ import annotations

import io
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from fontTools.ttLib import TTFont
from svgwrite import Drawing
from svgwrite.container import Group

from gbdraw.api import render as api_render
from gbdraw.render import export as export_module
from gbdraw.render.export import has_cairosvg
from gbdraw.svg.circular_ticks import (
    generate_circular_tick_labels,
    resolve_circular_tick_label_geometry,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FONT_FAMILY = "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
FONT_PATH = REPO_ROOT / "gbdraw" / "data" / "LiberationSans-Regular.ttf"
SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"

requires_cairosvg = pytest.mark.skipif(not has_cairosvg(), reason="CairoSVG is not installed")


def _tick_kwargs(font_size: float) -> dict:
    return {
        "total_len": 10_000,
        "size": "large",
        "font_size": font_size,
        "font_family": FONT_FAMILY,
        "track_type": "tuckin",
        "strandedness": True,
        "dpi": 96,
        "label_side": "outside",
        "tick_side": "inside",
        "tick_length_px": 10.0,
        "tick_width": 2.0,
    }


def _font_oracle(text: str) -> tuple[float, float, float, float]:
    """Return ascent, descent, ink top, and ink bottom of ``text`` in em units."""
    font = TTFont(FONT_PATH)
    units_per_em = font["head"].unitsPerEm
    cmap = font.getBestCmap()
    glyphs = [font["glyf"][cmap[ord(char)]] for char in text if not char.isspace()]
    ink_top = max(glyph.yMax for glyph in glyphs) / units_per_em
    ink_bottom = min(glyph.yMin for glyph in glyphs) / units_per_em
    return (
        font["hhea"].ascent / units_per_em,
        -font["hhea"].descent / units_per_em,
        ink_top,
        ink_bottom,
    )


def _tick_label_png(tick: int, *, radius: float, font_size: float, path_y: float) -> tuple[bytes, float]:
    """Rasterize one circular tick label with its path apex at ``path_y``."""
    kwargs = _tick_kwargs(font_size)
    elements = generate_circular_tick_labels(
        radius,
        kwargs["total_len"],
        kwargs["size"],
        [tick],
        "none",
        "black",
        font_size,
        "normal",
        FONT_FAMILY,
        kwargs["track_type"],
        kwargs["strandedness"],
        kwargs["dpi"],
        label_side=kwargs["label_side"],
        tick_side=kwargs["tick_side"],
        tick_length_px=kwargs["tick_length_px"],
        tick_width=kwargs["tick_width"],
    )
    geometry = resolve_circular_tick_label_geometry(
        center_radius_px=radius,
        tick=tick,
        label_text=f"{tick // 1000} kbp",
        **kwargs,
    )
    path_radius = geometry.path_radius_px
    drawing = Drawing(size=("400px", "240px"), viewBox="0 0 400 240")
    drawing.add(drawing.rect(insert=(0, 0), size=(400, 240), fill="white"))
    center_y = path_y + path_radius if tick == 0 else path_y - path_radius
    group = Group(transform=f"translate(200, {center_y})")
    for element in elements:
        group.add(element)
    drawing.add(group)
    return api_render.render_to_bytes(drawing, "png"), path_radius


def _ink_rows(png_bytes: bytes) -> tuple[int, int]:
    from PIL import Image

    image = Image.open(io.BytesIO(png_bytes)).convert("L")
    width, height = image.size
    pixels = image.load()
    rows = [y for y in range(height) if any(pixels[x, y] < 128 for x in range(width))]
    assert rows, "the tick label left no ink in the raster"
    return rows[0], rows[-1]


@requires_cairosvg
@pytest.mark.circular
def test_png_export_places_circular_tick_labels_on_browser_baselines() -> None:
    pytest.importorskip("PIL")
    font_size = 40.0
    path_y = 120.0
    tolerance = 3.0
    # Large radius keeps the label arc nearly flat over its width.
    radius = 2000.0

    # Lower half (180 degrees): text-before-edge on a reversed path. Browsers
    # put the font ascent on the path, so the label sits outside the axis.
    ascent, descent, ink_top, ink_bottom = _font_oracle("5 kbp")
    png_bytes, _ = _tick_label_png(5_000, radius=radius, font_size=font_size, path_y=path_y)
    top_row, bottom_row = _ink_rows(png_bytes)
    expected_baseline = path_y + ascent * font_size
    assert top_row == pytest.approx(expected_baseline - ink_top * font_size, abs=tolerance)
    assert bottom_row == pytest.approx(expected_baseline - ink_bottom * font_size, abs=tolerance)

    # Upper half (0 degrees): text-after-edge. Browsers put the font descent
    # on the path, so the baseline sits one descent outside the axis.
    ascent, descent, ink_top, ink_bottom = _font_oracle("0 kbp")
    png_bytes, _ = _tick_label_png(0, radius=radius, font_size=font_size, path_y=path_y)
    top_row, bottom_row = _ink_rows(png_bytes)
    expected_baseline = path_y - descent * font_size
    assert top_row == pytest.approx(expected_baseline - ink_top * font_size, abs=tolerance)
    assert bottom_row == pytest.approx(expected_baseline - ink_bottom * font_size, abs=tolerance)


def _svg(body: str) -> bytes:
    return (
        f'<?xml version="1.0" encoding="utf-8" ?>'
        f'<svg xmlns="{SVG_NS}" xmlns:xlink="{XLINK_NS}" width="200" height="100">{body}</svg>'
    ).encode("utf-8")


def _elements(svg_bytes: bytes, tag: str) -> list[ET.Element]:
    return list(ET.fromstring(svg_bytes).iter(f"{{{SVG_NS}}}{tag}"))


@pytest.mark.parametrize(
    ("baseline", "expected_em"),
    (
        ("text-before-edge", 1854 / 2048),
        ("text-after-edge", -434 / 2048),
        ("middle", 1082 / 2048 / 2),
        ("central", (1854 - 434) / 2048 / 2),
        ("hanging", 0.8 * 1854 / 2048),
    ),
)
def test_textpath_baseline_becomes_dy_for_cairosvg(baseline: str, expected_em: float) -> None:
    source = _svg(
        '<path id="p" d="M 0,50 L 200,50"/>'
        f'<text dominant-baseline="{baseline}"><textPath xlink:href="#p" font-size="20" '
        f'font-family="{FONT_FAMILY}" dominant-baseline="{baseline}">8 kbp</textPath></text>'
    )

    prepared = export_module.prepare_svg_for_cairosvg(source)

    (text_path,) = _elements(prepared, "textPath")
    assert text_path.get("dominant-baseline") == "alphabetic"
    assert float(text_path.get("dy")) == pytest.approx(expected_em * 20.0, abs=1e-3)
    assert text_path.get(f"{{{XLINK_NS}}}href") == "#p"


def test_cairosvg_rewrite_leaves_handled_text_unchanged() -> None:
    body = "".join(
        f'<text x="10" y="50" font-size="20" dominant-baseline="{baseline}">A</text>'
        for baseline in ("auto", "central", "text-before-edge", "text-after-edge")
    )
    body += (
        '<path id="p" d="M 0,50 L 200,50"/>'
        '<text><textPath xlink:href="#p" dominant-baseline="auto" font-size="20">label</textPath></text>'
    )
    source = _svg(body)

    assert export_module.prepare_svg_for_cairosvg(source) is source
    plain = _svg('<text x="10" y="50">A</text>')
    assert export_module.prepare_svg_for_cairosvg(plain) is plain


def test_cairosvg_rewrite_shifts_plain_hanging_and_middle_text() -> None:
    source = _svg(
        '<text x="10" y="50" font-size="20" dominant-baseline="hanging">Legend</text>'
        '<text x="10" y="80" font-size="10" dominant-baseline="middle"><tspan font-style="italic">Def</tspan></text>'
        '<text x="10" y="90" font-size="10" dominant-baseline="middle"><tspan y="95">Kept</tspan></text>'
    )

    hanging, middle, positioned = _elements(export_module.prepare_svg_for_cairosvg(source), "text")

    assert hanging.get("dominant-baseline") == "alphabetic"
    assert float(hanging.get("dy")) == pytest.approx(0.8 * 1854 / 2048 * 20.0, abs=1e-3)
    assert middle.get("dominant-baseline") == "alphabetic"
    assert float(middle.get("dy")) == pytest.approx(1082 / 2048 / 2 * 10.0, abs=1e-3)
    assert positioned.get("dominant-baseline") == "middle"
    assert positioned.get("dy") is None


def test_every_cairosvg_export_path_converts_the_rewritten_svg(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    converted: list[bytes] = []

    class CapturingCairoSvg:
        @staticmethod
        def svg2png(*, bytestring, write_to=None):
            converted.append(bytestring)
            if write_to is None:
                return b"png"
            if isinstance(write_to, str):
                Path(write_to).write_bytes(b"png")
            else:
                write_to.write(b"png")
            return None

    monkeypatch.setattr(export_module, "get_cairosvg", lambda: CapturingCairoSvg)

    def drawing(name: str) -> Drawing:
        canvas = Drawing(filename=str(tmp_path / f"{name}.svg"), size=("200px", "200px"))
        group = Group(transform="translate(100, 100)")
        for element in generate_circular_tick_labels(
            80.0, 10_000, "large", [5_000], "none", "black", 14.0, "normal", FONT_FAMILY,
            "tuckin", True, 96, label_side="outside", tick_side="inside",
            tick_length_px=10.0, tick_width=2.0,
        ):
            group.add(element)
        canvas.add(group)
        return canvas

    api_render.save_figure_to(drawing("strict"), ["svg", "png"])
    api_render.render_to_bytes(drawing("bytes"), "png")
    with pytest.warns(DeprecationWarning):
        export_module.save_figure(drawing("legacy"), ["svg", "png"])

    assert len(converted) == 3
    for source in converted:
        (text_path,) = _elements(source, "textPath")
        assert text_path.get("dominant-baseline") == "alphabetic"
        assert float(text_path.get("dy")) > 0
    for name in ("strict", "legacy"):
        saved = (tmp_path / f"{name}.svg").read_text(encoding="utf-8")
        assert 'dominant-baseline="text-before-edge"' in saved
        assert " dy=" not in saved


def test_only_the_export_owner_calls_cairosvg_converters() -> None:
    owner = REPO_ROOT / "gbdraw" / "render" / "export.py"
    pattern = re.compile(r"\.svg2(?:png|pdf|ps|eps|svg)\(|\bimport cairosvg\b")
    offenders = [
        str(path.relative_to(REPO_ROOT))
        for root in ("gbdraw", "tools")
        for path in sorted((REPO_ROOT / root).rglob("*.py"))
        if path != owner and pattern.search(path.read_text(encoding="utf-8"))
    ]
    assert offenders == []
