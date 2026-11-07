"""PNG and PDF size limits: a wrapped figure stops with a message naming what would fit; SVG never stops."""

from __future__ import annotations

import json
import struct
from pathlib import Path

import pytest
from svgwrite import Drawing

from gbdraw import render
from gbdraw.api import render_to_bytes, save_figure_to
from gbdraw.api.requests import GenBankInputSource, LinearDiagramOptions, LinearDiagramRequest, RecordInput
from gbdraw.api.request_render import plan_request
from gbdraw.exceptions import ExportError
from gbdraw.layout.linear_wrap import LinearWrapSummary, longer_row_lengths
from gbdraw.render import export_limits
from gbdraw.render.export import convert_svg_with_cairosvg
from gbdraw.render.export_limits import bind_export_wrap, check_export_limit
from gbdraw.web_support.export_limit import check_web_export_limit_json

SORANGIUM_BP = 13_033_779
# Measured: examples/NC_010162.gb at 50 kb per row, 261 rows.
SORANGIUM = LinearWrapSummary(50_000, 50_000, (SORANGIUM_BP,))
SORANGIUM_SIZE = (2_000.0, 45_967.0)
LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"


def _refusal(fmt: str, size: tuple[float, float], wrap: LinearWrapSummary | None, dpi: float | None = None) -> str:
    with pytest.raises(ExportError) as caught:
        check_export_limit(fmt, *size, wrap=wrap, dpi=dpi)
    return str(caught.value)


def _drawing(tmp_path: Path, width: int, height: int, wrap: LinearWrapSummary | None = None) -> Drawing:
    drawing = Drawing(filename=str(tmp_path / "figure.svg"), size=(f"{width}px", f"{height}px"))
    bind_export_wrap(drawing, wrap)
    return drawing


def test_a_sorangium_png_names_the_limit_the_next_row_length_that_fits_and_svg() -> None:
    message = _refusal("png", SORANGIUM_SIZE, SORANGIUM, dpi=96)
    assert "32,767 px" in message and "at 96 dpi." in message
    # 131 rows of the 261 drawn: about 23,072 px.
    assert "--wrap_length 100kb" in message and "23,072 px" in message
    assert "estimated from the drawn height" in message
    assert "SVG" in message
    # 68 dpi is the last that fits; the caller's DPI is named by its value only.
    assert "A DPI of 68 or lower would fit." in message
    assert "default" not in message


def test_a_png_without_a_dpi_is_checked_at_the_default_and_names_it() -> None:
    """A PNG without a DPI is drawn at 300 dpi (DEFAULT_PNG_DPI); its refusal says so and how to draw it smaller."""

    message = _refusal("png", SORANGIUM_SIZE, SORANGIUM)
    assert message == _refusal("png", SORANGIUM_SIZE, SORANGIUM, dpi=300).replace(
        "at 300 dpi.", "at 300 dpi (the default PNG DPI)."
    ).replace(
        "A DPI of 68 or lower would fit.",
        "A DPI of 68 or lower would fit; --dpi on the command line or dpi= in Python sets it.",
    )
    assert "143,647 px at 300 dpi (the default PNG DPI)." in message
    # 10,485 SVG px is 32,766 px at 300 dpi; 10,486 is over.
    check_export_limit("png", 10_485.0, 10.0)
    assert "A DPI of 299 or lower would fit" in _refusal("png", (10_486.0, 10.0), None)


def test_a_sorangium_pdf_skips_the_row_length_that_is_still_too_tall() -> None:
    message = _refusal("pdf", SORANGIUM_SIZE, SORANGIUM)
    assert "14,400 pt" in message
    # 100 kb is 131 rows, 17,303 pt; 200 kb is 66 rows, 8,718 pt.
    assert "34,475 pt" in message and "--wrap_length 200kb" in message and "8,718 pt" in message
    assert "SVG" in message
    assert "DPI" not in message


def test_a_dpi_from_the_caller_scales_the_pixels_and_names_a_lower_dpi() -> None:
    message = _refusal("png", SORANGIUM_SIZE, SORANGIUM, dpi=300)
    # 45,967 px at 300 dpi is 143,647 px; 68 dpi is the last that fits; 27 rows at 500 kb fit at 300 dpi.
    assert "143,647 px" in message
    assert "DPI of 68 or lower" in message
    assert "--wrap_length 500kb" in message
    # The same figure at 68 dpi fits: the caller's DPI decides.
    check_export_limit("png", *SORANGIUM_SIZE, wrap=SORANGIUM, dpi=68)


def test_a_pdf_never_names_a_dpi_a_pdf_page_does_not_depend_on_it() -> None:
    assert "DPI" not in _refusal("pdf", SORANGIUM_SIZE, SORANGIUM, dpi=300)


def test_no_row_length_fits_when_even_one_row_each_is_too_tall() -> None:
    many_records = LinearWrapSummary(50_000, 50_000, (60_000,) * 200)
    message = _refusal("png", (2_000.0, 80_000.0), many_records, dpi=96)
    assert "No longer row length would fit" in message
    assert "--wrap_length" not in message
    assert "SVG" in message


def test_no_row_length_fits_when_the_width_is_over_the_limit() -> None:
    message = _refusal("png", (40_000.0, 100.0), LinearWrapSummary(50_000, 50_000, (120_000,)), dpi=96)
    assert "No longer row length would fit" in message


def test_the_last_row_length_that_fits_may_draw_one_row() -> None:
    assert "--wrap_length 100kb" in _refusal(
        "png", (2_000.0, 40_000.0), LinearWrapSummary(50_000, 50_000, (60_000,)), dpi=96
    )


def test_row_lengths_run_up_the_1_2_5_ladder_to_the_first_that_holds_the_longest_record() -> None:
    assert longer_row_lengths(50_000, 300_000) == (100_000, 200_000, 500_000)
    assert longer_row_lengths(15_000, 19_000) == (20_000,)
    assert longer_row_lengths(10_000, 10_001) == (20_000,)


@pytest.mark.parametrize(
    ("fmt", "fits", "over"),
    [
        ("png", (32_767.0, 10.0), (32_768.0, 10.0)),
        ("png", (10.0, 32_767.0), (10.0, 32_768.0)),
        # CairoSVG rounds a PNG's size to whole pixels.
        ("png", (10.0, 32_767.4), (10.0, 32_767.6)),
        ("pdf", (19_200.0, 10.0), (19_201.0, 10.0)),
        ("pdf", (10.0, 19_200.0), (10.0, 19_201.0)),
    ],
)
def test_the_limit_is_on_each_side(fmt: str, fits: tuple[float, float], over: tuple[float, float]) -> None:
    """At 96 dpi a PNG has one pixel per SVG pixel."""

    wrap = LinearWrapSummary(50_000, 50_000, (120_000,))
    check_export_limit(fmt, *fits, wrap=wrap, dpi=96)
    _refusal(fmt, over, wrap, dpi=96)


@pytest.mark.parametrize("fmt", ["svg", "interactive_svg", "eps", "ps"])
def test_svg_and_the_formats_without_a_limit_are_never_refused(fmt: str) -> None:
    check_export_limit(fmt, 100_000.0, 100_000.0, wrap=SORANGIUM)


def test_an_unwrapped_png_over_the_limit_fails_with_the_limit_named_as_it_failed_before(tmp_path: Path) -> None:
    drawing = _drawing(tmp_path, 33_000, 100)
    with pytest.raises(ExportError, match=r"32,767 px.*SVG") as caught:
        save_figure_to(drawing, ["svg", "png"], output_prefix=str(tmp_path / "figure"))
    assert "--wrap_length" not in str(caught.value)
    assert not list(tmp_path.iterdir())


def test_an_unwrapped_pdf_over_the_limit_is_still_written(tmp_path: Path) -> None:
    pytest.importorskip("cairosvg")
    drawing = _drawing(tmp_path, 20_000, 100)
    paths = save_figure_to(drawing, ["pdf"], output_prefix=str(tmp_path / "figure"))
    assert any(path.endswith(".pdf") for path in paths)
    assert render_to_bytes(drawing, "pdf").startswith(b"%PDF")


def test_a_wrapped_refusal_stops_before_anything_is_converted_or_written(tmp_path: Path, monkeypatch) -> None:
    def no_conversion(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("converted")

    monkeypatch.setattr(render.export, "convert_svg_with_cairosvg", no_conversion)
    existing = tmp_path / "figure.svg"
    existing.write_text("before", encoding="utf-8")
    drawing = _drawing(tmp_path, 2_000, 45_967, SORANGIUM)
    with pytest.raises(ExportError, match="--wrap_length 100kb"):
        save_figure_to(drawing, ["svg", "png"], output_prefix=str(tmp_path / "figure"), overwrite=True, dpi=96)
    with pytest.raises(ExportError, match="--wrap_length 200kb"):
        render_to_bytes(drawing, "pdf")
    assert existing.read_text(encoding="utf-8") == "before"
    assert [path.name for path in tmp_path.iterdir()] == ["figure.svg"]


def test_a_wrapped_svg_over_every_limit_is_written(tmp_path: Path) -> None:
    drawing = _drawing(tmp_path, 40_000, 45_967, SORANGIUM)
    paths = save_figure_to(drawing, ["svg", "interactive_svg"], output_prefix=str(tmp_path / "figure"))
    assert len(paths) == 2 and all(Path(path).stat().st_size for path in paths)
    assert render_to_bytes(drawing, "svg").startswith(b"<?xml") or b"<svg" in render_to_bytes(drawing, "svg")


def test_a_drawn_wrapped_figure_carries_its_rows_to_the_check(tmp_path: Path, monkeypatch) -> None:
    built = plan_request(
        LinearDiagramRequest(
            records=(RecordInput(GenBankInputSource(str(LAMBDA)), record_key="r0"),),
            options=LinearDiagramOptions(config_overrides={"canvas.linear.wrap_length": 10_000}),
        )
    ).build()
    monkeypatch.setattr(export_limits, "PNG_MAX_PIXELS", 10)
    with pytest.raises(ExportError, match="No longer row length would fit"):
        render_to_bytes(built.drawing, "png")


@pytest.mark.parametrize(
    ("size", "dpi"),
    [
        # 1,236 x 3.125 = 3,862.5: CairoSVG rounds half to even, as Python does (a JavaScript Math.round gives 3,863).
        ((1_236.0, 100.0), 300),
        # 1,237 x 3.125 = 3,865.625: a canvas that truncated drew 3,865.
        ((1_237.0, 101.0), 300),
        ((10.0, 10.0), 72),
    ],
)
def test_the_web_png_check_returns_the_pixel_size_the_cli_png_has(size: tuple[float, float], dpi: int) -> None:
    # OV-2682: the Web sizes its PNG canvas by the size the check returns, so a Web PNG and a CLI PNG of the same
    # figure at the same DPI are the same size.
    pytest.importorskip("cairosvg")
    width, height = size
    reply = json.loads(check_web_export_limit_json(json.dumps(
        {"format": "png", "width": width, "height": height, "dpi": dpi, "linearWrap": None})))
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}"><rect width="{width}" height="{height}"/></svg>')
    png = convert_svg_with_cairosvg(svg, "png", dpi=dpi)
    cli_width, cli_height = struct.unpack(">II", png[16:24])
    assert reply == {"ok": True, "pixels": {"width": cli_width, "height": cli_height}}
    assert (cli_width, cli_height) == export_limits.png_pixel_size(width, height, dpi)


def test_a_pdf_check_names_no_pixel_size() -> None:
    reply = check_web_export_limit_json(json.dumps(
        {"format": "pdf", "width": 1_237.0, "height": 101.0, "dpi": None, "linearWrap": None}))
    assert json.loads(reply) == {"ok": True, "pixels": None}


# A Linear figure fitted to 180 mm wide whose height is converted to mm: its written root states 626.5757mm
# (export_root_size rounds to 4 decimals), 7,400.5004 px at 300 dpi, while the drawn 2,368.159945801535 px is
# 7,400.4998 px (OV-2683).
_SIZED_HEIGHT_PX = 2_368.159945801535


def _sized_drawing(tmp_path: Path) -> Drawing:
    from gbdraw.core.figure_length import RequestedFigureSize, parse_figure_length
    from gbdraw.render.document import bind_export_size

    width = parse_figure_length("180mm")
    drawing = Drawing(filename=str(tmp_path / "sized.svg"), size=(f"{width.px}px", f"{_SIZED_HEIGHT_PX}px"),
                      viewBox=f"0 0 {width.px} {_SIZED_HEIGHT_PX}")
    drawing.add(drawing.rect((10, 10), (100, 100), fill="#54bcf8"))
    bind_export_size(drawing, RequestedFigureSize(width, None))
    return drawing


def _ihdr(png: bytes) -> tuple[int, int]:
    return struct.unpack(">II", png[16:24])


def test_a_sized_png_has_the_pixel_size_of_the_drawn_px_on_every_python_path(tmp_path: Path) -> None:
    # OV-2683: the PNG is drawn from the drawing's px root, not from the root written in mm for the vector formats,
    # so its size is png_pixel_size of the drawn px (the size the check and the Web PNG take); the SVG keeps its mm.
    pytest.importorskip("cairosvg")
    drawing = _sized_drawing(tmp_path)
    width_px = float(drawing.attribs["width"].removesuffix("px"))
    expected = export_limits.png_pixel_size(width_px, _SIZED_HEIGHT_PX, None)
    assert expected[1] == 7_400
    assert _ihdr(render_to_bytes(drawing, "png")) == expected
    paths = save_figure_to(drawing, ["png"], output_prefix=str(tmp_path / "figure"))
    png = next(Path(path) for path in paths if path.endswith(".png"))
    assert _ihdr(png.read_bytes()) == expected
    assert 'height="626.5757mm"' in (tmp_path / "figure.svg").read_text(encoding="utf-8")
    (tmp_path / "deprecated_dir").mkdir()
    with pytest.warns(DeprecationWarning):
        render.export.save_figure(_sized_drawing(tmp_path / "deprecated_dir"), ["png"])
    assert _ihdr((tmp_path / "deprecated_dir" / "sized.png").read_bytes()) == expected


def test_vector_formats_convert_the_written_root_and_an_unsized_png_the_same_text(tmp_path: Path, monkeypatch) -> None:
    # The PDF, EPS and PS convert the written document (a requested size in its unit), as before; an unsized PNG
    # converts that same text, so its bytes do not change.
    from gbdraw.render.document import write_svg_document

    sized = _sized_drawing(tmp_path)
    unsized = _drawing(tmp_path, 1_237, 101)
    sources: list[tuple[str, bytes]] = []

    def recording(svg_source, fmt, *, write_to=None, **_options):
        sources.append((fmt, bytes(svg_source)))
        if write_to is not None:
            write_to.write(b"")
        return b""

    monkeypatch.setattr(render.export, "convert_svg_with_cairosvg", recording)
    monkeypatch.setattr(render.export, "get_cairosvg", lambda: object())
    monkeypatch.setattr(render.export, "require_cairosvg", lambda formats: object())
    for drawing, name in ((sized, "sized"), (unsized, "unsized")):
        sources.clear()
        written = write_svg_document(drawing).body.encode("utf-8")
        save_figure_to(drawing, ["png", "pdf", "eps", "ps"], output_prefix=str(tmp_path / name), overwrite=True)
        by_format = dict(sources)
        assert [by_format[fmt] for fmt in ("pdf", "eps", "ps")] == [written] * 3, name
        sources.clear()
        render_to_bytes(drawing, "pdf")
        assert sources == [("pdf", written)], name
        sources.clear()
        render_to_bytes(drawing, "png")
        if name == "unsized":
            assert by_format["png"] == written and sources == [("png", written)]
        else:
            assert b'height="626.5757mm"' in written
            assert b"626.5757mm" not in by_format["png"] and sources == [("png", by_format["png"])]


def test_an_unsized_png_is_the_png_of_the_written_document(tmp_path: Path) -> None:
    # The bytes an unsized PNG had before OV-2683: CairoSVG of the written document.
    pytest.importorskip("cairosvg")
    from gbdraw.render.document import write_svg_document

    built = plan_request(
        LinearDiagramRequest(records=(RecordInput(GenBankInputSource(str(LAMBDA)), record_key="r0"),))
    ).build()
    before = convert_svg_with_cairosvg(write_svg_document(built.drawing).body.encode("utf-8"), "png")
    assert render_to_bytes(built.drawing, "png") == before
    paths = save_figure_to(built.drawing, ["png"], output_prefix=str(tmp_path / "lambda"))
    assert next(Path(path) for path in paths if path.endswith(".png")).read_bytes() == before
