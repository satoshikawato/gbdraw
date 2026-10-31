"""PNG and PDF size limits: a wrapped figure stops with a message naming what would fit; SVG never stops."""

from __future__ import annotations

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
from gbdraw.render.export_limits import bind_export_wrap, check_export_limit

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
    message = _refusal("png", SORANGIUM_SIZE, SORANGIUM)
    assert "32,767 px" in message
    # 131 rows of the 261 drawn: about 23,072 px.
    assert "--wrap_length 100kb" in message and "23,072 px" in message
    assert "estimated from the drawn height" in message
    assert "SVG" in message
    assert "DPI" not in message


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
    message = _refusal("png", (2_000.0, 80_000.0), many_records)
    assert "No longer row length would fit" in message
    assert "--wrap_length" not in message
    assert "SVG" in message


def test_no_row_length_fits_when_the_width_is_over_the_limit() -> None:
    message = _refusal("png", (40_000.0, 100.0), LinearWrapSummary(50_000, 50_000, (120_000,)))
    assert "No longer row length would fit" in message


def test_the_last_row_length_that_fits_may_draw_one_row() -> None:
    assert "--wrap_length 100kb" in _refusal("png", (2_000.0, 40_000.0), LinearWrapSummary(50_000, 50_000, (60_000,)))


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
    wrap = LinearWrapSummary(50_000, 50_000, (120_000,))
    check_export_limit(fmt, *fits, wrap=wrap)
    _refusal(fmt, over, wrap)


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
        save_figure_to(drawing, ["svg", "png"], output_prefix=str(tmp_path / "figure"), overwrite=True)
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
