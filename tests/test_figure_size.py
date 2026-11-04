"""A requested figure size (issue #1427, PR-A): one parser, one config key, the composition planner's requested
canvas, the Linear width fit, its refusals, and export units and DPI."""

from __future__ import annotations

import re
import struct
import sys
import zlib
from dataclasses import asdict
from pathlib import Path

import pytest

import gbdraw
from gbdraw.api.config import apply_config_overrides
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.canvas import FIGURE_HEIGHT_PATH, FIGURE_WIDTH_PATH, figure_length_from_config
from gbdraw.config.toml import load_config_toml
from gbdraw.core.figure_length import FigureLength, argparse_figure_length, parse_figure_length
from gbdraw.exceptions import ValidationError
from gbdraw.layout.composition import CompositionItem, CompositionRequest, LegendPlacement, plan_composition
from gbdraw.layout.spatial import Aabb

INPUTS = Path(__file__).resolve().parent / "test_inputs"
MM_180_PX = 180 * 96 / 25.4


def _default_config() -> GbdrawConfig:
    return GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))


# --- one parser (D-02, D-08) ---


@pytest.mark.parametrize(
    ("text", "value", "unit", "canonical"),
    [
        ("180mm", 180.0, "mm", "180mm"),
        ("180 MM", 180.0, "mm", "180mm"),
        ("7.5in", 7.5, "in", "7.5in"),
        ("18cm", 18.0, "cm", "18cm"),
        (" 680px ", 680.0, "px", "680px"),
        (".5in", 0.5, "in", "0.5in"),
        ("180.50mm", 180.5, "mm", "180.5mm"),
    ],
)
def test_a_length_with_a_unit_parses_to_its_canonical_text(text, value, unit, canonical) -> None:
    length = parse_figure_length(text)
    assert length == FigureLength(value, unit)
    assert length is not None and length.text == canonical
    assert argparse_figure_length(text) == canonical


def test_lengths_convert_at_96_px_per_inch() -> None:
    assert FigureLength(180, "mm").px == pytest.approx(MM_180_PX)
    assert FigureLength(18, "cm").px == pytest.approx(MM_180_PX)
    assert FigureLength(1, "in").px == 96
    assert FigureLength(680, "px").px == 680


@pytest.mark.parametrize("text", ["auto", "AUTO", " Auto ", None])
def test_auto_is_no_requested_length(text) -> None:
    assert parse_figure_length(text) is None


@pytest.mark.parametrize("text", ["180", "0mm", "-5mm", "5pt", "none", "", "1e3mm", "nanmm", "infmm", "mm"])
def test_other_spellings_are_refused_with_the_units_and_an_example(text) -> None:
    with pytest.raises(ValidationError, match=r"mm, cm, in, or px\), for example 180mm"):
        parse_figure_length(text)


def test_the_config_leaf_stores_the_canonical_text() -> None:
    assert figure_length_from_config("180 MM", path=FIGURE_WIDTH_PATH) == "180mm"
    assert figure_length_from_config(None, path=FIGURE_WIDTH_PATH) == "auto"
    with pytest.raises(ValidationError) as caught:
        figure_length_from_config(180, path=FIGURE_HEIGHT_PATH)
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "configPath": FIGURE_HEIGHT_PATH}

    cfg = _default_config()
    assert (cfg.canvas.figure_width, cfg.canvas.figure_height) == ("auto", "auto")
    sized = apply_config_overrides(cfg, {FIGURE_WIDTH_PATH: "7 IN", FIGURE_HEIGHT_PATH: "12cm"})
    assert (sized.canvas.figure_width, sized.canvas.figure_height) == ("7in", "12cm")
    assert GbdrawConfig.from_dict(asdict(sized)).canvas.figure_width == "7in"
    assert sized.canvas.figure_size.width == FigureLength(7, "in")
    with pytest.raises(ValidationError):
        apply_config_overrides(cfg, {FIGURE_WIDTH_PATH: "180"})


# --- the planner's requested canvas (D-03) ---


def _request(**requested: float) -> CompositionRequest:
    return CompositionRequest(
        primary=CompositionItem("primary", Aabb(0.0, 0.0, 400.0, 200.0)),
        legend=CompositionItem("legend", Aabb(0.0, 0.0, 80.0, 60.0)),
        legend_placement=LegendPlacement.RIGHT,
        **requested,
    )


def test_no_requested_canvas_plans_as_before() -> None:
    plan = plan_composition(_request())
    assert plan.content_bounds == plan.canvas_bounds
    assert (plan.width, plan.height) == (400 + 24 + 80 + 32, 232)


def test_a_larger_requested_canvas_centres_the_content() -> None:
    plain = plan_composition(_request())
    plan = plan_composition(_request(requested_width_px=700.0, requested_height_px=300.0))
    assert (plan.width, plan.height) == (700.0, 300.0)
    dx, dy = (700.0 - plain.width) / 2, (300.0 - plain.height) / 2
    assert plan.content_bounds == plain.canvas_bounds.translated(dx, dy)
    for before, after in zip(plain.placements, plan.placements):
        assert after.final_bounds == before.final_bounds.translated(dx, dy)


def test_a_smaller_requested_canvas_keeps_the_content_size() -> None:
    plain = plan_composition(_request())
    plan = plan_composition(_request(requested_width_px=100.0, requested_height_px=50.0))
    assert plan.canvas_bounds == plain.canvas_bounds
    assert plan.placements == plain.placements


# --- the Linear width fit (D-01, D-09) ---


def _svg_size(svg: str) -> tuple[str, str, list[float]]:
    head = svg[svg.index("<svg") : svg.index(">", svg.index("<svg"))]
    width = re.search(r'\swidth="([^"]+)"', head)
    height = re.search(r'\sheight="([^"]+)"', head)
    view_box = re.search(r'viewBox="([^"]+)"', head)
    assert width and height and view_box
    return width[1], height[1], [float(value) for value in view_box[1].split()]


def _px(text: str) -> float:
    assert text.endswith("px"), text
    return float(text[:-2])


def _font_sizes(svg: str) -> set[str]:
    return set(re.findall(r'font-size="([^"]+)"', svg))


def _mito() -> list:
    return gbdraw.read_genbank(INPUTS / "HmmtDNA.gbk")


def test_a_requested_width_gives_that_width_with_the_same_font_sizes() -> None:
    auto = gbdraw.draw_linear(_mito()).to_svg()
    sized_diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm"))
    sized = sized_diagram.to_svg()
    width, _height, view_box = _svg_size(sized)
    assert width == "180mm"
    assert view_box[2] == pytest.approx(MM_180_PX, abs=0.01)
    # The drawing keeps px; only the written document carries the unit.
    assert _px(sized_diagram._drawing.attribs["width"]) == pytest.approx(MM_180_PX, abs=0.01)
    assert _font_sizes(sized) == _font_sizes(auto)
    assert _px(_svg_size(auto)[0]) > MM_180_PX + 100


@pytest.mark.parametrize("legend", ["right", "top", "left", "bottom", "none"])
def test_every_legend_place_fits_the_requested_width(legend) -> None:
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="250mm", legend=legend))
    assert _svg_size(diagram.to_svg())[2][2] == pytest.approx(250 * 96 / 25.4, abs=0.01)


def test_a_multi_record_figure_fits_the_requested_width() -> None:
    records = [*_mito(), *gbdraw.read_genbank(INPUTS / "BGC0000712.gbk")]
    diagram = gbdraw.draw_linear(records, options=gbdraw.LinearOptions(figure_width="7in"))
    width, _height, view_box = _svg_size(diagram.to_svg())
    assert width == "7in"
    assert view_box[2] == pytest.approx(7 * 96, abs=0.01)


def test_a_wrapped_record_draws_its_rows_at_the_fitted_width() -> None:
    records = gbdraw.read_genbank(INPUTS / "NC_005856.1.gb")
    options = gbdraw.LinearOptions(
        figure_width="180mm", config_overrides={"canvas.linear.wrap_length": "auto"}
    )
    diagram = gbdraw.draw_linear(records, options=options)
    assert _svg_size(diagram.to_svg())[2][2] == pytest.approx(MM_180_PX, abs=0.01)


def test_a_requested_height_centres_the_figure_and_keeps_the_width() -> None:
    auto_width, auto_height, _ = _svg_size(gbdraw.draw_linear(_mito()).to_svg())
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_height="10in"))
    width, height, view_box = _svg_size(diagram.to_svg())
    assert height == "10in"
    assert view_box[3] == pytest.approx(960, abs=0.01)
    # The width is written in the height's unit, at the auto figure's width.
    assert width.endswith("in")
    assert view_box[2] == pytest.approx(_px(auto_width), abs=0.01)
    assert _px(auto_height) < 960


# --- refusals (D-03) ---


def test_a_width_narrower_than_the_fixed_parts_names_them_and_the_legend_remedy() -> None:
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="40mm"))
    message = str(caught.value)
    assert message.startswith("figure_width 40mm is narrower than the parts of the figure that do not shrink: ")
    assert "record names and labels past the sequence ends" in message
    assert "legend" in message and "padding" in message
    assert "place the legend at the top or bottom, or leave figure_width auto." in message
    diagnostic = caught.value.diagnostic
    assert diagnostic["code"] == "INPUT_INVALID" and diagnostic["configPath"] == FIGURE_WIDTH_PATH
    minimum_px = diagnostic["minimumPx"]
    assert isinstance(minimum_px, int) and minimum_px > 40 * 96 / 25.4
    # The same width fits once the legend no longer sits beside the records.
    at_minimum = f"{minimum_px}px"
    with pytest.raises(ValidationError, match="narrower"):
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width=at_minimum))
    top = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width=at_minimum, legend="top"))
    assert _svg_size(top.to_svg())[2][2] == pytest.approx(minimum_px, abs=0.01)


def test_a_narrow_width_without_a_side_legend_offers_no_legend_remedy() -> None:
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="10mm", legend="top"))
    assert "place the legend" not in str(caught.value)


def test_a_height_shorter_than_the_figure_is_refused_with_its_parts() -> None:
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(
            _mito(), options=gbdraw.LinearOptions(figure_width="180mm", figure_height="10mm")
        )
    message = str(caught.value)
    assert message.startswith("figure_height 10mm is shorter than the figure at 180 mm wide: ")
    assert "records" in message and "padding" in message
    assert message.endswith("or leave figure_height auto.")
    assert caught.value.diagnostic["configPath"] == FIGURE_HEIGHT_PATH


# --- Circular waits for its fit (D-07) ---


@pytest.mark.parametrize(
    "options",
    [
        gbdraw.CircularOptions(figure_width="180mm"),
        gbdraw.CircularOptions(config_overrides={FIGURE_HEIGHT_PATH: "120mm"}),
    ],
)
def test_circular_refuses_a_requested_size(options) -> None:
    with pytest.raises(ValidationError, match="not available for Circular maps yet") as caught:
        gbdraw.draw_circular(_mito(), options=options)
    assert caught.value.diagnostic["code"] == "INPUT_INVALID"


def test_circular_multi_record_refuses_a_requested_size() -> None:
    records = [*_mito(), *gbdraw.read_genbank(INPUTS / "BGC0000712.gbk")]
    with pytest.raises(ValidationError, match="not available for Circular maps yet"):
        gbdraw.draw_circular(
            records, options=gbdraw.CircularOptions(figure_width="180mm"), layout=gbdraw.CircularLayout()
        )


# --- export units and DPI (D-10) ---


def _png_size(payload: bytes) -> tuple[int, int]:
    assert payload[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", payload[16:24])


def test_saved_svg_carries_the_requested_unit(tmp_path) -> None:
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm"))
    path = diagram.save(tmp_path / "sized.svg")
    width, height, view_box = _svg_size(path.read_text(encoding="utf-8"))
    assert width == "180mm"
    assert height.endswith("mm")
    assert float(height[:-2]) == pytest.approx(view_box[3] * 25.4 / 96, abs=0.001)


def test_png_dpi_scales_the_pixels() -> None:
    pytest.importorskip("cairosvg")
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm"))
    _w, _h, view_box = _svg_size(diagram.to_svg())
    assert _png_size(diagram.to_bytes("png")) == (round(view_box[2]), round(view_box[3]))
    assert _png_size(diagram.to_bytes("png", dpi=300)) == (
        round(view_box[2] * 300 / 96),
        round(view_box[3] * 300 / 96),
    )


def test_pdf_page_has_the_requested_width() -> None:
    pytest.importorskip("cairosvg")
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm"))
    payload = diagram.to_bytes("pdf")
    # cairo writes the page object into a compressed object stream.
    texts = [payload]
    for start in re.finditer(rb"stream\r?\n", payload):
        try:
            texts.append(zlib.decompressobj().decompress(payload[start.end():]))
        except zlib.error:
            continue
    boxes = [box for text in texts for box in re.findall(rb"/?MediaBox \[ ?([0-9.]+) ([0-9.]+) ([0-9.]+)", text)]
    assert boxes
    assert float(boxes[0][2]) == pytest.approx(180 / 25.4 * 72, abs=0.1)


# --- CLI ---


def _run_cli(tmp_path: Path, monkeypatch, mode: str, *args: str) -> Path:
    from gbdraw import cli

    output = tmp_path / "cli-out"
    monkeypatch.setattr(sys, "argv", ["gbdraw", mode, *args, "-o", str(output)])
    cli.main()
    return output


def test_cli_figure_width_and_dpi(tmp_path, monkeypatch) -> None:
    pytest.importorskip("cairosvg")
    output = _run_cli(
        tmp_path, monkeypatch, "linear", "--gbk", str(INPUTS / "HmmtDNA.gbk"),
        "--figure_width", "180 mm", "--dpi", "192", "-f", "svg,png",
    )
    width, _height, view_box = _svg_size(output.with_suffix(".svg").read_text(encoding="utf-8"))
    assert width == "180mm"
    assert _png_size(output.with_suffix(".png").read_bytes())[0] == round(view_box[2] * 2)


@pytest.mark.parametrize("args", [["--figure_width", "180"], ["--dpi", "0"]])
def test_cli_refuses_bad_values(tmp_path, monkeypatch, args) -> None:
    with pytest.raises(SystemExit) as caught:
        _run_cli(tmp_path, monkeypatch, "linear", "--gbk", str(INPUTS / "HmmtDNA.gbk"), *args)
    assert caught.value.code == 2


def test_cli_warns_when_dpi_has_no_png(tmp_path, monkeypatch, caplog) -> None:
    _run_cli(tmp_path, monkeypatch, "linear", "--gbk", str(INPUTS / "HmmtDNA.gbk"), "--dpi", "300", "-f", "svg")
    assert "--dpi applies to PNG output only" in caplog.text


def test_cli_circular_refuses_figure_width(tmp_path, monkeypatch, capsys) -> None:
    with pytest.raises(SystemExit) as caught:
        _run_cli(
            tmp_path, monkeypatch, "circular", "--gbk", str(INPUTS / "HmmtDNA.gbk"), "--figure_width", "180mm"
        )
    assert caught.value.code == 1
    assert "not available for Circular maps yet" in capsys.readouterr().err
    assert not (tmp_path / "cli-out.svg").exists()
