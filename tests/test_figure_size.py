"""A requested figure size (issue #1427, PR-A): one parser, one config key, the composition planner's requested
canvas, the Linear width fit, its refusals, and export units and DPI."""

from __future__ import annotations

import json
import math
import re
import struct
import sys
import warnings
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
from gbdraw.layout.figure_size import MeasuredBase, fit_base_length
from gbdraw.layout.composition import CompositionItem, CompositionRequest, LegendPlacement, plan_composition
from gbdraw.layout.spatial import Aabb
from gbdraw.render.composition import applied_composition_plan
from tests.utils.png_chunks import png_chunks, png_phys
from tests.utils.web_render_result import web_result_svg

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


# --- the base length resolver (D-09) ---


def _fit(size_of_base, target: float, *, first_base: int | None = None, slope: float = 1.0, left_out=None):
    """Fit a measured size to ``target`` as a mode does: the measure returns the excess over the target, and whether
    the pass left content out (``left_out(base)``, none by default)."""

    bases: list[int] = []

    def measure(base: int) -> MeasuredBase[int] | None:
        bases.append(base)
        # The walk measures at most 8 costly and 24 cheap bases, the first, and one middle base.
        assert len(bases) <= 40, f"the walk measured more than 40 bases: {bases[:12]} ..."
        size = size_of_base(base)
        if size is None:
            return None
        return MeasuredBase(size - target, base, (size,), left_out=left_out is not None and left_out(base))

    first = max(1, math.floor(target)) if first_base is None else first_base
    return fit_base_length(measure, first_base=first, slope=slope), bases


def test_a_fixed_part_that_stays_put_fits_in_two_passes() -> None:
    fit, bases = _fit(lambda base: base + 100.0, 500.0)
    assert fit.fits and -3.0 <= fit.excess_px <= 0.0 and fit.built == fit.base
    assert len(bases) == 2


def test_a_fixed_part_that_grows_as_the_base_shrinks_still_fits() -> None:
    fit, bases = _fit(lambda base: base + 100.0 + 40000.0 / base, 680.0)
    assert fit.fits and -3.0 <= fit.excess_px <= 0.0
    assert len(bases) <= 8


def test_a_size_flat_below_a_base_fits_on_the_flat_part() -> None:
    fit, _bases = _fit(lambda base: max(718.0, base + 150.0), 718.0)
    assert fit.fits and fit.excess_px == 0.0


def test_a_target_that_the_shortest_base_passes_is_refused_at_base_one() -> None:
    """Base 1 is over, so the middle of base 1 and the next pass is measured once before refusing (N1)."""

    fit, bases = _fit(lambda base: base + 600.0, 500.0)
    assert not fit.fits and fit.base == 1 and fit.built == 1
    assert bases == [500, 1, 250]


def test_the_first_base_and_slope_come_from_the_mode() -> None:
    """A size that grows 2 px per px of base (a Circular map) fits in two passes from the configured base."""

    fit, bases = _fit(lambda base: 2.0 * base + 180.0, 680.3, first_base=390, slope=2.0)
    assert bases[0] == 390 and len(bases) == 2
    assert fit.fits and -3.0 <= fit.excess_px <= 0.0


def test_a_bracket_is_searched_until_a_pass_lands_close_under() -> None:
    """N2: a secant rounded down onto the pass under steps one base up instead of halving the bracket, so a steep
    size lands within 3 px instead of returning a pass far under after 8 passes."""

    fit, bases = _fit(lambda base: 3.0 * base + 100.0, 680.3)
    assert fit.fits and -3.0 <= fit.excess_px <= 0.0 and fit.base == 193
    assert len(bases) <= 4


def test_a_middle_base_that_fits_is_found_before_refusing() -> None:
    """N1: base 1 is over, but a base between it and the next pass fits (a measure that is not monotonic). Outside
    the dip the size still shrinks with the base: two passes of the same size would end the walk before base 1."""

    fit, bases = _fit(lambda base: 679.0 if 30 <= base <= 40 else 700.0 + 0.1 * base, 680.0, first_base=100)
    assert fit.fits and 30 <= fit.base <= 40
    assert 1 in bases


def test_bases_too_small_to_draw_are_not_searched_below_the_floor() -> None:
    """Bases under 40 cannot be drawn; the refusal keeps the smallest drawn pass and no step goes below the floor
    twice. The size still shrinks a little below radius 100 (a multi-record grid: 0.1 px per px), so the walk goes
    down to the floor."""

    fit, bases = _fit(
        lambda base: None if base < 40 else max(2.0 * base, 200.0 + 0.1 * base) + 580.0, 680.0, first_base=390,
        slope=2.0,
    )
    assert not fit.fits and fit.base == 40
    undrawable = [base for base in bases if base < 40]
    assert undrawable and len(undrawable) <= 12
    assert len(bases) - len(undrawable) <= 8
    assert len(set(bases)) == len(bases)


def test_a_base_too_small_to_draw_above_a_pass_under_bounds_the_search() -> None:
    """A base that cannot be drawn above the widest pass under is a bound like a pass over: the search halves the
    way to it instead of measuring it again, and ends."""

    def size(base: int) -> float | None:
        if base == 299:
            return 630.0
        return None if 310 <= base <= 330 else 2.0 * base + 80.0

    fit, bases = _fit(size, 680.0, first_base=390, slope=2.0)
    assert fit.fits and -3.0 <= fit.excess_px <= 0.0 and fit.base == 300
    assert len(set(bases)) == len(bases)


def test_a_refusal_finds_the_smallest_drawable_base_with_few_drawn_passes() -> None:
    """Every base is over the box and bases under 53 cannot be drawn (they fail early). Below every pass the walk tries
    a quarter of the way up from the largest base too small to draw, so it reaches the smallest drawable base, which
    the refusal states, after four drawn passes; halving the way drew eight, 390, 195, 98, 73, 61, 55, 53 and 54, as
    HmmtDNA with labels out at 85 mm does, where a drawn pass takes 4 to 20 s. The middle check is skipped: no dip
    four times as steep as the steepest slope reaches the box between 53 and the next pass."""

    fit, bases = _fit(lambda base: None if base < 53 else 2.0 * base + 900.0, 680.0, first_base=390, slope=2.0)
    assert not fit.fits and fit.base == 53
    drawn = [base for base in bases if base >= 53]
    assert drawn == [390, 98, 56, 53]
    assert len(set(bases)) == len(bases) and len(bases) - len(drawn) <= 24


def test_the_middle_check_runs_only_when_a_dip_from_both_passes_reaches_the_box() -> None:
    """The size at a base between two passes over the box is at least what a line four times as steep as the
    steepest assumed slope (16 px per px of radius for a Circular map) gives from each pass, so the middle can fit
    only if both lines reach 0 px there."""

    from gbdraw.layout.figure_size import BaseLengthFit, _dip_reaches_box

    lowest = BaseLengthFit(53, 5.0, 53, fits=False)
    # From 53, 5 - 16 = -11 px; from 56 at +100 px, 100 - 2 * 16 = 68 px: the middle cannot fit.
    assert not _dip_reaches_box(lowest, BaseLengthFit(56, 100.0, 56, fits=False), 54, 2.0)
    # From 56 at +20 px, 20 - 2 * 16 = -12 px: both lines reach the box, so the middle is measured.
    assert _dip_reaches_box(lowest, BaseLengthFit(56, 20.0, 56, fits=False), 54, 2.0)


def test_bases_too_small_to_draw_end_the_search_after_their_bound() -> None:
    """Every base below the first cannot be drawn: the refusal keeps the first pass after at most 24 such passes."""

    fit, bases = _fit(lambda base: None if base < 390 else 2.0 * base + 180.0, 680.0, first_base=390, slope=2.0)
    assert not fit.fits and fit.base == 390
    assert 1 < len(bases) <= 25 and len(set(bases)) == len(bases)


# The MjeNMV map against a 120 mm box, measured per radius (CPython, 2026-11-06): the centre definition keeps the
# tracks out, so from radius 220 down to 49, the smallest radius whose tracks can be drawn, the content is 191 px over
# the box (the same to the last digits but one); above 220 it grows 2 px per px of radius.
def _mjenmv_without_labels(base: int) -> float | None:
    return None if base < 49 else max(2.0 * base - 249.0, 191.0)


def test_a_refusal_ends_once_the_two_smallest_passes_have_the_same_size() -> None:
    """When every pass is over the box and the two smallest have the same size, the content no longer shrinks with
    the base: no smaller base fits, so the walk ends there and the refusal states the figure at its smallest, the same
    size as at the smallest drawable radius. Before, the walk halved its way down to radius 49 in 11 passes."""

    fit, bases = _fit(_mjenmv_without_labels, 0.0, first_base=390, slope=2.0)
    assert not fit.fits and fit.excess_px == 191.0 and fit.base > 49
    assert bases == [390, 123, 1, 31, 54]


def test_sizes_equal_up_to_float_rounding_are_flat() -> None:
    """Parts that do not shrink are added in another order at another radius, so the same size can differ in its
    last digits (MjeNMV with the legend on the left: 191.0002571393764 and 191.00025713937652)."""

    fit, bases = _fit(
        lambda base: None if base < 49 else max(2.0 * base - 249.0, 191.0 + (1e-13 if base % 2 else 0.0)),
        0.0,
        first_base=390,
        slope=2.0,
    )
    assert not fit.fits and len(bases) == 5


def test_a_refusal_goes_on_while_a_part_shrinks_under_the_same_excess() -> None:
    """A legend beside the map and taller than it keeps the excess over a requested height the same while the map
    still shrinks. A refusal names the map at its smallest, so the walk goes on while the parts differ and keeps the
    smallest drawable radius."""

    def measure(base: int) -> MeasuredBase[int] | None:
        if base < 49:
            return None
        map_px, legend_px = 2.0 * base + 40.0, 400.0
        height = max(map_px, legend_px)
        return MeasuredBase(height - 189.0, base, (map_px, height - map_px))

    fit = fit_base_length(measure, first_base=390, slope=2.0)
    assert not fit.fits and fit.base == 49


def test_a_flat_stretch_that_shrinks_again_below_is_walked_while_it_differs() -> None:
    """Only the same size at the two smallest passes ends the walk: a stretch nearly flat (labels past the ends of a
    Linear record) that shrinks again further down is still walked, and a base there that fits is found."""

    fit, _bases = _fit(
        lambda base: 700.0 + 0.01 * base if base >= 60 else 700.6 - (60 - base), 680.0, first_base=500
    )
    assert fit.fits


def test_a_size_that_jumps_returns_the_widest_pass_under() -> None:
    """Labels not placed (CAPACITY-STOP) make a pass much narrower than its neighbour with labels: the measure jumps
    over the target, and the fit keeps the widest pass under it."""

    fit, bases = _fit(lambda base: base + 520.0 if base < 150 else base + 840.0, 680.0, first_base=300, slope=2.0)
    assert fit.fits and fit.base < 150 and fit.excess_px <= 0.0
    assert len(bases) <= 8


# The MjeNMV map with labels on both sides against a 180 mm box, measured per radius (CPython, 2026-11-06): from
# radius 112 the labels are placed and the content is over the box, the slope flattening near 112; below it the
# inner labels need more than 1.5 times the room, so the map's labels are not placed (CAPACITY-STOP) and the content
# stays 91.3 px under. A pass with labels near 112 took over a minute.
_MJENMV_EXCESS_WITH_LABELS = (
    (112, 30.8), (115, 34.2), (120, 40.0), (130, 51.7), (140, 62.9), (150, 85.8), (160, 108.3), (180, 155.6),
    (200, 202.7), (220, 249.5), (250, 320.0), (300, 436.7), (390, 646.6),
)


def _mjenmv_walk(shift_px: float):
    """Fit the measured MjeNMV radii (above) to a box ``shift_px`` narrower than 180 mm; returns the fit and the radii
    measured with labels."""

    with_labels: list[int] = []

    def measure(base: int) -> MeasuredBase[int]:
        if base < 112:
            return MeasuredBase(-91.3 + shift_px, base, (-91.3,), left_out=True)
        with_labels.append(base)
        points = _MJENMV_EXCESS_WITH_LABELS
        (b0, e0), (b1, e1) = next(
            ((low, high) for low, high in zip(points, points[1:]) if base <= high[0]), points[-2:]
        )
        excess = e0 + (e1 - e0) * (base - b0) / (b1 - b0)
        return MeasuredBase(excess + shift_px, base, (excess,))

    return fit_base_length(measure, first_base=390, slope=2.0), with_labels


def test_passes_whose_labels_were_not_placed_do_not_steer_the_walk() -> None:
    """A pass under the box only because its labels were left out (CAPACITY-STOP) says nothing about the size with
    labels: the step comes from the passes with labels, and the walk ends once even a line twice as steep as the
    first step reaches the box only where the labels are left out. It keeps the pass it kept before (the first one
    under) after three passes with labels instead of five."""

    fit, with_labels = _mjenmv_walk(0.0)
    assert fit.fits and fit.base == 65
    assert with_labels == [390, 138, 117]


def test_a_radius_with_labels_just_above_the_capacity_stop_is_still_found() -> None:
    """In a box 10 mm wider, radii from 112 with labels fit: the walk finds one within 3 px under, with four passes
    with labels instead of six."""

    fit, with_labels = _mjenmv_walk(-37.8)
    assert fit.fits and fit.base >= 112 and -3.0 <= fit.excess_px <= 0.0
    assert len(with_labels) <= 4


def _stop_walk(stop: int, without_labels, with_labels, *, first_base: int):
    """Fit a size that leaves the labels out below radius ``stop``; returns the fit and every radius measured."""

    def size(base: int) -> float:
        return without_labels(base) if base < stop else with_labels(base)

    return _fit(size, 0.0, first_base=first_base, slope=2.0, left_out=lambda base: base < stop)


def test_a_pass_without_labels_over_the_box_is_not_measured_again() -> None:
    """Without labels the map reaches the box at radius 190, below the stop at 200, so passes without labels can be
    over the box too. The walk keeps the widest of them under the box and measures no radius twice."""

    fit, bases = _stop_walk(
        200,
        lambda base: 2.0 * (base - 190),
        lambda base: 20 + 0.8 * (base - 200) if base <= 260 else 68 + 2.0 * (base - 260),
        first_base=300,
    )
    assert fit.fits and fit.base == 189
    assert len(set(bases)) == len(bases)


def test_a_radius_with_labels_less_than_2_px_under_the_box_is_found() -> None:
    """With labels the size reaches the box at radius 102, 0.5 px per px above the stop at 100: a step that aims
    2 px under lands where the labels are left out, but the walk still tries where the line reaches the box."""

    fit, bases = _stop_walk(
        100,
        lambda base: -50.0,
        lambda base: -1.0 + 0.5 * (base - 100) if base <= 140 else 19 + 2.0 * (base - 140),
        first_base=300,
    )
    assert fit.fits and fit.base >= 100 and -3.0 <= fit.excess_px <= 0.0


def test_a_size_that_shrinks_faster_near_the_stop_still_fits_with_labels() -> None:
    """With labels the size shrinks 1.2 px per px of radius above 200 but 2.5 near the stop at 150, so the line
    through the passes above reaches the box below the stop while radii 150 to 153 fit. A line twice as steep as the
    first step still reaches the box above the stop, so the walk goes on and finds one of them."""

    fit, _bases = _stop_walk(
        150,
        lambda base: -60 + 2.0 * (base - 150),
        lambda base: -8 + 2.5 * (base - 150) if base <= 200 else 117 + 1.2 * (base - 200),
        first_base=390,
    )
    assert fit.fits and fit.base >= 150 and -3.0 <= fit.excess_px <= 0.0


def test_only_the_kept_pass_warns_once() -> None:
    """D-16: warnings raised in a pass that is not kept do not reach the caller; the kept pass's arrive once."""

    def measure(base: int) -> MeasuredBase[int]:
        warnings.warn(f"pass {base}", UserWarning)
        warnings.warn(f"pass {base} again", UserWarning)
        return MeasuredBase(base + 100.0 - 500.0, base, (base + 100.0,))

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = fit_base_length(measure, first_base=500, slope=1.0)
    assert fit.base == 398
    assert [str(item.message) for item in caught] == ["pass 398", "pass 398 again"]


def test_a_refused_fit_emits_no_warning() -> None:
    def measure(base: int) -> MeasuredBase[int]:
        warnings.warn(f"pass {base}", UserWarning)
        return MeasuredBase(base + 600.0 - 500.0, base, (base + 600.0,))

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = fit_base_length(measure, first_base=500, slope=1.0)
    assert not fit.fits and caught == []


def test_the_width_fit_measures_the_figure_drawn_with_region_label_rows() -> None:
    """The figure kept is the one with the feature label rows a region label asked for (RGN-YIELD-ROW), so the
    fit measures that one, not the first assembly."""

    from svgwrite import Drawing
    from svgwrite.container import Group

    from gbdraw.api.diagram import _checked_linear_build, _linear_label_fit
    from gbdraw.core.figure_length import RequestedFigureSize
    from gbdraw.labels.linear import (
        LINEAR_LABEL_FIT_ATTR, FeatureLabelRow, FeatureLabelRows, RecordRegionLabels, linear_label_fit,
    )
    from gbdraw.render.composition import apply_composition_plan

    rows = FeatureLabelRows((FeatureLabelRow(0, 3, 4),), overlaps=1, region_labels=((0, 0),))
    asking = linear_label_fit([("rec", [])], 8.0, [RecordRegionLabels((10.0, 10.0), 2, conflicted=(0, 1))], rows)
    cleared = linear_label_fit([("rec", [])], 8.0, [RecordRegionLabels((10.0, 10.0), 1, conflicted=(1,))])

    def assemble(overrides, feature_label_rows):
        base = float(overrides["canvas.linear.width"])
        width = base + (100.0 if feature_label_rows is None else 150.0)
        plan = plan_composition(CompositionRequest(primary=CompositionItem("primary", Aabb(0.0, 0.0, width, 50.0))))
        drawing = Drawing(debug=False)
        apply_composition_plan(drawing, plan, primary_targets=[Group()])
        setattr(drawing, LINEAR_LABEL_FIT_ATTR, asking if feature_label_rows is None else cleared)
        return drawing

    size = RequestedFigureSize(FigureLength(500.0, "px"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawing = _checked_linear_build(assemble, features_separated=False, figure_size=size)
    assert _linear_label_fit(drawing) is cleared
    assert 497.0 <= applied_composition_plan(drawing).content_bounds.width <= 500.0


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


LABELS_ALL = {"labels.linear.scope": "all"}
TITLE = gbdraw.TitleOptions(
    text="Homo sapiens mitochondrion, complete genome (NC_012920.1)", position="top", font_size=24
)


def _count_assemblies(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """The sequence width of every Linear assembly the build makes."""

    import gbdraw.api.diagram as diagram_module

    widths: list[int] = []
    original = diagram_module.assemble_linear_diagram_from_records

    def counting(*args, **kwargs):
        widths.append(args[1].cfg.canvas.linear.width)
        return original(*args, **kwargs)

    monkeypatch.setattr(diagram_module, "assemble_linear_diagram_from_records", counting)
    return widths


def _content_width(diagram: gbdraw.Diagram) -> float:
    return applied_composition_plan(diagram._drawing).content_bounds.width


def _content_width_at_base_one(**options) -> float:
    overrides = {**options.pop("config_overrides", {}), "canvas.linear.width": 1}
    return _content_width(
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(config_overrides=overrides, **options))
    )


def test_a_figure_without_labels_fits_in_two_assemblies(monkeypatch) -> None:
    widths = _count_assemblies(monkeypatch)
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm"))
    assert len(widths) == 2
    assert MM_180_PX - 3 <= _content_width(diagram) <= MM_180_PX


def test_labels_that_widen_as_the_sequence_narrows_still_fit(monkeypatch) -> None:
    """The labels past the sequence ends grow as the sequence width shrinks, so a slope-1 step lands over."""

    widths = _count_assemblies(monkeypatch)
    diagram = gbdraw.draw_linear(
        _mito(), options=gbdraw.LinearOptions(figure_width="180mm", config_overrides=LABELS_ALL)
    )
    assert _svg_size(diagram.to_svg())[2][2] == pytest.approx(MM_180_PX, abs=0.01)
    assert MM_180_PX - 3 <= _content_width(diagram) <= MM_180_PX
    assert len(widths) <= 3


def test_too_narrow_with_labels_states_the_width_at_the_shortest_sequence() -> None:
    minimum = _content_width_at_base_one(config_overrides=LABELS_ALL)
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="40mm", config_overrides=LABELS_ALL))
    assert caught.value.diagnostic["minimumPx"] == math.ceil(minimum)
    assert "record names and labels past the sequence ends" in str(caught.value)
    assert "place the legend at the top or bottom" in str(caught.value)


@pytest.mark.parametrize("width", ["85mm", "120mm"])
def test_a_title_wider_than_the_records_sets_the_minimum_width(width) -> None:
    minimum = _content_width_at_base_one(title=TITLE, legend="top")
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(title=TITLE, legend="top", figure_width=width))
    message = str(caught.value)
    assert ": title " in message and "padding" in message
    assert "record names" not in message and "place the legend" not in message
    assert caught.value.diagnostic["minimumPx"] == math.ceil(minimum)
    at_minimum = f"{math.ceil(minimum)}px"
    diagram = gbdraw.draw_linear(
        _mito(), options=gbdraw.LinearOptions(title=TITLE, legend="top", figure_width=at_minimum)
    )
    assert _svg_size(diagram.to_svg())[2][2] == pytest.approx(math.ceil(minimum), abs=0.01)


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
    # The stated minimum is drawn; one px less is refused.
    assert minimum_px == math.ceil(_content_width_at_base_one())
    at_minimum = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width=f"{minimum_px}px"))
    assert _svg_size(at_minimum.to_svg())[2][2] == pytest.approx(minimum_px, abs=0.01)
    with pytest.raises(ValidationError, match="narrower"):
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width=f"{minimum_px - 1}px"))


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


# --- the Circular radius fit (D-01, D-03, D-12, D-13) ---


def _circular_plan(diagram: gbdraw.Diagram):
    return applied_composition_plan(diagram._drawing)


def _draw_mito_circular(**options) -> gbdraw.Diagram:
    return gbdraw.draw_circular(_mito(), options=gbdraw.CircularOptions(**options))


def test_a_circular_width_fits_the_map_and_keeps_the_fonts() -> None:
    auto = _draw_mito_circular()
    sized = _draw_mito_circular(figure_width="180mm")
    plan = _circular_plan(sized)
    assert plan.width == pytest.approx(MM_180_PX, abs=0.01)
    assert MM_180_PX - 3 <= plan.content_bounds.width <= MM_180_PX
    # A side left auto is the content's.
    assert plan.height == pytest.approx(plan.content_bounds.height, abs=0.01)
    assert plan.height < _circular_plan(auto).height
    assert _font_sizes(sized.to_svg()) == _font_sizes(auto.to_svg())
    width, height, _view_box = _svg_size(sized.to_svg())
    assert width == "180mm" and height.endswith("mm")


def test_a_circular_height_fits_the_map() -> None:
    plan = _circular_plan(_draw_mito_circular(figure_height="120mm"))
    target = 120 * 96 / 25.4
    assert plan.height == pytest.approx(target, abs=0.01)
    assert target - 3 <= plan.content_bounds.height <= target
    assert plan.width == pytest.approx(plan.content_bounds.width, abs=0.01)


def test_a_circular_width_and_height_take_the_largest_map_that_fits_both() -> None:
    plan = _circular_plan(_draw_mito_circular(figure_width="180mm", figure_height="200mm"))
    height = 200 * 96 / 25.4
    assert (plan.width, plan.height) == (pytest.approx(MM_180_PX, abs=0.01), pytest.approx(height, abs=0.01))
    # The width binds; the height has room left, centred.
    assert MM_180_PX - 3 <= plan.content_bounds.width <= MM_180_PX
    assert plan.content_bounds.height <= height
    assert plan.content_bounds.min_y == pytest.approx(height - plan.content_bounds.max_y, abs=0.01)


def test_a_circular_size_larger_than_auto_grows_the_map(monkeypatch) -> None:
    import gbdraw.api.diagram as diagram_module

    radii: list[float] = []
    original = diagram_module.assemble_circular_diagram_from_record

    def counting(*args, **kwargs):
        radii.append(float(args[1].cfg.canvas.circular.radius))
        return original(*args, **kwargs)

    monkeypatch.setattr(diagram_module, "assemble_circular_diagram_from_record", counting)
    auto = _draw_mito_circular()
    radii.clear()
    sized = _draw_mito_circular(figure_width="400mm")
    plan = _circular_plan(sized)
    target = 400 * 96 / 25.4
    assert target - 3 <= plan.content_bounds.width <= target
    assert radii[0] == pytest.approx(390.0) and max(radii) > 390.0
    assert _font_sizes(sized.to_svg()) == _font_sizes(auto.to_svg())


def test_label_fit_rows_redraw_the_circular_map_at_the_fitted_radius(monkeypatch) -> None:
    """MjeNMV at 180 mm: its centre definition reaches a track, so label fit rows are verified, and each redraws the
    map at the radius the fit kept, not the configured one."""

    import gbdraw.api.diagram as diagram_module
    from gbdraw.labels.circular_capacity import label_capacity_checked

    assemblies: list[tuple[float, bool]] = []
    original = diagram_module.assemble_circular_diagram_from_record

    def recording(*args, **kwargs):
        assemblies.append((float(args[1].cfg.canvas.circular.radius), label_capacity_checked()))
        return original(*args, **kwargs)

    monkeypatch.setattr(diagram_module, "assemble_circular_diagram_from_record", recording)
    records = gbdraw.read_genbank(Path(__file__).resolve().parents[1] / "examples" / "MjeNMV.gb")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gbdraw.draw_circular(records, options=gbdraw.CircularOptions(figure_width="180mm"))
    fitted = [radius for radius, checked in assemblies if checked]
    rows = [radius for radius, checked in assemblies if not checked]
    assert rows and fitted[-1] < 390.0
    assert set(rows) == {fitted[-1]}


def test_a_circular_grid_fits_the_request_and_its_records_keep_their_own_canvas(monkeypatch) -> None:
    import gbdraw.api.diagram as diagram_module

    plans: list = []
    original = diagram_module.assemble_circular_diagram_from_record

    def keeping(*args, **kwargs):
        drawing = original(*args, **kwargs)
        plans.append(getattr(drawing, "_gbdraw_circular_assembly_result").composition_plan)
        return drawing

    monkeypatch.setattr(diagram_module, "assemble_circular_diagram_from_record", keeping)
    records = [*_mito(), *gbdraw.read_genbank(INPUTS / "BGC0000712.gbk")]
    diagram = gbdraw.draw_circular(
        records, options=gbdraw.CircularOptions(figure_width="180mm"), layout=gbdraw.CircularLayout()
    )
    plan = _circular_plan(diagram)
    assert plan.width == pytest.approx(MM_180_PX, abs=0.01)
    assert MM_180_PX - 3 <= plan.content_bounds.width <= MM_180_PX
    assert plans
    for record_plan in plans:
        assert record_plan.width == pytest.approx(record_plan.content_bounds.width, abs=0.01)
        assert record_plan.height == pytest.approx(record_plan.content_bounds.height, abs=0.01)


def test_a_circular_width_below_the_smallest_map_is_refused_with_its_parts() -> None:
    with pytest.raises(ValidationError) as caught:
        _draw_mito_circular(figure_width="100mm")
    message = str(caught.value)
    assert message.startswith("figure_width 100mm is narrower than the parts of the figure that do not shrink: ")
    assert "map at its smallest" in message and "legend" in message and "padding" in message
    assert "place the legend at the top or bottom" in message
    diagnostic = caught.value.diagnostic
    assert diagnostic["code"] == "INPUT_INVALID" and diagnostic["configPath"] == FIGURE_WIDTH_PATH
    # The stated minimum draws.
    minimum = diagnostic["minimumPx"]
    plan = _circular_plan(_draw_mito_circular(figure_width=f"{minimum}px"))
    assert plan.width == pytest.approx(minimum, abs=0.01)


def test_a_circular_height_below_the_smallest_map_is_refused_with_its_parts() -> None:
    with pytest.raises(ValidationError) as caught:
        _draw_mito_circular(figure_height="50mm")
    message = str(caught.value)
    assert message.startswith("figure_height 50mm is shorter than the figure at ")
    assert "map at its smallest" in message and "padding" in message
    assert caught.value.diagnostic["configPath"] == FIGURE_HEIGHT_PATH
    assert caught.value.diagnostic["code"] == "INPUT_INVALID"


@pytest.mark.parametrize(
    ("options", "drawn_passes"),
    [
        ({"figure_width": "60mm"}, 3),
        # The legend beside the map is taller than the map, so the excess over 50 mm stays the same while the map
        # still shrinks: the walk goes on, and the refusal states the map at its smallest.
        ({"figure_height": "50mm", "legend": "right", "config_overrides": {"objects.definition.circular.show": False}}, 4),
    ],
)
def test_a_circular_refusal_ends_where_the_map_stops_shrinking(
    monkeypatch: pytest.MonkeyPatch, options: dict, drawn_passes: int
) -> None:
    """Below about radius 100 the centre definition keeps the tracks of the mitochondrial map out, so the map has the
    same size down to the smallest radius whose tracks can be drawn. The walk ends at the second pass of that size,
    after three drawn radii, and states the same sizes as the smallest drawable map. Without the definition the map
    shrinks down to that radius, and the walk goes there as before, in four drawn radii: below every pass it tries
    radii a quarter of the way up from the largest one too small to draw, which fail early."""

    import gbdraw.api.diagram as diagram
    from gbdraw.layout.figure_size import BaseLengthFit

    fit_base_length = diagram.fit_base_length
    measured: list[int] = []

    def counted(measure, **kwargs):
        def measured_base(base: int):
            result = measure(base)
            if result is not None:
                measured.append(base)
            return result

        return fit_base_length(measured_base, **kwargs)

    def at_the_floor(measure, **kwargs):
        base = fit_base_length(measure, **kwargs).base
        while measure(base - 1) is not None:
            base -= 1
        floor = measure(base)
        return BaseLengthFit(base, floor.excess_px, floor.built, fits=False)

    refusals = []
    for walk in (counted, at_the_floor):
        monkeypatch.setattr(diagram, "fit_base_length", walk)
        with pytest.raises(ValidationError) as caught:
            _draw_mito_circular(**options)
        refusals.append((str(caught.value), caught.value.diagnostic))
    assert refusals[0] == refusals[1]
    assert len(measured) == drawn_passes


def test_tracks_that_do_not_fit_the_configured_radius_raise_as_before() -> None:
    from gbdraw.diagrams.circular.radial_layout import track_layout_reason

    overrides = {"canvas.circular.radius": 20}
    with pytest.raises(ValidationError) as auto:
        _draw_mito_circular(config_overrides=overrides)
    with pytest.raises(ValidationError) as sized:
        _draw_mito_circular(config_overrides=overrides, figure_width="180mm")
    assert track_layout_reason(sized.value) is not None
    assert str(sized.value) == str(auto.value)


def _crowded_map(count: int = 400):
    """``count`` CDS labels on a 80 kb map. 400 at 16 pt need about 3.7 times the room, so they are not placed."""

    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ATGC" * 20000), id="SYN1", name="SYN1", description="crowded map")
    record.annotations.update(molecule_type="DNA", topology="circular")
    step = 80000 // count
    for index in range(count):
        start = index * step
        record.features.append(SeqFeature(
            FeatureLocation(start, start + step * 3 // 4, strand=1 if index % 2 else -1), type="CDS",
            qualifiers={"product": [f"DNA polymerase subunit {index:03d}"]},
        ))
    return [record]


CROWDED_LABELS = {"labels.circular.scope": "outer", "labels.font_size.short": 16, "labels.font_size.long": 16}


def test_labels_not_placed_at_a_requested_size_offer_a_larger_figure() -> None:
    from gbdraw.exceptions import LabelsDoNotFitError

    with pytest.raises(LabelsDoNotFitError) as auto:
        gbdraw.draw_circular(_crowded_map(), options=gbdraw.CircularOptions(config_overrides=CROWDED_LABELS))
    assert "figure" not in str(auto.value)
    assert "figure_size" not in auto.value.report.room_remedies

    with pytest.raises(LabelsDoNotFitError) as sized:
        gbdraw.draw_circular(
            _crowded_map(),
            options=gbdraw.CircularOptions(config_overrides=CROWDED_LABELS, figure_width="150mm"),
        )
    report = sized.value.report
    assert report.room_remedies == ("filter", "label_text", "font", "figure_size", "linear")
    assert report.larger_figure is not None and report.larger_figure.height is None
    lines = str(sized.value).splitlines()
    (line,) = [line for line in lines if "--figure_width" in line]
    assert re.fullmatch(r"  - make the figure larger: --figure_width about \d+(\.\d)?mm;", line), line
    assert lines.index(line) == lines.index("  - or draw a Linear map (gbdraw linear).") - 1
    about = float(re.search(r"about ([\d.]+)mm", line)[1])
    assert about > 150 * 3.0


def test_no_fitting_row_at_a_requested_size_offers_a_larger_figure_without_an_estimate() -> None:
    """NO-ROW-WAYS at a requested size (PD-OI-114 r3): after the row search, when no row fits, the ways to make room
    also name the requested size options, without a value, since the room the labels need is not known. Without a
    requested size the words stay as they were, byte for byte."""

    from gbdraw.core.figure_length import RequestedFigureSize
    from gbdraw.labels.fit import LabelFit, LabelFitRecord, LabelFitReport, LabelFitRow, LimitReason, format_label_fit_report
    from gbdraw.layout.figure_size import LargerFigure, larger_figure

    plan = plan_composition(CompositionRequest(primary=CompositionItem("primary", Aabb(0.0, 0.0, 300.0, 200.0))))
    width, height = FigureLength(180.0, "mm"), FigureLength(6.0, "in")
    record = LabelFitRecord("AP027078.1", "circular_outer", 562, 6.0, LabelFit(30, 2))
    rows = (
        LabelFitRow("outward_both", "not_found", reason=LimitReason("xy_offset_limit", 2.0)),
        LabelFitRow("smaller_font", "not_found", reason=LimitReason("font_limit", 6.0)),
    )
    auto = LabelFitReport(records=(record,), rows=rows, drawn=False)
    assert format_label_fit_report(auto) == "\n".join([
        "AP027078.1: 562 outer labels do not fit at 6 pt (30 label pairs overlap, 2 leader pairs cross). "
        "Nothing was drawn.",
        "No single change of the label settings fits this map.",
        "To make room:",
        "  - label fewer features with --label_blacklist or --label_whitelist;",
        "  - shorten the label text with --qualifier_priority (for example gene names only);",
        "  - or draw a Linear map (gbdraw linear).",
    ])

    def sized(size: RequestedFigureSize) -> LabelFitReport:
        return LabelFitReport(records=(record,), rows=rows, drawn=False, larger_figure=larger_figure(plan, size, None))

    for size, words in (
        (RequestedFigureSize(width), "--figure_width"),
        (RequestedFigureSize(None, height), "--figure_height"),
        (RequestedFigureSize(width, height), "--figure_width and --figure_height"),
    ):
        report = sized(size)
        assert report.room_remedies == ("filter", "label_text", "figure_size", "linear")
        assert format_label_fit_report(report).splitlines()[2:] == [
            "To make room:",
            "  - label fewer features with --label_blacklist or --label_whitelist;",
            "  - shorten the label text with --qualifier_priority (for example gene names only);",
            f"  - make the figure larger: increase {words};",
            "  - or draw a Linear map (gbdraw linear).",
        ]
    # The requested sides, with no estimate: a side left auto stays None.
    assert sized(RequestedFigureSize(width)).larger_figure == LargerFigure(True, None)
    # The sides are all estimated or none, and at least one is requested (as the Web reader requires).
    for width_side, height_side in ((True, "240.3 mm"), ("240.3 mm", True), (None, None)):
        with pytest.raises(ValueError):
            LargerFigure(width_side, height_side)
    assert LargerFigure(True, True).height is True and LargerFigure("1 mm", "2 mm").width == "1 mm"
    # A row that fits leaves no ways, so no larger figure either.
    fitting = LabelFitReport(
        records=(record,),
        rows=(*rows, LabelFitRow("outward_vertical", "found", {"labels.unified_adjustment.outer_labels.y_radius_offset": 1.5}, (1.0, 1.0))),
        drawn=False,
        larger_figure=larger_figure(plan, RequestedFigureSize(width), None),
    )
    assert fitting.room_remedies == () and "make the figure larger" not in format_label_fit_report(fitting)


def test_the_fit_steps_with_the_radii_whose_labels_were_placed(monkeypatch: pytest.MonkeyPatch) -> None:
    """120 labels at 10 pt are placed at the configured radius and not placed below about 128 px (CAPACITY-STOP). A
    radius drawn without them is marked as leaving them out, so the walk stops after two radii with labels instead of
    measuring radii just above 128."""

    import gbdraw.api.diagram as diagram
    from gbdraw.exceptions import LabelsDoNotFitError

    passes: list[tuple[bool, bool]] = []
    fit_base_length = diagram.fit_base_length

    def recorded(measure, **kwargs):
        def measured(base: int):
            result = measure(base)
            if result is not None:
                drawn_labels = "DNA polymerase subunit" in result.built.drawing.tostring()
                passes.append((result.left_out, drawn_labels))
            return result

        return fit_base_length(measured, **kwargs)

    monkeypatch.setattr(diagram, "fit_base_length", recorded)
    labels = {"labels.circular.scope": "outer", "labels.font_size.short": 10, "labels.font_size.long": 10}
    with pytest.raises(LabelsDoNotFitError) as caught:
        gbdraw.draw_circular(
            _crowded_map(120), options=gbdraw.CircularOptions(config_overrides=labels, figure_width="130mm")
        )
    assert caught.value.report.larger_figure is not None
    assert all(left_out != drawn_labels for left_out, drawn_labels in passes)
    assert len(passes) == 3 and [left_out for left_out, _drawn_labels in passes].count(False) == 2


def test_the_larger_figure_estimate_adds_the_missing_room_to_each_requested_side() -> None:
    from gbdraw.core.figure_length import RequestedFigureSize
    from gbdraw.layout.figure_size import larger_figure

    plan = plan_composition(CompositionRequest(primary=CompositionItem("primary", Aabb(0.0, 0.0, 300.0, 200.0))))
    content = plan.content_bounds
    size = RequestedFigureSize(FigureLength(100.0, "mm"), FigureLength(1000.0, "px"))
    estimate = larger_figure(plan, size, 2.5)
    width_px = content.width + 1.5 * 300.0
    assert estimate.width == f"{math.ceil(round(width_px * 25.4 / 96 * 10, 6)) / 10:g} mm"
    assert estimate.height == f"{math.ceil(content.height + 1.5 * 200.0)} px"
    assert larger_figure(plan, RequestedFigureSize(FigureLength(100.0, "mm")), 2.5).height is None


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


def test_the_export_root_size_writes_each_requested_side_as_written_and_converts_the_other() -> None:
    from gbdraw.core.figure_length import RequestedFigureSize
    from gbdraw.render.document import export_root_size

    width = FigureLength(180.0, "mm")
    assert export_root_size(RequestedFigureSize(), 700.0, 500.0) is None
    assert export_root_size(RequestedFigureSize(width=width), width.px, 96 / 25.4 * 100) == {
        "width": "180mm", "height": "100mm"}
    assert export_root_size(RequestedFigureSize(height=FigureLength(4.0, "in")), 480.0, 384.0) == {
        "width": "5in", "height": "4in"}
    both = RequestedFigureSize(width=width, height=FigureLength(15.0, "cm"))
    assert export_root_size(both, width.px + 0.005, FigureLength(15.0, "cm").px) == {"width": "180mm", "height": "15cm"}


def test_the_export_root_size_converts_a_requested_side_the_canvas_grew_past() -> None:
    # A Web-only addition after the fit (Legend edits, Canvas padding) can make the canvas larger than the box.
    from gbdraw.core.figure_length import RequestedFigureSize
    from gbdraw.render.document import export_root_size

    size = RequestedFigureSize(width=FigureLength(180.0, "mm"))
    assert export_root_size(size, 200 * 96 / 25.4, 96.0) == {"width": "200mm", "height": "25.4mm"}


def test_a_sized_web_render_keeps_the_page_root_in_px(tmp_path) -> None:
    """OV-2512: the Web reply read the written file, whose root states 300mm, beside px pieces; the page root is px."""

    from tests.test_web_label_fit import _sars_request, _web_render

    reply = _web_render(_sars_request({"canvas.figure_width": "300mm", "labels.circular.scope": "none"}), tmp_path)
    result = reply["results"][0]
    width, _height, view_box = _svg_size(web_result_svg(result))
    assert width.endswith("px")
    assert float(width.removesuffix("px")) == pytest.approx(300 * 96 / 25.4, abs=0.01) == view_box[2]
    assert web_result_svg(result).startswith(result["document"]["head"])
    # At Auto the page's text is the written file's, byte for byte.
    auto = _web_render(_sars_request({"labels.circular.scope": "none"}), tmp_path / "auto")["results"][0]
    assert '<?xml version="1.0" encoding="utf-8" ?>\n' + web_result_svg(auto) == (
        tmp_path / "auto" / "out" / auto["name"]
    ).read_text(encoding="utf-8")


def test_the_web_export_check_states_a_sized_svg_root_as_python_writes_it() -> None:
    import json

    from gbdraw.web_support.export_limit import check_web_export_limit_json

    request = {"format": "svg", "width": MM_180_PX, "height": 96.0, "dpi": None, "linearWrap": None,
               "figureSize": {"width": "180 MM", "height": None}}
    assert json.loads(check_web_export_limit_json(json.dumps(request))) == {
        "ok": True, "root": {"width": "180mm", "height": "25.4mm"}}
    with pytest.raises(ValueError):
        check_web_export_limit_json(json.dumps({**request, "figureSize": None}))


def _composition_metadata(svg: str) -> dict:
    import html
    import json

    head = svg[svg.index("<svg") : svg.index(">", svg.index("<svg"))]
    found = re.search(r'data-gbdraw-composition="([^"]+)"', head)
    assert found
    return json.loads(html.unescape(found[1]))


def test_only_a_sized_figure_carries_its_requested_canvas_in_the_composition_metadata() -> None:
    """Section 6: the Web replans a Legend edit within the requested canvas; Auto metadata keeps its keys and bytes."""

    auto = _composition_metadata(gbdraw.draw_linear(_mito()).to_svg())
    assert "requestedCanvas" not in auto
    assert sorted(auto) == ["legend", "legendReflow", "legendSide", "lengthBar", "overlayObstacles", "overlayPolicy",
                            "primary", "spacing", "title", "titleSide"]
    sized = _composition_metadata(
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm")).to_svg())
    assert sized["requestedCanvas"] == {"widthPx": pytest.approx(MM_180_PX), "heightPx": None}
    circular = _composition_metadata(gbdraw.draw_circular(_mito(), options=gbdraw.CircularOptions(
        figure_height="150mm")).to_svg())
    assert circular["requestedCanvas"] == {"widthPx": None, "heightPx": pytest.approx(150 * 96 / 25.4)}


def test_png_dpi_scales_the_pixels() -> None:
    """Without a DPI a PNG is drawn at the default, 300 dpi, as the Web app draws it (D-83 revised)."""

    pytest.importorskip("cairosvg")
    diagram = gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm"))
    _w, _h, view_box = _svg_size(diagram.to_svg())
    at_300 = (round(view_box[2] * 300 / 96), round(view_box[3] * 300 / 96))
    assert _png_size(diagram.to_bytes("png")) == at_300
    assert _png_size(diagram.to_bytes("png", dpi=300)) == at_300
    assert _png_size(diagram.to_bytes("png", dpi=96)) == (round(view_box[2]), round(view_box[3]))


def test_the_default_png_dpi_is_one_python_constant_among_the_session_choices() -> None:
    """The Web reads the default generated; it is one of the DPI menu's choices."""

    from gbdraw.render.export_limits import DEFAULT_PNG_DPI
    from gbdraw.session import SESSION_PNG_DPI_CHOICES

    assert DEFAULT_PNG_DPI == 300
    assert DEFAULT_PNG_DPI in SESSION_PNG_DPI_CHOICES
    generated = (Path(__file__).resolve().parents[1] / "gbdraw/web/js/contracts/schema.generated.js").read_text(
        encoding="utf-8"
    )
    assert f"export const DEFAULT_PNG_DPI = {DEFAULT_PNG_DPI};" in generated


def test_a_png_carries_the_resolution_it_was_drawn_at() -> None:
    """pHYs: pixels per metre of the requested DPI; the default 300 without one."""

    pytest.importorskip("cairosvg")
    diagram = gbdraw.draw_linear(_mito())
    assert png_phys(diagram.to_bytes("png")) == (11811, 11811, 1)
    assert png_phys(diagram.to_bytes("png", dpi=96)) == (3780, 3780, 1)
    assert png_phys(diagram.to_bytes("png", dpi=300)) == (11811, 11811, 1)
    assert png_phys(diagram.to_bytes("png", dpi=72.5)) == (2854, 2854, 1)


def test_every_python_png_writer_carries_the_resolution(tmp_path) -> None:
    pytest.importorskip("cairosvg")
    from gbdraw.api import render_to_bytes, save_figure_to

    canvas = gbdraw.draw_linear(_mito())._drawing
    assert png_phys(render_to_bytes(canvas, "png", dpi=192)) == (7559, 7559, 1)
    paths = save_figure_to(canvas, ["png"], output_dir=str(tmp_path), output_prefix="saved", dpi=300)
    assert png_phys((tmp_path / "saved.png").read_bytes()) == (11811, 11811, 1)
    assert [Path(path).name for path in paths] == ["saved.svg", "saved.png"]
    save_figure_to(canvas, ["png"], output_dir=str(tmp_path), output_prefix="plain")
    assert png_phys((tmp_path / "plain.png").read_bytes()) == (11811, 11811, 1)
    assert png_phys(render_to_bytes(canvas, "png")) == (11811, 11811, 1)
    gbdraw.draw_linear(_mito()).save(tmp_path / "diagram.png", dpi=96)
    assert png_phys((tmp_path / "diagram.png").read_bytes()) == (3780, 3780, 1)


# The pHYs writer and its Web twin, setDpiInPng, read the same vectors (tools/twins.json).
PNG_PHYS_VECTORS = json.loads((Path(__file__).parent / "fixtures" / "png_phys_vectors.json").read_text(encoding="utf-8"))


def test_png_resolution_writes_the_chunk_the_web_writes() -> None:
    """Pixels per metre round half up (D-84); a pHYs chunk is replaced, else inserted before the first IDAT."""

    from gbdraw.render.export import png_with_resolution

    bare = bytes.fromhex(PNG_PHYS_VECTORS["chunks"][0]["png"])
    for case in PNG_PHYS_VECTORS["pixelsPerMetre"]:
        per_metre = case["pixelsPerMetre"]
        assert png_phys(png_with_resolution(bare, case["dpi"])) == (per_metre, per_metre, 1), case
    for case in PNG_PHYS_VECTORS["chunks"]:
        assert png_with_resolution(bytes.fromhex(case["png"]), case["dpi"]).hex() == case["expected"], case["case"]
    inserted, replaced = (bytes.fromhex(case["expected"]) for case in PNG_PHYS_VECTORS["chunks"])
    assert [kind for kind, _data in png_chunks(inserted)] == [b"IHDR", b"pHYs", b"IDAT", b"IDAT", b"IEND"]
    assert [kind for kind, _data in png_chunks(replaced)] == [b"IHDR", b"pHYs", b"IDAT", b"IEND"]


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


def test_cli_png_without_dpi_is_drawn_at_300_dpi(tmp_path, monkeypatch) -> None:
    pytest.importorskip("cairosvg")
    output = _run_cli(tmp_path, monkeypatch, "linear", "--gbk", str(INPUTS / "HmmtDNA.gbk"), "-f", "svg,png")
    _width, _height, view_box = _svg_size(output.with_suffix(".svg").read_text(encoding="utf-8"))
    png = output.with_suffix(".png").read_bytes()
    assert _png_size(png) == (round(view_box[2] * 300 / 96), round(view_box[3] * 300 / 96))
    assert png_phys(png) == (11811, 11811, 1)


@pytest.mark.parametrize("args", [["--figure_width", "180"], ["--dpi", "0"]])
def test_cli_refuses_bad_values(tmp_path, monkeypatch, args) -> None:
    with pytest.raises(SystemExit) as caught:
        _run_cli(tmp_path, monkeypatch, "linear", "--gbk", str(INPUTS / "HmmtDNA.gbk"), *args)
    assert caught.value.code == 2


def test_cli_warns_when_dpi_has_no_png(tmp_path, monkeypatch, caplog) -> None:
    _run_cli(tmp_path, monkeypatch, "linear", "--gbk", str(INPUTS / "HmmtDNA.gbk"), "--dpi", "300", "-f", "svg")
    assert "--dpi applies to PNG output only" in caplog.text


def test_cli_circular_figure_width(tmp_path, monkeypatch) -> None:
    output = _run_cli(
        tmp_path, monkeypatch, "circular", "--gbk", str(INPUTS / "HmmtDNA.gbk"), "--figure_width", "180mm"
    )
    width, _height, view_box = _svg_size(output.with_suffix(".svg").read_text(encoding="utf-8"))
    assert width == "180mm"
    assert view_box[2] == pytest.approx(MM_180_PX, abs=0.01)


# --- what the Web reads of a refusal, the remedy estimate, and a form's replay (PR-C) ---


def _web_context(error: BaseException) -> dict:
    from gbdraw.web_support.error_adapter import serialize_web_error

    return serialize_web_error(error, operation="generate", stage="render")["context"]


def test_a_refused_width_carries_the_lengths_of_its_message_to_the_web() -> None:
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="40mm"))
    message = str(caught.value)
    context = _web_context(caught.value)
    assert context["reason"] == "FIGURE_TOO_NARROW" and context["configPath"] == FIGURE_WIDTH_PATH
    assert (context["requested"], context["primaryItem"], context["legendRemedy"]) == ("40 mm", "record_ends", True)
    # Each length is the one the message writes: the Web words the same numbers.
    assert f"({context['minimum']} in all)" in message
    assert f"record names and labels past the sequence ends {context['primary']}, legend {context['legend']}" in message
    assert f"padding {context['padding']}" in message and "title" not in context


def test_a_refused_width_without_a_side_legend_offers_no_legend_remedy_to_the_web() -> None:
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="10mm", legend="top"))
    context = _web_context(caught.value)
    assert context["reason"] == "FIGURE_TOO_NARROW" and "legendRemedy" not in context


def test_a_refused_height_carries_the_width_it_was_measured_at() -> None:
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(figure_width="180mm", figure_height="10mm"))
    message = str(caught.value)
    context = _web_context(caught.value)
    assert (context["reason"], context["requested"], context["atWidth"]) == ("FIGURE_TOO_SHORT", "10 mm", "180 mm")
    assert context["primaryItem"] == "records" and f"records {context['primary']}" in message
    assert f"at least {context['minimum']}" in message


def test_a_refused_circular_size_names_the_map() -> None:
    with pytest.raises(ValidationError) as narrow:
        _draw_mito_circular(figure_width="100mm")
    with pytest.raises(ValidationError) as short:
        _draw_mito_circular(figure_height="50mm")
    assert _web_context(narrow.value)["primaryItem"] == "map"
    assert (_web_context(short.value)["reason"], _web_context(short.value)["primaryItem"]) == ("FIGURE_TOO_SHORT", "map")


def test_the_web_reads_figure_size_lengths_only_from_a_figure_size_refusal() -> None:
    lengths = {"requested": "40 mm", "minimum": "112.3 mm", "primaryItem": "map", "legendRemedy": True}
    other = ValidationError("PRIVATE", diagnostic={"code": "INPUT_INVALID", "reason": "WRAP", **lengths})
    assert not set(lengths) & set(_web_context(other))
    unbounded = ValidationError(
        "PRIVATE",
        diagnostic={
            "code": "INPUT_INVALID", "reason": "FIGURE_TOO_NARROW", "requested": "40 mm <b>", "minimum": "1e9 mm",
            "primaryItem": "PRIVATE", "legendRemedy": "yes", "padding": "12.2 mm",
        },
    )
    context = _web_context(unbounded)
    assert context["reason"] == "FIGURE_TOO_NARROW" and context["padding"] == "12.2 mm"
    assert not {"requested", "minimum", "primaryItem", "legendRemedy"} & set(context)


@pytest.mark.parametrize(
    ("form", "expected"),
    [
        ({"figure_width": "180mm"}, ["--figure_width", "180mm"]),
        ({"figure_width": "180mm", "figure_height": "5in"}, ["--figure_width", "180mm", "--figure_height", "5in"]),
        ({"figure_height": "1200px"}, ["--figure_height", "1200px"]),
        ({"figure_width": None, "figure_height": "auto"}, []),
        ({"figure_width": "", "figure_height": "AUTO"}, []),
        ({}, []),
    ],
)
def test_a_drawing_form_size_replays_as_its_cli_options(form: dict, expected: list[str]) -> None:
    from gbdraw.session_io import _append_common_gui_args

    run_args: list[str] = []
    _append_common_gui_args(run_args, [], form=form, adv={})
    assert [arg for index, arg in enumerate(run_args) if arg.startswith("--figure_") or (
        index > 0 and run_args[index - 1].startswith("--figure_"))] == expected


def test_the_web_words_each_primary_item_and_part_python_names() -> None:
    from gbdraw.layout.figure_size import FIGURE_PRIMARY_ITEMS, FIGURE_SIZE_PARTS

    wording = (Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "js" / "utils" / "error-normalization.js").read_text(
        encoding="utf-8"
    )
    items = wording[wording.index("const FIGURE_PRIMARY_ITEMS") :]
    parts = wording[wording.index("const FIGURE_PARTS") :]
    assert set(re.findall(r"\b([a-z_]+): '", items[: items.index("});")])) == set(FIGURE_PRIMARY_ITEMS)
    assert {*re.findall(r"\b([a-z_]+): '", parts[: parts.index("});")]), "primary"} == set(FIGURE_SIZE_PARTS)
    assert f"@typedef {{{' | '.join(repr(item) for item in FIGURE_PRIMARY_ITEMS)}}} FigurePrimaryItem" in wording
