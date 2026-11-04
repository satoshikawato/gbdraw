"""A requested figure size: the base length that fills it, and the refusals when the figure cannot fit (D-03, D-09).

A mode keeps its fonts and line widths and changes one base length (the Linear sequence width, the Circular radius)
until its composed content fills the requested box. The composition planner centres the
content in a larger requested canvas and never shrinks it, so the checks here read the drawn plan.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from typing import Generic, Literal, TypeVar

from gbdraw.config.models.canvas import FIGURE_HEIGHT_PATH, FIGURE_WIDTH_PATH
from gbdraw.core.figure_length import FigureLength, RequestedFigureSize
from gbdraw.exceptions import ValidationError
from gbdraw.layout.composition import CompositionPlan
from gbdraw.layout.spatial import Aabb, union_aabbs

Built = TypeVar("Built")
# Each step aims this far under the target, and a pass this close under it is accepted: the planner centres the
# rest in the requested canvas, so the canvas stays exact. A pass over the target is never accepted.
_AIM_UNDER_PX = 2.0
_ACCEPT_UNDER_PX = 3.0
_MAX_PASSES = 8
# A base too small to draw fails early (a Circular map's track layout, about 0.05 s), so such passes have their own
# bound: halving the distance between the largest of them and the smallest drawn base reaches 1 from a few thousand.
_MAX_UNDRAWABLE_PASSES = 12
# The composed size grows with the base length, but labels past the ends can make a stretch nearly flat; a slope
# measured there is raised to this so one step cannot jump far.
_MIN_SLOPE = 0.25


@dataclass(frozen=True)
class BaseLengthFit(Generic[Built]):
    """A measured base length, how far its composed content is over the requested box, and what it assembled.

    ``excess_px`` is the largest amount by which a fitted side of the content exceeds the requested side (negative
    when under). ``fits`` is False only for a refusal; ``built`` is then the smallest figure that could be drawn,
    whose parts the refusal names. The caller keeps a fitting ``built`` as the drawn figure.
    """

    base: int
    excess_px: float
    built: Built
    fits: bool


@dataclass(frozen=True)
class _Pass(Generic[Built]):
    """A drawn pass and the Python warnings it raised, which reach the caller only if the pass is kept."""

    fit: BaseLengthFit[Built]
    caught: tuple[warnings.WarningMessage, ...]


def _measured(
    measure: Callable[[int], tuple[float, Built] | None], base: int
) -> tuple[tuple[float, Built] | None, tuple[warnings.WarningMessage, ...]]:
    with warnings.catch_warnings(record=True) as caught:
        result = measure(base)
    return result, tuple(caught)


def _kept(kept: _Pass[Built]) -> BaseLengthFit[Built]:
    """The kept pass, with its warnings emitted once, in order (D-16)."""

    for caught in kept.caught:
        warnings.warn_explicit(caught.message, caught.category, caught.filename, caught.lineno, source=caught.source)
    return kept.fit


def _step(start: BaseLengthFit[Built], slope: float) -> int:
    return start.base + math.floor((-_AIM_UNDER_PX - start.excess_px) / max(_MIN_SLOPE, slope))


def _slope(a: BaseLengthFit[Built], b: BaseLengthFit[Built]) -> float:
    return (b.excess_px - a.excess_px) / (b.base - a.base)


def _next_base(
    drawn: dict[int, _Pass[Built]], order: list[int], floor: int, undrawable: int, first_slope: float
) -> int | None:
    """The next base to measure, or None when no unmeasured base can come closer to the target."""

    if len(drawn) >= _MAX_PASSES:
        return None
    passes = [drawn_pass.fit for drawn_pass in drawn.values()]
    under = max((fit for fit in passes if fit.fits), key=lambda fit: fit.base, default=None)
    over = min(
        (fit for fit in passes if not fit.fits and (under is None or fit.base > under.base)),
        key=lambda fit: fit.base,
        default=None,
    )
    if under is not None and over is not None:
        if over.base - under.base <= 1:
            return None
        # The secant between the bracket's ends, kept strictly inside it: one rounded down onto the pass under still
        # has a base above it to try.
        return min(max(_step(under, _slope(under, over)), under.base + 1), over.base - 1)
    last = drawn[order[-1]].fit
    guess = _step(last, first_slope if len(order) < 2 else _slope(drawn[order[-2]].fit, last))
    if under is not None:
        return max(guess, under.base + 1)
    lowest = min(drawn)
    guess = max(1, min(guess, lowest - 1))
    if guess <= floor:
        # Bases at or below the floor cannot be drawn: halve the way down to it instead.
        if undrawable >= _MAX_UNDRAWABLE_PASSES:
            return None
        guess = (floor + lowest) // 2
    return None if guess <= floor or guess in drawn else guess


def fit_base_length(
    measure: Callable[[int], tuple[float, Built] | None], *, first_base: int, slope: float
) -> BaseLengthFit[Built]:
    """The base length whose composed content is within a few px under the requested box.

    ``measure(base)`` returns how far the content drawn at ``base`` is over the box (``excess_px``, negative when
    under) and what it assembled, or None when ``base`` is too small to draw. ``first_base`` must be drawable; the
    mode supplies it and ``slope``, the px of content per px of base of the first step. Each step aims 2 px under:
    later ones with the secant of the last two drawn passes. Once one pass is under and one over, the next base is
    the secant between them, kept strictly inside, until no base lies between them; then the widest pass under is
    kept. A base at or below the largest one too small to draw is not tried: the step halves the way to it instead.
    After 8 drawn passes the widest pass under is kept.

    With no pass under, the target is refused (``fits`` False) with the smallest drawn pass. Before refusing, the
    middle of that pass and the next drawn pass above it is measured once, if not yet measured, since a measure
    that is not monotonic may fit there. Only the kept pass's Python warnings reach the caller, once; a refusal
    emits none.
    """

    drawn: dict[int, _Pass[Built]] = {}
    order: list[int] = []
    floor = 0
    undrawable = 0
    checked_middle = False
    base: int | None = first_base
    while base is not None:
        result, caught = _measured(measure, base)
        if result is None:
            if not drawn:
                raise ValueError(f"the first base length {base} cannot be drawn")
            floor, undrawable = max(floor, base), undrawable + 1
        else:
            excess, built = result
            current = _Pass(BaseLengthFit(base, excess, built, fits=excess <= 0.0), caught)
            if current.fit.fits and excess >= -_ACCEPT_UNDER_PX:
                return _kept(current)
            drawn[base] = current
            order.append(base)
        base = _next_base(drawn, order, floor, undrawable, slope)
        if base is None and not any(drawn_pass.fit.fits for drawn_pass in drawn.values()):
            lowest = min(drawn)
            above = min((other for other in drawn if other > lowest), default=None)
            if not checked_middle and len(drawn) < _MAX_PASSES and above is not None and above - lowest > 1:
                checked_middle = True
                base = (lowest + above) // 2
    fitting = [drawn_pass for drawn_pass in drawn.values() if drawn_pass.fit.fits]
    if fitting:
        return _kept(max(fitting, key=lambda drawn_pass: drawn_pass.fit.excess_px))
    return drawn[min(drawn)].fit


@dataclass(frozen=True)
class _AxisParts:
    """What one axis of the padded content is made of, in px; the parts add up to the content."""

    primary: float
    legend: float
    title: float
    padding: float

    @property
    def total(self) -> float:
        return self.primary + self.legend + self.title + self.padding


def _span(boxes: list[Aabb], axis: Literal["x", "y"]) -> float:
    union = union_aabbs(boxes)
    if union is None:
        return 0.0
    return union.width if axis == "x" else union.height


def _axis_parts(plan: CompositionPlan, axis: Literal["x", "y"], *, legend_role: str, title_role: str) -> _AxisParts:
    legend = plan.placement_for(legend_role)
    title = plan.placement_for(title_role)
    boxes = [plan.primary_bounds]
    primary = _span(boxes, axis)
    if legend is not None:
        boxes.append(legend.final_bounds)
    with_legend = _span(boxes, axis)
    if title is not None:
        boxes.append(title.final_bounds)
    painted = _span(boxes, axis)
    content = plan.content_bounds.width if axis == "x" else plan.content_bounds.height
    return _AxisParts(primary, with_legend - primary, painted - with_legend, content - painted)


def _beside_primary(plan: CompositionPlan, legend_role: str) -> bool:
    """Whether the legend sits left or right of the primary item, so it adds to the width."""

    legend = plan.placement_for(legend_role)
    if legend is None:
        return False
    bounds, primary = legend.final_bounds, plan.primary_bounds
    return bounds.max_x <= primary.min_x or bounds.min_x >= primary.max_x


def _parts_text(length: FigureLength, parts: list[tuple[str, float]]) -> str:
    return ", ".join(f"{name} {length.message_text(px)}" for name, px in parts if px > 0.0)


def figure_too_short(
    plan: CompositionPlan,
    height: FigureLength,
    width: FigureLength | None,
    *,
    primary_name: str,
    legend_role: str = "legend",
    title_role: str = "title",
) -> ValidationError:
    """The refusal of a requested height shorter than the padded content of ``plan``.

    ``primary_name`` names the primary item (Linear records, a Circular map at its smallest); ``width`` is the
    requested width, if any.
    """

    content = plan.content_bounds
    shown_width = width if width is not None else height.of_px(plan.width)
    parts = _axis_parts(plan, "y", legend_role=legend_role, title_role=title_role)
    minimum = height.message_text(content.height, round_up=True)
    detail = _parts_text(
        height,
        [(primary_name, parts.primary), ("legend", parts.legend), ("title", parts.title), ("padding", parts.padding)],
    )
    return ValidationError(
        f"figure_height {height.text} is shorter than the figure at {shown_width.message_text(plan.width)} "
        f"wide: {minimum} ({detail}). Use a figure_height of at least {minimum}, or leave figure_height auto.",
        diagnostic={
            "code": "INPUT_INVALID",
            "configPath": FIGURE_HEIGHT_PATH,
            "minimumPx": math.ceil(content.height),
        },
    )


def check_figure_height(plan: CompositionPlan, size: RequestedFigureSize, *, primary_name: str) -> None:
    """Refuse a drawn figure whose padded content is taller than a requested height.

    A requested width is the base-length fit's: the figure drawn is the one that fit it.
    """

    if size.height is not None and plan.content_bounds.height > size.height.px:
        raise figure_too_short(plan, size.height, size.width, primary_name=primary_name)


def figure_size_excess(plan: CompositionPlan, size: RequestedFigureSize) -> float:
    """How far the padded content of ``plan`` is over the requested box: the largest excess over a requested side."""

    content = plan.content_bounds
    excesses = [
        drawn - length.px
        for drawn, length in ((content.width, size.width), (content.height, size.height))
        if length is not None
    ]
    if not excesses:
        raise ValueError("an excess needs a requested side")
    return max(excesses)


def figure_size_refusal(plan: CompositionPlan, size: RequestedFigureSize, *, primary_name: str) -> ValidationError:
    """The refusal of a requested box that ``plan``, the smallest figure the fit drew, does not fit: the width
    refusal when it is wider than a requested width, else the height refusal."""

    if size.height is None or (size.width is not None and plan.content_bounds.width > size.width.px):
        if size.width is None:
            raise ValueError("a refusal needs a requested side")
        return figure_too_narrow(plan, size.width, fixed_primary_name=primary_name)
    return figure_too_short(plan, size.height, size.width, primary_name=primary_name)


@dataclass(frozen=True)
class LargerFigure:
    """About how large a requested figure must be for labels that were not placed for lack of room (CAPACITY-STOP):
    each requested side as an argument value (``240.3mm``), None for a side left auto."""

    width: str | None
    height: str | None


def larger_figure(plan: CompositionPlan, size: RequestedFigureSize, ratio: float) -> LargerFigure:
    """The size the labels need, from the drawn figure ``plan`` and the room the labels need as a multiple of the
    room around the map (``ratio``): each requested side grows by ``ratio - 1`` times the map's side, rounded up."""

    def side(length: FigureLength | None, content: float, primary: float) -> str | None:
        if length is None:
            return None
        # The message's rounded value ("240.3 mm") written as an argument.
        return length.message_text(content + (ratio - 1.0) * primary, round_up=True).replace(" ", "")

    content, primary = plan.content_bounds, plan.primary_bounds
    return LargerFigure(
        side(size.width, content.width, primary.width), side(size.height, content.height, primary.height)
    )


def figure_too_narrow(
    plan: CompositionPlan,
    width: FigureLength,
    *,
    fixed_primary_name: str,
    legend_role: str = "legend",
    title_role: str = "title",
) -> ValidationError:
    """The refusal of a requested width narrower than ``plan``, the smallest figure the fit drew.

    The narrowest figure is that plan's padded content. When the primary row (the primary item with a legend
    beside it) is its widest item, the message names ``fixed_primary_name`` (the primary item at its smallest), the
    legend beside it, and the padding; otherwise the widest item above or below (the title or the legend) and the
    padding.
    """

    content = plan.content_bounds.width
    legend = plan.placement_for(legend_role)
    title = plan.placement_for(title_role)
    beside = _beside_primary(plan, legend_role)
    row_boxes = [plan.primary_bounds]
    if legend is not None and beside:
        row_boxes.append(legend.final_bounds)
    primary = _span([plan.primary_bounds], "x")
    row = _span(row_boxes, "x")
    above_or_below = [
        (name, placement.final_bounds)
        for name, placement in (("title", title), ("legend", None if beside else legend))
        if placement is not None
    ]
    painted = _span([*row_boxes, *(bounds for _name, bounds in above_or_below)], "x")
    row_binds = painted - row < 0.5
    if row_binds:
        parts = [(fixed_primary_name, primary), ("legend", row - primary)]
    else:
        widest = max(above_or_below, key=lambda item: item[1].width)[0]
        parts = [(widest, painted)]
    minimum = width.message_text(content, round_up=True)
    detail = _parts_text(width, [*parts, ("padding", content - painted)])
    legend_remedy = " place the legend at the top or bottom," if row_binds and beside else ""
    return ValidationError(
        f"figure_width {width.text} is narrower than the parts of the figure that do not shrink: {detail} "
        f"({minimum} in all). Use a figure_width of at least {minimum},{legend_remedy} or leave figure_width auto.",
        diagnostic={"code": "INPUT_INVALID", "configPath": FIGURE_WIDTH_PATH, "minimumPx": math.ceil(content)},
    )


__all__ = [
    "BaseLengthFit",
    "LargerFigure",
    "check_figure_height",
    "figure_size_excess",
    "figure_size_refusal",
    "figure_too_narrow",
    "figure_too_short",
    "fit_base_length",
    "larger_figure",
]
