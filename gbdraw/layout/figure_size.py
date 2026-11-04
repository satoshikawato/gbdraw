"""A requested figure size: the base length that fills it, and the refusals when the figure cannot fit (D-03, D-09).

A mode keeps its fonts and line widths and changes one base length (the Linear sequence width; the Circular radius
in a later change) until its composed content fills the requested side. The composition planner centres the
content in a larger requested canvas and never shrinks it, so the checks here read the drawn plan.
"""

from __future__ import annotations

import math
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
# The composed size grows with the base length, but labels past the ends can make a stretch nearly flat; a slope
# measured there is raised to this so one step cannot jump far.
_MIN_SLOPE = 0.25


@dataclass(frozen=True)
class BaseLengthFit(Generic[Built]):
    """A measured base length, its composed size, and what that measurement assembled.

    ``fits`` is False only for base length 1 still over the target; ``built`` is then the figure at base 1, whose
    parts the refusal names. The caller keeps a fitting ``built`` as the drawn figure.
    """

    base: int
    size_px: float
    built: Built
    fits: bool


def _step(start: BaseLengthFit[Built], aim: float, slope: float) -> int:
    return start.base + math.floor((aim - start.size_px) / max(_MIN_SLOPE, slope))


def _slope(a: BaseLengthFit[Built], b: BaseLengthFit[Built]) -> float:
    return (b.size_px - a.size_px) / (b.base - a.base)


def fit_base_length(measure: Callable[[int], tuple[float, Built]], target_px: float) -> BaseLengthFit[Built]:
    """The base length whose composed size (``measure``) is within a few px under ``target_px``.

    The first measurement uses the target as the base. Each step aims 2 px under the target: the first with
    slope 1, later ones with the slope between the last two passes. Once one pass is under and one over, the next
    base stays between them. After 8 passes the widest pass under the target is taken. The base never goes below
    1, and the target is refused (``fits`` False) only when base 1 is over it.
    """

    aim = target_px - _AIM_UNDER_PX
    measured: dict[int, BaseLengthFit[Built]] = {}
    under: BaseLengthFit[Built] | None = None
    over: BaseLengthFit[Built] | None = None
    last: BaseLengthFit[Built] | None = None
    base = max(1, math.floor(target_px))
    while True:
        size, built = measure(base)
        current = BaseLengthFit(base, size, built, fits=size <= target_px)
        if current.fits and size >= target_px - _ACCEPT_UNDER_PX:
            return current
        if not current.fits and base == 1:
            return current
        if current.fits and (under is None or base > under.base):
            under = current
        if not current.fits and (over is None or base < over.base):
            over = current
        previous, last = last, current
        measured[base] = current
        if len(measured) >= _MAX_PASSES:
            if under is not None:
                return under
            base = 1
            continue
        if under is not None and over is not None:
            if over.base - under.base <= 1:
                return under
            base = _step(under, aim, _slope(under, over))
            if not under.base < base < over.base or base in measured:
                base = (under.base + over.base) // 2
            if base in measured:
                return under
            continue
        base = max(1, _step(current, aim, 1.0 if previous is None else _slope(previous, current)))
        if base in measured:
            base = current.base + (1 if current.fits else -1)
        if base in measured or base < 1:
            # Only reachable when the size is not monotonic in the base.
            if under is not None:
                return under
            base = 1


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


def check_figure_height(
    plan: CompositionPlan,
    size: RequestedFigureSize,
    *,
    legend_role: str = "legend",
    title_role: str = "title",
) -> None:
    """Refuse a drawn figure whose padded content is taller than a requested height.

    A requested width is the base-length fit's: the figure drawn is the one that fit it.
    """

    content = plan.content_bounds
    if size.height is None or content.height <= size.height.px:
        return
    width = size.width if size.width is not None else size.height.of_px(plan.width)
    parts = _axis_parts(plan, "y", legend_role=legend_role, title_role=title_role)
    minimum = size.height.message_text(content.height, round_up=True)
    detail = _parts_text(
        size.height,
        [("records", parts.primary), ("legend", parts.legend), ("title", parts.title), ("padding", parts.padding)],
    )
    raise ValidationError(
        f"figure_height {size.height.text} is shorter than the figure at {width.message_text(plan.width)} "
        f"wide: {minimum} ({detail}). Use a figure_height of at least {minimum}, or leave figure_height auto.",
        diagnostic={
            "code": "INPUT_INVALID",
            "configPath": FIGURE_HEIGHT_PATH,
            "minimumPx": math.ceil(content.height),
        },
    )


def figure_too_narrow(
    plan: CompositionPlan,
    width: FigureLength,
    *,
    fixed_primary_name: str,
    legend_role: str = "legend",
    title_role: str = "title",
) -> ValidationError:
    """The refusal of a requested width narrower than ``plan``, the figure at base length 1.

    The narrowest figure is that plan's padded content. When the primary row (the primary item with a legend
    beside it) is its widest item, the message names ``fixed_primary_name`` (the primary item at base 1), the
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


__all__ = ["BaseLengthFit", "check_figure_height", "figure_too_narrow", "fit_base_length"]
