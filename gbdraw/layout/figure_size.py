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
# Measurements after the first: the first gives the fixed part, the next lands on the target when the composed size
# grows one px per px of base length; the rest correct a slope that is not exactly one.
_CORRECTION_PASSES = 3


@dataclass(frozen=True)
class BaseLengthFit(Generic[Built]):
    """The widest measured base length whose composed size fits the target, or the probe when none fits.

    ``built`` is what that measurement assembled, so the caller can keep it as the drawn figure.
    """

    base: int
    size_px: float
    built: Built
    fits: bool


def fit_base_length(measure: Callable[[int], tuple[float, Built]], target_px: float) -> BaseLengthFit[Built]:
    """The base length whose composed size (``measure``) comes closest to ``target_px`` without passing it.

    The first measurement uses the target itself as the base; the composed size less the base is the part that
    does not grow with it. Each next base adds what is left to the target (slope 1), down when the size passed it
    and up while at least 1 px is left, for at most a few passes.
    """

    base = max(1, math.floor(target_px))
    size, built = measure(base)
    probe = BaseLengthFit(base, size, built, fits=False)
    if size - base >= target_px:
        return probe
    measured = {base}
    best = BaseLengthFit(base, size, built, fits=True) if size <= target_px else None
    for _ in range(_CORRECTION_PASSES):
        if best is not None and target_px - best.size_px < 1.0:
            break
        base = math.floor(base + target_px - size)
        if base < 1 or base in measured:
            break
        measured.add(base)
        size, built = measure(base)
        if size <= target_px and (best is None or base > best.base):
            best = BaseLengthFit(base, size, built, fits=True)
    return best if best is not None else probe


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


def check_figure_size(
    plan: CompositionPlan,
    size: RequestedFigureSize,
    *,
    base_px: float,
    fixed_primary_name: str,
    legend_role: str = "legend",
    title_role: str = "title",
) -> None:
    """Refuse a drawn figure whose padded content is larger than a requested side.

    ``base_px`` is the base length of the primary item along the width (the Linear sequence width), and
    ``fixed_primary_name`` names the rest of the primary item's width in the message.
    """

    content = plan.content_bounds
    if size.width is not None and content.width > size.width.px:
        raise figure_too_narrow(
            plan, size.width, base_px=base_px, fixed_primary_name=fixed_primary_name,
            legend_role=legend_role, title_role=title_role,
        )
    if size.height is not None and content.height > size.height.px:
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
    base_px: float,
    fixed_primary_name: str,
    legend_role: str = "legend",
    title_role: str = "title",
) -> ValidationError:
    """The refusal of a requested width that the parts that do not shrink fill (``plan`` drew ``base_px``)."""

    parts = _axis_parts(plan, "x", legend_role=legend_role, title_role=title_role)
    fixed = parts.total - base_px
    minimum = width.message_text(fixed, round_up=True)
    detail = _parts_text(
        width,
        [
            (fixed_primary_name, parts.primary - base_px),
            ("legend", parts.legend),
            ("title", parts.title),
            ("padding", parts.padding),
        ],
    )
    legend_remedy = " place the legend at the top or bottom," if _beside_primary(plan, legend_role) else ""
    return ValidationError(
        f"figure_width {width.text} is narrower than the parts of the figure that do not shrink: {detail} "
        f"({minimum} in all). Use a figure_width wider than {minimum},{legend_remedy} or leave figure_width auto.",
        diagnostic={"code": "INPUT_INVALID", "configPath": FIGURE_WIDTH_PATH, "minimumPx": math.ceil(fixed)},
    )


__all__ = ["BaseLengthFit", "check_figure_size", "figure_too_narrow", "fit_base_length"]
