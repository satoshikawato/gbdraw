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
from typing import Generic, Literal, TypeAlias, TypeVar, get_args

from gbdraw.config.models.canvas import FIGURE_HEIGHT_PATH, FIGURE_WIDTH_PATH
from gbdraw.core.figure_length import FigureLength, RequestedFigureSize
from gbdraw.exceptions import ValidationError
from gbdraw.layout.composition import CompositionPlacement, CompositionPlan
from gbdraw.layout.spatial import Aabb, union_aabbs

Built = TypeVar("Built")
# Each step aims this far under the target, and a pass this close under it is accepted: the planner centres the
# rest in the requested canvas, so the canvas stays exact. A pass over the target is never accepted.
_AIM_UNDER_PX = 2.0
_ACCEPT_UNDER_PX = 3.0
_MAX_PASSES = 8
# A base too small to draw fails early (a Circular map's track layout, about 0.05 s, where a drawn pass with labels
# can take 4-20 s), and a pass that left content out skips its costly part (a Circular map whose labels were not
# placed, about 0.05 s), so such passes have their own bound, after which the search stops: stepping a quarter of the
# way from the largest of them below the smallest drawn base to that base reaches it from a few hundred.
_MAX_CHEAP_PASSES = 24
# Below every pass, where the size reaches the box only under the smallest base that can be drawn, the search tries
# this fraction of the way up from the largest base too small to draw: a try that cannot be drawn is cheap, so the
# smallest drawn pass, which a refusal states, is found with fewer costly passes than by halving.
_CHEAP_SIDE_SPLIT = 4
# The composed size grows with the base length, but labels past the ends can make a stretch nearly flat; a slope
# measured there is raised to this so one step cannot jump far.
_MIN_SLOPE = 0.25
# The search assumes that a figure's content shrinks at most this many times the mode's first slope per px of base
# (Circular labels sit at up to twice the radius at the outward rows' limit). A search past a pass that left content
# out ends only when even this slope reaches the box at or below that pass. Content that jumps (a legend that moves to
# another side) can shrink faster; the search can then end before a base that keeps the content and fits.
_STEEPEST_SLOPE_FACTOR = 2.0
# Before a refusal, a base between the smallest pass and the next one is measured, since a size that is not monotonic
# may fit there; not when even this many times the steepest slope from either pass cannot reach the box there.
_MIDDLE_DIP_FACTOR = 4.0
# Two passes whose parts (``MeasuredBase.parts_px``) differ by no more than this have the same size. The parts that do
# not shrink are added in another order at another base, so the same size can differ in its last digits (1e-13 px
# measured).
_SAME_SIZE_PX = 1e-6

# The primary item a refusal names, by the id its diagnostic carries (``primaryItem``): the Linear records at their
# smallest (``record_ends``: what stays of them as the sequence shrinks), the Linear records, the Circular map, or a
# multi-record grid of maps. The Web words each (utils/error-normalization.js).
FigurePrimaryItem: TypeAlias = Literal["record_ends", "records", "map", "maps"]
FIGURE_PRIMARY_ITEMS: tuple[FigurePrimaryItem, ...] = get_args(FigurePrimaryItem)
_PRIMARY_NAMES: dict[FigurePrimaryItem, str] = {
    "record_ends": "record names and labels past the sequence ends",
    "records": "records",
    "map": "map at its smallest",
    "maps": "maps at their smallest",
}
# A part a refusal names, by its diagnostic key: the primary item, the legend, the title, the padding.
FigureSizePart: TypeAlias = Literal["primary", "legend", "title", "padding"]
FIGURE_SIZE_PARTS: tuple[FigureSizePart, ...] = get_args(FigureSizePart)
# The refusals of a requested size, by their diagnostic reason (the Web words each).
FIGURE_SIZE_REASONS = ("FIGURE_TOO_NARROW", "FIGURE_TOO_SHORT")


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
class MeasuredBase(Generic[Built]):
    """What a mode's ``measure`` drew at one base length: how far its composed content is over the requested box
    (negative when under) and what it assembled.

    ``parts_px`` is what a refusal at this base would state (``figure_size_parts``): two passes with the same parts
    state the same refusal. ``left_out`` marks content drawn without a part that does not fit at this base, which no
    smaller base fits either (Circular labels not placed for lack of room, CAPACITY-STOP): such a pass is under the
    box only because of what it left out, so its excess says nothing about the size of the figure with that part.
    """

    excess_px: float
    built: Built
    parts_px: tuple[float, ...]
    left_out: bool = False


@dataclass(frozen=True)
class _Pass(Generic[Built]):
    """A drawn pass, whether it left content out, and the Python warnings it raised, which reach the caller only if
    the pass is kept."""

    fit: BaseLengthFit[Built]
    parts_px: tuple[float, ...]
    left_out: bool
    caught: tuple[warnings.WarningMessage, ...]


def _measured(
    measure: Callable[[int], MeasuredBase[Built] | None], base: int
) -> tuple[MeasuredBase[Built] | None, tuple[warnings.WarningMessage, ...]]:
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


def _flat_at_smallest(drawn: dict[int, _Pass[Built]]) -> bool:
    """Whether the two smallest drawn passes have the same parts, so a refusal at either states the same sizes: the
    content no longer shrinks with the base there, because the parts that do not shrink set it (a Circular map whose
    centre definition keeps the tracks out, D-13). The excess alone can stay the same while the map still shrinks,
    when another part (a legend beside the map) sets the side that is over.

    The search assumes that such content shrinks no further below: it is the larger of the parts that do not shrink
    and the parts that grow with the base. Content that shrinks again below a stretch of the same size (a legend that
    moves to another side) could fit at a smaller base that the search no longer tries.
    """

    smallest = [drawn[base].parts_px for base in sorted(drawn)[:2]]
    return (
        len(smallest) == 2
        and len(smallest[0]) == len(smallest[1])
        and all(abs(a - b) <= _SAME_SIZE_PX for a, b in zip(*smallest))
    )


def _dip_reaches_box(lower: BaseLengthFit[Built], upper: BaseLengthFit[Built], base: int, first_slope: float) -> bool:
    """Whether a size that dips between two passes over the box could reach it at ``base``: from each pass, a line
    ``_MIDDLE_DIP_FACTOR`` times as steep as the steepest assumed slope reaches 0 px at or before ``base``."""

    steepest = _MIDDLE_DIP_FACTOR * _STEEPEST_SLOPE_FACTOR * first_slope
    return max(lower.excess_px - steepest * (base - lower.base), upper.excess_px - steepest * (upper.base - base)) <= 0.0


def _costly_passes(drawn: dict[int, _Pass[Built]]) -> int:
    """The drawn passes that kept their content, which count against ``_MAX_PASSES``."""

    return sum(not drawn_pass.left_out for drawn_pass in drawn.values())


def _next_base(
    drawn: dict[int, _Pass[Built]], order: list[int], undrawable: set[int], first_slope: float
) -> int | None:
    """The next base to measure, or None when no unmeasured base can come closer to the target."""

    costly = _costly_passes(drawn)
    if costly >= _MAX_PASSES or len(undrawable) + len(drawn) - costly >= _MAX_CHEAP_PASSES:
        return None
    passes = [drawn_pass.fit for drawn_pass in drawn.values()]
    under = max((fit for fit in passes if fit.fits), key=lambda fit: fit.base, default=None)
    over = min(
        (fit for fit in passes if not fit.fits and (under is None or fit.base > under.base)),
        key=lambda fit: fit.base,
        default=None,
    )
    ceiling = min((base for base in undrawable if under is not None and base > under.base), default=None)
    if under is not None and ceiling is not None and (over is None or ceiling < over.base):
        # A base above the widest pass under that cannot be drawn bounds the search as a pass over the box would.
        return None if ceiling - under.base <= 1 else (under.base + ceiling) // 2
    if under is not None and over is not None:
        if over.base - under.base <= 1:
            return None
        above = [base for base in sorted(drawn) if base > over.base and not drawn[base].left_out]
        if drawn[under.base].left_out and not drawn[over.base].left_out and above:
            # The pass under left content out, and so does every base below it, so its excess says nothing about the
            # figure with that content: the line through the pass over and the next pass above that kept it does.
            # The search ends once even the steepest line from the pass over reaches the box (0 px) only at or below
            # the pass under. Otherwise it tries the line's step or where the line reaches the box, the first above
            # the pass under, else the middle of the pass under and where the steepest line reaches the box.
            slope = max(_MIN_SLOPE, _slope(over, drawn[above[0]].fit))
            steepest = max(slope, _STEEPEST_SLOPE_FACTOR * first_slope)
            reach, farthest = (over.base + math.floor(-over.excess_px / line) for line in (slope, steepest))
            if farthest <= under.base:
                return None
            middle = (under.base + farthest + 1) // 2
            return next((base for base in (_step(over, slope), reach) if base > under.base), middle)
        # The secant between the bracket's ends, kept strictly inside it: one rounded down onto the pass under still
        # has a base above it to try.
        return min(max(_step(under, _slope(under, over)), under.base + 1), over.base - 1)
    if under is None and _flat_at_smallest(drawn):
        # Every pass is over the box and the content stopped shrinking: no smaller base fits, and the smallest pass
        # has the size of the figure at its smallest.
        return None
    last = drawn[order[-1]].fit
    guess = _step(last, first_slope if len(order) < 2 else _slope(drawn[order[-2]].fit, last))
    if under is not None:
        return max(guess, under.base + 1)
    lowest = min(drawn)
    floor = max((base for base in undrawable if base < lowest), default=0)
    guess = max(1, min(guess, lowest - 1))
    if guess <= floor:
        # Bases at or below the floor cannot be drawn: step a quarter of the way up from it instead.
        guess = floor + max(1, (lowest - floor) // _CHEAP_SIDE_SPLIT)
    return None if guess <= floor or guess in drawn else guess


def fit_base_length(
    measure: Callable[[int], MeasuredBase[Built] | None], *, first_base: int, slope: float
) -> BaseLengthFit[Built]:
    """The base length whose composed content is within a few px under the requested box.

    ``measure(base)`` returns how far the content drawn at ``base`` is over the box (``excess_px``, negative when
    under), what it assembled, and whether it left content out (``MeasuredBase``), or None when ``base`` is too small
    to draw. ``first_base`` must be drawable; the mode supplies it and ``slope``, the px of content per px of base of
    the first step. Each step aims 2 px under: later ones with the secant of the last two drawn passes. Once one pass
    is under and one over, the next base is the secant between them, kept strictly inside, until no base lies
    between them; then the widest pass under is kept. When the widest pass under left content out (``left_out``) and
    the pass over and another pass above it kept it, the line through those two sets the step instead, and the
    search ends once a line from the pass over twice as steep as ``slope`` (or the measured one, if steeper), the
    steepest the content is assumed to shrink, reaches the box only at or below the pass under. A base too small to draw above the widest pass under bounds the search
    as a pass over would (its midpoint is tried); below the smallest drawn pass, a base at or below the largest one
    too small to draw is not tried: the step goes a quarter of the way up from it instead, since a base too small to
    draw is cheap to try. After 8 drawn passes that kept their content, or 24 that were too small to draw or left
    content out, the widest pass under is kept.

    With no pass under, the target is refused (``fits`` False) with the smallest drawn pass. The search then ends as
    soon as the two smallest drawn passes have the same parts (``parts_px``, within float rounding): the content
    stopped shrinking, so no smaller base fits, and the refusal states the sizes of the figure at its smallest. Otherwise, before
    refusing, the middle of the smallest pass and the next drawn pass above it is measured once, if not yet measured,
    since a measure that is not monotonic may fit there, unless even a dip four times as steep as the steepest assumed
    slope from either pass could not reach the box there. Only the kept pass's Python warnings reach the caller, once;
    a refusal emits none.
    """

    drawn: dict[int, _Pass[Built]] = {}
    order: list[int] = []
    undrawable: set[int] = set()
    checked_middle = False
    base: int | None = first_base
    while base is not None:
        result, caught = _measured(measure, base)
        if result is None:
            if not drawn:
                raise ValueError(f"the first base length {base} cannot be drawn")
            undrawable.add(base)
        else:
            excess = result.excess_px
            current = _Pass(
                BaseLengthFit(base, excess, result.built, fits=excess <= 0.0), result.parts_px, result.left_out, caught
            )
            if current.fit.fits and excess >= -_ACCEPT_UNDER_PX:
                return _kept(current)
            drawn[base] = current
            order.append(base)
        base = _next_base(drawn, order, undrawable, slope)
        if base is None and not any(drawn_pass.fit.fits for drawn_pass in drawn.values()):
            lowest = min(drawn)
            above = min((other for other in drawn if other > lowest), default=None)
            middle = None if above is None else (lowest + above) // 2
            if (
                not checked_middle
                and not _flat_at_smallest(drawn)
                and _costly_passes(drawn) < _MAX_PASSES
                and middle is not None
                and lowest < middle
                and middle not in undrawable
                and above is not None
                and _dip_reaches_box(drawn[lowest].fit, drawn[above].fit, middle, slope)
            ):
                checked_middle = True
                base = middle
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


def _parts_context(
    length: FigureLength, parts: list[tuple[FigureSizePart, float]], primary: FigurePrimaryItem
) -> tuple[str, dict[str, str]]:
    """The parts of a refusal as its message words them (``records 120.3 mm, legend 40 mm``) and as its diagnostic
    carries them (each part's length by its key), the parts of zero length left out."""

    shown = [(part, length.message_text(px)) for part, px in parts if px > 0.0]
    text = ", ".join(f"{_PRIMARY_NAMES[primary] if part == 'primary' else part} {value}" for part, value in shown)
    return text, dict(shown)


def figure_too_short(
    plan: CompositionPlan,
    height: FigureLength,
    width: FigureLength | None,
    *,
    primary: FigurePrimaryItem,
    legend_role: str = "legend",
    title_role: str = "title",
) -> ValidationError:
    """The refusal of a requested height shorter than the padded content of ``plan``.

    ``primary`` names the primary item (Linear records, a Circular map at its smallest); ``width`` is the
    requested width, if any. The diagnostic carries each length of the message as the message writes it
    (``FIGURE_TOO_SHORT``), so the Web words the same numbers.
    """

    content = plan.content_bounds
    shown_width = width if width is not None else height.of_px(plan.width)
    parts = _axis_parts(plan, "y", legend_role=legend_role, title_role=title_role)
    minimum = height.message_text(content.height, round_up=True)
    at_width = shown_width.message_text(plan.width)
    detail, lengths = _parts_context(
        height,
        [("primary", parts.primary), ("legend", parts.legend), ("title", parts.title), ("padding", parts.padding)],
        primary,
    )
    return ValidationError(
        f"figure_height {height.text} is shorter than the figure at {at_width} "
        f"wide: {minimum} ({detail}). Use a figure_height of at least {minimum}, or leave figure_height auto.",
        diagnostic={
            "code": "INPUT_INVALID",
            "reason": "FIGURE_TOO_SHORT",
            "configPath": FIGURE_HEIGHT_PATH,
            "minimumPx": math.ceil(content.height),
            "requested": height.spaced_text,
            "minimum": minimum,
            "atWidth": at_width,
            "primaryItem": primary,
            **lengths,
        },
    )


def check_figure_height(plan: CompositionPlan, size: RequestedFigureSize, *, primary: FigurePrimaryItem) -> None:
    """Refuse a drawn figure whose padded content is taller than a requested height.

    A requested width is the base-length fit's: the figure drawn is the one that fit it.
    """

    if size.height is not None and plan.content_bounds.height > size.height.px:
        raise figure_too_short(plan, size.height, size.width, primary=primary)


def figure_size_parts(plan: CompositionPlan) -> tuple[float, ...]:
    """What a refusal of a requested size states about ``plan``: the primary item, the legend, the title, and the
    padding along each axis (``_axis_parts``, as the refusals read them), so a mode's measure reports them
    (``MeasuredBase.parts_px``)."""

    return tuple(
        part
        for axis in ("x", "y")
        for parts in (_axis_parts(plan, axis, legend_role="legend", title_role="title"),)
        for part in (parts.primary, parts.legend, parts.title, parts.padding)
    )


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


def figure_size_refusal(
    plan: CompositionPlan, size: RequestedFigureSize, *, primary: FigurePrimaryItem
) -> ValidationError:
    """The refusal of a requested box that ``plan``, the smallest figure the fit drew, does not fit: the width
    refusal when it is wider than a requested width, else the height refusal."""

    if size.height is None or (size.width is not None and plan.content_bounds.width > size.width.px):
        if size.width is None:
            raise ValueError("a refusal needs a requested side")
        return figure_too_narrow(plan, size.width, primary=primary)
    return figure_too_short(plan, size.height, size.width, primary=primary)


# A requested side of a larger figure: about how long it must be, as a message writes it (``240.3 mm``), when the room
# the labels need is known (CAPACITY-STOP); True when it is not (NO-ROW-WAYS, PD-OI-114 r3).
LargerSide: TypeAlias = str | Literal[True]


@dataclass(frozen=True)
class LargerFigure:
    """The sides of a requested figure to make larger for labels that do not fit: each requested side as a
    ``LargerSide``, None for a side left auto. The CLI writes an estimate as an argument without the space
    (``240.3mm``). Both requested sides have an estimate, or neither has."""

    width: LargerSide | None
    height: LargerSide | None

    def __post_init__(self) -> None:
        sides = [side for side in (self.width, self.height) if side is not None]
        if not sides or len({side is True for side in sides}) > 1:
            raise ValueError("a larger figure names a requested side and estimates every requested side or none")


def larger_figure(plan: CompositionPlan, size: RequestedFigureSize, ratio: float | None) -> LargerFigure:
    """The requested sides to make larger, from the drawn figure ``plan`` and the room the labels need as a multiple
    of the room around the map (``ratio``): each requested side grows by ``ratio - 1`` times the map's side, rounded
    up. Without a ratio (no row of the label settings fits, so the room is not known) each requested side has no
    estimate."""

    def side(length: FigureLength | None, content: float, primary: float) -> LargerSide | None:
        if length is None:
            return None
        if ratio is None:
            return True
        return length.message_text(content + (ratio - 1.0) * primary, round_up=True)

    content, primary = plan.content_bounds, plan.primary_bounds
    return LargerFigure(
        side(size.width, content.width, primary.width), side(size.height, content.height, primary.height)
    )


def figure_too_narrow(
    plan: CompositionPlan,
    width: FigureLength,
    *,
    primary: FigurePrimaryItem,
    legend_role: str = "legend",
    title_role: str = "title",
) -> ValidationError:
    """The refusal of a requested width narrower than ``plan``, the smallest figure the fit drew.

    The narrowest figure is that plan's padded content. When the primary row (the primary item with a legend
    beside it) is its widest item, the message names ``primary`` (the primary item at its smallest), the
    legend beside it, and the padding; otherwise the widest item above or below (the title or the legend) and the
    padding. The diagnostic carries each length of the message as the message writes it (``FIGURE_TOO_NARROW``).
    """

    content = plan.content_bounds.width
    legend = plan.placement_for(legend_role)
    title = plan.placement_for(title_role)
    beside = _beside_primary(plan, legend_role)
    row_boxes = [plan.primary_bounds]
    if legend is not None and beside:
        row_boxes.append(legend.final_bounds)
    primary_px = _span([plan.primary_bounds], "x")
    row = _span(row_boxes, "x")
    candidates: tuple[tuple[FigureSizePart, CompositionPlacement | None], ...] = (
        ("title", title),
        ("legend", None if beside else legend),
    )
    above_or_below = [(name, placement.final_bounds) for name, placement in candidates if placement is not None]
    painted = _span([*row_boxes, *(bounds for _name, bounds in above_or_below)], "x")
    row_binds = painted - row < 0.5
    parts: list[tuple[FigureSizePart, float]]
    if row_binds:
        parts = [("primary", primary_px), ("legend", row - primary_px)]
    else:
        widest = max(above_or_below, key=lambda item: item[1].width)[0]
        parts = [(widest, painted)]
    minimum = width.message_text(content, round_up=True)
    detail, lengths = _parts_context(width, [*parts, ("padding", content - painted)], primary)
    legend_remedy = row_binds and beside
    return ValidationError(
        f"figure_width {width.text} is narrower than the parts of the figure that do not shrink: {detail} "
        f"({minimum} in all). Use a figure_width of at least {minimum},"
        f"{' place the legend at the top or bottom,' if legend_remedy else ''} or leave figure_width auto.",
        diagnostic={
            "code": "INPUT_INVALID",
            "reason": "FIGURE_TOO_NARROW",
            "configPath": FIGURE_WIDTH_PATH,
            "minimumPx": math.ceil(content),
            "requested": width.spaced_text,
            "minimum": minimum,
            "primaryItem": primary,
            **lengths,
            **({"legendRemedy": True} if legend_remedy else {}),
        },
    )


__all__ = [
    "FIGURE_PRIMARY_ITEMS",
    "FIGURE_SIZE_PARTS",
    "FIGURE_SIZE_REASONS",
    "BaseLengthFit",
    "FigurePrimaryItem",
    "LargerFigure",
    "LargerSide",
    "MeasuredBase",
    "check_figure_height",
    "figure_size_excess",
    "figure_size_parts",
    "figure_size_refusal",
    "figure_too_narrow",
    "figure_too_short",
    "fit_base_length",
    "larger_figure",
]
