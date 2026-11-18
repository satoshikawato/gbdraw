"""The paint of a numeric series track and the style of its axis, as values (internal).

They sit apart from the series geometry (``scalar_axis``) so the track configurators, which give them, and the
Legend, which reads the paint, import nothing that imports them back.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, get_args

if TYPE_CHECKING:
    from gbdraw.config.models.objects import ScalarAxisLabelStyle

    from .scalar_axis import ScalarAxisSource


# How a numeric series draws its values: the area between its line and its baseline, or the line alone.
SeriesStyle = Literal["area", "line"]
SERIES_STYLES: tuple[SeriesStyle, ...] = get_args(SeriesStyle)
# Where a series' area closes: ``auto`` (0 when the axis holds it, else the axis' lower end), the series' mean, or a
# value (held within the axis).
SeriesBaseline = Literal["auto", "mean"] | float


@dataclass(frozen=True)
class SeriesPaint:
    """How a numeric series track (GC content, GC skew, depth) is painted: the drawers of both modes paint its path
    with it and its Legend rows show it. ``below_fill`` and ``below_stroke`` are the fill and the stroke of the part
    below the baseline (GC skew; a depth track's negative colour), drawn by a copy of the path clipped to the part
    below the baseline; none for a series of one paint. Each track's configurator gives it (``series_paint``)."""

    fill: str
    stroke: str
    stroke_width: float
    fill_opacity: float
    below_fill: str | None
    below_stroke: str | None

    @property
    def writes_stroke_width(self) -> bool:
        """Whether the drawers of both modes write the path's ``stroke-width``: not for a stroke of width 0 in colour
        ``none``, which paints nothing either way."""

        return self.stroke_width > 0 or self.stroke.strip().lower() != "none"

    @property
    def low_fill(self) -> str:
        """The fill of the part below the baseline: ``below_fill``, else the one fill."""

        return self.fill if self.below_fill is None else self.below_fill

    @property
    def low_stroke(self) -> str:
        """The stroke of the part below the baseline: ``below_stroke``, else the one stroke."""

        return self.stroke if self.below_stroke is None else self.below_stroke

    @property
    def paints_below(self) -> bool:
        """Whether the part below the baseline has a paint of its own, drawn by a clipped copy of the path."""

        return self.below_fill is not None or self.below_stroke is not None


@dataclass(frozen=True)
class ScalarAxisStyle:
    """A scalar track's axis and ticks as a value: whether it draws them, their font size and their labels. The
    drawn axis carries it (``ScalarAxis.style``), and the room kept for the axis reads it (the Circular radial
    footprint, the Linear vertical footprint), so what is drawn and what is reserved follow one rule
    (``draws_numbered_labels``). A render stage that replans the radial layout takes this instead of the track's
    configurator, so a window or scale change does not rebuild it."""

    show_axis: bool
    show_ticks: bool
    tick_font_size: float | None
    axis_labels: ScalarAxisLabelStyle

    @classmethod
    def of(cls, source: ScalarAxisSource) -> ScalarAxisStyle:
        return cls(
            show_axis=bool(source.show_axis),
            show_ticks=bool(source.show_ticks),
            tick_font_size=source.tick_font_size,
            axis_labels=source.axis_labels,
        )

    @property
    def draws_numbered_labels(self) -> bool:
        """Whether the axis numbers its large ticks: their labels sit centred on the ticks, so the outermost reach
        half a label line past the band. Without them a range label sits inside the band."""

        return self.show_axis and self.show_ticks and self.axis_labels == "numbered"


__all__ = ["SERIES_STYLES", "ScalarAxisStyle", "SeriesBaseline", "SeriesPaint", "SeriesStyle"]
