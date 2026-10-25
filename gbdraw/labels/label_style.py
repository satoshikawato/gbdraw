"""How a feature label with a leader is drawn: its line, its protection, and the leader's width.

The one owner of the label style values (LP-D98, LP-D100). Each diagram mode
has a ``FeatureLabelStyle``, read from ``labels.linear`` / ``labels.circular``
(``line``, ``protection``) by the mode's render profile:

- ``line``: ``under`` draws a line under the letters of a horizontal label
  that continues the leader (``label_ink``); a label drawn turned gets none
  and its leader ends at its letters, as with ``none`` (Owner, LP-D112).
  ``FeatureLabelStyle.draws_underline`` is the one decision every drawer and
  obstacle reads, through the label's ink. ``none`` draws none.
- ``protection``: ``halo`` draws a white glyph-shaped outline behind the
  letters (``HALO_WIDTH_EM`` wide), so a leader passing behind them is cut
  away; ``box`` draws a white rounded box (``BOX_OPACITY``) behind them;
  ``none`` draws neither.

The halo width, the box opacity, the gaps and the width rule are constants,
not settings (Owner, LP-D100 Q-B). Leaders and underlines are drawn
``leader_width_px(font size)`` wide in both modes.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final, Literal, TypeAlias, cast, get_args

from ..exceptions import ValidationError

LabelStyleMode: TypeAlias = Literal["circular", "linear"]
LabelLine: TypeAlias = Literal["none", "under"]
LabelProtection: TypeAlias = Literal["none", "halo", "box"]
LABEL_LINES: tuple[LabelLine, ...] = get_args(LabelLine)
LABEL_PROTECTIONS: tuple[LabelProtection, ...] = get_args(LabelProtection)

# The packaged defaults (config.toml writes the same values; a test compares them): Linear an underline on a
# box, Circular neither a line nor a protection; a halo can be chosen (LP-D110).
DEFAULT_LABEL_LINE: Mapping[LabelStyleMode, LabelLine] = {"linear": "under", "circular": "none"}
DEFAULT_LABEL_PROTECTION: Mapping[LabelStyleMode, LabelProtection] = {"linear": "box", "circular": "none"}

# The halo: a white stroke of the letters this many font sizes wide (half of it beyond the glyphs).
HALO_WIDTH_EM = 0.08
HALO_COLOR = "#ffffff"
# The box: the region label backing's white rounded box, at this opacity (both modes, LP-D106).
BOX_OPACITY = 0.8
# The underline lies this many font sizes below the font's descent line, measured to its edge.
UNDERLINE_GAP_EM = 0.06
# Each end of the underline stops this many font sizes inside the letters' ink, at most a quarter of the
# line, so the lines of neighbours in one row stay apart (LP-D103).
UNDERLINE_TRIM_EM = 0.12
# A leader that would reach the underline through the letters stops this many font sizes above them.
GAP_STOP_EM = 0.20
# Neighbouring Linear feature labels in one label row keep at least this many font sizes between them
# (each label half of it on each side, of its own font size), so two names never read as one word.
ROW_LABEL_GAP_EM = 0.25
# Feature labels are laid out in the box their text is drawn in (``core.text.drawn_text_box``: the drawn
# width, the face's line box) plus these margins, in font sizes (OV-165).
# Circular: extra box width, the room kept between neighbours along the rim (the embed test reads it too).
CIRCULAR_LABEL_GAP_EM: Final = 0.30
# Linear: extra box height, the room kept between label rows stacked above or below a track.
LINEAR_LABEL_HEIGHT_MARGIN_EM: Final = 0.126
# Leader and underline width: max(LEADER_MIN_PX, LEADER_EM x font size).
LEADER_MIN_PX = 1.0
LEADER_EM = 0.07


@dataclass(frozen=True)
class FeatureLabelStyle:
    """The look of one mode's feature labels with leaders."""

    line: LabelLine
    protection: LabelProtection

    def draws_underline(self, turned: bool) -> bool:
        """Whether a feature label gets this style's underline: a turned one (a Linear label given a rotation,
        a Circular radial label, whatever its angle) does not (LP-D112)."""

        return self.line == "under" and not turned

    def draws_box(self, *, has_leader: bool, beside_feature: bool) -> bool:
        """Whether a feature label gets this style's white box (``protection = box``): a label with a leader,
        and a Linear ``above_feature`` label drawn on its feature without one (``beside_feature``), so that
        another label's leader passing through its letters is faded there (OV-1361, LP-D119). A label drawn
        inside its feature has neither and gets no box."""

        return self.protection == "box" and (has_leader or beside_feature)

    @property
    def draws_from_ink(self) -> bool:
        """Whether this style draws a part placed from a label's letters (an underline or a box)."""

        return self.line == "under" or self.protection == "box"


def leader_width_px(font_size: float) -> float:
    """The width of a leader and of its underline for a label of ``font_size`` px."""

    return max(LEADER_MIN_PX, LEADER_EM * float(font_size))


def label_leader_width_px(label: Mapping[str, object]) -> float:
    """The width of one prepared feature label's leader and underline, from its own ``font_size``: what the
    drawing strokes and what an obstacle or a bound keeps clear of its leader."""

    return leader_width_px(float(str(label["font_size"])))


def _normalized(value: object, choices: tuple[str, ...], path: str) -> str:
    text = value.strip().lower() if isinstance(value, str) else None
    if text not in choices:
        raise ValidationError(
            f"{path} must be one of: {', '.join(choices)}",
            diagnostic={"code": "INPUT_INVALID", "configPath": path},
        )
    return text


def normalize_label_line(value: object, *, path: str) -> LabelLine:
    """``value`` as a label line (``none`` or ``under``); ``path`` names the setting in the error."""

    return cast(LabelLine, _normalized(value, LABEL_LINES, path))


def normalize_label_protection(value: object, *, path: str) -> LabelProtection:
    """``value`` as a label protection (``none``, ``halo`` or ``box``); ``path`` names the setting."""

    return cast(LabelProtection, _normalized(value, LABEL_PROTECTIONS, path))


__all__ = [
    "BOX_OPACITY",
    "CIRCULAR_LABEL_GAP_EM",
    "DEFAULT_LABEL_LINE",
    "DEFAULT_LABEL_PROTECTION",
    "FeatureLabelStyle",
    "GAP_STOP_EM",
    "HALO_COLOR",
    "HALO_WIDTH_EM",
    "LABEL_LINES",
    "LABEL_PROTECTIONS",
    "LINEAR_LABEL_HEIGHT_MARGIN_EM",
    "LEADER_EM",
    "LEADER_MIN_PX",
    "LabelLine",
    "LabelProtection",
    "LabelStyleMode",
    "UNDERLINE_GAP_EM",
    "ROW_LABEL_GAP_EM",
    "UNDERLINE_TRIM_EM",
    "label_leader_width_px",
    "leader_width_px",
    "normalize_label_line",
    "normalize_label_protection",
]
