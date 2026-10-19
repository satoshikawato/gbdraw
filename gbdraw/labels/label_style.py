"""How a feature label with a leader is drawn: its line, its protection, and the leader's width.

The one owner of the label style values (LP-D98, LP-D100). Each diagram mode
has a ``FeatureLabelStyle``, read from ``labels.linear`` / ``labels.circular``
(``line``, ``protection``) by the mode's render profile:

- ``line``: ``under`` draws a line under the letters that continues the
  leader (``label_ink``); ``none`` draws none.
- ``protection``: ``halo`` draws a white glyph-shaped outline behind the
  letters (``HALO_WIDTH_EM`` wide), so a leader passing behind them is cut
  away; ``box`` draws a white rounded box (``BOX_OPACITY``) behind them;
  ``none`` draws neither.

The halo width, the box opacity, the gaps and the width rule are constants,
not settings (Owner, LP-D100 Q-B). Leaders and underlines are drawn
``leader_width_px(font size)`` wide in both modes.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Literal, TypeAlias, cast, get_args

from ..exceptions import ValidationError

LabelStyleMode: TypeAlias = Literal["circular", "linear"]
LabelLine: TypeAlias = Literal["none", "under"]
LabelProtection: TypeAlias = Literal["none", "halo", "box"]
LABEL_LINES: tuple[LabelLine, ...] = get_args(LabelLine)
LABEL_PROTECTIONS: tuple[LabelProtection, ...] = get_args(LabelProtection)

# The packaged defaults (config.toml writes the same values; a test compares them): Linear an underline on a
# box, Circular a halo (LP-D106).
DEFAULT_LABEL_LINE: Mapping[LabelStyleMode, LabelLine] = {"linear": "under", "circular": "none"}
DEFAULT_LABEL_PROTECTION: Mapping[LabelStyleMode, LabelProtection] = {"linear": "box", "circular": "halo"}

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
# Leader and underline width: max(LEADER_MIN_PX, LEADER_EM x font size).
LEADER_MIN_PX = 1.0
LEADER_EM = 0.07


@dataclass(frozen=True)
class FeatureLabelStyle:
    """The look of one mode's feature labels with leaders."""

    line: LabelLine
    protection: LabelProtection

    @property
    def underline(self) -> bool:
        return self.line == "under"


def leader_width_px(font_size: float) -> float:
    """The width of a leader and of its underline for a label of ``font_size`` px."""

    return max(LEADER_MIN_PX, LEADER_EM * float(font_size))


def widest_leader_px(labels: Iterable[Mapping[str, object]]) -> float:
    """The widest leader of ``labels`` (each drawn at its own ``font_size``), for readers that keep its
    half width clear; ``LEADER_MIN_PX`` without labels."""

    return max((leader_width_px(float(str(label["font_size"]))) for label in labels), default=LEADER_MIN_PX)


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
    "DEFAULT_LABEL_LINE",
    "DEFAULT_LABEL_PROTECTION",
    "FeatureLabelStyle",
    "GAP_STOP_EM",
    "HALO_COLOR",
    "HALO_WIDTH_EM",
    "LABEL_LINES",
    "LABEL_PROTECTIONS",
    "LEADER_EM",
    "LEADER_MIN_PX",
    "LabelLine",
    "LabelProtection",
    "LabelStyleMode",
    "UNDERLINE_GAP_EM",
    "ROW_LABEL_GAP_EM",
    "UNDERLINE_TRIM_EM",
    "leader_width_px",
    "normalize_label_line",
    "normalize_label_protection",
    "widest_leader_px",
]
