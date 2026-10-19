"""The translucent backing of a label: one geometry for region and feature labels.

A backed label is drawn on a white box with rounded corners at 70 % opacity
(LP-D82 ``RGN-FRONT``), ``BACKING_PADDING_EM`` font sizes beyond its letters
on every side, with the same corner radius. A box may stop at an edge on one
side (a region label's mark, a feature label's own leader end), so it never
covers what it belongs to.

Study prototype (LP-D85): feature labels take the same backing where
``GBDRAW_LABEL_BACKING`` asks for it: ``through`` backs a label that another
label's leader passes through, ``all`` backs every label with a leader.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Literal, TypeAlias

from gbdraw.layout.spatial import Aabb
from gbdraw.labels.label_ink import UNDERLINE_GAP_EM, DiagramMode, FeatureLabelInk, segment_through_box, study_switch
from gbdraw.svg.elements import Rect, Text, Use

BACKING_PADDING_EM = 0.18
BACKING_FILL = "#ffffff"
BACKING_OPACITY = 0.7

FeatureLabelBackingMode: TypeAlias = Literal["off", "through", "all"]
BACKING_MODE_ENV = "GBDRAW_LABEL_BACKING"
BACKING_OPACITY_ENV = "GBDRAW_LABEL_BACKING_OPACITY"
HALO_ENV = "GBDRAW_LABEL_HALO"
HALO_COLOR = "#ffffff"
LABEL_FEATURE_ID = "data-label-feature-id"


def feature_label_halo_em(mode: DiagramMode) -> float | None:
    """The prototype halo (``GBDRAW_LABEL_HALO`` or ``..._LINEAR`` / ``..._CIRCULAR``): the stroke width of the
    white outline behind each feature label's letters, in em (it reaches half of it beyond the glyphs); None
    when unset, zero or not a number."""

    try:
        width = float(study_switch(HALO_ENV, mode) or 0.0)
    except ValueError:
        return None
    return width if width > 0.0 else None


def label_halos(texts: Sequence[Text], width_em: float) -> list[Use]:
    """One halo per feature label text: a ``<use>`` of the text, stroked white with round joins, to be drawn
    behind it. It is a reference, not a second ``<text>``, so the Web label editor (which collects ``text``
    elements) never sees it, and it follows any edit of its text. CairoSVG has no ``paint-order``, so the
    halo cannot be a stroke on the text itself. Gives each text an ``id`` when it has none."""

    halos: list[Use] = []
    seen: dict[str, int] = {}
    for text in texts:
        feature_id = text.attribs.get(LABEL_FEATURE_ID)
        font_size = text.attribs.get("font-size")
        if not feature_id or font_size is None:
            continue
        if "id" not in text.attribs:
            count = seen.get(str(feature_id), 0)
            seen[str(feature_id)] = count + 1
            text.attribs["id"] = f"label_text_{feature_id}" + (f"_{count}" if count else "")
        halo = Use(
            href=f"#{text.attribs['id']}",
            fill="none",  # the clone draws no letters of its own (PDF text is not doubled)
            stroke=HALO_COLOR,
            stroke_width=width_em * float(font_size),
            stroke_linejoin="round",
            stroke_linecap="round",
        )
        halo.attribs["aria-hidden"] = "true"
        halo.attribs["pointer-events"] = "none"
        halos.append(halo)
    return halos


def backing_box(ink: Aabb, font_size: float, *, edge_x: float | None = None, edge_y: float | None = None) -> Aabb:
    """The box a backed label is drawn on: its ``ink`` padded on every side, stopped at ``edge_x`` /
    ``edge_y`` on the side of the letters where that edge lies (past the letters' far side when the
    edge is not before them)."""

    box = ink.expanded(BACKING_PADDING_EM * font_size)
    min_x, min_y, max_x, max_y = box.min_x, box.min_y, box.max_x, box.max_y
    if edge_y is not None:
        if ink.max_y <= edge_y:
            max_y = min(max_y, edge_y)
        else:
            min_y = max(min_y, edge_y)
    if edge_x is not None:
        if ink.max_x <= edge_x:
            max_x = min(max_x, edge_x)
        else:
            min_x = max(min_x, edge_x)
    return Aabb(min_x, min_y, max_x, max_y)


def backing_rect(box: Aabb, font_size: float, *, opacity: float = BACKING_OPACITY) -> Rect:
    """The drawn backing of ``box``."""

    radius = BACKING_PADDING_EM * font_size
    return Rect(
        insert=(box.min_x, box.min_y),
        size=(box.width, box.height),
        rx=radius,
        ry=radius,
        fill=BACKING_FILL,
        fill_opacity=opacity,
        stroke="none",
    )


def feature_label_backing_opacity() -> float:
    """The feature label backing's fill opacity (``GBDRAW_LABEL_BACKING_OPACITY``, LP-D88 study; the
    region label's ``BACKING_OPACITY`` when unset)."""

    try:
        return float(os.environ.get(BACKING_OPACITY_ENV, "") or BACKING_OPACITY)
    except ValueError:
        return BACKING_OPACITY


def feature_label_backing_mode(mode: DiagramMode) -> FeatureLabelBackingMode:
    """The prototype switch (``GBDRAW_LABEL_BACKING``, or ``..._LINEAR`` / ``..._CIRCULAR``); ``off`` when
    unset or unknown."""

    value = study_switch(BACKING_MODE_ENV, mode)
    if value == "through":
        return "through"
    if value == "all":
        return "all"
    return "off"


def feature_label_backing(label: FeatureLabelInk) -> Aabb:
    """A feature label's backing in the frame of its ``ink``, stopped at its own drawn leader's end
    on the side the leader arrives from, or short of its underline when it has one."""

    line = label.underline
    if line is not None:
        # Stop short of the underline by its gap, so the backing never paints over the line's stroke.
        line_y = label.to_label_frame(line[0])[1]
        gap = UNDERLINE_GAP_EM * label.font_size
        return backing_box(label.ink, label.font_size, edge_y=line_y - gap if line_y > label.ink.max_y else line_y + gap)
    drawn = label.drawn
    if not drawn:
        return backing_box(label.ink, label.font_size)
    end_x, end_y = label.to_label_frame(drawn[-1])
    ink = label.ink
    out_x = max(ink.min_x - end_x, end_x - ink.max_x)
    out_y = max(ink.min_y - end_y, end_y - ink.max_y)
    if out_x >= out_y:
        edge = max(end_x, ink.max_x) if end_x > 0.5 * (ink.min_x + ink.max_x) else min(end_x, ink.min_x)
        return backing_box(ink, label.font_size, edge_x=edge)
    edge = max(end_y, ink.max_y) if end_y > 0.5 * (ink.min_y + ink.max_y) else min(end_y, ink.min_y)
    return backing_box(ink, label.font_size, edge_y=edge)


def feature_label_backing_rect(label: FeatureLabelInk) -> Rect:
    rect = backing_rect(feature_label_backing(label), label.font_size, opacity=feature_label_backing_opacity())
    if label.rotation_deg != 0.0:
        rect.rotate(label.rotation_deg, center=label.center)
    return rect


def backing_rects(labels: Sequence[FeatureLabelInk], mode: FeatureLabelBackingMode) -> list[Rect]:
    """The backings of one record's feature labels, to draw before (under) all of their letters."""

    backed = backed_feature_labels(labels, mode)
    return [feature_label_backing_rect(label) for index, label in enumerate(labels) if index in backed]


def leaders_through_text(labels: Sequence[FeatureLabelInk]) -> list[tuple[int, int]]:
    """(leader of, text of) pairs: a label's leader drawn through another label's letters."""

    found: list[tuple[int, int]] = []
    for text_of, label in enumerate(labels):
        for leader_of, other in enumerate(labels):
            if leader_of == text_of or not other.drawn:
                continue
            points = [label.to_label_frame(point) for point in other.drawn]
            if any(segment_through_box(a, b, label.ink) for a, b in zip(points, points[1:])):
                found.append((leader_of, text_of))
    return found


def backed_feature_labels(labels: Sequence[FeatureLabelInk], mode: FeatureLabelBackingMode) -> frozenset[int]:
    """The labels drawn on a backing: none, those another leader passes through, or every one with a leader."""

    if mode == "all":
        return frozenset(index for index, label in enumerate(labels) if label.drawn)
    if mode == "through":
        return frozenset(text_of for _leader_of, text_of in leaders_through_text(labels))
    return frozenset()


__all__ = [
    "BACKING_FILL",
    "BACKING_MODE_ENV",
    "BACKING_OPACITY",
    "BACKING_PADDING_EM",
    "FeatureLabelBackingMode",
    "HALO_ENV",
    "feature_label_halo_em",
    "label_halos",
    "backed_feature_labels",
    "backing_box",
    "backing_rect",
    "backing_rects",
    "feature_label_backing",
    "feature_label_backing_opacity",
    "feature_label_backing_mode",
    "leaders_through_text",
]
