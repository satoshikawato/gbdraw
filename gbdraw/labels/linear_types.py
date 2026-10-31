"""The Linear feature label record that placement, the fit count, region labels and the renderer share."""

from __future__ import annotations

from typing import TypedDict


class _LinearLabelRequired(TypedDict):
    """Keys every Linear feature label has from ``linear_record_labels`` on."""

    label_text: str
    source_feature_index: int | None
    # The text anchor's x and the label box's x span (record coordinates); ``middle`` is where an
    # external label's leader starts.
    middle: float
    start: float
    end: float
    middle_x: float
    middle_y: float
    width_px: float
    height_px: float
    strand: str | None
    # A feature placed below the axis in a combined strand pool: its external labels stack below.
    secondary_below: bool
    feature_middle_y: float
    feature_top_y: float
    feature_bottom_y: float
    feature_start_x: float
    feature_end_x: float
    feature_anchor_x: float
    font_size: float
    font_family: str
    rotation_deg: float
    text_anchor: str
    is_embedded: bool


class LinearLabel(_LinearLabelRequired, total=False):
    """One Linear feature label, in its record's coordinates (x from the record's start, y from its
    features slot's origin).

    Embedded, external and ``above_feature`` labels travel in one list, so the keys below are
    optional for every kind:

    - ``track_id``: the label's row (``track_0`` on its feature); an external label gets it when its
      row is placed (``place_linear_label_row``);
    - ``above_feature`` labels only (``is_above_feature_label``): ``above_feature_place_above``,
      ``label_contact_x_offset``, ``label_contact_y_offset``;
    - a leader drawn from ``leader_start_*`` to ``leader_end_*`` when ``leader_line`` is true
      (``above_feature`` labels; an external label's leader comes from its row);
    - an external label's rows: ``row_origin_y`` and ``row_pitch`` (``linear_label_side``);
    - ``feature_id``: the rendered feature the drawn label binds to, added by the record renderer;
    - ``wrap_row``: the row of a wrapped record the label stands on (``display_label_segment``).
    """

    track_id: str
    above_feature_place_above: bool
    label_contact_x_offset: float
    label_contact_y_offset: float
    leader_line: bool
    leader_start_x: float
    leader_start_y: float
    leader_end_x: float
    leader_end_y: float
    row_origin_y: float
    row_pitch: float
    feature_id: str
    wrap_row: int


def is_above_feature_label(label: LinearLabel) -> bool:
    """Whether ``label`` was placed by ``label_placement = above_feature`` (beside its feature, in that
    mode's rows) rather than by the auto placement."""

    return "above_feature_place_above" in label


def is_level_rotation(rotation_deg: float) -> bool:
    """Whether a Linear label turned by ``rotation_deg`` reads level: no turn or whole turns (360, -360,
    720 ...), so it is placed, sized and drawn as a horizontal label, underline included."""

    return rotation_deg % 360.0 == 0.0


def is_horizontal_label(label: LinearLabel) -> bool:
    """Whether ``label`` is drawn level (``is_level_rotation``): the box, row pitch and ink tests that read
    only its width and height apply."""

    return is_level_rotation(label["rotation_deg"])


__all__ = ["LinearLabel", "is_above_feature_label", "is_horizontal_label", "is_level_rotation"]
