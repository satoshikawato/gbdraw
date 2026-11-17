"""Horizontal inner labels and the scale's numbers.

With a ticks slot between the feature track and the inner labels (``--separate_strands --track_type spreadout``),
inner labels were drawn over the numbers ("SimC" and "PmgC-like head morphogenesis" over "50 kbp" of NC_005856.1).
The figures are drawn by the command line; the test reads the labels the label drawer is given and the letters of
each number from the tick group the drawing builds, and meets them as the fit check does (each label's fit box
against the letters' rings, ``box_meets_ink``).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from gbdraw import circular as circular_cli_module
from gbdraw.diagrams.circular import assemble
from gbdraw.labels.circular import _label_fit_box
from gbdraw.layout.circular_ink import CircularInkRing, box_meets_ink
from gbdraw.layout.spatial import Aabb
from gbdraw.render.groups.circular.ticks import TickGroup

INPUTS = Path(__file__).resolve().parent / "test_inputs"
P1 = INPUTS / "NC_005856.1.gb"
BGC = INPUTS / "BGC0000708.gbk"
P1_ARGS = ("--separate_strands", "--track_type", "spreadout", "--labels", "both")
BGC_ARGS = ("--definition_font_size", "30", "--species", "<i>Genus species</i>", "--strain", "strain X-1", "--labels", "both")


class _Drawn:
    """The inner labels the drawer was given and the letters of the numbers the tick groups drew."""

    def __init__(self) -> None:
        self.numbers: list[CircularInkRing] = []
        self.inner_labels: list[tuple[dict[str, Any], int]] = []

    def over_numbers(self) -> list[str]:
        """The inner labels whose fit box meets a number's letters."""

        over: list[str] = []
        for label, length in self.inner_labels:
            (min_x, min_y), _, (max_x, max_y), _ = _label_fit_box(label, length)
            if box_meets_ink(Aabb(min_x, min_y, max_x, max_y), self.numbers, threshold_px=0.0):
                over.append(str(label["label_text"]))
        return over


def _draw(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, gbk: Path, args: tuple[str, ...]) -> tuple[str, _Drawn]:
    drawn = _Drawn()
    draw_ticks = assemble.add_tick_group_on_canvas
    draw_labels = assemble.add_labels_group_on_canvas

    def ticks(canvas: Any, gb_record: Any, canvas_config: Any, **kwargs: Any) -> Any:
        options = {key: value for key, value in kwargs.items() if key not in {"group_key", "feature_leaders"}}
        radius = options.pop("radius_override", None)
        group = TickGroup(gb_record, canvas_config.geometry(), profile=canvas_config.profile, radius=radius, **options)
        drawn.numbers.extend(group._label_rings())
        return draw_ticks(canvas, gb_record, canvas_config, **kwargs)

    def labels(canvas: Any, gb_record: Any, canvas_config: Any, **kwargs: Any) -> Any:
        if kwargs["phase"] == "text":
            length = kwargs.get("display_length") or len(gb_record.seq)
            drawn.inner_labels.extend(
                (label, length) for label in kwargs["precalculated_labels"]
                if label.get("is_inner") and not label.get("is_embedded")
            )
        return draw_labels(canvas, gb_record, canvas_config, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(assemble, "add_tick_group_on_canvas", ticks)
        patch.setattr(assemble, "add_labels_group_on_canvas", labels)
        out = tmp_path / "figure"
        circular_cli_module.circular_main(["--gbk", str(gbk), *args, "-o", str(out), "-f", "svg"])
    return (tmp_path / "figure.svg").read_text(encoding="utf-8"), drawn


def _without_tick_numbers(monkeypatch: pytest.MonkeyPatch) -> None:
    """Place the labels as if the label stage were given no scale numbers."""

    radii: Callable[..., Any] = assemble.circular_label_radii
    monkeypatch.setattr(assemble, "circular_label_radii", lambda layout, *_numbers: radii(layout))


@pytest.mark.circular
def test_inner_labels_keep_off_the_scale_numbers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """No inner label of NC_005856.1 (ticks inside the strands) is drawn over a scale number."""

    _, drawn = _draw(tmp_path, monkeypatch, P1, P1_ARGS)

    assert drawn.numbers and len(drawn.inner_labels) > 20
    assert drawn.over_numbers() == []


@pytest.mark.circular
def test_a_record_whose_inner_labels_reach_no_number_is_placed_as_before(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """BGC0000708 with a large definition: its inner labels come near the scale but over no number, so the numbers
    hold none of them and the figure is the one drawn without them; NC_005856.1, whose labels reach the numbers,
    is not."""

    svg, drawn = _draw(tmp_path / "with", monkeypatch, BGC, BGC_ARGS)
    assert drawn.numbers and drawn.inner_labels and drawn.over_numbers() == []
    held_p1, _ = _draw(tmp_path / "p1", monkeypatch, P1, P1_ARGS)

    _without_tick_numbers(monkeypatch)
    unheld, unheld_drawn = _draw(tmp_path / "without", monkeypatch, BGC, BGC_ARGS)
    unheld_p1, unheld_p1_drawn = _draw(tmp_path / "p1-without", monkeypatch, P1, P1_ARGS)

    assert unheld == svg
    assert unheld_drawn.over_numbers() == []
    assert unheld_p1 != held_p1 and unheld_p1_drawn.over_numbers()
