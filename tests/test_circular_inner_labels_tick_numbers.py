"""Horizontal inner labels and the scale's numbers.

With a ticks slot between the feature track and the inner labels (``--separate_strands --track_type spreadout``),
inner labels were drawn over the numbers ("SimC" and "PmgC-like head morphogenesis" over "50 kbp" of NC_005856.1).
The figures are drawn by the command line; the tests read the labels the label drawer is given, the letters of the
numbers the drawn tick groups have (``TickSlotLayout.numbers``) and what the tracks paint (``_record_track_ink``, as
the definition-over-track check reads it). Whether a label is over a number is the placement's own test
(``labels_over_tick_numbers``); whether a label is over a track's paint is checked here on its own, each label's fit
box against each piece of paint (``box_meets_ink``).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import pytest

from gbdraw import circular as circular_cli_module
from gbdraw.diagrams.circular import assemble, builders
from gbdraw.labels.circular import _label_fit_box, labels_over_tick_numbers
from gbdraw.labels.circular_types import CircularLabel
from gbdraw.layout.circular_ink import CircularInk, CircularInkRing, box_meets_ink
from gbdraw.layout.spatial import Aabb

INPUTS = Path(__file__).resolve().parent / "test_inputs"
P1 = INPUTS / "NC_005856.1.gb"
BGC = INPUTS / "BGC0000708.gbk"
P1_ARGS = ("--separate_strands", "--track_type", "spreadout", "--labels", "both")
BGC_ARGS = ("--definition_font_size", "30", "--species", "<i>Genus species</i>", "--strain", "strain X-1", "--labels", "both")
# GC content and GC skew inside the scale, under the inner labels.
INSIDE_TRACKS = (
    "--circular_track_slot", "features:features",
    "--circular_track_slot", "ticks:ticks@side=inside",
    "--circular_track_slot", "gc_content:dinucleotide_content@side=inside",
    "--circular_track_slot", "gc_skew:dinucleotide_skew@side=inside",
)


class _Drawn:
    """The inner labels the drawer was given, the letters of the numbers the tick groups drew, and what the
    record's tracks paint."""

    def __init__(self) -> None:
        self.numbers: list[CircularInkRing] = []
        self.inner_labels: list[CircularLabel] = []
        self.total_length = 0
        self.track_ink: tuple[CircularInk, ...] = ()

    def over_numbers(self) -> list[str]:
        """The inner labels over a number's letters, by the placement's own test."""

        return [str(label["label_text"]) for label in labels_over_tick_numbers(self.inner_labels, self.numbers, self.total_length)]

    def track_paint_met(self, label: CircularLabel) -> set[int]:
        """The pieces of the tracks' paint the label's fit box meets."""

        (min_x, min_y), _, (max_x, max_y), _ = _label_fit_box(label, self.total_length)
        box = Aabb(min_x, min_y, max_x, max_y)
        return {index for index, ink in enumerate(self.track_ink) if box_meets_ink(box, (ink,), threshold_px=0.0)}

    def newly_over_tracks(self, before: _Drawn) -> list[str]:
        """The inner labels over a track's paint that the same label (by stable id) of ``before`` is not over."""

        earlier = {str(label["stable_id"]): label for label in before.inner_labels}
        newly: list[str] = []
        for label in self.inner_labels:
            met = self.track_paint_met(label)
            had = earlier.get(str(label["stable_id"]))
            if met - (before.track_paint_met(had) if had is not None else set()):
                newly.append(str(label["label_text"]))
        return newly


def _draw(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, gbk: Path, args: tuple[str, ...]) -> tuple[str, _Drawn]:
    drawn = _Drawn()
    tick_layout = builders.tick_slot_layout
    track_ink = assemble._record_track_ink
    draw_labels = assemble.add_labels_group_on_canvas

    def ticks(*args: Any, **kwargs: Any) -> Any:
        layout = tick_layout(*args, **kwargs)
        drawn.numbers.extend(number.letters for number in layout.numbers)
        return layout

    def paint(axis: Any, slots: Sequence[Any], sources: Any, annotation_labels: Any) -> tuple[CircularInk, ...]:
        inks = track_ink(axis, slots, sources, annotation_labels)
        drawn.track_ink = inks
        return inks

    def labels(canvas: Any, gb_record: Any, canvas_config: Any, **kwargs: Any) -> Any:
        if kwargs["phase"] == "text":
            drawn.total_length = kwargs.get("display_length") or len(gb_record.seq)
            drawn.inner_labels.extend(
                label for label in kwargs["precalculated_labels"]
                if label.get("is_inner") and not label.get("is_embedded")
            )
        return draw_labels(canvas, gb_record, canvas_config, **kwargs)

    with monkeypatch.context() as patch:
        # The drawing builds its tick groups through ``builders.tick_slot_layout``; the label stage reads its own.
        patch.setattr(builders, "tick_slot_layout", ticks)
        patch.setattr(assemble, "_record_track_ink", paint)
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


@pytest.mark.circular
def test_held_inner_labels_go_onto_no_track_paint_they_were_not_on(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """NC_005856.1 with GC content and GC skew inside the scale: the inner labels held off the numbers would sit on
    the GC tracks where the figure drawn without the numbers had them on white, so no inner label is drawn over a
    track's paint that it was not over in that figure."""

    _, held = _draw(tmp_path / "held", monkeypatch, P1, P1_ARGS + INSIDE_TRACKS)
    _without_tick_numbers(monkeypatch)
    _, unheld = _draw(tmp_path / "unheld", monkeypatch, P1, P1_ARGS + INSIDE_TRACKS)

    assert held.track_ink and len(held.inner_labels) > 20 and unheld.over_numbers()
    assert held.newly_over_tracks(unheld) == []
