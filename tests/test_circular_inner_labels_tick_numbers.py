"""Horizontal inner labels and the scale's numbers.

With a ticks slot between the feature track and the inner labels (``--separate_strands --track_type spreadout``),
inner labels were drawn over the numbers ("SimC" and "PmgC-like head morphogenesis" over "50 kbp" of NC_005856.1).
The figures are drawn by the command line; the tests read the labels the label drawer is given, the letters of the
numbers the drawn tick groups have (``TickSlotLayout.numbers``) and what the data tracks paint (``_record_track_ink``
of ``DATA_TRACK_RENDERERS``). Whether a label is over a number is the placement's own test
(``labels_over_tick_numbers``). The inner labels are held inside the numbers whatever paint they land on, and under
the default protection (``auto``) every inner label whose letters lie on a data track's paint gets the white box of
``--label_protection box``; the tests find those labels on their own (each label's letters' box against the paint,
``box_meets_ink``) and the boxes in the SVG by their label. The paint is read after the ``circular_labels`` render
stage, so a Worker transaction still reuses the labels.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import pytest

from gbdraw import circular as circular_cli_module
from gbdraw.diagrams.circular import assemble, builders
from gbdraw.labels.circular import labels_over_tick_numbers
from gbdraw.labels.circular_types import CircularLabel
from gbdraw.labels.label_ink import circular_label_ink
from gbdraw.labels.label_style import FeatureLabelStyle
from gbdraw.layout.circular_ink import CircularInk, CircularInkRing, box_meets_ink
from gbdraw.layout.spatial import bounding_aabb
from tests.utils.render_stage_equivalence import EquivalenceHarness

INPUTS = Path(__file__).resolve().parent / "test_inputs"
P1 = INPUTS / "NC_005856.1.gb"
BGC = INPUTS / "BGC0000708.gbk"
HMM = INPUTS / "HmmtDNA.gbk"
P1_ARGS = ("--separate_strands", "--track_type", "spreadout", "--labels", "both")
BGC_ARGS = ("--definition_font_size", "30", "--species", "<i>Genus species</i>", "--strain", "strain X-1", "--labels", "both")
# GC content and GC skew inside the scale, under the inner labels.
INSIDE_TRACKS = (
    "--circular_track_slot", "features:features",
    "--circular_track_slot", "ticks:ticks@side=inside",
    "--circular_track_slot", "gc_content:dinucleotide_content@side=inside",
    "--circular_track_slot", "gc_skew:dinucleotide_skew@side=inside",
)
# GC content and GC skew inside the features, the scale outside: inner labels on the GC paint, none near a number.
INSIDE_GC_ONLY = (
    "--circular_track_slot", "features:features",
    "--circular_track_slot", "gc_content:dinucleotide_content@side=inside",
    "--circular_track_slot", "gc_skew:dinucleotide_skew@side=inside",
)
# A box behind a label's letters, bound to its label (``feature_label_box``).
BOX = re.compile(r'<rect aria-hidden="true" data-label-feature-id="([^"]*)" fill="#ffffff" fill-opacity="0.8"[^>]*/>')


class _Drawn:
    """The inner labels the drawer was given, the letters of the numbers the tick groups drew, and what the
    record's data tracks paint."""

    def __init__(self) -> None:
        self.numbers: list[CircularInkRing] = []
        self.inner_labels: list[CircularLabel] = []
        self.total_length = 0
        self.data_paint: tuple[CircularInk, ...] = ()

    def over_numbers(self) -> list[str]:
        """The inner labels over a number's letters, by the placement's own test."""

        return [str(label["label_text"]) for label in labels_over_tick_numbers(self.inner_labels, self.numbers, self.total_length)]

    def on_data_paint(self) -> set[str]:
        """The feature ids of the inner labels whose letters' box (turned with the letters) meets the data paint."""

        style = FeatureLabelStyle("none", "none", False)
        on_paint: set[str] = set()
        for label in self.inner_labels:
            ink = circular_label_ink(label, self.total_length, style)
            box = bounding_aabb(tuple(ink.from_label_frame(corner) for corner in (
                (ink.ink.min_x, ink.ink.min_y), (ink.ink.max_x, ink.ink.min_y),
                (ink.ink.max_x, ink.ink.max_y), (ink.ink.min_x, ink.ink.max_y),
            )))
            if box is not None and box_meets_ink(box, self.data_paint, threshold_px=0.0):
                on_paint.add(str(label["feature_id"]))
        return on_paint


def _draw(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, gbk: Path, args: tuple[str, ...]) -> tuple[str, _Drawn]:
    drawn = _Drawn()
    tick_layout = builders.tick_slot_layout
    track_ink = assemble._record_track_ink
    draw_labels = assemble.add_labels_group_on_canvas

    def ticks(*args: Any, **kwargs: Any) -> Any:
        layout = tick_layout(*args, **kwargs)
        drawn.numbers.extend(number.letters for number in layout.numbers)
        return layout

    def paint(axis: Any, slots: Sequence[Any], sources: Any, annotation_labels: Any, *renderers: Any) -> tuple[CircularInk, ...]:
        inks = track_ink(axis, slots, sources, annotation_labels, *renderers)
        if renderers and renderers[0] is not None:
            drawn.data_paint += inks
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


def _boxed(svg: str) -> set[str]:
    """The feature ids of the labels drawn with a white box."""

    return set(BOX.findall(svg))


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
def test_inner_labels_are_held_inside_the_numbers_over_data_tracks_and_boxed_there(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """NC_005856.1 with GC content and GC skew inside the scale: the inner labels are held inside the numbers although
    the held places lie on the GC paint, and each inner label whose letters lie on that paint, and no other label,
    gets a white box."""

    svg, drawn = _draw(tmp_path / "held", monkeypatch, P1, P1_ARGS + INSIDE_TRACKS)
    _without_tick_numbers(monkeypatch)
    _, unheld = _draw(tmp_path / "unheld", monkeypatch, P1, P1_ARGS + INSIDE_TRACKS)

    assert unheld.over_numbers() and len(drawn.inner_labels) > 20
    assert drawn.over_numbers() == []
    assert drawn.data_paint and len(drawn.on_data_paint()) > 10
    assert _boxed(svg) == drawn.on_data_paint()


@pytest.mark.circular
def test_inner_labels_on_data_tracks_gain_a_box_and_nothing_else(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """HmmtDNA with GC content and GC skew inside the features: no inner label is near a number, some lie on the GC
    paint; the default draws a box behind each of those, and the figure without the boxes is the one an explicit
    ``--label_protection none`` draws."""

    args = ("--labels", "both", *INSIDE_GC_ONLY)
    svg, drawn = _draw(tmp_path / "auto", monkeypatch, HMM, args)
    plain, _ = _draw(tmp_path / "none", monkeypatch, HMM, (*args, "--label_protection", "none"))

    assert drawn.on_data_paint() and _boxed(svg) == drawn.on_data_paint()
    assert BOX.sub("", svg) == plain


@pytest.mark.circular
@pytest.mark.parametrize("protection", ["none", "halo", "box"])
def test_an_explicit_label_protection_reads_no_data_paint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, protection: str
) -> None:
    """An explicit ``--label_protection`` draws as it drew before ``auto``: the data paint is not built, and the
    boxes are those of ``box`` on every label with a leader, or none."""

    args = ("--labels", "both", *INSIDE_GC_ONLY, "--label_protection", protection)
    svg, drawn = _draw(tmp_path, monkeypatch, HMM, args)

    assert drawn.inner_labels and drawn.data_paint == ()
    assert bool(_boxed(svg)) is (protection == "box")


@pytest.mark.circular
@pytest.mark.parametrize("tracks", [(), INSIDE_TRACKS], ids=["held", "held-boxed-on-gc-paint"])
def test_a_warm_render_reuses_the_labels_held_off_the_numbers(tmp_path: Path, tracks: tuple[str, ...]) -> None:
    """In a Worker transaction the ``circular_labels`` stage of NC_005856.1 (inner labels held off the numbers;
    with GC tracks inside, also boxed on their paint after the stage) is built once: the same render again and a GC
    window change, which the labels do not read, reuse it."""

    harness = EquivalenceHarness("circular", tmp_path)
    args = ["--gbk", str(P1), *P1_ARGS, *tracks]
    first, again, window = (
        harness.warm(arguments).stages for arguments in (args, args, [*args, "--window", "500", "--step", "100"])
    )

    counts = [
        (stages.builds.get("circular_labels", 0), stages.hits.get("circular_labels", 0),
         stages.uncacheable.get("circular_labels", 0))
        for stages in (first, again, window)
    ]
    assert counts == [(1, 0, 0), (0, 1, 0), (0, 1, 0)]


V0145_WEB_CIRCULAR_SESSION = (
    Path(__file__).resolve().parent / "fixtures" / "sessions" / "tobacco-chloroplast.v50-web0145.gbdraw-session.json.gz"
)


@pytest.mark.circular
def test_a_session_saved_before_auto_replays_without_boxes(tmp_path: Path) -> None:
    """The v0.14.5 Web saved this Circular drawing (inner labels on, the default protection) with no
    ``labels.circular.protection`` in its request: it drew no box. Replay reads the unstated protection as that
    writer's ``none`` and draws no box; the same request stating ``auto`` boxes inner labels on the GC paint."""

    from dataclasses import replace

    from gbdraw.api import render_request
    from gbdraw.session import materialize_session, session_to_request
    from gbdraw.session_legacy_values import CIRCULAR_LABEL_PROTECTION_PATH

    circular_cli_module.circular_main(
        ["--session", str(V0145_WEB_CIRCULAR_SESSION), "-o", str(tmp_path / "replay"), "-f", "svg"]
    )
    with materialize_session(V0145_WEB_CIRCULAR_SESSION, output_directory=tmp_path / "session") as materialized:
        request = session_to_request(materialized)
        overrides = dict(request.options.config_overrides or {})
        auto = render_request(replace(request, options=replace(
            request.options, config_overrides={**overrides, CIRCULAR_LABEL_PROTECTION_PATH: "auto"}
        )))

    assert overrides[CIRCULAR_LABEL_PROTECTION_PATH] == "none"
    assert _boxed((tmp_path / "replay.svg").read_text(encoding="utf-8")) == set()
    assert any(_boxed(item.drawing.tostring()) for item in auto.items)
