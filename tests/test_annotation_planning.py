from dataclasses import dataclass, field
from typing import Any, Mapping

import pytest

from gbdraw.annotations.models import ResolvedAnnotationBundle
from gbdraw.annotations.planning import annotation_sets_left_out, prepare_annotation_track_slots
from gbdraw.exceptions import ValidationError


@dataclass(frozen=True)
class _Slot:
    id: str
    renderer: str
    height: object | None = None
    params: Mapping[str, Any] = field(default_factory=dict)


def test_unknown_annotation_set_uses_shared_validation_error() -> None:
    slots = [
        _Slot(
            id="annotations",
            renderer="annotations",
            params={"set_id": "missing"},
        )
    ]
    bundle = ResolvedAnnotationBundle(annotations=(), set_ids=("known",))

    with pytest.raises(ValidationError, match="unknown set_id.*missing"):
        prepare_annotation_track_slots(
            bundle,
            [],
            slots,
            mode="linear",
            default_slots=list,
            slot_factory=_Slot,
        )


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_empty_resolved_set_has_no_native_automatic_slot_but_preserves_explicit(mode):
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    factory = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    bundle = ResolvedAnnotationBundle(annotations=(), set_ids=("empty",))
    slots, returned, _ = prepare_annotation_track_slots(
        bundle, [], None, mode=mode, default_slots=list, slot_factory=factory
    )
    assert slots == [] and returned is bundle
    explicit = [
        factory(
            id="empty",
            renderer="annotations",
            side="outside" if mode == "circular" else "above",
            params={"set_id": "empty"},
        )
    ]
    slots, returned, _ = prepare_annotation_track_slots(
        bundle, [], explicit, mode=mode, default_slots=list, slot_factory=factory
    )
    assert slots is explicit and returned is bundle


def _search_bundle(*, keep: bool = True) -> ResolvedAnnotationBundle:
    from gbdraw.annotations.models import AnnotationSet, CoordinateSpan, RegionAnnotation
    from gbdraw.annotations.resolve import resolve_annotations
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ATGC" * 250), id="rec")
    record.annotations["molecule_type"] = "DNA"
    other = AnnotationSet(id="other", annotations=(
        RegionAnnotation(id="o1", target=CoordinateSpan(record=None, start=10, end=50), mark="highlight"),
    ))
    search = AnnotationSet(id="search_1", annotations=(
        RegionAnnotation(id="q#1", target=CoordinateSpan(record=None, start=100, end=200), mark="highlight"),
    ), keep_highlight_overlay=keep)
    return resolve_annotations((other, search), [record], mode="circular")


def _stack(mode: str, *extra: object) -> list:
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    if mode == "circular":
        return [CircularTrackSlot(id="ticks", renderer="ticks"), CircularTrackSlot(id="feat", renderer="features", z=2),
                *extra]
    return [LinearTrackSlot(id="feat", renderer="features", z=2), *extra]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_explicit_stack_draws_a_kept_set_in_its_highlight_overlay(mode: str) -> None:
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    factory = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    explicit = _stack(mode)
    slots, _, auto_ids = prepare_annotation_track_slots(
        _search_bundle(), [], explicit, mode=mode, default_slots=list, slot_factory=factory
    )
    assert slots is not None and slots[: len(explicit)] == explicit
    # Only the kept set gets a slot, after the stack's own (its axis index stays), in front of its feature track.
    assert slots[len(explicit):] == [factory(
        id="annotations_2", renderer="annotations", side="overlay", z=3,
        params={"set_id": "search_1", "marks": ("highlight",), "anchor_slot": "feat", "layer": "foreground",
                "cover_anchor": True, "padding_px": 0.0},
    )]
    assert "annotations_2" in auto_ids


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_automatic_stack_draws_a_kept_set_in_front_and_others_behind(mode: str) -> None:
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    factory = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    slots, _, _ = prepare_annotation_track_slots(
        _search_bundle(), [], None, mode=mode, default_slots=list, slot_factory=factory
    )
    assert slots is not None
    layers = {slot.params["set_id"]: (slot.params["layer"], slot.z, slot.params["anchor_slot"]) for slot in slots}
    assert layers == {"other": ("underlay", -1, "features"), "search_1": ("foreground", 1, "features")}


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_slot_that_names_the_kept_set_wins(mode: str) -> None:
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    factory = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    named = factory(
        id="annotations_2", renderer="annotations", side="outside" if mode == "circular" else "above",
        params={"set_id": "search_1"},
    )
    explicit = _stack(mode, named)
    slots, _, _ = prepare_annotation_track_slots(
        _search_bundle(), [], explicit, mode=mode, default_slots=list, slot_factory=factory
    )
    assert slots is explicit


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_set_without_the_marker_follows_the_stack(mode: str) -> None:
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    factory = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    explicit = _stack(mode)
    slots, _, _ = prepare_annotation_track_slots(
        _search_bundle(keep=False), [], explicit, mode=mode, default_slots=list, slot_factory=factory
    )
    assert slots is explicit


def test_the_kept_overlay_takes_an_unused_slot_id() -> None:
    from gbdraw.tracks import CircularTrackSlot

    taken = CircularTrackSlot(id="annotations_2", renderer="annotations", side="outside", params={"set_id": "other"})
    slots, _, _ = prepare_annotation_track_slots(
        _search_bundle(), [], _stack("circular", taken), mode="circular", default_slots=list,
        slot_factory=CircularTrackSlot,
    )
    assert slots is not None and [slot.id for slot in slots][-1] == "annotations_2_2"


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_sets_an_explicit_stack_leaves_out_are_the_ones_it_draws_nothing_of(mode: str) -> None:
    """OV-3007: the input summary names a set no annotations slot names; a kept set drawn as its overlay, a set a
    slot names (on or off), a set without annotations, and every set without a stack are not left out."""
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    factory = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    side = "outside" if mode == "circular" else "above"
    bundle = _search_bundle()
    assert annotation_sets_left_out(bundle, None) == ()
    # The kept search set is drawn as its overlay in front of the feature track; "other" is left out.
    assert annotation_sets_left_out(bundle, _stack(mode)) == ("other",)
    named_off = factory(id="o", renderer="annotations", side=side, enabled=False, params={"set_id": "other"})
    assert annotation_sets_left_out(bundle, _stack(mode, named_off)) == ()
    # Without a feature track the kept set has no overlay to be drawn in.
    assert annotation_sets_left_out(bundle, [factory(id="gc", renderer="gc_content")]) == ("other", "search_1")
    empty = ResolvedAnnotationBundle(annotations=(), set_ids=("empty",))
    assert annotation_sets_left_out(empty, _stack(mode)) == ()
