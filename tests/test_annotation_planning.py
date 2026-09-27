from dataclasses import dataclass, field
from typing import Any, Mapping

import pytest

from gbdraw.annotations.models import ResolvedAnnotationBundle
from gbdraw.annotations.planning import prepare_annotation_track_slots
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
