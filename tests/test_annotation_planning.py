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
