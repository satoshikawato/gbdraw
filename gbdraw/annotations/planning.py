"""Shared semantic planning for circular and linear annotation tracks."""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Callable, Literal, Sequence, TypeVar

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from ..exceptions import ValidationError
from ..layout.record_coordinates import RecordDisplayTransform, SourceInterval
from .models import AnnotationOptions, ResolvedAnnotationBundle
from .resolve import resolve_annotations


logger = logging.getLogger(__name__)

SlotT = TypeVar("SlotT")


def prepare_annotation_track_slots(
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None,
    records: Sequence[SeqRecord],
    slots: list[SlotT] | None,
    *,
    mode: Literal["circular", "linear"],
    default_slots: Callable[[], list[SlotT]],
    slot_factory: Callable[..., SlotT],
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    record_indices: Sequence[int] | None = None,
) -> tuple[list[SlotT] | None, ResolvedAnnotationBundle, frozenset[str]]:
    """Resolve annotation inputs and bind their sets to track slots."""

    bundle = (
        annotations
        if isinstance(annotations, ResolvedAnnotationBundle)
        else resolve_annotations(annotations, records, mode=mode, record_transforms=record_transforms)
    )
    if not bundle.set_ids and not bundle.annotations:
        return slots, bundle, frozenset()

    if record_transforms is not None:
        transforms = dict(zip(record_indices or range(len(records)), record_transforms, strict=True))
        projected = []
        for item in bundle.annotations:
            transform = transforms.get(item.record_index)
            if transform is not None and transform.start_coordinate is not None:
                parts = transform.project_local_parts(tuple(SourceInterval(start, end) for start, end in item.segments))
                item = replace(item, display_parts=tuple(parts))
            projected.append(item)
        bundle = replace(bundle, annotations=tuple(projected))

    set_ids = bundle.set_ids or tuple(dict.fromkeys(item.set_id for item in bundle.annotations))
    if slots is None:
        generated: list[SlotT] = []
        for index, set_id in enumerate(set_ids, start=1):
            marks = tuple(
                dict.fromkeys(item.mark for item in bundle.annotations if item.set_id == set_id)
            )
            lane_marks = tuple(mark for mark in marks if mark != "highlight")
            if lane_marks:
                generated.append(
                    slot_factory(
                        id=f"annotations_{index}",
                        renderer="annotations",
                        side="outside" if mode == "circular" else "above",
                        params={"set_id": set_id, "marks": lane_marks},
                    )
                )
            if "highlight" in marks:
                generated.append(
                    slot_factory(
                        id=f"annotations_{index}{'_highlight' if lane_marks else ''}",
                        renderer="annotations",
                        side="overlay",
                        z=-1,
                        params={
                            "set_id": set_id,
                            "marks": ("highlight",),
                            "anchor_slot": "features",
                            "layer": "underlay",
                            "cover_anchor": True,
                            "padding_px": 0.0,
                        },
                    )
                )
        slots = [*generated, *default_slots()]

    requested = {
        str(getattr(slot, "params", {}).get("set_id", "")).strip()
        for slot in slots
        if str(getattr(slot, "renderer", "")).strip().lower() == "annotations"
    }
    unknown = requested - set(set_ids)
    if unknown:
        raise ValidationError(
            f"Annotation track references unknown set_id(s): {', '.join(sorted(unknown))}"
        )
    for set_id in set_ids:
        if set_id not in requested:
            logger.warning(
                "Annotation set %s is not referenced by a %s track slot.", set_id, mode
            )

    extent_name = "width" if mode == "circular" else "height"
    auto_slot_ids = frozenset(
        str(getattr(slot, "id"))
        for slot in slots
        if str(getattr(slot, "renderer", "")).strip().lower() == "annotations"
        and getattr(slot, extent_name) is None
    )
    return slots, bundle, auto_slot_ids
