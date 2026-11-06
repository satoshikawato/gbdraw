"""Shared semantic planning for circular and linear annotation tracks."""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Callable, Literal, Sequence, TypeVar, cast

from Bio.SeqRecord import SeqRecord

from ..exceptions import ValidationError
from ..layout.record_coordinates import DisplayFragment, RecordDisplayTransform, SourceInterval
from .models import AnnotationOptions, ResolvedAnnotationBundle
from .resolve import resolve_annotations


logger = logging.getLogger(__name__)

SlotT = TypeVar("SlotT")


def _renderer(slot: object) -> str:
    return str(getattr(slot, "renderer", "")).strip().lower()


def _highlight_overlay_slot(
    slot_factory: Callable[..., SlotT],
    slot_id: str,
    set_id: str,
    *,
    in_front: bool,
    anchor_slot: str = "features",
    anchor_z: int = 0,
) -> SlotT:
    """The automatic slot of a set's highlights, an overlay over the feature track: behind its features, or, for a
    set that keeps its highlight overlay (a sequence search's, PD-OI-116), in front of them, covering them."""

    return slot_factory(
        id=slot_id,
        renderer="annotations",
        side="overlay",
        z=anchor_z + 1 if in_front else anchor_z - 1,
        params={
            "set_id": set_id,
            "marks": ("highlight",),
            "anchor_slot": anchor_slot,
            "layer": "foreground" if in_front else "underlay",
            "cover_anchor": True,
            "padding_px": 0.0,
        },
    )


def _with_kept_highlight_overlays(
    slots: list[SlotT],
    bundle: ResolvedAnnotationBundle,
    set_ids: Sequence[str],
    slot_factory: Callable[..., SlotT],
) -> list[SlotT]:
    """An explicit stack with the automatic highlight overlay of each set that keeps it
    (``AnnotationSet.keep_highlight_overlay``, a sequence search's set, PD-OI-116) and that no enabled slot names:
    after the stack's own slots, so its axis index stays, in front of the stack's feature track. A slot that names
    the set wins; a stack without a feature track draws no overlay. Other sets follow the stack."""

    features = next(
        (slot for slot in slots if _renderer(slot) == "features" and getattr(slot, "enabled", True)), None
    )
    if features is None or not bundle.highlight_overlay_set_ids:
        return slots
    named = {
        str(getattr(slot, "params", {}).get("set_id", "")).strip()
        for slot in slots
        if _renderer(slot) == "annotations" and getattr(slot, "enabled", True)
    }
    highlighted = {item.set_id for item in bundle.annotations if item.mark == "highlight"}
    used = {str(getattr(slot, "id")) for slot in slots}
    added: list[SlotT] = []
    for index, set_id in enumerate(set_ids, start=1):
        if set_id not in bundle.highlight_overlay_set_ids or set_id in named or set_id not in highlighted:
            continue
        slot_id, suffix = f"annotations_{index}", 1
        while slot_id in used:
            suffix += 1
            slot_id = f"annotations_{index}_{suffix}"
        used.add(slot_id)
        added.append(_highlight_overlay_slot(
            slot_factory, slot_id, set_id, in_front=True, anchor_slot=str(getattr(features, "id")),
            anchor_z=int(getattr(features, "z", 0)),
        ))
    return [*slots, *added] if added else slots


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
            if transform is not None and transform.projects:
                parts = transform.project_local_parts(tuple(SourceInterval(start, end) for start, end in item.segments))
                # A projecting transform makes project_local_parts return DisplayFragments.
                item = replace(
                    item, display_parts=tuple(cast("Sequence[DisplayFragment]", parts))
                )
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
                    _highlight_overlay_slot(
                        slot_factory, f"annotations_{index}{'_highlight' if lane_marks else ''}", set_id,
                        in_front=set_id in bundle.highlight_overlay_set_ids,
                    )
                )
        slots = [*generated, *default_slots()]
    else:
        slots = _with_kept_highlight_overlays(slots, bundle, set_ids, slot_factory)

    requested = {
        str(getattr(slot, "params", {}).get("set_id", "")).strip()
        for slot in slots
        if _renderer(slot) == "annotations"
    }
    unknown = requested - set(set_ids)
    if unknown:
        raise ValidationError(
            f"Annotation track references unknown set_id(s): {', '.join(sorted(unknown))}"
        )
    for set_id in set_ids:
        if set_id not in requested and any(item.set_id == set_id for item in bundle.annotations):
            logger.warning(
                "Annotation set %s is not referenced by a %s track slot.", set_id, mode
            )

    extent_name = "width" if mode == "circular" else "height"
    auto_slot_ids = frozenset(
        str(getattr(slot, "id"))
        for slot in slots
        if _renderer(slot) == "annotations"
        and getattr(slot, extent_name) is None
    )
    return slots, bundle, auto_slot_ids
