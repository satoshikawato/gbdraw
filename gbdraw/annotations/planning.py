"""Shared semantic planning for circular and linear annotation tracks."""

from __future__ import annotations

from dataclasses import replace
from typing import Callable, Literal, Sequence, TypeVar, cast

from Bio.SeqRecord import SeqRecord

from ..exceptions import ValidationError
from ..layout.record_coordinates import DisplayFragment, RecordDisplayTransform, SourceInterval
from .models import AnnotationOptions, ResolvedAnnotationBundle
from .resolve import resolve_annotations


SlotT = TypeVar("SlotT")

# The side of the lane slot that a set without a stack gets, by mode. The Web's stack builders give a set the same row
# (generated into the Web contracts), so turning Use custom stack on does not move it.
AUTOMATIC_ANNOTATION_LANE_SIDES: dict[str, str] = {"circular": "outside", "linear": "above"}


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
    set that keeps its highlight overlay (a sequence search's, PD-OI-117), in front of them, covering them."""

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


def _bundle_set_ids(bundle: ResolvedAnnotationBundle) -> tuple[str, ...]:
    return bundle.set_ids or tuple(dict.fromkeys(item.set_id for item in bundle.annotations))


def _named_set_ids(slots: Sequence[object], *, enabled_only: bool) -> set[str]:
    return {
        str(getattr(slot, "params", {}).get("set_id", "")).strip()
        for slot in slots
        if _renderer(slot) == "annotations" and (not enabled_only or getattr(slot, "enabled", True))
    }


def _features_slot(slots: Sequence[SlotT]) -> SlotT | None:
    return next(
        (slot for slot in slots if _renderer(slot) == "features" and getattr(slot, "enabled", True)), None
    )


def _kept_overlay_set_ids(
    slots: Sequence[object], bundle: ResolvedAnnotationBundle, set_ids: Sequence[str]
) -> tuple[str, ...]:
    """The sets an explicit stack draws as their automatic highlight overlay: each set that keeps it
    (``AnnotationSet.keep_highlight_overlay``, a sequence search's set, PD-OI-117), has highlights, and that no
    enabled slot names, while the stack has an enabled feature track."""

    if _features_slot(slots) is None or not bundle.highlight_overlay_set_ids:
        return ()
    named = _named_set_ids(slots, enabled_only=True)
    highlighted = {item.set_id for item in bundle.annotations if item.mark == "highlight"}
    return tuple(
        set_id for set_id in set_ids
        if set_id in bundle.highlight_overlay_set_ids and set_id not in named and set_id in highlighted
    )


def _with_kept_highlight_overlays(
    slots: list[SlotT],
    bundle: ResolvedAnnotationBundle,
    set_ids: Sequence[str],
    slot_factory: Callable[..., SlotT],
) -> list[SlotT]:
    """An explicit stack with the automatic highlight overlay of each set that keeps it (``_kept_overlay_set_ids``):
    after the stack's own slots, so its axis index stays, in front of the stack's feature track. A slot that names
    the set wins; a stack without a feature track draws no overlay. Other sets follow the stack."""

    kept = set(_kept_overlay_set_ids(slots, bundle, set_ids))
    features = _features_slot(slots)
    if not kept or features is None:
        return slots
    used = {str(getattr(slot, "id")) for slot in slots}
    added: list[SlotT] = []
    for index, set_id in enumerate(set_ids, start=1):
        if set_id not in kept:
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
                    item,
                    display_parts=tuple(cast("Sequence[DisplayFragment]", parts)),
                    shares_circle=transform.circle is not None,
                )
            projected.append(item)
        bundle = replace(bundle, annotations=tuple(projected))

    set_ids = _bundle_set_ids(bundle)
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
                        side=AUTOMATIC_ANNOTATION_LANE_SIDES[mode],
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

    unknown = _named_set_ids(slots, enabled_only=False) - set(set_ids)
    if unknown:
        raise ValidationError(
            f"Annotation track references unknown set_id(s): {', '.join(sorted(unknown))}"
        )

    extent_name = "width" if mode == "circular" else "height"
    auto_slot_ids = frozenset(
        str(getattr(slot, "id"))
        for slot in slots
        if _renderer(slot) == "annotations"
        and getattr(slot, extent_name) is None
    )
    return slots, bundle, auto_slot_ids


def annotation_sets_left_out(
    bundle: ResolvedAnnotationBundle, slots: Sequence[object] | None
) -> tuple[str, ...]:
    """The sets with annotations that the explicit track stack ``slots`` draws nothing of: no annotations slot names
    them, and they get no kept highlight overlay (``_kept_overlay_set_ids``). Without a stack (``None``) every set
    gets its automatic slots, so none is left out. The input summary reports them (OV-3007)."""

    if slots is None:
        return ()
    set_ids = _bundle_set_ids(bundle)
    drawn = _named_set_ids(slots, enabled_only=False) | set(_kept_overlay_set_ids(slots, bundle, set_ids))
    with_annotations = {item.set_id for item in bundle.annotations}
    return tuple(set_id for set_id in set_ids if set_id not in drawn and set_id in with_annotations)
