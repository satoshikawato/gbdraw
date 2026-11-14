"""Worker helper behind the Web's New drawing from a region (PD-OI-090 〔作成〕(4)) and Copy to another drawing.

The page sends the source drawing's current draft projection (its canonical
request, as the next Generate would send it, with its staged resources), the
regions and whether sizes adapt, and, for Copy to another drawing, the copy's
choices: a new drawing's settings and mode, or the existing drawing's draft
projection with what comes along and what a record already there does. The
derivations are the ones the CLI and the Python API use
(:func:`gbdraw.api.region_drawing.derive_region_request`,
:func:`gbdraw.api.region_drawing.paste_region_request`); this module only
decodes the requests, encodes the result with the resource IDs of the drawings
whose files it reads, and types the reply. It renders nothing.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Literal, TypedDict

from gbdraw.api.region_drawing import (
    CopyRefusalReason,
    CopyRefused,
    DuplicatesChoice,
    FeatureRegionSelection,
    NotCarried,
    NotCarriedKind,
    PasteEdits,
    PasteNotCarriedKind,
    RegionChoice,
    RegionSelection,
    SettingsChoice,
    StartAfterEndReason,
    StartAfterEndRegion,
    derive_region_request,
    paste_region_request,
)
from gbdraw.api.requests import DiagramRequest
from gbdraw.auto_sizes import AutoValue, SettingAdaptation
from gbdraw.exceptions import ValidationError
from gbdraw.session_drawings import DrawingMode
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request
from gbdraw.session_resources import SessionResourceTable, canonical_resource_ids


class CoordinateRegionInput(TypedDict):
    """One region: the source drawing's record key and source coordinates, 1-based and inclusive."""

    recordKey: str
    start: int
    end: int


class FeatureRegionInput(TypedDict):
    """One region: the span of a feature, by its record key and biological feature ID."""

    recordKey: str
    featureId: str


RegionInput = CoordinateRegionInput | FeatureRegionInput


class ReturnedSize(TypedDict):
    """A hand-set size that returned to Auto: the setting, the hand-set value, and the Auto value it takes."""

    setting: str
    explicit: AutoValue
    auto: AutoValue | None


class NotCarriedItem(TypedDict):
    """A part of the source drawing the new drawing does not carry, by kind, with its locators."""

    kind: NotCarriedKind
    count: int
    items: list[str]


class DerivedRegionDrawing(TypedDict):
    """The derived drawing: its mode, its canonical request (as a Session stores a committed one), and the lists.

    ``resourceSources`` names, for each resource the request names, the
    source drawing's resource whose file it reads. Each input keeps a
    resource of its own, so two records of one file name two resources of
    the same bytes.
    """

    status: Literal["derived"]
    mode: DrawingMode
    renderRequest: dict[str, object]
    resourceSources: dict[str, str]
    returnedToAuto: list[ReturnedSize]
    notCarried: list[NotCarriedItem]


class RefusedRegionDrawing(TypedDict):
    """A refused region: the reason (worded for the record's topology) and the region it names."""

    status: Literal["refused"]
    reason: StartAfterEndReason
    recordKey: str
    start: int
    end: int


class NewDrawingCopy(TypedDict):
    """Copy to a new drawing: its settings (the source drawing's, or the mode's defaults) and its mode (``None``:
    the source's mode when every region is a whole record, else Linear)."""

    into: Literal["new"]
    settings: SettingsChoice
    mode: DrawingMode | None


class PasteTarget(TypedDict):
    """The existing drawing's draft projection: its canonical request and its staged resources by resource ID."""

    renderRequest: dict[str, object]
    resourcePaths: dict[str, str]


class ExistingDrawingCopy(TypedDict):
    """Copy into an existing drawing: what comes along with the records, and what a selection the drawing draws
    already does."""

    into: Literal["existing"]
    target: PasteTarget
    edits: PasteEdits
    duplicates: DuplicatesChoice


CopyChoice = NewDrawingCopy | ExistingDrawingCopy


class AlreadyPresentItem(TypedDict):
    """A selected span the existing drawing draws already: the source record key, the span in source coordinates
    (the whole record when ``whole``), and the existing drawing's record that draws it."""

    recordKey: str
    start: int
    end: int
    whole: bool
    targetRecordKey: str


class PasteNotCarriedItem(TypedDict):
    """A part a paste does not bring, or the pasted records lack, by kind, with its count."""

    kind: PasteNotCarriedKind
    count: int


class PastedRecords(TypedDict):
    """The existing drawing with the pasted records: its mode, its canonical request, and the lists.

    Each resource the request names reads either a resource of the existing
    drawing (``targetResourceSources``) or of the source drawing
    (``resourceSources``), by that drawing's resource ID. ``pastedRecordKeys``
    are the pasted records' keys in the request; empty when every selection
    was already there and left out.
    """

    status: Literal["pasted"]
    mode: DrawingMode
    renderRequest: dict[str, object]
    resourceSources: dict[str, str]
    targetResourceSources: dict[str, str]
    pastedRecordKeys: list[str]
    alreadyPresent: list[AlreadyPresentItem]
    notCarried: list[PasteNotCarriedItem]


class RefusedCopy(TypedDict):
    """A copy the drawing cannot take, by reason."""

    status: Literal["refused"]
    reason: CopyRefusalReason


RegionDrawingReply = DerivedRegionDrawing | RefusedRegionDrawing | PastedRecords | RefusedCopy


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL"})


def _regions(value: object) -> list[RegionChoice]:
    if not isinstance(value, list):
        raise _protocol_error("A region drawing needs a list of regions.")
    regions: list[RegionChoice] = []
    for item in value:
        if not (isinstance(item, Mapping) and isinstance(item.get("recordKey"), str)):
            raise _protocol_error("A region needs a record key.")
        if isinstance(item.get("featureId"), str) and "start" not in item and "end" not in item:
            regions.append(FeatureRegionSelection(item["recordKey"], item["featureId"]))
        elif "featureId" not in item and all(type(item.get(key)) is int for key in ("start", "end")):
            regions.append(RegionSelection(item["recordKey"], item["start"], item["end"]))
        else:
            raise _protocol_error("A region needs a start and an end, or a feature ID.")
    return regions


def _copy_choice(value: object) -> CopyChoice | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise _protocol_error("A copy needs its choices.")
    if value.get("into") == "new" and set(value) == {"into", "settings", "mode"} and value["settings"] in {
        "keep", "defaults"
    } and value["mode"] in {None, "circular", "linear"}:
        return {"into": "new", "settings": value["settings"], "mode": value["mode"]}
    target = value.get("target")
    if (
        value.get("into") == "existing"
        and set(value) == {"into", "target", "edits", "duplicates"}
        and value["edits"] in {"edits", "records"}
        and value["duplicates"] in {"skip", "add"}
        and isinstance(target, Mapping)
        and set(target) == {"renderRequest", "resourcePaths"}
        and isinstance(target["renderRequest"], Mapping)
        and isinstance(target["resourcePaths"], Mapping)
        and all(isinstance(item, str) for pair in target["resourcePaths"].items() for item in pair)
    ):
        return {
            "into": "existing",
            "target": {"renderRequest": dict(target["renderRequest"]), "resourcePaths": dict(target["resourcePaths"])},
            "edits": value["edits"],
            "duplicates": value["duplicates"],
        }
    raise _protocol_error("A copy goes into a new drawing with its settings and mode, or into an existing drawing.")


def derive_web_region_drawing(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    regions: Sequence[RegionInput],
    *,
    adapt_sizes: bool,
    copy: object = None,
) -> RegionDrawingReply:
    """The region drawing of the source drawing's canonical request ``request_payload``, or the paste of its
    regions into an existing drawing.

    ``resource_paths`` are its staged resources by resource ID. A derived
    request names the same resources: the source's files are its only inputs.
    ``copy`` is a :data:`CopyChoice`; without it (New drawing from a region) a
    Circular drawing takes one record, as before Copy to another drawing.
    """

    if not isinstance(adapt_sizes, bool):
        raise _protocol_error("adaptSizes must be true or false.")
    choice = _copy_choice(copy)
    source = decode_canonical_request(request_payload, resource_paths=resource_paths, output_directory=output_directory)
    selections = _regions(regions)
    try:
        if choice is not None and choice["into"] == "existing":
            return _pasted_reply(source, selections, choice, resource_paths, output_directory)
        derived = derive_region_request(
            source,
            selections,
            mode=None if choice is None else choice["mode"],
            adapt_sizes=adapt_sizes,
            settings="keep" if choice is None else choice["settings"],
            circular_canvas=choice is not None,
        )
    except StartAfterEndRegion as refusal:
        return {
            "status": "refused", "reason": refusal.reason,
            "recordKey": refusal.record_key, "start": refusal.start, "end": refusal.end,
        }
    except CopyRefused as refusal:
        return {"status": "refused", "reason": refusal.reason}
    return derived_drawing_reply(
        derived.request, derived.mode, derived.adaptation, derived.not_carried, resource_paths
    )


def _pasted_reply(
    source: DiagramRequest,
    selections: Sequence[RegionChoice],
    choice: ExistingDrawingCopy,
    resource_paths: Mapping[str, str],
    output_directory: str,
) -> PastedRecords:
    """The paste of ``selections`` of ``source`` into the existing drawing of ``choice``, encoded with the resources
    of the drawing whose file each record reads (the existing drawing's first)."""

    target_paths = choice["target"]["resourcePaths"]
    target = decode_canonical_request(
        choice["target"]["renderRequest"], resource_paths=target_paths, output_directory=output_directory
    )
    pasted = paste_region_request(source, selections, target, edits=choice["edits"], duplicates=choice["duplicates"])
    payload, named = _encoded(pasted.request, target_paths)
    target_sources, source_sources = _resource_sources(named, target_paths, resource_paths)
    return {
        "status": "pasted",
        "mode": pasted.mode,
        "renderRequest": payload,
        "resourceSources": source_sources,
        "targetResourceSources": target_sources,
        "pastedRecordKeys": list(pasted.pasted_record_keys),
        "alreadyPresent": [
            {
                "recordKey": item.record_key, "start": item.start, "end": item.end, "whole": item.whole,
                "targetRecordKey": item.target_record_key,
            }
            for item in pasted.already_present
        ],
        "notCarried": [{"kind": item.kind, "count": item.count} for item in pasted.not_carried],
    }


def _encoded(request: DiagramRequest, held: Mapping[str, str]) -> tuple[dict[str, object], dict[str, str]]:
    """``request`` encoded with the staged files ``held`` under their own IDs, and the file of each resource it
    names: the encoder takes a held resource with equal bytes, and adds one for a second input of a file."""

    table = SessionResourceTable()
    for resource_id, path in held.items():
        table.request().add_path(resource_id, kind="web-file", value=path)
    encoded = encode_canonical_request(request, table=table)
    added = {resource.resource_id: str(resource.source_path) for resource in encoded.resources}
    return encoded.payload, {
        resource_id: held[resource_id] if resource_id in held else added.get(resource_id, "")
        for resource_id in canonical_resource_ids(encoded.payload)
    }


def _resource_sources(named: Mapping[str, str], *drawings: Mapping[str, str]) -> tuple[dict[str, str], ...]:
    """For each drawing (its staged resources by ID), the resources of ``named`` that read its files, by its
    resource ID: a held resource reads itself, an added one the resource of its path. The first drawing whose
    files a resource reads has it."""

    sources: tuple[dict[str, str], ...] = tuple({} for _ in drawings)
    for resource_id, path in named.items():
        for paths, found in zip(drawings, sources, strict=True):
            source_id = resource_id if paths.get(resource_id) == path else next(
                (held_id for held_id, held_path in paths.items() if str(held_path) == path), None
            )
            if source_id is not None:
                found[resource_id] = source_id
                break
        else:
            raise _protocol_error("The derived drawing names a resource its source drawing does not have.")
    return sources


def derived_drawing_reply(
    request: DiagramRequest,
    mode: DrawingMode,
    adaptation: SettingAdaptation,
    not_carried: Sequence[NotCarried],
    resource_paths: Mapping[str, str],
) -> DerivedRegionDrawing:
    """A drawing derived from a source drawing (New drawing from a region, Find in other files), encoded with the
    source's resources: the request names the source's resource of each file it reads."""

    payload, named = _encoded(request, resource_paths)
    (resource_sources,) = _resource_sources(named, resource_paths)
    return {
        "status": "derived",
        "mode": mode,
        "renderRequest": payload,
        "resourceSources": resource_sources,
        "returnedToAuto": [
            {"setting": reset.setting, "explicit": reset.explicit, "auto": reset.target_auto}
            for reset in adaptation.reset
        ],
        "notCarried": [
            {"kind": item.kind, "count": item.count, "items": list(item.items)} for item in not_carried
        ],
    }


def derive_web_region_drawing_json(
    request_json: str,
    resource_paths_json: str,
    output_directory: str,
    regions_json: str,
    adapt_sizes_json: str,
    copy_json: str = "null",
) -> str:
    """JSON boundary of the diagram Worker; ``copy_json`` is the copy's choices, ``null`` for New drawing from a
    region."""

    return json.dumps(
        derive_web_region_drawing(
            json.loads(str(request_json)),
            json.loads(str(resource_paths_json)),
            str(output_directory),
            json.loads(str(regions_json)),
            adapt_sizes=json.loads(str(adapt_sizes_json)),
            copy=json.loads(str(copy_json)),
        )
    )


__all__ = [
    "AlreadyPresentItem",
    "CopyChoice",
    "DerivedRegionDrawing",
    "ExistingDrawingCopy",
    "NewDrawingCopy",
    "PasteNotCarriedItem",
    "PasteTarget",
    "PastedRecords",
    "RefusedCopy",
    "RefusedRegionDrawing",
    "RegionDrawingReply",
    "derive_web_region_drawing",
    "derive_web_region_drawing_json",
    "derived_drawing_reply",
]
