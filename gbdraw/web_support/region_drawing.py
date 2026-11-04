"""Worker helper behind the Web's New drawing from a region (PD-OI-090 〔作成〕(4)).

The page sends the source drawing's current draft projection (its canonical
request, as the next Generate would send it, with its staged resources), the
regions and whether sizes adapt. The derivation is the one the CLI and the
Python API use (:func:`gbdraw.api.region_drawing.derive_region_request`); this
module only decodes the request, encodes the derived one with the source's
resource IDs, and types the reply. It renders nothing.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Literal, TypedDict

from gbdraw.api.region_drawing import (
    FeatureRegionSelection,
    NotCarriedKind,
    RegionChoice,
    RegionSelection,
    StartAfterEndReason,
    StartAfterEndRegion,
    derive_region_request,
)
from gbdraw.auto_sizes import AutoValue
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


RegionDrawingReply = DerivedRegionDrawing | RefusedRegionDrawing


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


def derive_web_region_drawing(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    regions: Sequence[RegionInput],
    *,
    adapt_sizes: bool,
) -> RegionDrawingReply:
    """The region drawing of the source drawing's canonical request ``request_payload``.

    ``resource_paths`` are its staged resources by resource ID. The derived
    request names the same resources: the source's files are its only inputs.
    """

    if not isinstance(adapt_sizes, bool):
        raise _protocol_error("adaptSizes must be true or false.")
    source = decode_canonical_request(request_payload, resource_paths=resource_paths, output_directory=output_directory)
    try:
        derived = derive_region_request(source, _regions(regions), adapt_sizes=adapt_sizes)
    except StartAfterEndRegion as refusal:
        return {
            "status": "refused", "reason": refusal.reason,
            "recordKey": refusal.record_key, "start": refusal.start, "end": refusal.end,
        }
    # The source's staged files, held under their own IDs: the encoder takes
    # a held resource with equal bytes, so the derived request names them.
    table = SessionResourceTable()
    for resource_id, path in resource_paths.items():
        table.request().add_path(resource_id, kind="web-file", value=path)
    encoded = encode_canonical_request(derived.request, table=table)
    # Only the resources the derived request names: a held one reads itself,
    # an added one (a second input of one file) the held one of its path.
    source_of_path = {str(path): resource_id for resource_id, path in resource_paths.items()}
    added = {resource.resource_id: str(resource.source_path) for resource in encoded.resources}
    resource_sources: dict[str, str] = {}
    for resource_id in canonical_resource_ids(encoded.payload):
        source_id = resource_id if resource_id in resource_paths else source_of_path.get(added.get(resource_id, ""))
        if source_id is None:
            raise _protocol_error("The region drawing names a resource its source drawing does not have.")
        resource_sources[resource_id] = source_id
    return {
        "status": "derived",
        "mode": derived.mode,
        "renderRequest": encoded.payload,
        "resourceSources": resource_sources,
        "returnedToAuto": [
            {"setting": reset.setting, "explicit": reset.explicit, "auto": reset.target_auto}
            for reset in derived.adaptation.reset
        ],
        "notCarried": [
            {"kind": item.kind, "count": item.count, "items": list(item.items)} for item in derived.not_carried
        ],
    }


def derive_web_region_drawing_json(
    request_json: str,
    resource_paths_json: str,
    output_directory: str,
    regions_json: str,
    adapt_sizes_json: str,
) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(
        derive_web_region_drawing(
            json.loads(str(request_json)),
            json.loads(str(resource_paths_json)),
            str(output_directory),
            json.loads(str(regions_json)),
            adapt_sizes=json.loads(str(adapt_sizes_json)),
        )
    )


__all__ = [
    "DerivedRegionDrawing",
    "RefusedRegionDrawing",
    "RegionDrawingReply",
    "derive_web_region_drawing",
    "derive_web_region_drawing_json",
]
