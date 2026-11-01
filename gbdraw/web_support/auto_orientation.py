"""Worker helper behind the Web's "Orient records from comparisons" button.

The page sends the request of the last Generate, which holds the comparisons
in the search frame as the Session keeps them. The helper plans it with
``auto_orient`` as ``--auto_orient`` does on the command line, so the rule
(:func:`gbdraw.linear_comparison.auto_orient_reversed_records`) decides in its
one owner, and returns the orientation each record is then drawn with, by
record key. The page writes the records whose orientation changes into their
Reverse complement checkboxes as one History step; the next Generate draws
them.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import replace
from typing import TypedDict

from gbdraw.api.record_planning import record_is_reverse_complemented
from gbdraw.api.request_render import plan_linear_request
from gbdraw.api.requests import LinearDiagramRequest
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import decode_canonical_request


class RecordOrientation(TypedDict):
    """The orientation a record is drawn with after the rule: its request record key and its reverse."""

    recordKey: str
    reverseComplement: bool


class AutoOrientationReply(TypedDict):
    """Every record of the request, in drawn order."""

    records: list[RecordOrientation]


def decide_auto_orientation(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
) -> AutoOrientationReply:
    """Plan the Linear request with ``auto_orient`` and return each record's drawn orientation.

    A record the request already reverses keeps its orientation; the rule
    reverses others from the nucleotide comparisons of the request.
    """

    request = decode_canonical_request(
        request_payload, resource_paths=resource_paths, output_directory=output_directory,
    )
    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError(
            "Record orientation from comparisons requires a Linear request.",
            diagnostic={"code": "HELPER_PROTOCOL"},
        )
    plan = plan_linear_request(replace(request, options=replace(request.options, auto_orient=True)))
    return {"records": [
        {"recordKey": provenance.record_key, "reverseComplement": record_is_reverse_complemented(record)}
        for record, provenance in zip(plan.records, plan.provenance, strict=True)
    ]}


def decide_auto_orientation_json(
    request_json: str,
    resource_paths_json: str,
    output_directory: str,
) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(decide_auto_orientation(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory),
    ))


__all__ = [
    "AutoOrientationReply",
    "RecordOrientation",
    "decide_auto_orientation",
    "decide_auto_orientation_json",
]
