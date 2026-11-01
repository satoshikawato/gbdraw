"""Worker helper behind the Web's "Orient records from comparisons" button.

The page sends the request of the last Generate, which holds the comparisons
in the search frame as the Session keeps them. The button sets the Reverse
complement of every record below the top row from the comparisons: the
helper plans the request with ``auto_orient`` as ``--auto_orient`` does on the
command line, with the reverse flags of the records below the top row cleared
in its input, so the rule
(:func:`gbdraw.linear_comparison.auto_orient_reversed_records`) decides them
in its one owner against the top row as the user drew it. It returns the
orientation each record is then drawn with, by record key. (On the command
line an explicit reverse still wins; the rule is unchanged, only its input
differs.) The page writes the records whose orientation changes into their
Reverse complement checkboxes as one History step; the next Generate draws
them.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import replace
from typing import TypedDict

from gbdraw.api.record_planning import ResolvedRecordProvenance, record_is_reverse_complemented
from gbdraw.api.request_render import LinearRequestPlan, plan_linear_request
from gbdraw.api.requests import LinearDiagramRequest, RecordCardinality, RecordInput
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.layout.record_placement import resolve_record_row_positions
from gbdraw.session_request_codec import decode_canonical_request


class RecordOrientation(TypedDict):
    """The orientation a record is drawn with after the rule: its request record key and its reverse."""

    recordKey: str
    reverseComplement: bool


class AutoOrientationReply(TypedDict):
    """Every record of the request, in drawn order."""

    records: list[RecordOrientation]


def _forward(item: RecordInput) -> RecordInput:
    """The record input with its reverse flags (presentation and region) cleared."""

    return replace(
        item,
        presentation=replace(item.presentation, reverse_complement=False),
        region=None if item.region is None else replace(item.region, reverse_complement=False),
    )


def _is_reversed(item: RecordInput) -> bool:
    return item.presentation.reverse_complement or bool(item.region and item.region.reverse_complement)


def _record_input_of(item: RecordInput, provenance: ResolvedRecordProvenance) -> RecordInput:
    """The input that draws only this record: an input that draws all its file's records, narrowed to this one."""

    if item.cardinality is not RecordCardinality.ALL:
        return item
    index = provenance.source_record_index
    return replace(
        item,
        selector=RecordSelector(raw=f"#{index + 1}", record_id=None, record_index=index),
        cardinality=RecordCardinality.EXACTLY_ONE,
        record_key=provenance.record_key,
    )


def _with_rows_below_top_unreversed(request: LinearDiagramRequest, plan: LinearRequestPlan) -> LinearDiagramRequest | None:
    """One input per drawn record, the reverse flags cleared below the top row; ``None`` when nothing is cleared.

    The rows are the planner's (``resolve_record_row_positions`` on the plan's
    records and layout). An input that draws several records is narrowed to
    one input per record, so a record below the top row is cleared also when
    another record of its file is in the top row.
    """

    rows = resolve_record_row_positions(
        plan.records, plan.layout.multi_record_positions if plan.layout is not None else None
    )[1]
    below = [row != 0 and _is_reversed(request.records[provenance.input_index])
             for provenance, row in zip(plan.provenance, rows, strict=True)]
    if not any(below):
        return None
    # Every other input draws exactly one record, so the drawn records name the inputs in order.
    narrowed = (_record_input_of(request.records[provenance.input_index], provenance) for provenance in plan.provenance)
    return replace(request, records=tuple(
        _forward(item) if clear else item for item, clear in zip(narrowed, below, strict=True)
    ))


def decide_auto_orientation(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
) -> AutoOrientationReply:
    """Decide the orientation of every record below the top row from the request's comparisons.

    The records of the top row keep their orientation, as the rule keeps them.
    A record below it is asked with its reverse flags cleared, so the rule
    decides it against the row above as drawn, also when the request reverses
    it; a record that the comparisons do not decide is drawn forward. The
    comparisons are in the search frame, which no reverse flag changes, so
    they are read as given.
    """

    request = decode_canonical_request(
        request_payload, resource_paths=resource_paths, output_directory=output_directory,
    )
    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError(
            "Record orientation from comparisons requires a Linear request.",
            diagnostic={"code": "HELPER_PROTOCOL"},
        )
    oriented = replace(request, options=replace(request.options, auto_orient=True))
    plan = plan_linear_request(oriented)
    unreversed = _with_rows_below_top_unreversed(oriented, plan)
    if unreversed is not None:
        # Only a reversed record below the top row changes the rule's input; plan again without its flags.
        plan = plan_linear_request(unreversed)
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
