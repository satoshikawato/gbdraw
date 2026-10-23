"""Worker helper that puts the Web's ring LOSAT rows in the search frame (PD-OI-073, Owner D-28).

The Web searches the whole reference file for its Circular similarity rings
(Owner D-25) and keeps that output in its raw LOSAT cache. Before Generate
renders the request, the page converts the rows here, with the owner the CLI
planner uses (``ring_rows_in_search_frame``), against the records the request
draws. The request and its Session then store ring tables in the search frame,
as every ring table is.
"""

from __future__ import annotations

import base64
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

from gbdraw.api.request_render import resolve_request_records
from gbdraw.comparisons.circular_losat import ring_rows_in_search_frame
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import decode_canonical_request


def circular_ring_tables_in_search_frame(
    resource_ids: Sequence[str],
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
) -> dict[str, bytes]:
    """The search-frame bytes, by resource ID, of the ring tables ``resource_ids`` of the request.

    ``request_payload`` is the canonical request whose ring tables are the
    whole-file search output, and ``resource_paths`` its staged files.
    """

    if isinstance(resource_ids, str) or not all(
        isinstance(resource_id, str) and resource_id in resource_paths for resource_id in resource_ids
    ):
        raise ValidationError(
            "Ring table conversion requires staged ring table resources.",
            diagnostic={"code": "HELPER_PROTOCOL"},
        )
    request = decode_canonical_request(
        request_payload, resource_paths=resource_paths, output_directory=output_directory,
    )
    drawn_records = resolve_request_records(request).records
    return {
        resource_id: ring_rows_in_search_frame(
            Path(resource_paths[resource_id]).read_text(encoding="utf-8"), drawn_records=drawn_records,
        ).encode("utf-8")
        for resource_id in resource_ids
    }


def circular_ring_tables_in_search_frame_json(
    resource_ids_json: str,
    request_json: str,
    resource_paths_json: str,
    output_directory: str,
) -> str:
    """JSON boundary of the diagram Worker: base64 bytes and size by resource ID."""

    converted = circular_ring_tables_in_search_frame(
        json.loads(str(resource_ids_json)),
        json.loads(str(request_json)),
        json.loads(str(resource_paths_json)),
        str(output_directory),
    )
    return json.dumps({"resources": {
        resource_id: {"data": base64.b64encode(content).decode("ascii"), "size": len(content)}
        for resource_id, content in converted.items()
    }})


__all__ = ["circular_ring_tables_in_search_frame", "circular_ring_tables_in_search_frame_json"]
