"""Web boundary of the Session 47 protein result migrator (OV-1145).

A Session up to 47 stored the protein results of a reverse-complemented Linear
record in the drawn frame. The Web converts one drawing's stored request
through this helper, the same Python owner the CLI uses, and replaces the
returned resource bytes under their stored IDs.
"""

from __future__ import annotations

import base64
import json
from typing import Any, Mapping

from gbdraw.api.request_render import plan_request
from gbdraw.api.session_compat import (
    DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION,
    decode_session_request,
    describe_unmatched_frame_item,
    drawn_frame_protein_results_to_search_frame,
    search_frame_protein_resources,
)
from gbdraw.exceptions import ValidationError
from gbdraw.session import SessionConversionError


def stored_protein_results_to_search_frame(
    payload: Mapping[str, Any],
    resource_paths: Mapping[str, str],
    output_directory: str,
) -> dict[str, bytes]:
    """The new bytes, by stored resource ID, of the drawing's protein results in the search frame.

    ``payload`` is ``{"version", "renderRequest"}`` of one drawing and
    ``resource_paths`` the staged files of its request. A result matching
    neither strand of its record raises SessionConversionError; nothing is
    returned for a drawing already in the search frame or a Session after 47.
    """

    version = payload.get("version")
    request_payload = payload.get("renderRequest")
    if isinstance(version, bool) or not isinstance(version, int) or not isinstance(request_payload, Mapping):
        raise ValidationError("Stored protein result conversion requires a Session version and request.")
    if version > DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION:
        return {}
    request = decode_session_request(
        version,
        request_payload,
        resource_paths=resource_paths,
        output_directory=output_directory,
    )
    converted, unmatched = drawn_frame_protein_results_to_search_frame(plan_request(request), request)
    if unmatched:
        raise SessionConversionError(
            f"Session protein results cannot be converted: {describe_unmatched_frame_item(unmatched[0])}."
        )
    if converted is request:
        return {}
    return search_frame_protein_resources(request_payload, request, converted)


def stored_protein_results_to_search_frame_json(
    payload_json: object,
    resource_paths_json: object = "{}",
    output_directory: str = "/tmp",
) -> str:
    """JSON boundary used by the existing lazy diagram Worker: base64 bytes and size by resource ID."""

    converted = stored_protein_results_to_search_frame(
        json.loads(str(payload_json)), json.loads(str(resource_paths_json)), output_directory,
    )
    return json.dumps({"resources": {
        resource_id: {"data": base64.b64encode(content).decode("ascii"), "size": len(content)}
        for resource_id, content in converted.items()
    }})


__all__ = ["stored_protein_results_to_search_frame", "stored_protein_results_to_search_frame_json"]
