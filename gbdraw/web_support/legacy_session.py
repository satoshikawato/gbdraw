"""Worker entries for Web Load of an older Session.

A Session 27-44: one call returns its Session 51.

The Web shows the stored figure and makes this one call; the migrations
are Python's (:func:`gbdraw.session.drawings_document_from_legacy`). The
stored Results a Session 51 cannot hold (Sessions 27-39) come back beside
the document, by drawing ID; the Web shows them in memory and, as before,
saves none of them until one Generate (PD-OI-079). The reply also names the
file Result set each drawing took (``resultSetSources``), so Load keeps a
set's repair paths the file's own (D-15-130), and how many ``hash=``
annotation targets moved to their source feature, which Load's notice names
(``movedAnnotationTargets``, OV-2072).

A Session 46-50, which the Web reads itself: a Linear drawing whose
committed LOSATP search its saved LOSAT rows answer takes the results of
those rows, as a Session 27-44's does in its move (OV-3545), so the browser,
which has no LOSAT, redraws it before the first Generate.
"""

from __future__ import annotations

import base64
import json
from typing import TypedDict, cast

from gbdraw.session import (
    LegacyResultSetSources,
    StoredResults,
    drawings_document_from_legacy,
    mode_scoped_saved_protein_results,
)
from gbdraw.session_contract import CanonicalRequestJson, GbdrawSession51
from gbdraw.session_io import SessionResources
from gbdraw.session_migration import DraftMigrationRepair
from gbdraw.web_support.error_adapter import serialize_web_error

# The Worker operation of this entry. The error adapter reports it as
# ``unknown`` until the Web registers it with its Worker protocol (B2d-4):
# the Python and Web operation lists are checked equal.
UPGRADE_LEGACY_SESSION_OPERATION = "upgradeLegacySession"
STORE_SAVED_PROTEIN_RESULTS_OPERATION = "storeSavedProteinResults"


class LegacySessionDrawings(TypedDict):
    """The Session 51 of a Session 27-44, the stored Results it cannot hold,
    the file Result set each drawing took, the Load repairs, and the count of
    moved annotation targets."""

    document: GbdrawSession51
    storedResults: StoredResults
    resultSetSources: LegacyResultSetSources
    warnings: list[DraftMigrationRepair]
    movedAnnotationTargets: int


def legacy_session_to_drawings_json(document_json: str, stored_request_schema_json: str = "null") -> str:
    """``{"document": <Session 51>, "storedResults": {<drawing ID>: [...]},
    "resultSetSources": {<drawing ID>: ""}, "warnings": [...],
    "movedAnnotationTargets": <count>}`` for a Session 27-44 as JSON.

    ``stored_request_schema_json`` is the schema of the request the file
    stored (an integer, or ``null`` to read it from the document's request):
    Web Load's reader read the request it sends.

    A Session this step cannot move returns ``{"error": ...}`` in the shape
    of :func:`gbdraw.web_support.error_adapter.serialize_web_error`.
    """

    try:
        stored_request_schema = json.loads(stored_request_schema_json)
        if stored_request_schema is not None and type(stored_request_schema) is not int:
            raise ValueError("The stored request schema must be an integer or null.")
        legacy = drawings_document_from_legacy(
            json.loads(document_json), stored_request_schema=stored_request_schema
        )
    except Exception as exc:
        return json.dumps(
            {"error": serialize_web_error(exc, operation=UPGRADE_LEGACY_SESSION_OPERATION, stage="helper")}
        )
    document = cast(GbdrawSession51, {**legacy.document, "resources": _embedded(legacy.resources)})
    payload: LegacySessionDrawings = {
        "document": document,
        "storedResults": legacy.stored_results,
        "resultSetSources": legacy.result_set_sources,
        "warnings": list(legacy.repairs),
        "movedAnnotationTargets": legacy.moved_annotation_targets,
    }
    return json.dumps(payload)


def _embedded(resources: SessionResources) -> dict[str, dict[str, object]]:
    # The Web reads the resources of a reply as it reads a JSON Session's:
    # each with its bytes in base64 (``services/session-resource-backing.js``).
    return {
        resource_id: {**descriptor, "encoding": "base64",
                      "data": base64.b64encode(resources.read(resource_id)).decode("ascii")}
        for resource_id, descriptor in resources.descriptors.items()
    }


class SavedProteinResultRequests(TypedDict):
    """The committed requests of a Session 46-50 that took the results of their
    saved LOSAT rows, by drawing ID (the drawing's mode), and the resources they
    add, by ID, with their bytes."""

    requests: dict[str, CanonicalRequestJson]
    resources: dict[str, dict[str, object]]


def saved_protein_results_json(document_json: str) -> str:
    """``{"requests": {<drawing ID>: <request>}, "resources": {<ID>: ...}}``
    for a Session 46-50 as JSON (:func:`gbdraw.session.mode_scoped_saved_protein_results`):
    both empty when no request changes.

    A Session this step cannot read returns ``{"error": ...}`` in the shape
    of :func:`gbdraw.web_support.error_adapter.serialize_web_error`.
    """

    try:
        saved = mode_scoped_saved_protein_results(json.loads(document_json))
    except Exception as exc:
        return json.dumps(
            {"error": serialize_web_error(exc, operation=STORE_SAVED_PROTEIN_RESULTS_OPERATION, stage="helper")}
        )
    payload: SavedProteinResultRequests = {
        "requests": saved.requests,
        "resources": _embedded(saved.resources),
    }
    return json.dumps(payload)


__all__ = [
    "LegacySessionDrawings",
    "STORE_SAVED_PROTEIN_RESULTS_OPERATION",
    "SavedProteinResultRequests",
    "UPGRADE_LEGACY_SESSION_OPERATION",
    "legacy_session_to_drawings_json",
    "saved_protein_results_json",
]
