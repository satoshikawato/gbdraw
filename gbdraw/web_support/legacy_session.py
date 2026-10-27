"""Worker entry for Web Load of a Session 27-44: one call returns its Session 49.

The Web shows the stored figure and makes this one call; the migrations
are Python's (:func:`gbdraw.session.drawings_document_from_legacy`). The
stored Results a Session 49 cannot hold (Sessions 27-39) come back beside
the document, by drawing ID; the Web shows them in memory and, as before,
saves none of them until one Generate (PD-OI-079). The reply also names the
file Result set each drawing took (``resultSetSources``), so Load keeps a
set's repair paths the file's own (D-15-130).
"""

from __future__ import annotations

import json
from typing import TypedDict

from gbdraw.session import LegacyResultSetSources, StoredResults, drawings_document_from_legacy
from gbdraw.session_contract import GbdrawSession49
from gbdraw.session_migration import DraftMigrationRepair
from gbdraw.web_support.error_adapter import serialize_web_error

# The Worker operation of this entry. The error adapter reports it as
# ``unknown`` until the Web registers it with its Worker protocol (B2d-4):
# the Python and Web operation lists are checked equal.
UPGRADE_LEGACY_SESSION_OPERATION = "upgradeLegacySession"


class LegacySessionDrawings(TypedDict):
    """The Session 49 of a Session 27-44, the stored Results it cannot hold,
    the file Result set each drawing took, and the Load repairs."""

    document: GbdrawSession49
    storedResults: StoredResults
    resultSetSources: LegacyResultSetSources
    warnings: list[DraftMigrationRepair]


def legacy_session_to_drawings_json(document_json: str, stored_request_schema_json: str = "null") -> str:
    """``{"document": <Session 49>, "storedResults": {<drawing ID>: [...]},
    "resultSetSources": {<drawing ID>: ""}, "warnings": [...]}`` for a Session
    27-44 as JSON.

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
    payload: LegacySessionDrawings = {
        "document": legacy.document,
        "storedResults": legacy.stored_results,
        "resultSetSources": legacy.result_set_sources,
        "warnings": list(legacy.repairs),
    }
    return json.dumps(payload)


__all__ = ["LegacySessionDrawings", "UPGRADE_LEGACY_SESSION_OPERATION", "legacy_session_to_drawings_json"]
