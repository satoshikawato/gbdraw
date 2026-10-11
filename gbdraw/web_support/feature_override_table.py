"""Worker helper of the Web's Load Feature Edits TSV (design Q4 6.4, R4).

Python owns the ``--feature_override_table`` reader and its resolution, so the
Web reads a table against the records of its committed request here.
"""

from __future__ import annotations

import json

from gbdraw.api.request_render import read_request_feature_override_table
from gbdraw.session_request_codec import decode_canonical_request


def read_feature_override_table_json(
    table_path: str, request_json: str, resource_paths_json: str, output_directory: str,
) -> str:
    """Return the table's request rows and the row numbers that name no record or feature."""
    request = decode_canonical_request(
        json.loads(str(request_json)),
        resource_paths=json.loads(str(resource_paths_json)),
        output_directory=output_directory,
    )
    unmatched: list[int] = []
    rows = read_request_feature_override_table(request, table_path, unmatched=unmatched)
    return json.dumps({"rows": [row.to_mapping() for row in rows], "unmatchedRows": unmatched})
