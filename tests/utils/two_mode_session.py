"""A Session 46 that keeps a Result of each diagram mode (``otherModeResult``).

The Circular Gallery Session supplies the top-level set; the Linear Gallery
Session supplies the other set. Both name their GenBank source
``record-1-genbank`` with different bytes, so the Linear set's resources join
the one table under new IDs and names, as the Web writer allocates them.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

from tests.utils.session_drawing_view import catalogs_at, drawing_view, view_at
from tests.utils.session_files import read_session_payload

REPO_ROOT = Path(__file__).resolve().parents[2]
GALLERY_SESSIONS = REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions"
CIRCULAR_SESSION = GALLERY_SESSIONS / "HmmtDNA_basic_circular.gbdraw"
LINEAR_SESSION = GALLERY_SESSIONS / "lambda_basic_linear.gbdraw"
REFERENCE_FIELDS = ("resourceId", "gffResourceId", "fastaResourceId")


def _rename_references(value: Any, aliases: dict[str, str]) -> Any:
    if isinstance(value, dict):
        return {
            key: aliases.get(item, item) if key in REFERENCE_FIELDS and isinstance(item, str)
            else _rename_references(item, aliases)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_rename_references(item, aliases) for item in value]
    return value


def two_mode_session(linear: dict[str, Any] | None = None) -> dict[str, Any]:
    """The Session; ``linear`` (a one-drawing view) replaces the Linear Gallery Session as the other set."""

    circular = drawing_view(read_session_payload(CIRCULAR_SESSION))
    linear = drawing_view(read_session_payload(LINEAR_SESSION)) if linear is None else linear
    aliases = {resource_id: f"linear-{resource_id}" for resource_id in linear["resources"]}
    resources = dict(circular["resources"])
    for resource_id, descriptor in linear["resources"].items():
        resources[aliases[resource_id]] = {**descriptor, "name": f"linear-{descriptor['name']}"}
    # The 46-48 layout this helper builds; a 49 holds no otherModeResult.
    session = view_at(copy.deepcopy(circular), 48)
    session["resources"] = resources
    session["otherModeResult"] = {
        "renderRequest": _rename_references(linear["renderRequest"], aliases),
        "results": linear["results"],
        "editorState": {
            "featureCatalog": linear["editorState"]["featureCatalog"],
            "alignmentResetReceipt": None,
            "legend": {
                "originalOrder": linear["editorState"].get("legend", {}).get("originalOrder", []),
                "originalColors": linear["editorState"].get("legend", {}).get("originalColors", {}),
            },
            "originalSvgStroke": linear["editorState"].get("originalSvgStroke", {"color": None, "width": None}),
        },
        "ui": {"selectedResultIndex": 0, "generatedLegendPosition": "bottom"},
        **({"runMetadata": linear["runMetadata"]} if "runMetadata" in linear else {}),
    }
    return catalogs_at(session, 48)
