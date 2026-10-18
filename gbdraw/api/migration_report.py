"""The compatibility evidence that a Session replay returns with its render.

``gbdraw.api.session_compat`` fills it; a fresh render leaves it empty.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class SessionMigrationReport:
    """Compatibility evidence retained when a session enters the current core."""

    protein_raw_candidates: tuple[Mapping[str, Any], ...] = ()
    protein_derived_evidence: tuple[Mapping[str, Any], ...] = ()
    protein_id_map: Mapping[str, str] | None = None
    warnings: tuple[str, ...] = ()


__all__ = ["SessionMigrationReport"]
