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
    # Stored protein results of a reversed record that matched neither frame
    # (OV-1145): a replay draws them as stored; writing them as Session 48 fails.
    unconverted_protein_results: tuple[str, ...] = ()


__all__ = ["SessionMigrationReport"]
