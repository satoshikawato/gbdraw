"""Values that older Sessions and schema 1-2 requests wrote, and what they read as now.

The canonical-request decoder (:mod:`gbdraw.session_request_codec`) and the
Session migrations of a draft or a recorded command line
(:mod:`gbdraw.session_io`) read the same renames from here.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Mapping

# The feature types an older writer drew when it named none.
LEGACY_FEATURE_TYPES: tuple[str, ...] = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "repeat_region",
)
# The shapes an older writer drew these feature types as when it left the shape unset.
LEGACY_FEATURE_SHAPES: Mapping[str, str] = MappingProxyType({"repeat_region": "rectangle"})
# The Circular multi-record size modes an older writer wrote under another name.
LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES: Mapping[str, str] = MappingProxyType({"sqrt": "auto"})

__all__ = [
    "LEGACY_FEATURE_SHAPES",
    "LEGACY_FEATURE_TYPES",
    "LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES",
]
