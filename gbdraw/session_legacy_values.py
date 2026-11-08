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
# The depth and GC content axes of a config written before they had a label style (v0.14.5 and earlier) drew a
# number on each large tick: the config sections, by name, their override paths, and the style they read as
# (D-16, D-24, 2026-11-08). A stored full config and stored config overrides without the style read so; every
# writer of config overrides since states it.
LEGACY_SCALAR_AXIS_LABEL_SECTIONS: tuple[str, ...] = ("depth", "gc_content")
SCALAR_AXIS_LABEL_PATHS: tuple[str, ...] = tuple(
    f"objects.{section}.axis_labels" for section in LEGACY_SCALAR_AXIS_LABEL_SECTIONS
)
LEGACY_SCALAR_AXIS_LABELS = "numbered"

__all__ = [
    "LEGACY_SCALAR_AXIS_LABELS",
    "LEGACY_SCALAR_AXIS_LABEL_SECTIONS",
    "SCALAR_AXIS_LABEL_PATHS",
    "LEGACY_FEATURE_SHAPES",
    "LEGACY_FEATURE_TYPES",
    "LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES",
]
